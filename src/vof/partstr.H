#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef struct {
  int Np; // number of particles per string
  int Ns; // number of strings (cross sections) per cell
  Real Hp; // length of string relative to cell size
  Real eps; // threshold for convergence
  int itermax; // maximum number of iterations
  Real eta; // relaxation factor
  bool nohf; // disable height functions, always use particles
  int stat_sh; // number of cells with height functions
  int stat_sc; // number of cells with particles
} Partstr;

// Transformation of coordinates.
//    |n
//    |
//    |
//    o-----t
//   /
//  /u
typedef struct {
  XDim3 o; // origin
  XDim3 t, n, u; // orthonormal positive-oriented base
} Trans;

const int kMaxSection = 125 * 2; // maximum number of endpoints
                                 // from neighbor cells
                                 // (5x5x5 stencil and 2 endpoints in each)
const int kMaxFacet = 12; // maximum number of vertices per facet

const int kMaxNp = 31; // maximum value of partstr_conf.Np


#if AMREX_SPACEDIM==2
#define kNs 1
#else
#define kNs 3
#endif
static Partstr partstr_conf = {7, kNs, 4., 1e-5, 20, 0.5, false};
#undef kNs

// Unit vector at angle ph
static XDim3 Unit(Real ph) {
  XDim3 p;
  p.x = cos(ph);
  p.y = sin(ph);
  p.z = 0;
  return p;
}

// Rotate planar vector e by planar unit vector de
static XDim3 Rotate(XDim3 e, XDim3 de) {
  XDim3 p;
  p.x = e.x * de.x - e.y * de.y;
  p.y = e.x * de.y + e.y * de.x;
  p.z = 0;
  return p;
}

// Rotate planar vector e by planar unit vector de in the negative direction
static XDim3 Rotatem(XDim3 e, XDim3 de) {
  XDim3 p;
  p.x = e.x * de.x + e.y * de.y;
  p.y = -e.x * de.y + e.y * de.x;
  p.z = 0;
  return p;
}

// Third component of cross product
static double Cross3(XDim3 a, XDim3 b) {
  return a.x * b.y - a.y * b.x;
}

static XDim3 Cross(XDim3 a, XDim3 b) {
  XDim3 r;
  r.x = a.y * b.z - a.z * b.y;
  r.y = a.z * b.x - a.x * b.z;
  r.z = a.x * b.y - a.y * b.x;
  return r;
}

static XDim3 Add(XDim3 a, XDim3 b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

static XDim3 Sub(XDim3 a, XDim3 b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

static XDim3 Mul(XDim3 a, Real k) {
  a.x *= k;
  a.y *= k;
  a.z *= k;
  return a;
}

static XDim3 Div(XDim3 a, Real k) {
  a.x /= k;
  a.y /= k;
  a.z /= k;
  return a;
}

static XDim3 Coord(Real x, Real y, Real z) {
  XDim3 r = {x, y, z};
  return r;
}

static XDim3 Zero() {
  XDim3 r;
  r.x = 0;
  r.y = 0;
  r.z = 0;
  return r;
}

static Real Dot(XDim3 a, XDim3 b) {
  Real s = 0;
  s += a.x * b.x;
  s += a.y * b.y;
  s += a.z * b.z;
  return s;
}

static Real Sqnorm(XDim3 a) {
  return Dot(a, a);
}

static Real Norm(XDim3 a) {
  return sqrt(Sqnorm(a));
}

static Real Norm1(XDim3 a) {
  return fabs(a.x) + fabs(a.y) + fabs(a.z);
}

static Real NormMax(XDim3 a) {
  Real r = fabs(a.x);
  if (fabs(a.y) > r) {
    r = fabs(a.y);
  }
  if (fabs(a.z) > r) {
    r = fabs(a.z);
  }
  return r;
}

static Real Sqdist(XDim3 a, XDim3 b) {
  return Sqnorm(Sub(a, b));
}

static Real Dist(XDim3 a, XDim3 b) {
  return sqrt(Sqdist(a, b));
}

static Real Dotv(int np, const XDim3* aa, const XDim3* bb) {
  Real s = 0;
  int i;
  for (i = 0; i < np; ++i) {
    s += Dot(aa[i], bb[i]);
  }
  return s;
}

// Positions of particles.
// p: central particle
// ph: orientation angle (phi)
// th: bending angel (theta)
// np: number of particles
// hp: distance between particles
// Output:
// xx: array of length np
static void XX(XDim3 p, Real ph, Real th, int np, Real hp, XDim3* xx) {
  int c = np / 2, j;
  xx[c] = p;
  XDim3 ep = Mul(Unit(ph + 0.5 * th), hp);
  XDim3 em = Mul(Unit(ph - 0.5 * th), -hp);
  XDim3 de = Unit(th);
  for (j = 0; j < c; ++j) {
    xx[c + j + 1] = Add(xx[c + j], ep);
    xx[c - j - 1] = Add(xx[c - j], em);
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
  }
}

// Returns nearest to x point on line segment [a, b].
static XDim3 Nearest(XDim3 a, XDim3 b, XDim3 x) {
  b = Sub(b, a);
  x = Sub(x, a);
  Real q = Dot(b, x) / Dot(b, b);
  q = CLAMP(q, 0., 1.);
  a.x += b.x * q;
  a.y += b.y * q;
  a.z = 0;
  return a;
}

// Derivative dX/dph
static void DxDph(XDim3 p, Real ph, Real th, int np, Real hp, XDim3* xx) {
  (void)p;
  int c = np / 2, j;
  xx[c] = Zero();
  XDim3 ep = Mul(Unit(ph + 0.5 * th + PI * 0.5), hp);
  XDim3 em = Mul(Unit(ph - 0.5 * th + PI * 0.5), -hp);
  XDim3 de = Unit(th);
  for (j = 0; j < c; ++j) {
    xx[c + j + 1] = Add(xx[c + j], ep);
    xx[c - j - 1] = Add(xx[c - j], em);
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
  }
}

// Derivative dX/dth
static void DxDth(XDim3 p, Real ph, Real th, int np, Real hp, XDim3* xx) {
  (void)p;
  int c = np / 2;
  xx[c] = Zero();
  XDim3 ep = Mul(Unit(ph + 0.5 * th + PI * 0.5), hp);
  XDim3 em = Mul(Unit(ph - 0.5 * th + PI * 0.5), hp);
  XDim3 de = Unit(th);
  for (int j = 0; j < c; ++j) {
    Real jp = j + 0.5;
    xx[c + j + 1] = Add(xx[c + j], Mul(ep, jp));
    xx[c - j - 1] = Add(xx[c - j], Mul(em, jp));
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
  }
}

// Forces on particles.
// np: number of particles
// xx: positions
// nl: number of points in ll
// ll: flat array of endpoints of line segments
// eta: relaxation factor
// k: curvature
// Output:
// ff: forces
static void FF(int np, const XDim3* xx, int nl, 
    const XDim3 * ll, Real eta, Real k,XDim3* ff){     
  if (nl == 0) {
    for (int i = 0; i < np; ++i) {
      ff[i] = Zero();
    }
    return;
  }
  for (int i = 0; i < np; ++i) {
    const XDim3 x = xx[i];
    XDim3 pm;
    Real distm = -1;

    for (int l = 0; l < nl; l += 2) {
      const XDim3 p = Nearest(ll[l], ll[l + 1], x);
      const Real dist = Sqdist(x, p);
      if (distm < 0 || dist < distm) {
        pm = p;
        distm = dist;
      }
    }
    ff[i] = Mul(Sub(pm, x), eta);
  }
}

// Iteration of evolution of particles.
// p_,ph_,th_: current configuration
// np: number of particles
// hp: distance between particles
// ff: forces (output of F())
// xx: buffer of size np
// Output:
// p_,ph_,th_: new configuration
// ff, xx: modified
// Returns maximum absolute difference.
static double Iter(XDim3& p_, Real& ph_, Real& th_, 
                  int np, Real hp, XDim3* ff, XDim3 * xx) 
{    
  XDim3 p = p_;
  Real ph = ph_;
  Real th = th_;

  XDim3 tt[kMaxNp];
  int c = np / 2;

  XX(p, ph, th, np, hp, xx);

  // Skip correction of p, assume that p is initialized on a line segment
  // // correct p
  // p = Add(p, ff[c]);
  // X(p, ph, th, np, hp, tt);
  // for (int i = 0; i < np; ++i) {
  //   XDim3 dx = Sub(tt[i], xx[i]);
  //   xx[i] = Add(xx[i], dx);
  //   ff[i] = Sub(ff[i], dx);
  // }

  // correct phi
  DxDph(p, ph, th, np, hp, tt);
  ph += Dotv(np, tt, ff) / Dotv(np, tt, tt);
  XX(p, ph, th, np, hp, tt);
  for (int i = 0; i < np; ++i) {
    XDim3 dx = Sub(tt[i], xx[i]);
    xx[i] = Add(xx[i], dx);
    ff[i] = Sub(ff[i], dx);
  }

  // correct theta
  DxDth(p, ph, th, np, hp, tt);
  th += Dotv(np, tt, ff) / Dotv(np, tt, tt);
  XX(p, ph, th, np, hp, tt);

  Real r = 0;
  for (int i = 0; i < np; ++i) {
    r = MAX(r, NormMax(Sub(xx[i], tt[i])));
  }
  p_ = p;
  ph_ = ph;
  th_ = th;
  return r;
}

static Real Curv(Real hp, Real th) {
  return 2. * sin(th * 0.5) / hp;
}

static void Swap(XDim3& a, XDim3& b) {
  XDim3 t;
  t = a;
  a = b;
  b = t;
}

static XDim3 Abs(XDim3 p) {
  p.x = fabs(p.x);
  p.y = fabs(p.y);
  p.z = fabs(p.z);
  return p;
}

static Real Get(XDim3 p, int i) {
  return i == 0 ? p.x : i == 1 ? p.y : p.z;
}

// Index of minimal component.
static int Argmin(XDim3 p) {
  int im = 0;
  for (int i = 1; i < 3; ++i) {
    if (Get(p, i) < Get(p, im)) {
      im = i;
    }
  }
  return im;
}

// Sets component i of point p to a
static void Set(XDim3& p, int i, Real a) {
  if (i == 0) {
    p.x = a;
  } else if (i == 1) {
    p.y = a;
  } else {
    p.z = a;
  }
}

// Returns 3D base aligned with mesh directions.
// n: unit vector
// Output:
// *t, *u: vectors to form orthonormal base <t,u,n>
static void GetBase(XDim3 n, XDim3& t, XDim3& u) {
  int i = Argmin(Abs(n));
  XDim3 e = Zero();
  Set(e, i, 1);
  n = Div(n, Norm(n));
  t = Cross(n, e);
  t = Div(t, Norm(t));
  u = Cross(t, n);
}

// Transformation from global to local coordinates.
// p: global coordinates
// o: origin
// t,n,u: orthonormal base
// Output:
// l: local coordinates of p in <o,t,n,u>
static XDim3 GlbToLoc(XDim3 p, const Trans& w) {
  p = Sub(p, w.o);
  XDim3 l = {Dot(w.t, p), Dot(w.n, p), Dot(w.u, p)};
  return l;
}

// Transformation from local to global coordinates.
// l: local coordinates in <t,n,u,o>
// t,n,u: orthonormal base
// o: origin
// Output:
// *p: global coordinates
static XDim3 LocToGlb(XDim3 l, const Trans& w) {
  XDim3 p = w.o;
  p = Add(p, Mul(w.t, l.x));
  p = Add(p, Mul(w.n, l.y));
  p = Add(p, Mul(w.u, l.z));
  return p;
}

// Set z-component to zero if 2D.
static void SetZeroZ(XDim3 & p) {
#if AMREX_SPACEDIM==2 
  p.z = 0;
#endif
}
#define EPS 1e-4 
/**
 * vof_facet:
 * @dx: grid size in x-,y-, z-direction.
 * @q: cell center. 
 * @p: return the vertices coordinates (of size 2 in 2D and 6 in 3D)
 * @m:  the normal to the interface.
 
 * Fills @p with the coordinates of points defining the
 * VOF-reconstructed interface facet defined by @t.
 *
 * Returns: the number of points defining the facet.
 */
int vof_facet (GpuArray<Real, AMREX_SPACEDIM> dx,
               XDim3 q, XDim3 m, Real alpha, XDim3 * p)
         
{
 int c,n =0;
 XDim3 o=q; 
#if AMREX_SPACEDIM==2
  Real x, y;
  
  if (fabs (m.y) > EPS) {
    y = (alpha - m.x)/m.y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x + dx[0]/2.; p[n].y = q.y + dx[1]*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m.x) > EPS) {
    x = (alpha - m.y)/m.x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + dx[0]*(x - 0.5); p[n].y = q.y + dx[1]/2.; p[n++].z = 0.;
    }
  }
  if (fabs (m.y) > EPS) {
    y = alpha/m.y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x - dx[0]/2.; p[n].y = q.y + dx[1]*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m.x) > EPS) {
    x = alpha/m.x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + dx[0]*(x - 0.5); p[n].y = q.y - dx[1]/2.; p[n++].z = 0.;
    }
  }
  AMREX_ASSERT (n <= 2);  
#else /* 3D */ 
  Real max = fabs (m.x);
  c = 0;
  if (fabs (m.y) > max) {
    max = fabs (m.y);
    c = 1;
  }
  if (fabs (m.z) > max)
    c = 2;
  q.x -= dx[0]/2.; q.y -= dx[1]/2.; q.z -= dx[2]/2.;
  (&q.x)[c] += dx[c]*alpha/(&m.x)[c];

  int d[12];
  /* array of node coordinates for a cut face */
  Array<XDim3, 12>  nodecutface;
  n = cut_cube_vertices (o, dx, q, m, nodecutface, d);
  for (c=0; c<n; ++c)
    p[c]=nodecutface[c];      
    
  AMREX_ASSERT (n <= 6);
#endif /* 3D */
  return n; 
}

// Cross-section of interface from neighbor cells in plane coordinates.
// point: center cell
// w: transformation
// ll: buffer for at least kMaxSection more points
// *nl: current size of ll
// Output:
// ll: appended by local coordinates of endpoints,
//     p = o + t*l.x  + t*l.y ,  [pa,pb] is one line segment
// *nl: new size of ll
static void Section(Array<int,3> & cell, const Trans & w, XDim3* ll, int& nl,
                    GpuArray<Real, AMREX_SPACEDIM> dx,
                    GpuArray<Real, AMREX_SPACEDIM> problo,                   
                    Array4<Real const> const & vof,
                    Array4<Real const> const & mv,
                    Array4<Real const> const & alpha) 
{
#if AMREX_SPACEDIM==3 /* 3D */
 for (int dk = -2; dk <= 2; ++dk) {
#else /* 2D */
 int dk=0;    
#endif     
 for (int dj = -2; dj <= 2; ++dj) {
 for (int di = -2; di <= 2; ++di) {
   if (dk==0&&dj==0&&di==0)
       continue;
   Array<int,3> ncell={cell[0]+di,cell[1]+dj, cell[2]+dk}; 
   Real fvol = vof(ncell[0],ncell[1],ncell[2]);                          
   if (!CELL_IS_FULL (fvol)) {
      int q = 0; // number of intersections found
      XDim3 rn;// neighbor cell center
      for (int c = 0; c < AMREX_SPACEDIM; c++)
        (&rn.x)[c] = problo[c] + dx[c]*(ncell[c] + .5); 
#if AMREX_SPACEDIM==3    
      // if cell intersects plane
      if (fabs(Dot(w.u, Sub(w.o, rn))) <= dx[0] * Norm1(w.u) * 0.5) {
#endif          
        // normal to facet
        XDim3 m={mv(ncell[0],ncell[1],ncell[2],0),
                 mv(ncell[0],ncell[1],ncell[2],1),
                 mv(ncell[0],ncell[1],ncell[2],2)};
        Real al = alpha(ncell[0],ncell[1],ncell[2],0);         
        XDim3  pp[kMaxFacet];
        int nf = vof_facet(dx, rn, m, al, pp);
        /*Print()<<"vof_facet"<<nf<<cell[0]<<" "<<cell[1]<<" "<<cell[2]<<", "
                                 << ncell[0]<<" "<<ncell[1]<<" "<<ncell[2]<<", "
                                 <<fvol<<"\n";
                                  exit(1);*/
#if AMREX_SPACEDIM==2 /* 2D */
        if (nf == 2) {
          XDim3 l = GlbToLoc(pp[0], w);
          XDim3 lb = GlbToLoc(pp[1], w);  
          XDim3 dl = Sub(l, lb);
          Real mt = Dot(w.t, m);
          Real mn = Dot(w.n, m);
          Real dt = dl.x;
          Real dn = dl.y;
          if (Cross3(Coord(mt, mn, 0), Coord(dt, dn, 0)) < 0) 
            Swap(l, lb);
          if (Sqnorm(dl) > 0) {
            ll[nl] = l;
            ll[nl + 1] = lb;
            nl += 2;
          }
        }
#else   /* 3D */             
        AMREX_ASSERT(nf <= kMaxFacet);
        // find two intersections with facet edge and plane
        for (int i = 0; i < nf && q < 2; ++i) {
          int ib = (i + 1) % nf;
          XDim3 l = GlbToLoc(pp[i], w);
          XDim3 lb = GlbToLoc(pp[ib], w);
          // intersection with l.z=0:
          //   0 = l.z * s + lb.z * (1 - s)
          if (l.z * lb.z < 0) { // line crosses the plane
            Real s = lb.z / (lb.z - l.z);
            ll[nl + q].x = l.x * s + lb.x * (1 - s);
            ll[nl + q].y = l.y * s + lb.y * (1 - s);
            ll[nl + q].z = 0.;
            if (++q == 2) {
              XDim3 dl = Sub(ll[nl + 1], ll[nl]);
              Real mt = Dot(w.t, m); // local coordinates of normal
              Real mn = Dot(w.n, m);
              Real dt = dl.x; // local coordinates of line segment
              Real dn = dl.y;
              // ensure positive orientation with normal m
              if (Cross3(Coord(dt, dn, 0), Coord(mt, mn, 0)) < 0)
                Swap(ll[nl], ll[nl + 1]);
              if (Sqnorm(dl) > 0) 
                nl += 2;
            }
          }
        }    
      }
#endif      
   }    
 }
 }
#if AMREX_SPACEDIM==3     
 }
#endif 
}

// Curvature of a set line segments.
// ll: flat array endpoints of line segments
// nl: size of ll
// delta: cell size
// Output:
// res_: difference at last iteration
// it_: number of iterations
static double GetLinesCurv(XDim3 *ll, int nl, Real delta, const Partstr* conf,
                           Real& res_, int& it_, const Trans& w)
{
  if (nl >= 4) { // require at least two segments
    const int Np = conf->Np;
    const Real eta = conf->eta;
    const int itermax = conf->itermax;
    const Real hp = conf->Hp * delta / (Np - 1);

    XDim3 xx[kMaxNp]; // positions
    XDim3 ff[kMaxNp]; // forces
    XDim3 p = {0., 0, 0.};
    Real ph = 0.;
    Real th = 0.;
    Real kappa = 0;
    Real res = 0;
    int it = 0;
    for (it = 0; it < itermax; ++it) {
      XX(p, ph, th, Np, hp, xx);
      FF(Np, xx, nl, ll, eta, kappa, ff);

      kappa = Curv(hp, th);
      res = Iter(p, ph, th, Np, hp, ff, xx);
      if (res / (eta * delta) < conf->eps) {
        break;
      }
    }
    res_ = res;
    it_ = it;
    // write the file of particle coordinates
    FILE* o;
    o = fopen("particle.txt", "w");
    fprintf(o, " x,y,z\n");
    
    for (int i = 0; i < Np; ++i){ 
        xx[i] = LocToGlb(xx[i], w);
    }

    for (int i = 0; i < Np; ++i) {
      fprintf(o, "%g,%g,%g", xx[i].x, xx[i].y, xx[i].z);
      fprintf(o, "\n");
    }
    fclose(o); 
    return kappa;
  }
  return 0;
}

// Curvature of cross section by plane through a.n,a.t and point a.o
// point: target cell
// nn: normal
// w: transformation to local coordinates
static Real GetCrossCurv(Array<int,3> & cell, const Trans & w, const Partstr* conf,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         GpuArray<Real, AMREX_SPACEDIM> problo,
                         Array4<Real const> const & vof,
                         Array4<Real const> const & mv,
                         Array4<Real const> const & alpha) 
{
  XDim3 ll[kMaxSection];
  int nl = 0;
  Section(cell, w, ll, nl, dx, problo, vof, mv, alpha);
  Real res;
  int it;
  //fix me: grid size
  return GetLinesCurv(ll, nl, dx[0], conf, res, it, w);

}


// Transformation b rotated at angle.
// s: index of cross section
// Ns: number of cross sections
static Trans GetSectionTrans(int s, int Ns, const Trans & b) {
  const Real g = PI * s / Ns;
  Trans w = b;
  w.t = Add(Mul(b.t, cos(g)), Mul(b.u, sin(g)));
  w.u = Cross(w.t, w.n);
  return w;
}

// Mean curvature over multiple cross sections by planes rotated around b.n
static Real GetMeanCurv(Array<int,3> &cell, const Trans & b, const Partstr* conf,
                        GpuArray<Real, AMREX_SPACEDIM> dx,
                        GpuArray<Real, AMREX_SPACEDIM> problo,
                        Array4<Real const> const & vof,
                        Array4<Real const> const & mv,
                        Array4<Real const> const & alpha) 
{
  Real ksum = 0;
  const int Ns = conf->Ns;
  for (int s = 0; s < Ns; ++s) {
    const Trans w = GetSectionTrans(s, Ns, b);
    ksum += GetCrossCurv(cell, w, conf, dx, problo, vof, mv, alpha);
  }
  return ksum / Ns;
}

void CheckConf(const Partstr* conf) {
  if (!(conf->Np <= kMaxNp)) {
    fprintf(
        stderr, "error: Too many particles per string Np=%d > %d\n", conf->Np,
        kMaxNp);
    exit(1);
  }
  if (!(conf->Np % 2 == 1)) {
    fprintf(stderr, "error: Np=%d is not an odd number\n", conf->Np);
    exit(1);
  }
}

// Method of connected particles for curvature estimation (January 13, 2020)
// based on paper
// Karnakov P., Litvinov S. and Koumoutsakos P.,
// A hybrid particle volume-of-fluid method for curvature estimation in
// multiphase flows. International Journal of Multiphase Flow. 2020
// https://doi.org/10.1016/j.ijmultiphaseflow.2020.103209 
Real partstr_curvature (int i,int j,int k, 
                        GpuArray<Real, AMREX_SPACEDIM> dx,
                        GpuArray<Real, AMREX_SPACEDIM> problo, 
                        Array4<Real const> const & vof,
                        Array4<Real const> const & mv,
                        Array4<Real const> const & alpha)
{

  Partstr* conf = &partstr_conf;
  CheckConf(conf);
  XDim3 m={mv(i,j,k,0),mv(i,j,k,1),mv(i,j,k,2)},p;
  //find transformation to local coordinates at the interface.
  Trans b;   
  plane_area_center (m, alpha(i,j,k,0),p);
  for (int c = 0; c < AMREX_SPACEDIM; c++)
    (&p.x)[c] = problo[c] + dx[c]*((&p.x)[c] + (c<1?i:c<2?j:k)); 
  SetZeroZ(p);

  b.o = p;
  b.n = Div(m, Norm(m));
  GetBase(b.n, b.t, b.u);   
  Array<int,3> cell={i,j,k};
  Real kappa = GetMeanCurv(cell, b, conf, dx, problo, vof,mv, alpha);
  return kappa;     
}


#define PARABOLA_FIT_CENTER_WEIGHT .1
#define PARABOLA_SIMPLER 0
// Function to calculate the norm of the vector
Real vector_norm(const VofVector& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Function to normalize the vector
void vector_normalize(VofVector& v) {
    Real norm = vector_norm(v);
    if (norm > 0.0) {
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }
}
// Function to calculate the cross product of two vectors
void vector_cross(VofVector& C, const VofVector& A, const VofVector& B) {
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
}

/**
 * matrix_inverse:
 * @m: a square matrix.
 * @n: size of the matrix.
 * @pivmin: the minimum value of the pivoting coefficient.
 *
 * Replaces @m with its inverse.
 *
 * Returns: 0. if the inversion encounters a pivot coefficient smaller
 * than or equal to @pivmin (i.e. @m is non-invertible), the minimum
 * absolute value of the pivoting coefficient otherwise.
 */
Real matrix_inverse (Real ** m, int n, Real pivmin)
{
  std::vector<int> indxc(n), indxr(n), ipiv(n, -1);
  int i, icol = 0, irow = 0, j, k, l, ll;
  Real big, dum, pivinv, temp, minpiv = VOF_NODATA;
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
 if (n>3){
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
    for (k = 0; k < n; k++) {
      if (ipiv[k] == -1) {
        if (fabs (m[j][k]) >= big) {
          big = fabs (m[j][k]);
          irow = j;
          icol = k;
        }
      }
    }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
    SWAP (m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin) {
      return 0.;
    }
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
    dum = m[ll][icol];
    m[ll][icol] = 0.0;
    for (l = 0; l < n; l++)
      m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
        SWAP (m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
 }
 else if (n==3){

  Real det = (m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) - 
     m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) + 
     m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]));
  if (det == 0.0)
    return 0.;

  Array<Array<Real,3>,3> mi;

  mi[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det; 
  mi[0][1] = (m[2][1]*m[0][2] - m[0][1]*m[2][2])/det;
  mi[0][2] = (m[0][1]*m[1][2] - m[1][1]*m[0][2])/det; 
  mi[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det; 
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])/det; 
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])/det; 
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])/det; 
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])/det; 
  mi[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det; 
  for (i=0;i<3;i++)
      for (j=0;j<3;j++)
          m[i][j]=mi[i][j];
  return 1.;          
 }
 return 0.; 
}


typedef struct {
  VofVector o;  
  Real ** M;
#if AMREX_SPACEDIM==2 /* y = a[0]*x^2 + a[0]*x + a[1] */
  VofVector m, rhs, a;
#else /* 3D */  
# if PARABOLA_SIMPLER /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y */ 
  VofVector rhs, a;
# else /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y + a[3]*x + a[4]*y + a[5] */
  Real rhs[6], a[6];
# endif
  VofVector t[3];
#endif  /* 3D */ 
} ParabolaFit;

static void parabola_fit_init (ParabolaFit & p, XDim3 & o, XDim3 & m)
{
  p.o[0] = o.x; p.o[1] = o.y; p.o[2] = o.z;
  int nn = 3;
#if AMREX_SPACEDIM==2 /* 2D */ 
  p.m[0] = m.x; p.m[1] = m.y; p.m[2] = 0.;
  vector_normalize (p.m);
#else /* 3D */   
  Real max;
  VofVector nx = {0., 0., 0.}, ny, nz;
  int d = 0;
  nz[0] = m.x; nz[1] = m.y; nz[2] = m.z;
  vector_normalize (nz);
  max = nz[0]*nz[0];
  /* build a vector orthogonal to nz */
  if (nz[1]*nz[1] > max) { max = nz[1]*nz[1]; d = 1; }
  if (nz[2]*nz[2] > max) d = 2;
  switch (d) {
   case 0: nx[0] = - nz[2]/nz[0]; nx[2] = 1.0; break;
   case 1: nx[1] = - nz[2]/nz[1]; nx[2] = 1.0; break;
   case 2: nx[2] = - nz[0]/nz[2]; nx[0] = 1.0; break;
  }
  vector_normalize (nx);

  /* build a second vector orthogonal to nx and nz */
  vector_cross (ny, nz, nx);

  /* transformation matrix from (i,j,k) to (nx, ny, nz) */
  p.t[0][0] = nx[0]; p.t[0][1] = nx[1]; p.t[0][2] = nx[2];
  p.t[1][0] = ny[0]; p.t[1][1] = ny[1]; p.t[1][2] = ny[2];
  p.t[2][0] = nz[0]; p.t[2][1] = nz[1]; p.t[2][2] = nz[2];
  
# if !PARABOLA_SIMPLER
    nn=6;
# endif
#endif /* 3D */
  p.M = new Real*[nn];
  for (int i = 0; i < nn; ++i) {
   p.M[i] = new Real[nn]();
   p.rhs[i] =0.;
  }  

}

static void parabola_fit_add (ParabolaFit & p, VofVector const & m, Real w)
{
#if AMREX_SPACEDIM==2
  Real x1 = m[0] - p.o[0];
  Real y1 = m[1] - p.o[1];
  Real x = p.m[1]*x1 - p.m[0]*y1;
  Real y = p.m[0]*x1 + p.m[1]*y1;
  Real x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p.M[0][0] += x4;
  p.M[1][0] += x3; p.M[1][1] += x2;
  p.M[2][1] += w*x; p.M[2][2] += w;
  p.rhs[0] += x2*y;
  p.rhs[1] += w*x*y;
  p.rhs[2] += w*y;
#else /* 3D */
  Real x1 = m[0] - p.o[0];
  Real y1 = m[1] - p.o[1];
  Real z1 = m[2] - p.o[2];
  Real x = p.t[0][0]*x1 + p.t[0][1]*y1 + p.t[0][2]*z1;
  Real y = p.t[1][0]*x1 + p.t[1][1]*y1 + p.t[1][2]*z1;
  Real z = p.t[2][0]*x1 + p.t[2][1]*y1 + p.t[2][2]*z1;
  Real x2 = x*x, x3 = x2*x, x4 = x3*x;
  Real y2 = y*y, y3 = y2*y, y4 = y3*y;
# if PARABOLA_SIMPLER
  p.M[0][0] += w*x4;
  p.M[1][0] += w*x2*y2; p.M[1][1] += w*y4;
  p.M[2][0] += w*x3*y;  p.M[2][1] += w*x*y3;
  p.rhs[0] += w*z*x2;   p.rhs[1] += w*z*y2;  p.rhs[2] += w*z*x*y;
# else
  p.M[0][0] += w*x4; p.M[1][1] += w*y4; p.M[2][2] += w*x2*y2; 
  p.M[3][3] += w*x2; p.M[4][4] += w*y2; p.M[5][5] += w;
  p.M[0][2] += w*x3*y; p.M[0][3] += w*x3; p.M[0][4] += w*x2*y;
  p.M[1][2] += w*x*y3; p.M[1][3] += w*x*y2; p.M[1][4] += w*y3;
  p.M[2][5] += w*x*y;
  p.M[3][5] += w*x;
  p.M[4][5] += w*y;
  p.rhs[0] += w*x2*z; p.rhs[1] += w*y2*z; p.rhs[2] += w*x*y*z;
  p.rhs[3] += w*x*z; p.rhs[4] += w*y*z; p.rhs[5] += w*z;
# endif
#endif /* 3D */
}

static void parabola_fit_solve (ParabolaFit & p)
{
#if  AMREX_SPACEDIM==2
  p.M[0][1] = p.M[1][0];
  p.M[0][2] = p.M[2][0] = p.M[1][1];
  p.M[1][2] = p.M[2][1];
  if (matrix_inverse (p.M, 3, 1e-10)) {
    p.a[0] = p.M[0][0]*p.rhs[0] + p.M[0][1]*p.rhs[1] + p.M[0][2]*p.rhs[2];
    p.a[1] = p.M[1][0]*p.rhs[0] + p.M[1][1]*p.rhs[1] + p.M[1][2]*p.rhs[2];

  }
  else /* this may be a degenerate/isolated interface fragment */
    p.a[0] = p.a[1] = 0.;
#else /* 3D */
# if PARABOLA_SIMPLER
  p.M[0][1] = p.M[1][0]; p.M[0][2] = p.M[2][0];
  p.M[1][2] = p.M[2][1]; p.M[2][2] = p.M[1][0];
  if (matrix_inverse (p.M, 3, 1e-10)) {
    p.a[0] = p.M[0][0]*p.rhs[0] + p.M[0][1]*p.rhs[1] + p.M[0][2]*p.rhs[2];
    p.a[1] = p.M[1][0]*p.rhs[0] + p.M[1][1]*p.rhs[1] + p.M[1][2]*p.rhs[2];
    p.a[2] = p.M[2][0]*p.rhs[0] + p.M[2][1]*p.rhs[1] + p.M[2][2]*p.rhs[2];
  }
  else /* this may be a degenerate/isolated interface fragment */
    p.a[0] = p.a[1] = p.a[2] = 0.;
# else
  p.M[0][1] = p.M[2][2]; p.M[0][5] = p.M[3][3];
  p.M[1][5] = p.M[4][4];
  p.M[2][3] = p.M[0][4]; p.M[2][4] = p.M[1][3];
  p.M[3][4] = p.M[2][5];
  int i, j;
  for (i = 1; i < 6; i++)
    for (j = 0; j < i; j++)
      p.M[i][j] = p.M[j][i];
  if (matrix_inverse (p.M, 6, 1e-10)) {
    for (i = 0; i < 6; i++) {
      p.a[i] = 0.;
      for (j = 0; j < 6; j++)
        p.a[i] += p.M[i][j]*p.rhs[j];
    }
  }
  else { /* this may be a degenerate/isolated interface fragment */
    AllPrint()<<"warning: singular matrix in parabola_fit_solve"<<"\n";
    //amrex_error_module::amrex_warning("warning: singular matrix in parabola_fit_solve");    
    p.a[0] = p.a[1] = p.a[2] = 0.;
  }
# endif
#endif /* 3D */
}

static Real parabola_fit_curvature (ParabolaFit & p, Real kappamax)
{
  Real kappa;
#if AMREX_SPACEDIM==2
  Real dnm = 1. + p.a[1]*p.a[1];
  kappa = 2.*p.a[0]/sqrt (dnm*dnm*dnm);
#else /* 3D */
  Real hxx = 2.*p.a[0];
  Real hyy = 2.*p.a[1];
  Real hxy = p.a[2];
  Real hx, hy;
# if PARABOLA_SIMPLER
  hx = hy = 0.;
# else
  hx = p.a[3];
  hy = p.a[4];
# endif
  Real dnm = 1. + hx*hx + hy*hy;
  kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/sqrt (dnm*dnm*dnm);
#endif /* 3D */
  if (fabs (kappa) > kappamax) 
    return kappa > 0. ? kappamax : - kappamax;
  return kappa;
}