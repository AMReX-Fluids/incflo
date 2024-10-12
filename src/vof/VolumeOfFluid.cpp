#include <incflo.H>
#include <VolumeOfFluid.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#endif

using namespace amrex;
#define EPS 1e-4
#define THRESHOLD(c) {if ((c) < 0.) c = 0.; else if ((c) > 1.) c = 1.;}
#define CELL_IS_FULL(f)             ((f) == 0. || (f) == 1.)
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define CLAMP(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define VOF_NODATA std::numeric_limits<Real>::max()
#if AMREX_SPACEDIM==2
#define CELL_IS_BOUNDARY(c,min,max) (c[0]< min[0]||c[1]< min[1]||c[0]>max[0]||c[1]> max[1])
#else
#define CELL_IS_BOUNDARY(c,min,max) (c[0]< min[0]||c[1]< min[1]||c[0]>max[0]||c[1]> max[1]||c[2]< min[2]||c[2]> max[2])
#endif
// Define VofVector as an array of 3 doubles
using VofVector = Array<Real, 3>;
static_assert(sizeof(XDim3)==3*sizeof(Real));
typedef struct {
  Real min, max, sum, sum2, mean, stddev;
  int n;
} VofRange;
//
// range_init:
// @r: a #VofRange.
//
// Initializes a #VofRange.
//
void range_init (VofRange & r)
{
  r.max = std::numeric_limits<Real>::lowest();
  r.min = std::numeric_limits<Real>::max();
  r.sum = r.sum2 = 0.0;
  r.n = 0;
}
/**
 * range_add_value:
 * @r: a #VofRange.
 * @val: a value to add to @r.
 *
 * Adds @val to @r.
 */
void range_add_value (VofRange & r, Real val)
{
  if (val < r.min) r.min = val;
  if (val > r.max) r.max = val;
  r.sum += val;
  r.sum2 += val*val;
  r.n++;
}
/**
 * range_update:
 * @r: a #VofRange.
 *
 * Updates the fields of @r.
 */
void range_update (VofRange & r)
{
  if (r.n > 0) {
    if (r.sum2 - r.sum*r.sum/(Real) r.n >= 0.)
      r.stddev = sqrt ((r.sum2 - r.sum*r.sum/(Real) r.n)
            /(Real) r.n);
    else
      r.stddev = 0.;
    r.mean = r.sum/(Real) r.n;
  }
  else
    r.min = r.max = r.mean = r.stddev = 0.;
}

static void range_reduce (void * i, void * o, int * len,
                          MPI_Datatype * type)
{
  Real * in = (Real *) i;
  Real * inout = (Real *) o;
  if (in[0] < inout[0]) /* min */
    inout[0] = in[0];
  if (in[1] > inout[1]) /* max */
    inout[1] = in[1];
  inout[2] += in[2];    /* sum */
  inout[3] += in[3];    /* sum2 */
  inout[4] += in[4];    /* n */
}

static void domain_range_reduce ( VofRange & s)
{

  double in[5];
  double out[5] = { std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::lowest(), 0., 0., 0. };
  MPI_Op op;

  MPI_Op_create (range_reduce, true, &op);
  in[0] = s.min; in[1] = s.max; in[2] = s.sum; in[3] = s.sum2;
  in[4] = s.n;
  MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
  MPI_Op_free (&op);
  s.min = out[0]; s.max = out[1]; s.sum = out[2]; s.sum2 = out[3];
  s.n = out[4];

}


//////////////////////////////////////////////////////////////////////////////
////
////
/////////////////////////////////////////////////////////////////////////////

VolumeOfFluid::VolumeOfFluid (incflo* a_incflo) : v_incflo(a_incflo)
{
    finest_level = v_incflo->finestLevel();

    // *************************************************************************************
    // Allocate space for the information of the interface segments
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; ++lev){
        normal.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], AMREX_SPACEDIM, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev));
         alpha.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev));
           tag.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev));
        Array<MultiFab, 2> new_height={
         MultiFab(v_incflo->grids[lev], v_incflo->dmap[lev], AMREX_SPACEDIM, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev)),
         MultiFab(v_incflo->grids[lev], v_incflo->dmap[lev], AMREX_SPACEDIM, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev))
        };
        height.emplace_back(std::move(new_height));
        kappa.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev));
        force.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], AMREX_SPACEDIM, v_incflo->nghost_state(), MFInfo(), v_incflo->Factory(lev));
        //fixme
        force[lev].setVal(0,0,AMREX_SPACEDIM,v_incflo->nghost_state());
    }

    ParmParse pp("incflo");
    pp.query("output_drop_frequence", output_drop_frequence);

}
static XDim3 edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},{{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},{{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},{{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};
/* first index is the edge number, second index is the edge orientation
   (0 or 1), third index are the edges which this edge may connect to
   in order and the corresponding face direction */
static int connect[12][2][4] = {
  {{9, 1, 8, 3}, {4, 3, 7, 5}},     /* 0 */
  {{6, 2, 5, 4},  {8, 0, 9, 3}},    /* 1 */
  {{10, 3, 11, 2},  {5, 1, 6, 4}},  /* 2 */
  {{7, 0, 4, 5},   {11, 2, 10, 2}}, /* 3 */
  {{3, 7, 0, 5},   {8, 5, 11, 1}},  /* 4 */
  {{11, 4, 8, 1},  {1, 6, 2, 4}},   /* 5 */
  {{2, 5, 1, 4},  {9, 7, 10, 1}},   /* 6 */
  {{10, 6, 9, 1}, {0, 4, 3, 5}},    /* 7 */
  {{5, 11, 4, 1},  {0, 9, 1, 3}},   /* 8 */
  {{1, 8, 0, 3}, {7, 10, 6, 1}},    /* 9 */
  {{6, 9, 7, 1},  {3, 11, 2, 2}},   /* 10 */
  {{2, 10, 3, 2},   {4, 8, 5, 1}}   /* 11 */
};

static void cube_plane_intersection (XDim3 center, GpuArray <Real,AMREX_SPACEDIM> dx,
                                     XDim3 const & O, XDim3 const & n, Array<XDim3, 12> & p, int orient[12])
{
  XDim3 o;
  int i;
  /* use the corner of the cube as origin of the coordinate sys*/
  for (i=0; i<AMREX_SPACEDIM; i++)
    (&o.x)[i] = (&center.x)[i]-Real(0.5)*dx[i];
  for (i = 0; i < 12; i++) {
    XDim3 e, d;
    Real h= i<4? dx[0]:(i<8?dx[1]:dx[2]);
    /* d and e represent the coordinates of the starting and ending
      point of each edge of the cube, respectively */
    d.x = o.x + h*edge[i][0].x; d.y = o.y + h*edge[i][0].y; d.z = o.z + h*edge[i][0].z;
    e.x = o.x + h*edge[i][1].x; e.y = o.y + h*edge[i][1].y; e.z = o.z + h*edge[i][1].z;
    /* calculate the product of the edge vector and the intersecting plane normal vector */
    Real den = n.x*(e.x - d.x) + n.y*(e.y - d.y) + n.z*(e.z - d.z);
    orient[i] = -1;
    /* only when the edge is not parallel with the plane */
    if (fabs (den) > 1e-10) {
      Real t = (n.x*(O.x - d.x) + n.y*(O.y - d.y) + n.z*(O.z - d.z))/den;
      if (t >= 0. && t <= 1.) {
    p[i].x = d.x + t*(e.x - d.x); p[i].y = d.y + t*(e.y - d.y); p[i].z = d.z + t*(e.z - d.z);
    orient[i] = (n.x*(e.x - O.x) + n.y*(e.y - O.y) + n.z*(e.z - O.z) >= 0.);
      }
    }
  }

 /* for (i = 0; i < 12; i++)
   if(orient[i]>=0)
    for (int j=i+1; j<12; j++)
      if(orient[j]>=0) {
       Real dist = (p[i].x-p[j].x)*(p[i].x-p[j].x)+
                   (p[i].y-p[j].y)*(p[i].y-p[j].y)+
                   (p[i].z-p[j].z)*(p[i].z-p[j].z);
       if (dist < EPS*(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]))
             orient[j]=-1;
    }*/


}

/**
 * cut_cube_vertices:
 * @center: the coordinates of the cell center
 * @p: a point on the plane.
 * @n: the normal to the plane.
 * @v: where to return the vertices coordinates.
 * @d: where to return the direction.
 *
 * Fills @v, @d and @val with the coordinates/values of the vertices,
 * intersections of @cell with the plane defined by @p and @n.
 *
 * The vertices are ordered consistently to define a consistent,
 * oriented polygon.
 *
 * Returns: the number of vertices (0 if the plane does not cut the cell).
 */
static int cut_cube_vertices (XDim3 center, GpuArray <Real,AMREX_SPACEDIM> dx,
                           XDim3 const & p,  XDim3 const &  n,
                           Array<XDim3, 12> & v, int d[12])
{
  Array<XDim3, 12> a;
  int orient[12];
  int i;

  cube_plane_intersection (center, dx, p, n, a, orient);
  for (i = 0; i < 12; i++) {
    int nv = 0, e = i;
    bool duplicate = false;
    while (orient[e] >= 0) {
      int m = 0, * ne = connect[e][orient[e]];
      d[nv] = ne[3];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
        e = ne[m++];
    }
    if (nv > 2){

     if (nv > 4){
  /* there may exists duplicate points stored in v[] */
      int fnv=0;
      for (i = 0; i < nv; i++) {
        if (d[i]>=0)
          for (int j=i+1; j<nv; j++)
            if (d[j]>=0) {
               Real dist = (v[i].x-v[j].x)*(v[i].x-v[j].x)+
                             (v[i].y-v[j].y)*(v[i].y-v[j].y)+
                             (v[i].z-v[j].z)*(v[i].z-v[j].z);
                if (dist < EPS*(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])){
                     d[j]=-1;
                     fnv++;
                }
            }
      }
      if (fnv > 0){
        i=1;
       while (i<nv){
        if (d[i] < 0){
             for (int j=i; j<nv-1; j++){
              d[j]=d[j+1];
              v[j]=v[j+1];
             }
           nv--;
        }
        else
         i++;
       }
      }
     }
      return nv;
    }
  }
  return 0;
}



struct Segment{
  int nnodes;             /* number of nodes (2, 3 or 4) */
#if AMREX_SPACEDIM==2     /* 2D */
  XDim3 node[2];          /* node coordinates */
#else                     /* 3D */
  XDim3 node[4];
#endif
  XDim3 mv;
  Real alpha;
  //vof, tag
  Array<Real,3> vars;
 // Constructor to initialize the Segment
  Segment(int n,
#if AMREX_SPACEDIM==2  /* 2D */
  Array<XDim3, 4> const& nodes,
#else                  /* 3D */
  Array<XDim3, 12> const& nodes,
#endif
  XDim3 m, Real a, Array<Real,3> v, int ns=0)
   : nnodes(n), mv (m), alpha(a), vars(v) {
   for (int i = 0; i < n; ++i)
      node[i]= nodes[ns==0?i:(i + 3)%(n + 2)];
  }
};

static void add_segment (XDim3 const & center, GpuArray <Real,AMREX_SPACEDIM> const & dx,
                         Real alpha, XDim3 const & o, XDim3 const & m,
                         Vector<Segment> & segments, int & nt, Array<Real,3> vars)
{

   /*Print() <<" add_segment "<< *o<<"  "<<" vector "<<*m
                     <<"vof"<<"  "<<vof<<"  "<<"alpha " <<alpha<<"\n";  */

#if AMREX_SPACEDIM==2 /* 2D*/
  /* array of node coordinates for a cut face */
  Array<XDim3, 4>  nodecutface;
  Real x, y, h=dx[0];
  int n=0, nnodecutface;
  if (fabs (m.y) > EPS) {
    y = (alpha - m.x)/m.y;
    if (y >= 0. && y <= 1.) {
      nodecutface[n].x = center.x + h/2.; nodecutface[n].y = center.y + h*(y - 0.5); nodecutface[n++].z = 0.;
    }
  }
  if (fabs (m.x) > EPS) {
    x = (alpha - m.y)/m.x;
    if (x >= 0. && x <= 1.) {
      nodecutface[n].x = center.x + h*(x - 0.5); nodecutface[n].y = center.y + h/2.; nodecutface[n++].z = 0.;
    }
  }
  if (fabs (m.y) > EPS) {
    y = alpha/m.y;
    if (y >= 0. && y <= 1.) {
      nodecutface[n].x = center.x - h/2.; nodecutface[n].y = center.y + h*(y - 0.5); nodecutface[n++].z = 0.;
    }
  }
  if (fabs (m.x) > EPS) {
    x = alpha/m.x;
    if (x >= 0. && x <= 1.) {
      nodecutface[n].x = center.x + h*(x - 0.5); nodecutface[n].y = center.y - h/2.; nodecutface[n++].z = 0.;
    }
  }
  nnodecutface = n;
  if (n > 2) {
  /*check if there are duplicated points*/
   int i,j;
   bool ok[n];
   for (i=0; i<n; i++)
    ok[i]=true;
   XDim3 diff;
   for (i=0; i<n-1; i++)
    if (ok[i]){
     for (j=i+1; j<n; j++){
       diff.x = nodecutface[i].x - nodecutface[j].x;
       diff.y = nodecutface[i].y - nodecutface[j].y;
       if (diff.x*diff.x+diff.y*diff.y<1e-10)
         ok[j]=false;
     }
    }
   for (i=1; i<n; i++)
    if (!ok[i]){
     if (i!=n-1)
      for (j=i+1; j<n; j++){
          nodecutface[j-1].x = nodecutface[j].x;
          nodecutface[j-1].y = nodecutface[j].y;
      }
     nnodecutface--;
    }
  }
  AMREX_ASSERT (nnodecutface <= 2);
  /* assign data to nodeinfo array, increment number of wall faces and number of nodes */
  if (nnodecutface == 2) {
    nt += nnodecutface;
    segments.emplace_back(nnodecutface, nodecutface, m, alpha, vars);
  }
#else     /* 3D */
  int d[12];
  /* array of node coordinates for a cut face */
  Array<XDim3, 12>  nodecutface;
  int inode, inode2, jnode_max_sintheta = 0,
  nnodecutface = cut_cube_vertices (center, dx, o, m, nodecutface, d);
  AMREX_ASSERT (nnodecutface <= 6);
/*   Print()<<" add_segment "<<nnodecutface <<"\n";
   for (int i=0;i<nnodecutface;++i)
    Print()<<"node "<<i<<" "<<nodecutface[i]<<"\n";  */
  if (nnodecutface > 3) {   /* reorder faces if necessary */
   /* Tecplot can think that opposite vertices of the quadrilateral surface
      element are connected. This may result in ugly X-shaped surface elements.
      Reorder the array of node if necessary, using bubble sort,
      to maximize the sine of the angle between the vector from each
      intersection point to the cell center and the vector from this point
      to the next point in the array */
    int i_switchnodes = 0;        /* counter to avoid infinite loop */
    bool switchnodes = false;   /* logical variable used to reorder cut face nodes */

    do {
      i_switchnodes++;

      for (inode = 0; inode < nnodecutface; inode++) {
         XDim3 node = nodecutface[inode];     /* face node coordinates */
         XDim3 diff1 = {o.x - node.x, o.y - node.y, o.z - node.z};
         Real length_diff1 = sqrt(diff1.x*diff1.x+diff1.y*diff1.y+diff1.z*diff1.z);
         if (length_diff1 < 1e-20) /*degenerated case*/
              return;
         Real max_sintheta = 0.;
     /*  cycle through all other nodes (jnode) where cut face intersects cell edges */
         for (inode2 = 1; inode2 < nnodecutface; inode2++) {
           int jnode = (inode + inode2)%nnodecutface;
           XDim3 diff2 = {nodecutface[jnode].x - node.x,
                            nodecutface[jnode].y - node.y,
                            nodecutface[jnode].z - node.z};
           Real length_diff2 = sqrt(diff2.x*diff2.x+diff2.y*diff2.y+diff2.z*diff2.z);
           if (length_diff2 < 1e-20)
              return;
           Real sintheta = ((diff1.y*diff2.z - diff1.z*diff2.y)*m.x +
                              (diff1.z*diff2.x - diff1.x*diff2.z)*m.y +
                              (diff1.x*diff2.y - diff1.y*diff2.x)*m.z)/
                            (length_diff1*length_diff2);

           if (sintheta > max_sintheta) {
             max_sintheta = sintheta;
             jnode_max_sintheta = jnode;
           }
          }
         /* terminate if cannot find positive angle between cut face nodes */
         if (max_sintheta == 0.)
             return;
         AMREX_ASSERT (max_sintheta != 0.);
         inode2 = (inode + 1)%nnodecutface;
         if (jnode_max_sintheta != inode2) {
           node = nodecutface[jnode_max_sintheta];
           nodecutface[jnode_max_sintheta] = nodecutface[inode2];
           nodecutface[inode2] = node;
           switchnodes = true;
         }
      } /* inode-loop */
    } while (switchnodes && i_switchnodes < 1000);   /* avoid infinite loop */
  } /* reorder faces if necessary */
 // Print()<<"inside add_segment"<<nnodecutface<<"\n";
  if (nnodecutface < 3)
   return;
 /* If there are more than 4 points in the array,
    divide the cut face into a quadrilateral face and
    either a quadrilateral or triangular face */

 /* assign data to nodeinfo array, increment number of wall faces and number of nodes */
 if (nnodecutface <= 4) {
   nt += nnodecutface;
   segments.emplace_back(nnodecutface, nodecutface, m, alpha, vars);
   //Print() << " normal direction " <<nnodecutface<<" "<<*m<<""<<vof<<"\n";
 }
 else { /* cut face must be divided into 2 quadrilateral/triangular faces */
   /* first face is quadrilateral */
   nt += 4;
   segments.emplace_back(4, nodecutface, m, alpha, vars);
   /* second face is triangular if nnodecutface=5; otherwise quadrilateral */
//   WallFace * face2 = g_malloc (sizeof (WallFace));

   nnodecutface -= 2;

   nt += nnodecutface;
   segments.emplace_back(nnodecutface, nodecutface, m, alpha, vars, 3);
 } /* cut face must be divided into 2 quadrilateral/triangular faces */

#endif

}
/**
 * line_area:
 * @m: normal to the line.
 * @alpha: line constant.
 *
 * Returns: the area of the fraction of a cell lying under the line
 * (@m,@alpha).
 */
Real line_area (Array <Real, AMREX_SPACEDIM> &m, Real alpha)
{
  XDim3 n;
  Real alpha1, a, v, area;

  n.x = m[0], n.y = m[1];
  alpha1 = alpha;
  if (n.x < 0.) {
    alpha1 -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha1 -= n.y;
    n.y = - n.y;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y)
    return 1.;

  if (n.x == 0.)
    area = alpha1/n.y;
  else if (n.y == 0.)
    area = alpha1/n.x;
  else {
    v = alpha1*alpha1;

    a = alpha1 - n.x;
    if (a > 0.)
      v -= a*a;

    a = alpha1 - n.y;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*n.x*n.y);
  }

  return CLAMP (area, 0., 1.);
}

/**
 * line_alpha:
 * @m: a #FttVector.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the area of a square cell
 * lying under the line defined by @m.@x = @alpha is equal to @c.
 */
Real line_alpha (XDim3 & m, Real c)
{
  Real alpha, m1, m2, v1;

  m1 = fabs (m.x); m2 = fabs (m.y);
  if (m1 > m2) {
    v1 = m1; m1 = m2; m2 = v1;
  }

  v1 = m1/2.;
  if (c <= v1/m2)
    alpha = sqrt (2.*c*m1*m2);
  else if (c <= 1. - v1/m2)
    alpha = c*m2 + v1;
  else
    alpha = m1 + m2 - sqrt (2.*m1*m2*(1. - c));

  if (m.x < 0.)
    alpha += m.x;
  if (m.y < 0.)
    alpha += m.y;

  return alpha;
}

/**
 * line_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @a: area of cell fraction.
 * @p: a #coordinates.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a square cell lying under the line (@m,@alpha).
 */
void line_center (XDim3 const & m, Real alpha, Real a, XDim3 & p)
{
  XDim3 n;
  Real b;
  n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p.z = 0.;
  if (alpha <= 0.) {
    p.x = p.y = 0.;
    return;
  }

  if (alpha >= n.x + n.y) {
    p.x = p.y = 0.5;
    return;
  }


  if (n.x < EPS) {
    p.x = 0.5;
    p.y = m.y < 0. ? 1. - a/2. : a/2.;
    return;
  }

  if (n.y < EPS) {
    p.y = 0.5;
    p.x = m.x < 0. ? 1. - a/2. : a/2.;
    return;
  }

  p.x = p.y = alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p.x -= b*b*(alpha + 2.*n.x);
    p.y -= b*b*b;
  }

  b = alpha - n.y;
  if (b > 0.) {
    p.y -= b*b*(alpha + 2.*n.y);
    p.x -= b*b*b;
  }

  p.x /= 6.*n.x*n.x*n.y*a;
  p.y /= 6.*n.x*n.y*n.y*a;

  if (m.x < 0.)
    p.x = 1. - p.x;
  if (m.y < 0.)
    p.y = 1. - p.y;
}

/**
 * line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #amrex::XDim3.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
Real line_area_center (XDim3 const & m, Real alpha, XDim3 & p)
{
  XDim3 n;

  n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p.z = 0.;
  //made change from "alpha>=n.x+n.y" to consider
  //the extreme cases
  if (alpha <= 0. || alpha > n.x + n.y) {
    p.x = p.y = 0.;
    return 0.;
  }

  if (n.x < EPS) {
    p.x = 0.5;
    p.y = m.y < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  if (n.y < EPS) {
    p.y = 0.5;
    p.x = m.x < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  p.x = p.y = 0.;

  if (alpha >= n.x) {
    p.x += 1.;
    p.y += (alpha - n.x)/n.y;
  }
  else
    p.x += alpha/n.x;

  Real ax = p.x, ay = p.y;
  if (alpha >= n.y) {
    p.y += 1.;
    ay -= 1.;
    p.x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p.y += alpha/n.y;
    ay -= alpha/n.y;
  }

  p.x /= 2.;
  p.y /= 2.;

  THRESHOLD (p.x);
  THRESHOLD (p.y);

  if (m.x < 0.)
    p.x = 1. - p.x;
  if (m.y < 0.)
    p.y = 1. - p.y;

  return sqrt (ax*ax + ay*ay);

}

#if AMREX_SPACEDIM == 2 /* 3D */
#  define plane_volume         line_area
#  define plane_alpha          line_alpha
#  define plane_center         line_center
#  define plane_area_center    line_area_center
#else /* 3D */

/**
 * plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #amrex::XDim3.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 * note: the coordinate reference is (0.,0.,0.), the cube size is (1,1,1)
 * Returns: the area of the facet.
 */
Real plane_area_center (XDim3 const & m, Real alpha, XDim3 & p)
{

  if (fabs (m.x) < EPS) {
    XDim3 n, q;
    n.x = m.y;
    n.y = m.z;
    Real area = line_area_center (n, alpha, q);
    p.x = 0.5;
    p.y = q.x;
    p.z = q.y;
    return area;
  }
  if (fabs (m.y) < EPS) {
    XDim3 n, q;
    n.x = m.z;
    n.y = m.x;
    Real area = line_area_center (n, alpha, q);
    p.x = q.y;
    p.y = 0.5;
    p.z = q.x;
    return area;
  }
  if (fabs (m.z) < EPS) {
    Real area = line_area_center (m, alpha, p);
    p.z = 0.5;
    return area;
  }

  XDim3 n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }
  if (n.z < 0.) {
    alpha -= n.z;
    n.z = - n.z;
  }

  Real amax = n.x + n.y + n.z;
//  Print() << " plane_area_center " <<"("<<n.x<<","<<n.y<<","<<n.z<<")"<<" "<<amax<<"  "<< alpha<<"\n";
  if (alpha <= 0. || alpha > amax) {
    p.x = p.y = p.z = 0.;
    return 0.;
  }

  Real area = alpha*alpha;
  p.x = p.y = p.z = area*alpha;

  Real b = alpha - n.x;
  if (b > 0.) {
    area -= b*b;
    p.x -= b*b*(2.*n.x + alpha);
    p.y -= b*b*b;
    p.z -= b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    area -= b*b;
    p.y -= b*b*(2.*n.y + alpha);
    p.x -= b*b*b;
    p.z -= b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    area -= b*b;
    p.z -= b*b*(2.*n.z + alpha);
    p.x -= b*b*b;
    p.y -= b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    area += b*b;
    p.y += b*b*(2.*n.y + alpha - n.z);
    p.z += b*b*(2.*n.z + alpha - n.y);
    p.x += b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    area += b*b;
    p.x += b*b*(2.*n.x + alpha - n.z);
    p.z += b*b*(2.*n.z + alpha - n.x);
    p.y += b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    area += b*b;
    p.x += b*b*(2.*n.x + alpha - n.y);
    p.y += b*b*(2.*n.y + alpha - n.x);
    p.z += b*b*b;
  }

  area *= 3.;
  p.x /= area*n.x;
  p.y /= area*n.y;
  p.z /= area*n.z;

  THRESHOLD (p.x);
  THRESHOLD (p.y);
  THRESHOLD (p.z);

  if (m.x < 0.) p.x = 1. - p.x;
  if (m.y < 0.) p.y = 1. - p.y;
  if (m.z < 0.) p.z = 1. - p.z;

  return area*sqrt (1./(n.x*n.x*n.y*n.y) + 1./(n.x*n.x*n.z*n.z) + 1./(n.z*n.z*n.y*n.y))/6.;
}


/**
 * plane_volume:
 * @m: normal to the plane.
 * @alpha: plane constant.
 *
 * Returns: the volume of a cell lying under the plane (@m,@alpha).
 */
Real plane_volume (Array <Real, AMREX_SPACEDIM> &m, Real alpha)
{

  Real al = alpha + MAX(0., -m[0]) + MAX(0., -m[1]) + MAX(0., -m[2]);
  if (al <= 0.)
    return 0.;
  Real tmp = fabs(m[0]) + fabs(m[1]) + fabs(m[2]);
  if (al >= tmp)
    return 1.;
  AMREX_ASSERT (tmp > 0.);
  Real n1 = fabs(m[0])/tmp;
  Real n2 = fabs(m[1])/tmp;
  Real n3 = fabs(m[2])/tmp;
  al = MAX(0., MIN(1., al/tmp));
  Real al0 = MIN(al, 1. - al);
  Real b1 = MIN(n1*1, n2);
  Real b3 = MAX(n1*1, n2);
  Real b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  Real b12 = b1 + b2;
  Real bm = MIN(b12, b3);
  Real pr = MAX(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) + b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
       b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  Real volume = al <= 0.5 ? tmp : 1. - tmp;
  return CLAMP (volume, 0., 1.);
}
/**
 * plane_alpha:
 *
 * Returns: the value @alpha such that the volume of a cubic cell
 * lying under the plane defined by @m.@x = @alpha is equal to @c.
 */
Real plane_alpha (XDim3 & m, Real c)
{

  AMREX_ASSERT(c >= 0. && c <= 1.);
  // made change to avoid the numerical issue for
  // full cell that has empty cell in the face neighbor.
  if (c == 1. )
   c =1.-EPS;

  Real alpha;
  XDim3 n;
//  m.x =0., m.y=-0.5, m.z=0.5, c=.5;
//  Print()<<"vector"<<m.x<<" "<<m.y<<" "<<m.z<<"  "<<"vof "<<c<<"\n";
  n.x = fabs (m.x); n.y = fabs (m.y); n.z = fabs (m.z);

  Real m1, m2, m3;
  m1 = MIN(n.x, n.y);
  m3 = MAX(n.x, n.y);
  m2 = n.z;
//  Print()<<"vector 1 "<<m1<<" "<<m2<<" "<<m3<<"  "<<"vof "<<c<<"\n";
  if (m2 < m1) {
    Real tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    Real tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
//  Print()<<"vector 2 "<<m1<<" "<<m2<<" "<<m3<<"  "<<"vof "<<c<<"\n";
  Real m12 = m1 + m2;
  Real pr = MAX(6.*m1*m2*m3, 1e-50);
  Real V1 = m1*m1*m1/pr;
  Real V2 = V1 + (m2 - m1)/(2.*m3), V3;
  Real mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  Real ch = MIN(c, 1. - c);
//  Print()<<"vector 3 "<<mm<<" "<<V1<<" "<<V2<<"  "<< V3<<"\n";
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    Real p = 2.*m1*m2;
    Real q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    Real p12 = sqrt (p);
    Real teta = acos(q/(p*p12))/3.;
    Real cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 <= m3)
    alpha = m3*ch + mm/2.;
  else {
    Real p = m1*(m2 + m3) + m2*m3 - 1./4.;
    Real q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    Real p12 = sqrt(p);
//    Print()<<"p q p12  "<<p<<" "<<q<<" "<<p12<<"  "<<"vof"<<c<<"\n";
    Real teta = acos(q/(p*p12))/3.;
    Real cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (m.x < 0.)
    alpha += m.x;
  if (m.y < 0.)
    alpha += m.y;
  if (m.z < 0.)
    alpha += m.z;
//  Print()<<"alpha----  "<<alpha<<"\n";
  return alpha;
}

#endif

#if AMREX_SPACEDIM == 2
#include "myc2d.h"
#else
#include "myc.h"
#endif

#include <partstr.H>

void stencil (int const i, int const j, int const k,
              Array4<Real const> const & v,
              AMREX_D_PICK(      ,
                    Real fv[3][3],
                    Real fv[3][3][3]))
{
  int x, y, z = 0;
  AMREX_D_PICK(0,fv[1][1],fv[1][1][1]) = v (i,j,k);
#if AMREX_SPACEDIM == 3
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
    if (x != 0 || y != 0 || z != 0)
      AMREX_D_PICK(,fv[x + 1][y+1],fv[x + 1][y+1][z+1])=v(i+x,j+y,k+z);


  /* boundary conditions (symmetry) */
#if AMREX_SPACEDIM == 2
  for (x = 0; x <= 2; x++) {
    if (fv[x][0] < 0.) fv[x][0] = fv[x][1];
    if (fv[x][2] < 0.) fv[x][2] = fv[x][1];
  }
  for (y = 0; y <= 2; y++) {
    if (fv[0][y] < 0.) fv[0][y] = fv[1][y];
    if (fv[2][y] < 0.) fv[2][y] = fv[1][y];
  }
#else /* 3D */
  for (x = 0; x <= 2; x++)
    for (y = 0; y <= 2; y++) {
      if (fv[x][y][0] < 0.) fv[x][y][0] = fv[x][y][1];
      if (fv[x][y][2] < 0.) fv[x][y][2] = fv[x][y][1];
    }
  for (x = 0; x <= 2; x++)
    for (z = 0; z <= 2; z++) {
      if (fv[x][0][z] < 0.) fv[x][0][z] = fv[x][1][z];
      if (fv[x][2][z] < 0.) fv[x][2][z] = fv[x][1][z];
    }
  for (z = 0; z <= 2; z++)
    for (y = 0; y <= 2; y++) {
      if (fv[0][y][z] < 0.) fv[0][y][z] = fv[1][y][z];
      if (fv[2][y][z] < 0.) fv[2][y][z] = fv[1][y][z];
    }
#endif /* 3D */
}

bool interface_cell (int const i, int const j, int const k,
                     Array4<Real const> const & v, Real fc)
{
  if (fc == 1.){
//when a full cell has an empty cell in its face neighbor (i.e., +-x,+-y,+-z)
// we also need to calculate the normal and alpha of the cut plane.
  for (int dim=0; dim<AMREX_SPACEDIM; dim++)
    for (int step = -1; step <= 1; step+=2)
    if (dim==0?v(i+step,j,k) == 0.:
        dim==1?v(i,j+step,k) == 0.:
               v(i,j,k+step) == 0.){
          return true;
        }
  }
  else if (fc == 0.)
// empty cell
    return false;
  else
// 0.< vof value of cell < 1.
    return true;
  return false;
}

#define HMAX 5
#define SIGN(x) ((x) > 0. ? 1. : -1.)
#define BOUNDARY_HIT (2.*HMAX)


static int half_height (Array <int, 3> cell, Array4<Real const> const & fv, int d,
                        Real & H, int & n, Array<int,2> range)
{
  int s = 0, dim=d/2;
  n = 0;
  cell[dim]+=d%2?-1:1;
  while (n < HMAX && !s) {
    Real f = fv (cell[0],cell[1],cell[2],0);
    if (!CELL_IS_FULL(f)) { /* interfacial cell */
  //  if (f > EPS && f < 1. - EPS) { /* interfacial cell */
     //hit the boundary
     if (cell[dim]<range[0]||cell[dim]>range[1])
       return 2;
     H += f;
     n++;
    }
    else /* full or empty cell */
      s = (f - 0.5)>0.? 1.: -1;
    cell[dim]+=d%2?-1:1;
  }
  return s;
}

#define DMAX 3.5

static void height_propagation (Array <int, 3> cell, int dim, Array4<Real const> const & fv,
                                Array4<Real > const & hght, Array<int,2> range, Real orientation)
{
  for (int d = 1; d >= -1; d-=2, orientation = - orientation) {
    Array <int, 3> neighbor=cell;
    Real H = hght(cell[0],cell[1],cell[2],dim);
    neighbor[dim]+=d;
    bool interface = !CELL_IS_FULL(fv(neighbor[0],neighbor[1],neighbor[2],0));//false;
    while (fabs (H) < DMAX - 1.&& !interface &&
           neighbor[dim]>=range[0]&& neighbor[dim]<=range[1]) {
      H -= orientation;
      hght(neighbor[0],neighbor[1],neighbor[2],dim) = H;
      auto fvol = fv(neighbor[0],neighbor[1],neighbor[2],0);
      interface = !CELL_IS_FULL(fvol);
      neighbor[dim]+=d;
    }
  }
}

void calculate_height(int i, int j, int k, int dim, Array4<Real const> const & vof,
                      Array4<Real > const & hb, Array4<Real > const & ht, Array<int,2> range)
{
    Real H = vof(i,j,k,0);
    Array <int, 3> cell={i,j,k};
    // top part of the column
    int nt, st = half_height (cell, vof, 2*dim, H, nt, range);
    if (!st) /* still an interfacial cell */
      return;
    // bottom part of the column
    int nb, sb = half_height (cell, vof, 2*dim + 1, H, nb, range);
    if (!sb) /* still an interfacial cell */
      return;
    if (sb != 2 && st != 2) {
      if (st*sb > 0) /* the column does not cross the interface */
        return;
    }
    else { /* column hit a boundary */
      if (sb == 2 && st == 2) /* cannot hit a boundary on both sides */
        return;
      if (sb == 2)
        sb = st > 0.? -1.: 1;
      H += BOUNDARY_HIT;
    }
    if (sb > 0) {
      hb(i,j,k,dim) = H - 0.5 - nb;
      height_propagation (cell, dim, vof, hb, range, 1.);
    }
    else {
      ht(i,j,k,dim) = H - 0.5 - nt;
      height_propagation (cell, dim, vof, ht, range, -1.);
    }
}

static Array4<Real > const * boundary_hit (int i,int j,int k, int d, Array4<Real > const & hb,
                                           Array4<Real > const & ht)
{
  if (hb(i,j,k,d)!= VOF_NODATA && hb(i,j,k,d)> BOUNDARY_HIT/2.)
    return &hb;
  if (ht(i,j,k,d)!= VOF_NODATA && ht(i,j,k,d)> BOUNDARY_HIT/2.)
    return &ht;
  return nullptr;
}


static void height_propagation_from_boundary (Array <int, 3> cell, int dim, int d, Array4<Real const> const & fv,
                                              Array4<Real > const & hght, Array<int,2> range, int hb)
{
  Real orientation = (d % 2 ? -1 : 1)*hb;
  Real H = hght(cell[0],cell[1],cell[2],dim);;
  cell[dim]+=(d % 2 ? 1 : -1);
  Real H0=hght(cell[0],cell[1],cell[2],dim);
  while ( H0!=VOF_NODATA && H0 > BOUNDARY_HIT/2. &&
          cell[dim]>=range[0]&&cell[dim]<=range[1]) {
    H += orientation;
    hght(cell[0],cell[1],cell[2],dim) = H;
    cell[dim]+=(d % 2 ? 1 : -1);
    H0=hght(cell[0],cell[1],cell[2],dim);
  }
  /* propagate to non-interfacial cells up to DMAX */
  auto fvol = fv(cell[0],cell[1],cell[2],0);
  bool interface = !CELL_IS_FULL(fvol);
  while (fabs (H) < DMAX - 1. && !interface &&
         cell[dim]>=range[0]&&cell[dim]<=range[1]) {
    H += orientation;
    hght(cell[0],cell[1],cell[2],dim) = H;
    cell[dim]+=(d % 2 ? 1 : -1);
  }
}

Array4<Real const> const * closest_height (int i,int j,int k, int d, Array4<Real const> const & hb,
                           Array4<Real const> const & ht, Real * orientation)
{
  Array4<Real const> const * hv = nullptr;
  Real o = 0.;
  if (hb(i,j,k,d)!=VOF_NODATA) {
    hv = &hb; o = 1.;
    if (ht(i,j,k,d)!=VOF_NODATA &&
        fabs (ht(i,j,k,d)) < fabs (hb(i,j,k,d))) {
      hv = & ht; o = -1.;
    }
  }
  else if (ht(i,j,k,d)!=VOF_NODATA) {
    hv = & ht; o = -1.;
  }
  if (orientation) *orientation = o;
  return hv;
}
/* Returns: the height @h of the neighboring column in direction @c or
   VOF_NODATA if it is undefined. Also fills @x with the coordinates
   of the cell  */
static Real neighboring_column (int i,int j,int k, int c,
                      Array4<Real const> const * h , int d, Real * x)
{
  Array<int,3> neighbor={i,j,k};
  neighbor[d/2]+=d%2?-1:1;
  Real height=(*h)(neighbor[0],neighbor[1],neighbor[2],c);
   if (height!=VOF_NODATA) {
     *x = 1.;
     return height;
    }
  return VOF_NODATA;
}
/*
 The function is similar to neighboring_column().
 The difference is that neighboring_column_corner() returns height @h of the neighboring column in
 direction @(d[0], d[1]). kind of corner neighbors
 */

static Real neighboring_column_corner (int i,int j,int k, int c,
                                          Array4<Real const> const * h, int * d, Real (*x)[2])
{
  Array<int,3> neighbor={i,j,k};
  neighbor[d[0]/2]+=d[0]%2?-1:1;
  neighbor[d[1]/2]+=d[1]%2?-1:1;
  Real height=(*h)(neighbor[0],neighbor[1],neighbor[2],c);
  if (height!=VOF_NODATA) {
      (*x)[0] = d[0] % 2 ? -1. : 1.;
      (*x)[1] = d[1] % 2 ? -1. : 1.;
     return height;
  }
  else
   return VOF_NODATA;

}

static bool height_normal (int i,int j,int k, Array4<Real const> const & hb,
                           Array4<Real const> const & ht, XDim3 & m )
{
  Real slope = VOF_NODATA;
  static int oc[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
  for (int d = 0; d < AMREX_SPACEDIM; d++){
    Real orientation;
    Array4<Real const> const * hv = closest_height (i,j,k,d,hb,ht,&orientation);
    if (hv != nullptr && fabs ((*hv)(i,j,k,d) <= 1.)) {
      Real H = (*hv)(i,j,k,d);
      Real x[2], h[2][2]={0.}, hd[2]={0.};
      for (int nd = 0; nd < (AMREX_SPACEDIM==3?2:1); nd++) {
        h[nd][0] = neighboring_column (i,j,k, d, hv, 2*oc[d][nd], &x[0]);
        if (h[nd][0] == VOF_NODATA)
          break;
        h[nd][1] = neighboring_column (i,j,k, d, hv, 2*oc[d][nd] + 1, &x[1]);
        if (h[nd][1] == VOF_NODATA)
          break;
        x[1] = - x[1];
        Real det = x[0]*x[1]*(x[0] - x[1]), a = x[1]*(h[nd][0] - H), b = x[0]*(h[nd][1] - H);
        hd[nd] = (x[0]*b - x[1]*a)/det;
      }
      if (h[0][0] == VOF_NODATA || h[0][1] == VOF_NODATA ||
          h[1][0] == VOF_NODATA || h[1][1] == VOF_NODATA)
          continue;
      if (hd[0]*hd[0] + hd[1]*hd[1] < slope) {
          slope = hd[0]*hd[0] + hd[1]*hd[1];
          (&m.x)[d] = orientation;
          (&m.x)[oc[d][0]] = - hd[0];
#if AMREX_SPACEDIM==3
          (&m.x)[oc[d][1]] = - hd[1];
#endif
      }

    }
  }
  //Print()<<"-------slope---"<<slope<<"  "<<(slope < VOF_NODATA)<<"\n";
  return slope < VOF_NODATA;

}

static Real curvature_from_h (Real *x, Real *h, int c)
{
   Real det, a, b, hxx,hx,hyy,hy,hxy, kappa;
#if AMREX_SPACEDIM == 2
  det = x[0]*x[1]*(x[0] - x[1]), a = x[1]*(h[0] - h[2]), b = x[0]*(h[1] - h[2]);
  hxx = 2.*(a - b)/det;
  hx = (x[0]*b - x[1]*a)/det;
  Real dnm = 1. + hx*hx;
  kappa = hxx/sqrt(dnm*dnm*dnm);
#else /* 3D*/
   /*det = x[0]*x[1]*(x[0] - x[1]); a = x[1]*(h[0] - h[8]); b = x[0]*(h[1] - h[8]);
   Real hxx = 2.*(a - b)/det;
   Real hx = (x[0]*b - x[1]*a)/det;

   det = x[2]*x[3]*(x[2] - x[3]); a = x[3]*(h[2] - h[8]); b = x[2]*(h[3] - h[8]);
   Real hyy = 2.*(a - b)/det;
   Real hy = (x[2]*b - x[3]*a)/det;

   det = x[4]*x[5]*(x[4] - x[5]); a = x[5]*(h[4] - h[2]); b = x[4]*(h[5] - h[2]);
   Real hx0 = (x[4]*b - x[5]*a)/det;
   det = x[6]*x[7]*(x[6] - x[7]); a = x[7]*(h[6] - h[3]); b = x[6]*(h[7] - h[3]);
   Real hx1 = (x[6]*b - x[7]*a)/det;
   det = x[2]*x[3]*(x[2] - x[3]); a = x[3]*(hx0 - hx); b = x[2]*(hx1 - hx);
   Real hxy = (x[2]*b - x[3]*a)/det;*/

   hxx = h[0] - 2.*h[8] + h[1];
   hyy = h[2] - 2.*h[8] + h[3];
   hx = (h[0] - h[1])/2.;
   hy = (h[2] - h[3])/2.;
   hxy = (h[4] + h[7] - h[5] - h[6])/4.;

   /*hxx = (h[4] - 2.*h[8] + h[7])/2.;
   hyy = (h[5] - 2.*h[8] + h[6])/2.;
    hx = (h[4] - h[7])/2./sqrt(2.);
    hy = (h[5] - h[6])/2./sqrt(2.);
    hxy = (h[2] + h[3] - h[1] - h[0])/2.;*/

   Real dnm = 1. + hx*hx + hy*hy;
   kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/sqrt (dnm*dnm*dnm);
#endif
   return kappa;
}


/**
 * curvature_along_direction:
 * @(i,j,k): current cell.
 * @d: x, y or z.
 * @kappa: the curvature.
 *
 * Tries to compute an interface curvature for @cell using
 * height-functions on equally-spaced columns in direction @d.
 *
 * Returns: %true if the curvature was successfully computed, %false
 * otherwise.
 */

bool curvature_along_direction (int i,int j,int k, int d, GpuArray<Real, AMREX_SPACEDIM> dx,
                                Array4<Real const> const & hb,
                                Array4<Real const> const & ht,
                                Real & kappa)
{
  Real x[AMREX_SPACEDIM==3?9:3], h[AMREX_SPACEDIM==3?9:3];
  Real orientation;
  static int oc[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
  Array4<Real const> const * hv = closest_height (i,j,k,d,hb,ht,&orientation);

  if (!hv) {
    bool loop=true;
     /* no data for either directions, look four neighbors to collect potential interface positions */
    for (int nd = 0; nd < (AMREX_SPACEDIM==3?2:1) && loop; nd++)
     for (int ndd = 0; ndd < 2; ndd++) {
       Array<int,3> neighbor={i,j,k};
       neighbor[oc[d][nd]]+=ndd%2?-1:1;
       hv = closest_height (neighbor[0],neighbor[1],neighbor[2],d,hb,ht,&orientation);
       if (hv){
         loop = false;
         break;
       }
    }
    if (!hv) /* give up */
      return false;
  }
  else if (fabs((*hv)(i,j,k,d))>1.)
     return false;
  int n=0;
  for (int nd = 0; nd < (AMREX_SPACEDIM==3?2:1); nd++) {
    h[n] = neighboring_column (i, j, k, d, hv, 2*oc[d][nd], &x[n]);
    if (h[n] == VOF_NODATA )
      break;
    n++;
    h[n] = neighboring_column (i, j, k, d, hv, 2*oc[d][nd]+1, &x[n]);
    if (h[n] == VOF_NODATA )
      break;
    x[n] = - x[n];
    n++;
  }
#if AMREX_SPACEDIM==2 /* 2D */
  if (n == 2) {
    h[n] = (*hv)(i,j,k,d); x[n] = 0.;
    if(h[n] != VOF_NODATA) {
      kappa = curvature_from_h (x, h, d)/dx[0];
      return true;
    }
  }
#else   /* 3D */
  if (n == 4) {
   int od[2];
   Real xd[2];
   bool search = true;
   for (int nd = 0; nd < 2 && search; nd++)
     for (int nj = 0; nj < 2 && search; nj++) {
       od[0] = 2*oc[d][0] + nj;
       od[1] = 2*oc[d][1] + nd;
       h[n] = neighboring_column_corner (i, j, k, d, hv, od, &xd);
       x[n] = xd[0];
       if (h[n] == VOF_NODATA || fabs (x[n]) != 1.) {
         search = false;
         break;
       }
       n++;
     }
    if (n == 8) {
    /* all nine height function are found
     * use 9-point stensil to calculate curvature
     */
      h[n] = (*hv)(i,j,k,d); x[n] = 0.;
      if(h[n] != VOF_NODATA) {
        kappa = curvature_from_h (x, h, d)/dx[0];
        return true;
      }
    }
    else {
      h[4] = (*hv)(i,j,k,d); x[4] = 1.;
      int ni;
      for (ni = 0; ni < 5; ni++)
        if (h[ni] == VOF_NODATA || fabs(x[ni])<1e-20)
           break;
      if (ni == 5 ) {
        Real s[4];
        for (int nk=0; nk<4; nk++){
          s[nk] = (h[nk]-h[4])/x[nk];
          s[nk] = s[nk] /sqrt(1 + s[nk]*s[nk]);
        }
        kappa = (s[0] - s[1] + s[2] - s[3]) /dx[0];
        return true;
      }
    }
  }
#endif
  return false;
}


bool curvature_along_direction_new (int i,int j,int k, int d, GpuArray<Real, AMREX_SPACEDIM> dx,
                                Array4<Real const> const & hb,
                                Array4<Real const> const & ht,
                                Real & kappa, Vector<VofVector> &interface)
{
  Real x[AMREX_SPACEDIM==3?9:3], h[AMREX_SPACEDIM==3?9:3];
  Real orientation;
  static int oc[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
  Array4<Real const> const * hv = closest_height (i,j,k,d,hb,ht,&orientation);

  if (!hv) {
    bool loop=true;
    /* no data for either directions, look four neighbors to collect potential interface positions */
    for (int nd = 0; nd < (AMREX_SPACEDIM==3?2:1) && loop; nd++)
      for (int ndd = 0; ndd < 2; ndd++) {
        Array<int,3> neighbor={i,j,k};
        neighbor[oc[d][nd]]+=ndd%2?-1:1;
        hv = closest_height (neighbor[0],neighbor[1],neighbor[2],d,hb,ht,&orientation);
        if (hv){
         loop = false;
         break;
        }
    }
    if (!hv) /* give up */
      return false;
  }
  else if (fabs((*hv)(i,j,k,d))>1.)
     return false;
  int n=0;
  for (int nd = 0; nd < (AMREX_SPACEDIM==3?2:1); nd++) {
    h[n] = neighboring_column (i, j, k, d, hv, 2*oc[d][nd], &x[n]);
    n++;
    h[n] = neighboring_column (i, j, k, d, hv, 2*oc[d][nd]+1, &x[n]);
    x[n] = - x[n];
    n++;
  }

#if AMREX_SPACEDIM==2 /* 2D */
  h[n] = (*hv)(i,j,k,d); x[n] = 0.;
  if (h[2] != VOF_NODATA && h[0] != VOF_NODATA && h[1] != VOF_NODATA) {
    kappa = curvature_from_h (x, h, d)/dx[0];
    return true;
  }
  else { /* h[2] == VOF_NODATA || h[0] == VOF_NODATA || h[1] == VOF_NODATA */
    /* collect interface positions (based on height function) */
    VofVector pos;
    for (n = 0; n < 3; n++)
     if (h[n] != VOF_NODATA) {
       pos[oc[d][0]] = x[n];
       pos[d] = orientation*h[n];
       interface.emplace_back(std::move(pos));
     }
    return false;
  }
#else   /* 3D */
  int od[2],m=0;
  Real xd[4][2];
  for (int nd = 0; nd < 2 ; nd++)
    for (int nj = 0; nj < 2 ; nj++) {
      od[0] = 2*oc[d][0] + nj;
      od[1] = 2*oc[d][1] + nd;
      h[n] = neighboring_column_corner (i, j, k, d, hv, od, &xd[m]);
      x[n] = xd[m][0];
      n++,m++;
    }
  h[n] = (*hv)(i,j,k,d); x[n] = 0.;
  for (n = 0; n < 9; n++)
    if (h[n] == VOF_NODATA)
       break;
  if (n == 9) {
    /* all nine height functions are found */
    kappa = curvature_from_h (x, h, d)/dx[0];
    return true;
  }
  else {
      /* collect interface positions (based on height function) */
    VofVector pos;
    for (n = 0; n < 9; n++)
     if (h[n] != VOF_NODATA) {
      if (n < 2) {
         pos[oc[d][0]] = x[n];
         pos[oc[d][1]] = 0.;
       }
       else if (n < 4) {
         pos[oc[d][1]] = x[n];
         pos[oc[d][0]] = 0.;
       }
       else if (n == 8) {
         pos[oc[d][0]] = 0.;
         pos[oc[d][1]] = 0.;
       }
       else {
         pos[oc[d][0]] = xd[n-4][0];
         pos[oc[d][1]] = xd[n-4][1];
       }
       pos[d] = orientation*h[n];
       interface.emplace_back(std::move(pos));
     }
    return false;
  }
#endif   /* 3D */
  return false;
}

static void orientation (VofVector m, Array<int,3> &c)
{
  int i, j;
  for (i = 0; i < AMREX_SPACEDIM; i++)
    c[i] = i;
  for (i = 0; i < AMREX_SPACEDIM - 1; i++)
   for (j = 0; j < AMREX_SPACEDIM - 1 - i; j++)
    if (fabs (m[c[j + 1]]) > fabs (m[c[j]])) {
      int tmp = c[j];
      c[j] = c[j + 1];
      c[j + 1] = tmp;
    }
}

static int independent_positions (Vector<VofVector> &interface)
{
  if (interface.size() < 2)
    return interface.size();

  int j, ni = 1;
  for (j = 1; j < interface.size() ; j++) {
    int i;
    bool depends = false;
    for (i = 0; i < j && !depends; i++) {
      Real d2 = 0.;
      for (int c = 0; c < AMREX_SPACEDIM; c++)
       d2 += (interface[i][c] - interface[j][c])*(interface[i][c] - interface[j][c]);
      depends = d2 < 0.5*0.5;
    }
    ni += !depends;
  }
  return ni;
}



/**
 * height_curvature_combined:
 * @(i,j,k): current cell containing an interface.
 * @dx: the grid size of the box.
 * @hb: stores the height values in the axis positive direction.
 * @ht: stores the height values in the axis negative direction.
 * @mv: stores the normal direction of interface.
 * @alpha: stores the plane constant of reconstructed segment.
 *
 * Tries to estimate the curvature of an interface using
 * height-functions, either on equally-spaced columns,
 * or using parabola fits of interface positions
 * defined using the height-functions in all directions.
 *
 *
 * Returns: (double in 3D) the mean curvature of the interface
 * contained in @(i,j,k), or %VOF_NODATA if the HF method could not
 * compute a consistent curvature.
 */
Real height_curvature_combined (int i,int j,int k, GpuArray<Real, AMREX_SPACEDIM> dx,
                                Array4<Real const> const & hb,
                                Array4<Real const> const & ht,
                                Array4<Real const> const & mv,
                                Array4<Real const> const & alpha)
{
   VofVector m;
   Array<int, 3> try_dir;
   for (int d = 0; d < AMREX_SPACEDIM; d++)
     m[d] = mv(i,j,k,d);

   orientation (m, try_dir); /* sort directions according to normal */

   Real kappa = 0.;
   Vector<VofVector> interface;
   for (int d = 0; d < AMREX_SPACEDIM; d++) /* try each direction */
     if (curvature_along_direction_new (i, j, k, try_dir[d], dx, hb, ht, kappa, interface))
      return kappa;
  /* Could not compute curvature from the simple algorithm along any direction:
   * Try parabola fitting of the collected interface positions */

   if (independent_positions (interface) < 3*(AMREX_SPACEDIM - 1))
    return VOF_NODATA;

   ParabolaFit fit;
   XDim3 mx={AMREX_D_DECL(m[0],m[1],m[2])},p;

   Real area=plane_area_center (mx, alpha(i,j,k,0),p);
   //shift the coordinates of the center of the interfacial segment
   //by using the center of the cube. plane_area_center() gives the
   //coordinates of area center with the coordinate origin as (0.,0.,0.)
   //After shifting, the origin becomes cell center.
   for (int c = 0; c < AMREX_SPACEDIM; c++)
    (&p.x)[c] -= 0.5;
   // initialize the parameters for parabola fit
   parabola_fit_init (fit, p, mx);

////#if AMREX_SPACEDIM==2
////  parabola_fit_add (&fit, &fc.x, PARABOLA_FIT_CENTER_WEIGHT);
////#elif !PARABOLA_SIMPLER
   parabola_fit_add (fit, {AMREX_D_DECL(p.x,p.y,p.z)}, area*100.);
////#endif
   for (int c = 0; c < interface.size(); c++)
     parabola_fit_add (fit, interface[c], 1.);
   parabola_fit_solve (fit);
   kappa = parabola_fit_curvature (fit, 2.)/dx[0];
# if PARABOLA_SIMPLER || AMREX_SPACEDIM==2
  int nn=3;
# else
  int nn=6;
# endif
   for (int c = 0; c < nn; c++)
        delete[] fit.M[c]; // Delete each row
   delete[] fit.M; // Delete the array of pointers
  return kappa;
}

/**
 * curvature_fit:
 * @(i,j,k): current cell containing an interface.
 * @dx: the grid size of the box
 * @vof: stores the volume fraction value of the cell.
 * @mv: stores the normal direction of interface.
 * @alpha: stores the plane constant of reconstructed segment.
 *
 * Computes an approximation of the curvature of the interface
 * contained in @cell(i,j,k) using paraboloid fitting of the centroids of the
 * reconstructed interface segments.
 *
 *
 * Returns: (double in 3D) the mean curvature of the interface contained in @cell.
 */
Real curvature_fit (int i,int j,int k, GpuArray<Real, AMREX_SPACEDIM> dx,
                    Array4<Real const> const & vof,
                    Array4<Real const> const & mv,
                    Array4<Real const> const & alpha)
{
   VofVector m;

   for (int d = 0; d < AMREX_SPACEDIM; d++)
     m[d] = mv(i,j,k,d);

   ParabolaFit fit;
   XDim3 mx={AMREX_D_DECL(m[0],m[1],m[2])},p;

   Real area=plane_area_center (mx, alpha(i,j,k,0),p);
   //shift the coordinates of the center of the interfacial segment
   //by using the center of the cube. plane_area_center() gives the
   //coordinates of area center with the coordinate origin as (0.,0.,0.)
   //After shifting, the origin becomes cell center.
   for (int c = 0; c < AMREX_SPACEDIM; c++)
    (&p.x)[c] -=  0.5;
   // initialize the parameters for parabola fit
   parabola_fit_init (fit, p, mx);
   // add the center of the segment with the area of the segment as weight
   parabola_fit_add (fit, {AMREX_D_DECL(p.x,p.y,p.z)}, area);
   int di=0,dj=0,dk=0;
#if AMREX_SPACEDIM==3
   for (dk = -2; dk <= 2; dk++)
#endif
    for (dj = -2; dj <= 2; dj++)
     for (di = -2; di <= 2; di++)
       if (di != 0|| dj != 0|| dk != 0) {
         int ni=i+di,nj=j+dj,nk=k+dk;
         Real fvol=vof(ni,nj,nk,0);
         if (!CELL_IS_FULL(fvol)){
            mx={AMREX_D_DECL(mv(ni,nj,nk,0),mv(ni,nj,nk,1),mv(ni,nj,nk,2))};
            area=plane_area_center (mx, alpha(ni,nj,nk,0),p);
            for (int c = 0; c < AMREX_SPACEDIM; c++)
              (&p.x)[c] += (c==0?di:c==1?dj:dk) - 0.5;
            parabola_fit_add (fit, {p.x,p.y,p.z}, area);
         }
       }
   parabola_fit_solve (fit);
   Real kappa = parabola_fit_curvature (fit, 2.)/dx[0];
# if PARABOLA_SIMPLER  || AMREX_SPACEDIM==2
  int nn=3;
# else
  int nn=6;
# endif
   for (int c = 0; c < nn; c++)
        delete[] fit.M[c]; // Delete each row
   delete[] fit.M; // Delete the array of pointers
  return kappa;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
///////
///////   Update VOF properties including height values and normal direction
///////
/////////////////////////////////////////////////////////////////////////////////////////////////

void
VolumeOfFluid::tracer_vof_update (int lev, MultiFab & vof_mf, Array<MultiFab,2> & height)
{
  Geometry const& geom =v_incflo->geom[lev];
  auto const& dx = geom.CellSizeArray();
  auto const& problo = geom.ProbLoArray();
  auto const& probhi = geom.ProbHiArray();
///////////////////////////////////////////////////
//          update height using vof field
///////////////////////////////////////////////////
  for (int dim = 0; dim < AMREX_SPACEDIM; dim++){

    height[0].setVal(VOF_NODATA,dim,1,v_incflo->nghost_state());
    height[1].setVal(VOF_NODATA,dim,1,v_incflo->nghost_state());
    //fix me: have not thought of a way to deal with the MFIter with tiling
    //an option is to use similar way as MPI's implementation.
    for (MFIter mfi(vof_mf); mfi.isValid(); ++mfi) {
       Box const& bx = mfi.validbox();
       Array<int,2> range ={bx.smallEnd()[dim], bx.bigEnd()[dim]};
       Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
       Array4<Real > const& hb_arr = height[0].array(mfi);
       Array4<Real > const& ht_arr = height[1].array(mfi);
       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         auto fvol = vof_arr(i,j,k,0);
         if (!CELL_IS_FULL(fvol)){
            calculate_height(i, j, k, dim, vof_arr, hb_arr, ht_arr, range);
         }// end if
       }); //end ParallelFor

    } //end MFIter
    //fix me: temporary solution for MPI boundaries
    height[0].FillBoundary(geom.periodicity());
    height[1].FillBoundary(geom.periodicity());

//deal with the situation where interface goes across the MPI or periodic boundaries.
if(1){
    for (MFIter mfi(vof_mf); mfi.isValid(); ++mfi) { /*fix me: no titling*/
       Box const& bx = mfi.validbox();
       Array<IntVect, 2> face_min_max;
       Array4<Real > const& hb_arr = height[0].array(mfi);
       Array4<Real > const& ht_arr = height[1].array(mfi);
       Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
       //search the cells on each boundary of the validbox
       //we do it by creating a new indexing space (i.e., bbx) with a constant
       //value for one coordinate direction. i.e., for +X face of the box, we can
       // set i=imax and just vary j and k index.
         Array<int,2> range = {bx.smallEnd()[dim],bx.bigEnd()[dim]};
         auto ijk_min= bx.smallEnd();
         auto ijk_max= bx.bigEnd();
//only loop through cells on two faces in the axis (defined by 'dim')
         for (int nn = 0; nn < 2; nn++){
//Note: we use the notation of Gerris for the direction of the Box (i.e.,FttDirection)
// FACE direction = 0,1,2,3,4,5 in 3D
// X+ (Right):0, X- (Left):1, Y+ (Top): 2, Y- (Bottom): 3, Z+ (Front): 4, Z- (Back):5
// direction%2=0 means the positive direction of a given axis direction (i.e.,int direction/2)
// direction%2=1 means the negative direction of a given axis direction    (i.e.,int direction/2)
// Axis direction = 0 (X-axis), 1(Y-axis), 2(Z-axis)
// therefore, 'nn=0' here means the positive direction.
           ijk_min[dim]= range[nn?0:1];
           ijk_max[dim]= range[nn?0:1];
           Box bbx(ijk_min, ijk_max);
// loop through the cells on the face of the box ('bbx')
           ParallelFor(bbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
             Array<int, 3> cell={i,j,k}, ghost=cell;
             ghost[dim]+=nn%2?-1:1;
             Array4<Real > const * h=boundary_hit (i,j,k, dim, hb_arr,ht_arr);
    /*if (i==7 && j==0 && k==5){
        AllPrint()<<"test_height_function   "<<"hb  "<<hb_arr(i,j,k,1)<<" ht " <<ht_arr(i,j,k,1)<<"\n";
        AllPrint()<<"----------"<<"\n";
    }*/
             // column hit boundary
             if(h){
               // column hit boundary
               Array4<Real > const *hn=boundary_hit (ghost[0],ghost[1],ghost[2], dim, hb_arr,ht_arr);
               if(h==hn){
                // the column crosses the interface
                // propagate column height correction from one side (or PE) to the other
                  Real orientation = (nn%2 ? -1:1)*(h == &hb_arr ? 1 : -1);
                Real h_ghost=(*h)(ghost[0],ghost[1],ghost[2],dim);
                Real Hn = h_ghost + 0.5 + (orientation - 1.)/2. - 2.*BOUNDARY_HIT;
                (*h)(i,j,k,dim) += Hn;
                height_propagation_from_boundary (cell, dim, 2*dim+nn, vof_arr, *h, range, h == &hb_arr ? 1 : -1);
               }
               else{
                // the column does not cross the interface
                Real hgh=(*h)(cell[0],cell[1],cell[2],dim);
                while (!CELL_IS_BOUNDARY(cell,bx.smallEnd(),bx.bigEnd()) &&
                        hgh!= VOF_NODATA && hgh> BOUNDARY_HIT/2.) {
                  (*h)(cell[0],cell[1],cell[2],dim) = VOF_NODATA;
                    cell[dim]+=nn%2?1:-1;
                }
               }
             }
             else{
              // column did not hit a boundary, propagate height across PE boundary */
              if (hb_arr(ghost[0],ghost[1],ghost[2],dim)!= VOF_NODATA)
                 height_propagation (ghost, dim, vof_arr, hb_arr, range, 1.);
              if (ht_arr(ghost[0],ghost[1],ghost[2],dim)!= VOF_NODATA)
                 height_propagation (ghost, dim, vof_arr, ht_arr, range, -1.);
             }
            //Print()<<"face_loop   "<<"i  "<<i<<" j " <<j<<" k "<<k<<"\n";
           });//end ParallelFor

         }// end loop of 6 faces of the valid box
    }// end MFIter
}
    for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real > const& hb_arr = height[0].array(mfi);
        Array4<Real > const& ht_arr = height[1].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           if (hb_arr(i,j,k,dim)!= VOF_NODATA && hb_arr(i,j,k,dim)> BOUNDARY_HIT/2)
             hb_arr(i,j,k,dim)= VOF_NODATA;
           if (ht_arr(i,j,k,dim)!= VOF_NODATA && ht_arr(i,j,k,dim)> BOUNDARY_HIT/2)
             ht_arr(i,j,k,dim)= VOF_NODATA;
        });
    } // end MFIter
    //fix me: temporary solution for MPI boundaries
    height[0].FillBoundary(geom.periodicity());
    height[1].FillBoundary(geom.periodicity());
  }//end for dim


/////////////////////////////////////////////////////////////////////////////////////////
//update the normal and alpha
/////////////////////////////////////////////////////////////////////////////////////////
    for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
        Array4<Real> const& mv = normal[lev].array(mfi);
        Array4<Real> const& al = alpha[lev].array(mfi);
        Array4<Real const > const& hb_arr = height[0].const_array(mfi);
        Array4<Real const > const& ht_arr = height[1].const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          XDim3 m={0.,0.,0.};
          auto fvol = vof_arr(i,j,k,0);
          THRESHOLD(fvol);
   /*if (i==5&&j==6&&k==8){
       int dddd;
       Print()<<"------------"<<"\n";
   }*/
          if (!height_normal (i,j,k, hb_arr, ht_arr, m)){
//          if(1){
            if (!interface_cell (i,j,k, vof_arr, fvol)) {
                AMREX_D_TERM(mv(i,j,k,0) = Real(0.);,
                             mv(i,j,k,1) = Real(0.);,
                             mv(i,j,k,2) = Real(0.););
                al(i,j,k) = fvol;
            }
            else {
               AMREX_D_PICK( ,Real f[3][3];, Real f[3][3][3];)
               stencil (i,j,k, vof_arr, f);
               mycs (f, &m.x);
            }
          }
          Real n = 0.;
          for (int d = 0; d < AMREX_SPACEDIM; d++)
               n += fabs ((&m.x)[d]);
          if (n > 0.)
             for (int d = 0; d < AMREX_SPACEDIM; d++)
                mv(i,j,k,d)= (&m.x)[d]/n;
          else {/* fixme: this is a small fragment */
             AMREX_D_TERM(mv(i,j,k,0) = Real(1.);,
                          mv(i,j,k,1) = Real(0.);,
                          mv(i,j,k,2) = Real(0.););
          }
          for (int d = 0; d < AMREX_SPACEDIM; d++)
            (&m.x)[d]= mv(i,j,k,d);
       /*    Print() <<" normal direction "<< m.x<<" "<<m.y<<" "<<m.z<<"("<<i<<","<<j<<","<<k<<")"
                       <<"vof"<<"  "<<fvol<<"\n";*/
               //if (i==22&&j==20&&k==19)
              al(i,j,k)= plane_alpha (m, fvol);

                  /*    if (i==12&&j==7&&k==7)
              Print() <<" normal direction "<<"("<<i<<","<<j<<","<<k<<") "
                     <<"vof"<<"  "<<fvol<<"  "<<"alpha " <<al(i,j,k)<<" "<<
                     interface_cell (i,j,k, vof, fvol)<<"\n";*/
        }); //  ParallelFor

    } // end MFIter

    //!!!!!!!!fix me: a temporary solution for the normal and alpha!!!!!!!!!!!
    // fill value of ghost cells (BCs, MPI info.)
    normal[lev].FillBoundary(geom.periodicity());
    alpha[lev].FillBoundary(geom.periodicity());
}

///////////////////////////////////////////////////////////////////////////////////////////////
/////
/////                          curvature calculation
/////
///////////////////////////////////////////////////////////////////////////////////////////////
void
VolumeOfFluid::curvature_calculation (int lev, MultiFab & vof_mf, Array<MultiFab,2> & height, MultiFab & kappa)
{

  Geometry const& geom =v_incflo->geom[lev];
  auto const& dx = geom.CellSizeArray();
  auto const& problo = geom.ProbLoArray();
  auto const& probhi = geom.ProbHiArray();
  MultiFab n_max(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(),
                 MFInfo(), v_incflo->Factory(lev));
  //fixme: need to change for BCs
  kappa.setVal(VOF_NODATA,0,1,v_incflo->nghost_state());
  n_max.setVal(-1.0);
// use height function method to calculate curvature
  for (int dim = 0; dim < AMREX_SPACEDIM; dim++){
    for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       Box const& bx = mfi.tilebox();
       Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
       Array4<Real const> const& mv = normal[lev].const_array(mfi);
       Array4<Real const> const& hb_arr = height[0].const_array(mfi);
       Array4<Real const> const& ht_arr = height[1].const_array(mfi);
       Array4<Real > const& kappa_arr = kappa.array(mfi);
       Array4<Real > const& nmax_arr = n_max.array(mfi);
       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         auto fvol = vof_arr(i,j,k,0);
         Real kappa0;
         if (!CELL_IS_FULL(fvol)){
  /* if ((i==4||i==11)&&j==9&&k==0){
       int dddd;
       Print()<<"------------"<<"\n";
   }*/
          if (curvature_along_direction(i, j, k, dim, dx, hb_arr, ht_arr, kappa0)) {
           if (fabs (mv(i,j,k,dim)) > nmax_arr(i,j,k,0)) {
             kappa_arr(i,j,k,0) = kappa0;
             nmax_arr(i,j,k,0) = fabs (mv(i,j,k,dim));
           }
           //propagate the curvature
           Real orientation;
           Array4<Real const> const * hv = closest_height (i,j,k,dim,hb_arr,ht_arr,&orientation);
           for (int d = 0; d <= 1; d++) {
             Array<int,3> neighbor={i,j,k};
             neighbor[dim]+=d?-1:1;
             int *np=&neighbor[0];
             while (!CELL_IS_BOUNDARY(neighbor,bx.smallEnd(),bx.bigEnd()) &&
                    !CELL_IS_FULL(vof_arr(*np,*(np+1),*(np+2),0)) &&
                    closest_height (*np,*(np+1),*(np+2),dim,hb_arr,ht_arr,&orientation) == hv) {
               if (fabs (mv(*np,*(np+1),*(np+2),dim)) > nmax_arr(*np,*(np+1),*(np+2),0)) {
                 kappa_arr(*np,*(np+1),*(np+2),0) = kappa0;
                 nmax_arr(*np,*(np+1),*(np+2),0) = fabs (mv(*np,*(np+1),*(np+2),dim));
               }
               neighbor[dim]+=d?-1:1;
             }
           }
          }
         }
       }); //  ParallelFor
    }//end MFIter
  }

//remaining_curvatures
  for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Box const& bx = mfi.tilebox();
    Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
    Array4<Real const> const& mv = normal[lev].const_array(mfi);
    Array4<Real const> const& alpha_arr = alpha[lev].const_array(mfi);
    Array4<Real const> const& hb_arr = height[0].const_array(mfi);
    Array4<Real const> const& ht_arr = height[1].const_array(mfi);
    Array4<Real > const& kappa_arr = kappa.array(mfi);
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      auto fvol = vof_arr(i,j,k,0);

   /*if ((i==4||i==11)&&j==11&&k==0){
       int dddd;
       Print()<<"------------"<<"\n";
   }*/
      if (!CELL_IS_FULL(fvol)){
        if (kappa_arr(i,j,k,0)==VOF_NODATA){
          // try height function and paraboloid fitting
          Real kappa0= height_curvature_combined (i,j,k, dx,hb_arr,ht_arr,mv,alpha_arr);
          if (kappa0!=VOF_NODATA)
            kappa_arr(i,j,k,0)=kappa0;
          //else
            // try particle method (defined in partstr.H)
            //kappa_arr(i,j,k,0)= partstr_curvature (i,j,k,dx,problo,vof_arr,mv,alpha_arr);
        }
      }
    }); //  ParallelFor
  }//end MFIter
  //!!!!!!!!fixme: a temporary solution for the curvature!!!!!!!!!!!
  // fill value of ghost cells (BCs, MPI info.)
    kappa.FillBoundary(geom.periodicity());

// diffuse curvatures
  int iter = 0;
  if (iter >0){
    MultiFab temp_K(kappa.boxArray(), kappa.DistributionMap(),1,kappa.nGrow());
    //fixme: need to change for BCs
    temp_K.setVal(VOF_NODATA,0,1,kappa.nGrow());
    while (iter--){
     for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       Box const& bx = mfi.tilebox();
       Array4<Real > const& temp_arr = temp_K.array(mfi);
       Array4<Real > const& kappa_arr = kappa.array(mfi);
       ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
        if (kappa_arr(i,j,k,0)!=VOF_NODATA)
         temp_arr(i,j,k,0)=kappa_arr(i,j,k,0);
        else{
        Real sa=0., s=0.;
/*#if AMREX_SPACEDIM==3
      for (int dk = -1; dk <= 1; dk++)
#endif
      for (int dj = -1; dj <= 1; dj++)
      for (int di = -1; di <= 1; di++)
       if (di != 0|| dj != 0|| dk != 0) {
         int ni=i+di,nj=j+dj,nk=k+dk;
         if (kappa_arr(ni,nj,nk,0)!=VOF_NODATA){
           s += kappa_arr(ni,nj,nk,0);
           sa += 1.;
         }
       }*/
        Array<int,3>nei;
        for (int c = 0; c < AMREX_SPACEDIM; c++){
         nei[0]=i,nei[1]=j,nei[2]=k;
         for (int di = -1; di <= 1; di+=2){
           nei[c]=(c==0?i:c==1?j:k)+di;
           if (kappa_arr(nei[0],nei[1],nei[2],0)!=VOF_NODATA){
             s += kappa_arr(nei[0],nei[1],nei[2],0);
             sa += 1.;
           }
         }
        }
        if (sa > 0.)
          temp_arr(i,j,k,0)=s/sa;
        else
          temp_arr(i,j,k,0)=VOF_NODATA;
        }
       }); //  ParallelFor
     }//end MFIter
    }
    MultiFab::Copy(kappa, temp_K, 0, 0, 1, kappa.nGrow());
  }
//fit_curvatures using paraboloid fitting of the centroids of the
//reconstructed interface segments
  for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Box const& bx = mfi.tilebox();
    Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
    Array4<Real const> const& mv = normal[lev].const_array(mfi);
    Array4<Real const> const& alpha_arr = alpha[lev].const_array(mfi);
    Array4<Real const> const& hb_arr = height[0].const_array(mfi);
    Array4<Real const> const& ht_arr = height[1].const_array(mfi);
    Array4<Real > const& kappa_arr = kappa.array(mfi);
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
   /*if ((i==4||i==11)&&j==11&&k==0){
       int dddd;
       Print()<<"------------"<<kappa_arr(-1,0,0,0)<<"\n";
   }*/
      auto fvol = vof_arr(i,j,k,0);
      if (!CELL_IS_FULL(fvol)){
        if (kappa_arr(i,j,k,0)==VOF_NODATA){
          // paraboloid fitting    of the centroids of the reconstructed interface segments
         kappa_arr(i,j,k,0) = curvature_fit (i,j,k,dx,vof_arr,mv,alpha_arr);
        }
      }
    }); //  ParallelFor
  }//end MFIter

  //!!!!!!!!fix me: a temporary solution for the curvature!!!!!!!!!!!
  // fill value of ghost cells (BCs, MPI info.)
    kappa.FillBoundary(geom.periodicity());

/////////////////////////////////////////////////////////////////////////
///   The following is used to do statistics of calculated curvature
///   for numerical tests.
/////////////////////////////////////////////////////////////////////////
if (0){
  VofRange kappa_range;
  range_init(kappa_range);

  struct KappaPrint{
    Real kappa, angle;
    XDim3 center;
    int i,j,k;
    // Default constructor
    KappaPrint() : kappa(0), angle(0), center(), i(0), j(0), k(0) {}
    // Constructor to initialize the KappaPrint
    KappaPrint(Real ka, Real a, XDim3 o,int i, int j, int k): kappa(ka),angle(a),center(o),
              i(i),j(j),k(k){}
  // Copy assignment operator
    KappaPrint& operator=(const KappaPrint& other) {
        if (this != &other) { // self-assignment check
         kappa = other.kappa;
         angle = other.angle;
         i = other.i;
         j = other.j;
         k = other.k;
           for (int c = 0; c < AMREX_SPACEDIM; c++)
            (&center.x)[c]=(&other.center.x)[c];
         }
        return *this;
    }
  };
  Vector<KappaPrint> kout,removed_elements;
  Box const& domain  = geom.Domain();
  IntVect half= (domain.smallEnd()+domain.bigEnd())/2;
  Array<Real,AMREX_SPACEDIM> center{AMREX_D_DECL(0.5*(problo[0]+probhi[0]),
                                                 0.5*(problo[1]+probhi[1]),
                                                 0.5*(problo[2]+probhi[2]))};
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    Box const& bx = mfi.tilebox();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    Array4<Real const> const& vof_arr = vof_mf.const_array(mfi);
    Array4<Real const> const& mv = normal[lev].const_array(mfi);
    Array4<Real const> const& alpha_arr = alpha[lev].const_array(mfi);
    Array4<Real const> const& kappa_arr = kappa.const_array(mfi);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
     /* if(i==6&&j==4&&k==7)
         Print() <<" ---- "<<"("<<i<<","<<j<<","<<k<<") "
                 <<"Kappa "<<kappa_arr(i,j,k,0) <<"vof "<<vof_arr(i,j,k,0) <<"\n";    */
      if (kappa_arr(i,j,k,0)!=VOF_NODATA){
       range_add_value(kappa_range,kappa_arr(i,j,k,0));

       //if(k==half[2]){

       if (i==j){
         XDim3 m={AMREX_D_DECL(mv(i,j,k,0),mv(i,j,k,1),mv(i,j,k,2))},p;
         plane_area_center (m, alpha_arr(i,j,k,0),p);
         Real nn=0.;
         for (int c = 0; c < AMREX_SPACEDIM; c++){
           (&p.x)[c] = problo[c]+((c==0?i:c==1?j:k)+(&p.x)[c])*dx[c]-center[c];
           nn+=((&p.x)[c])*((&p.x)[c]);
         }
         nn=sqrt(nn);
         //Real angle =acos(p.x/nn)*180./PI+(p.y>0? 0.:180.);
         Real nnxy=(p.x>0.?1.:-1.)*sqrt(p.x*p.x+p.y*p.y);
         Real angle =acos(nnxy/nn)*180./PI+(p.z>0? 0.:180.);
        /* Print() <<" ---- "<<"("<<i<<","<<j<<","<<k<<") "
                 <<"Kappa "<<kappa_arr(i,j,k,0) <<"  "<<p.x<<"  "<<p.y<<" "<<p.z<<"\n";*/
         kout.emplace_back(kappa_arr(i,j,k,0),angle,p,i,j,k);
       }
      }

    }
    }
    }


    }
    int nn=kout.size();
    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();
    if (nprocs > 1){
     MPI_Datatype mpi_kappa_type;
    // Define the lengths of each block
     int lengths[6] = {1, 1, 3, 1, 1, 1};

    // Define the displacements of each block
     const MPI_Aint base = offsetof(KappaPrint, kappa);
     MPI_Aint disp[6] = {
        offsetof(KappaPrint, kappa) - base,
        offsetof(KappaPrint, angle) - base,
        offsetof(KappaPrint, center) - base,
        offsetof(KappaPrint, i) - base,
        offsetof(KappaPrint, j) - base,
        offsetof(KappaPrint, k) - base
     };

    // Define the types of each block
     MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT};

    // Create the MPI datatype
     MPI_Type_create_struct(6, lengths, disp, types, &mpi_kappa_type);
     MPI_Type_commit(&mpi_kappa_type);

// Gather data from all processes
     Vector<int> recvcounts(nprocs);
     Vector<int> displs(nprocs, 0);
     MPI_Gather(&nn, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
     for (int i = 1; i < nprocs; ++i)
       displs[i] = displs[i-1] + recvcounts[i-1];

     Vector<KappaPrint> all_data(displs[nprocs-1] + recvcounts[nprocs-1]);
     int kk=sizeof(KappaPrint),tt=sizeof(XDim3), dd=sizeof(all_data);
     MPI_Gatherv(kout.data(), kout.size(), mpi_kappa_type, all_data.data(),
                recvcounts.data(), displs.data(), mpi_kappa_type, 0, MPI_COMM_WORLD);
      kout=all_data;
    }

    if (myproc==0){
// Sort the vector by center.x from high to low
      std::sort(kout.begin(), kout.end(), [](const KappaPrint& a, const KappaPrint& b) {
        return a.center.x > b.center.x;
      });

    // Use remove_if and copy elements that match center.y<0. to removed_elements
      auto it = std::remove_if(kout.begin(), kout.end(), [&](const KappaPrint& kp) {
        if (kp.center.z < 0.) {
          removed_elements.push_back(kp);
          return true;
        }
        return false;
      });
    // Erase the removed elements from the original vector
      kout.erase(it, kout.end());
   // Append the removed elements back to the original vector
      kout.insert(kout.end(), removed_elements.begin(), removed_elements.end());


      Print()<<"# of interfacial cells"<<kout.size()<<"\n";
      FILE* o;
      o = fopen("curvature.dat", "w");
      fprintf(o, " angle,k,x,y,z\n");

      for (int i = 0; i < kout.size(); ++i) {
       fprintf(o, "%g,%g,%g,%g,%g,%d,%d,%d",kout[i].angle,fabs(kout[i].kappa),kout[i].center.x,
             kout[i].center.y, kout[i].center.z, kout[i].i,kout[i].j,kout[i].k);
       fprintf(o, "\n");
      }
      fclose(o);
    }

    domain_range_reduce(kappa_range);
    range_update (kappa_range);
    Print()<<"curvature min: "<<kappa_range.min<<"  max: "<<kappa_range.max<<"  mean: "<<kappa_range.mean
           <<"  stddev: "<<kappa_range.stddev<<"\n";
  }


}

/////////////////////////////////////////////////////////////////////////////////////
///////
///////              Advect the VOF tracer uing geometry-based scheme
///////
////////////////////////////////////////////////////////////////////////////////////

void
VolumeOfFluid::tracer_vof_advection(Vector<MultiFab*> const& tracer,
                                    AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                                 Vector<MultiFab const*> const& v_mac,
                                                 Vector<MultiFab const*> const& w_mac),
                                    Real dt)
{
    static int start = 0;
    //amrex::Print() << " VOF Level#" << finest_level<<"\n";


    // ***********************************************************************
    // Allocate space for the fluxes for vof advection
    // ***********************************************************************
    Vector<MultiFab> m_total_flux, vof_total_flux;
    //auto& ld = *v_incflo->m_leveldata[lev];
    for (int lev = 0; lev <= finest_level; ++lev) {
      m_total_flux.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(),
                                MFInfo(), v_incflo->Factory(lev));
      vof_total_flux.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(),
                                MFInfo(), v_incflo->Factory(lev));
    }
    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();

//The vof advection is to scheme is to use dimension-splitting i.e. advect the vof tracer
//along each dimension successively using a one-dimensional scheme.
    for (int lev = 0; lev <= finest_level; ++lev) {
      Geometry const& geom = v_incflo->Geom(lev);
      auto const& dx = geom.CellSizeArray();

//define the effective volume 'vol_eff' to consider the non-zero for the MAC velocity at
//the cell face in the sweep direction due to dimension splitting scheme (see Lörstad &Fuchs
// JCP, 200(2004),pp153-176).
      MultiFab vol_eff(v_incflo->grids[lev], v_incflo->dmap[lev], 1, v_incflo->nghost_state(),
                 MFInfo(), v_incflo->Factory(lev));
// set its initial value to be 1. when the MAC velocity is divergence free, 'vol_eff' will
// be still one after the sweep of all dimensions.
      vol_eff.setVal(1.0);
      for (int d = 0; d < AMREX_SPACEDIM; d++){
    // the starting direction of the sweep for i,j,k direction for vof advection is alternated
    // during the solution to minimize the errors associated with the sweep direction.
       int dir = (start+d)%AMREX_SPACEDIM;
       m_total_flux[lev].setVal(0.0);
       vof_total_flux[lev].setVal(0.0);
       MultiFab const * U_MF = dir < 1? u_mac[lev]:
#if AMREX_SPACEDIM == 3
                               dir >= 2? w_mac[lev]:
#endif
                               v_mac[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*U_MF,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
//The following note is not relevant.It helps me remember some key concepts of
//indexing space implemented in AMREX.
//
//note 'tracer' is a cell-centered MultiFab, so the index of the titlebox()
//(i.e., object 'bx' defined in the following) is cell-centered indexing space.
//In order to calculate the vof flux, we need to sweep along the cell faces,
//i.e. the node-centered indexing space.It is why we grow the titlebox() along the
//sweep direction by 1. So the grown box size (i.e., object 'bxg') can equivalently cover the same
//node-centered indexing space of the computational grid as required by calculation
//of flux over cell faces in the sweep direction. (a side issue is that the sweep can
// start with left/bottom/back face of the ghost cells of the computational domain, which is
//not necessary for the flux calculation. As we will apply BCs for ghost cells, it does not
//really matter.)

// use for calculating cell-centered MultiFabs
           Box const& bx = mfi.tilebox();
           //auto const& ijk_min= bx.smallEnd();
           //auto const& ijk_max= bx.bigEnd();
           // use for calculating the vof flux over the cell faces
           //Box const& bxg = amrex::grow(bx,IntVect::TheDimensionVector(dir));
           Array4<Real> const& vof = tracer[lev]->array(mfi);
           Array4<Real const> const& mv = normal[lev].const_array(mfi);
           Array4<Real const> const& al =  alpha[lev].const_array(mfi);
           Array4<Real> const& m_flux_arr =  m_total_flux[lev].array(mfi);
           Array4<Real> const& vof_flux_arr =  vof_total_flux[lev].array(mfi);
           Array4<Real> const& vof_eff_arr = vol_eff.array(mfi);
           Array4<Real const> const& vel_mac_arr =  U_MF->const_array(mfi);
           // calculate the vof flux by doing the scanning of the cell faces
           // i.e., loop through the node-centered MultiFab.
           ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
            /* if (i==4 &&j==8 &&k==4&& v_incflo->m_cur_time>1.99) {
            Print() <<" vof_advection---dir "<<dir<<"  "<<m_flux_arr(i,j,k)<<"  "
                      <<"("<<i<<","<<j<<","<<k<<")"<<"vof"<<"  "<<vof(i,j,k)<<"  "
                      <<"vof_flux"<<"  "<<vof_flux_arr(i,j,k)<< "\n";
            Print() <<" vof_advection---dir "<<dir<<"  "<<m_flux_arr(i,j-1,k)<<"  "
                      <<"("<<i<<","<<(j-1)<<","<<k<<")"<<"vof"<<"  "<<vof(i,j-1,k)<<"  "
                      <<"vof_flux"<<"  "<<vof_flux_arr(i,j-1,k)<< "\n";
             }*/

             Real un=vel_mac_arr(i,j,k)*dt/dx[dir], s = (un < 0) ? -1 : (un > 0);
             /* if (fabs (un) > 0.51) {
                Real x = Real(i+0.5)*dx[0];
                Real y = Real(j+0.5)*dx[1];
                Real z = Real(k+0.5)*dx[2];
                Print()<< "Warning: CFL "<<un<<" at ("<<x<<", "<<y<<", "<<z<<") is larger than 0.51!"<<"\n";
              }*/


             int det_u = -(s + 1.)/2., det_d=(s - 1.)/2.;
             Array <int, 3> index={i,j,k}, // store upwinding index
                          index_d={i,j,k}; // store downwinding index
             // index for upwinding/downwinding cell
             index[dir]+=det_u;index_d[dir]+=det_d;
             Real fvol = vof(index[0],index[1],index[2]), cf;
             if (fvol <= 0. || fvol >= 1.)
                cf = fvol;
             else{
             // the normal vector and alpha of the interface in upwinding cell
               Array <Real, AMREX_SPACEDIM> m_v ={AMREX_D_DECL(
                                               mv(index[0],index[1],index[2],0),
                                               mv(index[0],index[1],index[2],1),
                                               mv(index[0],index[1],index[2],2)
                                               )};
               Real alpha_v = al(index[0],index[1],index[2]);
               if (un < 0.) {
                   m_v[dir]=-m_v[dir];
                   alpha_v+=m_v[dir];
               }
               Array<Real, AMREX_SPACEDIM> q0={AMREX_D_DECL(0.,0.,0.)},q1={AMREX_D_DECL(1.,1.,1.)};
               q0[dir]=1.-fabs(un);
               for (int dd = 0; dd < AMREX_SPACEDIM; dd++) {
                  alpha_v -= m_v[dd]*q0[dd];
                  m_v[dd] *= q1[dd] - q0[dd];
               }
               cf = plane_volume (m_v, alpha_v);
             }
             //Make sure we just update the cells in the valid box
             //upwinding cells
             //if (AMREX_D_TERM(index[0]>=ijk_min[0] && index[0]<=ijk_max[0],
             //              && index[1]>=ijk_min[1] && index[1]<=ijk_max[1],
             //              && index[2]>=ijk_min[2] && index[2]<=ijk_max[2])){
                m_flux_arr(index[0],index[1],index[2]) -=fabs(un);
                vof_flux_arr(index[0],index[1],index[2]) -=fabs(un)*cf;
            // }
             //downstream cells
             //if (AMREX_D_TERM(index_d[0]>=ijk_min[0] && index_d[0]<=ijk_max[0],
             //              && index_d[1]>=ijk_min[1] && index_d[1]<=ijk_max[1],
             //              && index_d[2]>=ijk_min[2] && index_d[2]<=ijk_max[2])){
               m_flux_arr(index_d[0],index_d[1],index_d[2]) +=fabs(un);
               vof_flux_arr(index_d[0],index_d[1],index_d[2]) +=fabs(un)*cf;
             //}
            }); //  end ParallelFor

            //loop through cell-centered MultiFab to update their value
            //Box const& bxc = mfi.tilebox(IntVect::TheZeroVector());
            //ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           // {
            //  vof(AMREX_D_DECL(i,j,k))*=vof_eff_arr(AMREX_D_DECL(i,j,k));
            //  vof(AMREX_D_DECL(i,j,k))+=vof_flux_arr(AMREX_D_DECL(i,j,k));
            //  vof_eff_arr(AMREX_D_DECL(i,j,k))+= m_flux_arr(AMREX_D_DECL(i,j,k));
            //  Real f = vof(AMREX_D_DECL(i,j,k))/vof_eff_arr(AMREX_D_DECL(i,j,k));
            //  vof(AMREX_D_DECL(i,j,k))= f< 1e-10? 0.:f>1.-1e-10? 1.:f;
           //}); //  end ParallelFor
        }// end MFIter
        //fix me: temporary solution for MPI boundary
        m_total_flux[lev].FillBoundary(geom.periodicity());
        vof_total_flux[lev].FillBoundary(geom.periodicity());


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*tracer[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //loop through cell-centered MultiFab to update their value
            //Box const& bxc = mfi.tilebox(IntVect::TheZeroVector());
           Box const& bx = mfi.tilebox();
           Array4<Real> const& vof = tracer[lev]->array(mfi);
           Array4<Real> const& m_flux_arr =  m_total_flux[lev].array(mfi);
           Array4<Real> const& vof_flux_arr =  vof_total_flux[lev].array(mfi);
           Array4<Real> const& vof_eff_arr = vol_eff.array(mfi);
           ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
              vof(i,j,k)*=vof_eff_arr(i,j,k);
              vof(i,j,k)+=vof_flux_arr(i,j,k);
              vof_eff_arr(i,j,k)+= m_flux_arr(i,j,k);
              Real f = vof(i,j,k)/vof_eff_arr(i,j,k);
              vof(i,j,k)= f< 1e-10? 0.:f>1.-1e-10? 1.:f;
              /*if (f > 0. && f < 1.)
              Print() <<" vof_advection---dir "<<dir<<"  "<<vof_eff_arr(i,j,k)<<"  "
                      <<"("<<i<<","<<j<<","<<k<<")"<<"vof"<<"  "<<f<<"  "
                      <<"vof_flux"<<"  "<<vof_flux_arr(i,j,k)<< "\n";*/
           }); //  end ParallelFor
        }// end MFIter
        //fix me: temporary solution for MPI boundary
        tracer[lev]->FillBoundary(geom.periodicity());



        // update the normal and alpha of the plane in each interface cell after each sweep
        tracer_vof_update (lev, *tracer[lev], height[lev]);
     }// end i-,j-,k-sweep: calculation of vof advection
        curvature_calculation (lev, *tracer[lev], height[lev], kappa[lev]);
  }// end lev
  start = (start + 1) % AMREX_SPACEDIM;

}
// add surface tension force (F) to the MAC velocity at the center of cell faces.
// F = dt*sigma*kappa*grad(VOF)/rho
// Umac <- Uma-F, note minus sign before F because of the way which curvature being calculated.
//   kappa and rho need to be estimated at the cell face center. The simple average of the cell-centered
//      values of two neighboring cells dilimited by the face is used to calculate the face-centered value.
//   grad(VOF) is also estimated at the face center using the center-difference method for two cells
//   i.e., in x-dir, grad(VOF)= (VOFcell[1]-VOFcell[0])/dx[0].
// u_mac/v_mac/w_mac stores the face-centered velocity (MAC).
//
// gu_mac/gv_mac/gw_mac stores the face-centered value of F/dt (i.e., sigma*kappa*grad(VOF)/rho)
// gu_mac/gv_mac/gw_mac will be averaged to the cell center when correcting the cell-centered velocity
// after final cell-centered projection.
void
VolumeOfFluid:: velocity_face_source (int lev, Real dt, AMREX_D_DECL(MultiFab& u_mac, MultiFab& v_mac, MultiFab& w_mac),
                                      AMREX_D_DECL(MultiFab* gu_mac, MultiFab* gv_mac, MultiFab* gw_mac))
{
   auto& ld = *v_incflo->m_leveldata[lev];
   Geometry const& geom = v_incflo->Geom(lev);
   auto const& dx = geom.CellSizeArray();
   Real sigma = v_incflo->m_sigma[0];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(kappa[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
       Box const& bx = mfi.tilebox();
       Box const& xbx = mfi.nodaltilebox(0);
       Box const& ybx = mfi.nodaltilebox(1);
       Box const& zbx = mfi.nodaltilebox(2);
       Array4<Real const> const& rho   = ld.density.const_array(mfi);
       Array4<Real const> const& tra   = ld.tracer.const_array(mfi);
       Array4<Real const> const& kap   = kappa[lev].const_array(mfi);
       AMREX_D_TERM(Array4<Real > const& umac = u_mac.array(mfi);,
                    Array4<Real > const& vmac = v_mac.array(mfi);,
                    Array4<Real > const& wmac = w_mac.array(mfi););
       AMREX_D_TERM(Array4<Real > gumac;, Array4<Real > gvmac;,Array4<Real > gwmac;);
       AMREX_D_TERM(if(gu_mac) gumac = gu_mac->array(mfi);,
                    if(gv_mac) gvmac = gv_mac->array(mfi);,
                    if(gw_mac) gwmac = gw_mac->array(mfi););
       ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i-1,j,k,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i-1,j,k,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i-1,j,k,0)!=VOF_NODATA)
            kaf=kap(i-1,j,k);
         else
            kaf=0.;
        // density estimated at the face center
          //Print()<<rho(i,j,k)<<"\n";
          //Print()<<tra(i-1,j-1,k-1)<<"\n";
          //Print()<<umac(i,j,k)<<"\n";
          umac(i,j,k) -= dt*kaf*sigma*2./(rho(i,j,k)+rho(i-1,j,k))*(tra(i,j,k,0)-tra(i-1,j,k,0))/dx[0];
          if (gu_mac)
             gumac(i,j,k)=kaf*sigma*2./(rho(i,j,k)+rho(i-1,j,k))*(tra(i,j,k,0)-tra(i-1,j,k,0))/dx[0];
       });

       ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j-1,k,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j-1,k,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i,j-1,k,0)!=VOF_NODATA)
            kaf=kap(i,j-1,k);
         else
            kaf=0.;
          vmac(i,j,k) -= dt*kaf*sigma*2./(rho(i,j,k)+rho(i,j-1,k))*(tra(i,j,k,0)-tra(i,j-1,k,0))/dx[1];
          if(gv_mac)
              gvmac(i,j,k) = kaf*sigma*2./(rho(i,j,k)+rho(i,j-1,k))*(tra(i,j,k,0)-tra(i,j-1,k,0))/dx[1];
       });
#if AMREX_SPACEDIM==3 /* 3D */
       ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         Real kaf;
         if(kap(i,j,k,0)!=VOF_NODATA && kap(i,j,k-1,0)!=VOF_NODATA)
            kaf=Real(0.5)*(kap(i,j,k,0)+kap(i,j,k-1,0));
         else if (kap(i,j,k,0)!=VOF_NODATA)
            kaf=kap(i,j,k);
         else if (kap(i,j,k-1,0)!=VOF_NODATA)
            kaf=kap(i,j,k-1);
         else
            kaf=0.;
         wmac(i,j,k) -= dt*kaf*sigma*2./(rho(i,j,k)+rho(i,j,k-1))*(tra(i,j,k,0)-tra(i,j,k-1,0))/dx[2];
         if(gw_mac)
            gwmac(i,j,k) = kaf*sigma*2./(rho(i,j,k)+rho(i,j,k-1))*(tra(i,j,k,0)-tra(i,j,k-1,0))/dx[2];
       });
#endif
    }

}

void VolumeOfFluid:: variable_filtered(MultiFab & variable)
{
    const auto& ba = variable.boxArray();
    const auto& dm = variable.DistributionMap();
    const auto& fact = variable.Factory();
    MultiFab node_val(amrex::convert(ba,IntVect::TheNodeVector()),dm, 1, 0 , MFInfo(), fact);
    MultiFab center_val(ba,dm,1,0,MFInfo(), fact);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(variable,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
         Box const& nbx = surroundingNodes(mfi.tilebox());
         Array4<Real > const& nv = node_val.array(mfi);
         Array4<Real const> const& var   = variable.const_array(mfi);
         ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
           nv(i,j,k,0)=0.;
           int nt=0, nrho=0, detk;
#if AMREX_SPACEDIM==3
           for (detk = 0; detk > -2; --detk)
#endif
            for (int detj = 0; detj > -2; --detj)
             for (int deti = 0; deti > -2; --deti) {
               Array<int,3> in{i+deti,j+detj,k+detk};
               //averaging to nodes
                 nv(i,j,k,0)+= var(in[0],in[1],in[2],0);
                 nrho++;
             }
           nv(i,j,k,0)/= Real(nrho);
      //    Print()<<i<<","<<j<<","<<k<<" "<<nv(i,j,k,0)<<"\n";
         });
       }
       average_node_to_cellcenter(center_val, 0, node_val, 0, 1,0);
       Copy(variable, center_val , 0, 0, 1, 0);
}


//////////////////////////////////////////////////////////////////////////////////
///////
///////  Initialize the VOF value using the EB implicit surface function
///////
/////////////////////////////////////////////////////////////////////////////////
Real myFunction(AMREX_D_DECL(Real x, Real y, Real z))
{
    return y - 0.01*cos (2.*3.14159265357*x);
}

void
VolumeOfFluid::tracer_vof_init_fraction (int lev, MultiFab& a_tracer)
{
    a_tracer.setVal(0.0);
    int vof_init_with_eb = 1;
    ParmParse pp("incflo");
    pp.query("vof_init_with_eb", vof_init_with_eb);

    Geometry const& geom = v_incflo->Geom(lev);
    auto const& dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
    auto const& probhi = geom.ProbHiArray();

#ifdef AMREX_USE_EB
    if (vof_init_with_eb) {
        if (lev == 0) {
           /* Array<Real,AMREX_SPACEDIM> center{AMREX_D_DECL((problo[0]+.45),
                                                           (problo[1]+.1),
                                                           (problo[2]+.45))};*/
            Array<Real,AMREX_SPACEDIM> center1{AMREX_D_DECL((problo[0]+.5),
                                                           (problo[1]+.75),
                                                           (problo[2]+.35))};
            Array<Real,AMREX_SPACEDIM> center{AMREX_D_DECL(0.5*(problo[0]+probhi[0]),
                                                           0.5*(problo[1]+probhi[1]),
                                                           0.5*(problo[2]+probhi[2]))};
            Real radius = .2; //5.0*dx[0];
            bool fluid_is_inside = true;

            EB2::SphereIF my_sphere(radius, center, fluid_is_inside);
            EB2::SphereIF my_sphere1(radius, center1, fluid_is_inside);

            Array<Real,AMREX_SPACEDIM> radii{AMREX_D_DECL(.2, .3, .25)};
            EB2::EllipsoidIF my_ellipsoid(radii, center, fluid_is_inside);

    // Initialise cylinder parameters
    int direction = 2;
    Real height = 1.6;


    center[0]=0.5*(problo[0]+probhi[0]);
    center[1]=0.5*(problo[1]+probhi[1]);
    center[2]=0.5*(problo[1]+probhi[1]);
    // Build the Cylinder implficit function representing the curved walls
    //EB2::CylinderIF my_cyl(radius, height, direction, center, true);
    //radius = 8.0*dx[0];
    //EB2::CylinderIF my_cyl_1(radius, height, direction, center, fluid_is_inside);

    //box
   /* Array<Real,AMREX_SPACEDIM> low{AMREX_D_DECL((problo[0]+10.3*dx[0]),
                                                  (problo[1]+10.3*dx[1]),
                                                (problo[2]+10.3*dx[2]))};
    Array<Real,AMREX_SPACEDIM> high{AMREX_D_DECL((probhi[0]-11.2*dx[0]),
                                                 (probhi[1]-11.2*dx[1]),
                                                 (probhi[2]-11.2*dx[2]))};    */
   Array<Real,AMREX_SPACEDIM> low{AMREX_D_DECL( (0.5*(problo[0]+probhi[0])-.2),
                                                (0.5*(problo[1]+probhi[1])-.2),
                                                (0.5*(problo[2]+probhi[2])-.2))};
    Array<Real,AMREX_SPACEDIM> high{AMREX_D_DECL((0.5*(problo[0]+probhi[0])+.2),
                                                 (0.5*(problo[1]+probhi[1])+.2),
                                                 (0.5*(problo[2]+probhi[2])+.2))};
    auto my_box= EB2::BoxIF( low,  high, fluid_is_inside);
    //auto my_box=  EB2::rotate(EB2::BoxIF( low,  high, fluid_is_inside), .3, 1);
    //auto my_box1=  EB2::rotate(my_box, .3, 0);
    //auto my_box2=  EB2::rotate(my_box1, .2, 2);
    //auto two = EB2::makeIntersection(my_sphere, my_box);
   //auto two = EB2::makeComplement(EB2::makeUnion(my_cyl_1, my_cyl));
   auto my_sin = EB2::DevicePtrIF(&myFunction);

    // Generate GeometryShop
    //auto gshop = EB2::makeShop(two);
    //auto gshop = EB2::makeShop(my_box);
    auto gshop = EB2::makeShop(my_sin);
   //auto gshop = EB2::makeShop(my_cyl);
            int max_level = v_incflo->maxLevel();
            EB2::Build(gshop, v_incflo->Geom(max_level), max_level, max_level);
        }

        auto fact = amrex::makeEBFabFactory(geom, a_tracer.boxArray(), a_tracer.DistributionMap(),
                                            {1,1,0}, EBSupport::volume);
        auto const& volfrac = fact->getVolFrac();
        MultiFab::Copy(a_tracer, volfrac, 0, 0, 1, 1);

        if (lev == v_incflo->finestLevel()) {
            EB2::IndexSpace::pop();
        }
    }
    else
#endif
    {

   struct VOFPrint{
  Real vof;
  int i,j,k;
  // Default constructor
    VOFPrint() : vof(0), i(0), j(0), k(0) {}
  // Constructor to initialize the VOFPrint
  VOFPrint(Real ka,int i, int j, int k): vof(ka),i(i),j(j),k(k){}
// Copy assignment operator
  VOFPrint& operator=(const VOFPrint& other) {
        if (this != &other) { // self-assignment check
            vof = other.vof;
            i = other.i;
            j = other.j;
            k = other.k;
        }
        return *this;
    }
  };
    Vector<VOFPrint> vout;
   // Define the file name
    std::string filename = "vof_value-32.dat";
    // Open the file
    std::ifstream infile(filename);

    if (!infile) {
        std::cerr << "Unable to open file " << filename << std::endl;
        exit;
    }
    // Read the file line by line
    std::string line;
    while (std::getline(infile, line)) {
       std::istringstream iss(line);
       int i, j, k;
       Real value;
       if (!(iss >> i >> j >> k >> value)) {
          std::cerr << "Error reading line: " << line << std::endl;
          continue;
       }
       vout.emplace_back(value,i,j,k);
    }
    infile.close();
#ifdef AMRE_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a_tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //Box const& vbx = mfi.validbox();
            Box const& vbx = amrex::grow(mfi.tilebox(),a_tracer.nGrow());
            auto const& tracer = a_tracer.array(mfi);

            for(int n=0;n<vout.size();++n)
              if(vout[n].i>=vbx.smallEnd()[0]&&vout[n].i<=vbx.bigEnd()[0]&&
                 vout[n].j>=vbx.smallEnd()[1]&&vout[n].j<=vbx.bigEnd()[1]
#if AMREX_SPACEDIM==2
                 &&vout[n].k==7)
#else
                 &&vout[n].k>=vbx.smallEnd()[2]&&vout[n].k<=vbx.bigEnd()[2])
#endif
              {
             //if(vbx.contains(vout[n].i,vout[n].j,vout[n].k)){
#if AMREX_SPACEDIM==2
                tracer(vout[n].i,vout[n].j,0,0)=vout[n].vof;
#else
                tracer(vout[n].i,vout[n].j,vout[n].k,0)=vout[n].vof;
#endif
              }

            /*amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = problo[0] + Real(i+0.5)*dx[0];
                Real y = problo[1] + Real(j+0.5)*dx[1];
                Real z = problo[2] + Real(k+0.5)*dx[2];

                Real r = std::sqrt(x*x + y*y + z*z);

                Real rad = 5.0*dx[0];
                Real dia = 10.0*dx[0];

                Real cenx = 0.5*(problo[0] + probhi[0]);
                Real ceny = 0.5*(problo[1] + probhi[1]);
                Real cenz = 0.5*(problo[2] + probhi[2]);

                Real xs = x - cenx;
                Real ys = y - ceny;
                Real zs = z - cenz;

                Real rs = (std::sqrt(xs*xs + ys*ys + zs*zs) - rad)/std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

                if (rs > 0.5) tracer(i,j,k) = 0.0;
                else if (rs < -.5) tracer(i,j,k) = 1.0;
                else tracer(i,j,k) = 0.5-rs;
            });*/
        }

    }

    // Once vof tracer is initialized, we calculate the normal direction and alpha of the plane segment
    // intersecting each interface cell.
    tracer_vof_update(lev, a_tracer, height[lev]);
    curvature_calculation(lev, a_tracer, height[lev], kappa[lev]);
    int n_tag=domain_tag_droplets (finest_level, v_incflo->grids,v_incflo->geom,
                                  v_incflo->get_tracer_new (),GetVecOfPtrs(tag));

}


///////////////////////////////////////////////////////////////////////////////////////////////
/////////
////////     output results in tecplot format
////////
///////////////////////////////////////////////////////////////////////////////////////////////
void
VolumeOfFluid::write_tecplot_surface(Real time, int nstep)
{

    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();
//    amrex::AllPrint() << " Output surface file at process#" << myproc<<"  " << nprocs << " at time " << time << std::endl;
    const std::string& tecplotfilename = amrex::Concatenate("tecplot_surface_", nstep)+"_";//+std::to_string(myproc);

    const int nfiles = 1;

    for (NFilesIter nfi(nfiles, tecplotfilename, false, true); nfi.ReadyToWrite(); ++nfi)
    {
      auto& TecplotFile = (std::ofstream&) nfi.Stream();
//    Print()<<"surface_plot"<< nstep<<"  "<< tecplotfilename<<"\n";

        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *v_incflo->m_leveldata[lev];
            Geometry const& geom = v_incflo->Geom(lev);
            Box const& domain  = geom.Domain();
            auto const& dx     = geom.CellSizeArray();
            auto const& problo = geom.ProbLoArray();
            auto const& probhi = geom.ProbHiArray();
            const BoxArray& ba = ld.tracer.boxArray();  //cell-centered multifab
            int nb = ba.size();
            const DistributionMapping& dm = ld.tracer.DistributionMap();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.validbox();
                const auto lo = lbound(bx);
                const auto hi = ubound(bx);
                Array4<Real const> const& vof = ld.tracer.const_array(mfi);
                Array4<Real const> const& mv = normal[lev].const_array(mfi);
                Array4<Real const> const& al =  alpha[lev].const_array(mfi);
                Array4<Real const> const& tag_arr = tag[lev].const_array(mfi);
                Array4<Real const> const& kappa_arr = kappa[lev].const_array(mfi);
                Vector<Segment> segments;
                int totalnodes = 0, k=0;
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    auto fvol = vof(i,j,k,0);

                    if (interface_cell (i,j,k, vof, fvol)){
                        Real   alpha;
                        XDim3    m, p, center;
                        for (int d = 0; d < AMREX_SPACEDIM; d++) {
                            (&m.x)[d]= mv(i,j,k,d);
                        }
                        alpha= al(i,j,k,0);
                        //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "<<" alpha "<<alpha<<" "<<fvol<<"\n";
                        plane_area_center(m, alpha, p);
                        /* convert the coord. of the center to the global sys*/
                        for (int dim = 0; dim < AMREX_SPACEDIM; dim++){
                            (&p.x)[dim] = problo[dim] + dx[dim]*((dim<1?i:dim<2?j:k)+(&p.x)[dim]);
                            (&center.x)[dim] = problo[dim] + dx[dim]*((dim<1?i:dim<2?j:k)+Real(0.5));
                        }
                        Array<Real, 3> vars={fvol,tag_arr(i,j,k),kappa_arr(i,j,k)};
                        /* Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "
                          " center "<< p <<" alpha "<<alpha<<" "<<fvol<<"\n";  */
                        //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" "<<fvol<<"\n";
                      //if (i==16 &&j==28&&k==3)
                        add_segment (center, dx, alpha, p, m, segments, totalnodes, vars);
                    }
                }
                }
                //           Print() << " Outside add_interface " << segments.size()<<" "<<totalnodes<<"\n";
                if (totalnodes > 0) {
                 // std::ofstream TecplotFile;
                //  TecplotFile.open(tecplotfilename, std::ios_base::trunc);
                   TecplotFile << "TITLE = \"incflow simulation from processor# " << myproc << "\" "<< "\n";
                   //spatial coordinates
                   TecplotFile << (AMREX_SPACEDIM== 2 ? "VARIABLES = \"X\", \"Y\"":"VARIABLES = \"X\", \"Y\", \"Z\"");
                   //output variables
                   TecplotFile <<", \"F\""<<", \"m_x\""<<", \"m_y\""
#if AMREX_SPACEDIM==3
                               <<", \"m_z\""
#endif
                               <<", \"alpha\""<<", \"tag\""<<", \"kappa\""<<"\n";
                   std::string zonetitle=("Level_"+std::to_string(lev)+
                                          "_Box_" +std::to_string(mfi.index())+
                                          "_Proc_"+std::to_string(myproc)+"_step_"+std::to_string(nstep));
                   TecplotFile <<(std::string("ZONE T=")+zonetitle);
                   TecplotFile <<", DATAPACKING=POINT"<<", NODES="<<totalnodes<<", ELEMENTS="<<segments.size()
                               << (AMREX_SPACEDIM== 2 ? ", ZONETYPE=FELINESEG":", ZONETYPE=FEQUADRILATERAL")
                               <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//      AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";
                   int nn=0;
                   for (const auto& seg : segments) {
                       for (int in = 0; in < seg.nnodes; in++)
                           TecplotFile <<seg.node[in].x<<" "<<seg.node[in].y<<" "
#if AMREX_SPACEDIM==3
                                       <<seg.node[in].z<<" "
#endif
                                       <<seg.vars[0]<<" "<<seg.mv.x<<"  "<<seg.mv.y<<"  "
#if AMREX_SPACEDIM==3
                                       <<seg.mv.z<<"  "
#endif
                                       <<seg.alpha<<"  "<<seg.vars[1]<<"  "<<seg.vars[2]<<"\n";

                   }
                   int inode = 1;
                   for (const auto& seg : segments) {
                       if (seg.nnodes == 4)
                           TecplotFile <<inode<<" "<<inode + 1<<" "<<inode + 2<<" "<<inode + 3<<"\n";
                       else if (seg.nnodes == 3)    /* triangular face */
                           TecplotFile <<inode<<" "<<inode + 1<<" "<<inode + 2<<" "<<inode + 2<<"\n";
                       else
                           TecplotFile <<inode<<" "<<inode + 1<<"\n";
                       inode += seg.nnodes;
                   }
                }
            } // end MFIter
        } // end lev
    } // NFiles
}

void VolumeOfFluid::WriteTecPlotFile(Real time, int nstep)
{
    BL_PROFILE("incflo::WriteTecPlotFile()");
    std::string m_tecplot_file{"tecplot_"};
    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();
    //amrex::AllPrint() << " Output file at process#" << myproc<<"  " << nprocs << " at time " << time << std::endl;
    const std::string& tecplotfilename = amrex::Concatenate(m_tecplot_file, nstep)+"_";

    const int nfiles = 1;

    for (NFilesIter nfi(nfiles, tecplotfilename, false, true); nfi.ReadyToWrite(); ++nfi)
    {
        auto& TecplotFile = (std::ofstream&) nfi.Stream();

        TecplotFile << "TITLE = \"incflow simulation from processor# " << myproc << "\" "<< "\n";
        //spatial coordinates
        TecplotFile << (AMREX_SPACEDIM== 2 ? "VARIABLES = \"X\", \"Y\"":"VARIABLES = \"X\", \"Y\", \"Z\"");
        //output variables
        TecplotFile <<", \"P\""<<", \"F\""<<", \"v_x\""<<", \"v_y\""<<
#if AMREX_SPACEDIM==3
                      ", \"v_z\""<<
#endif
                      ", \"m_x\""<<", \"m_y\""<<
#if AMREX_SPACEDIM==3
                      ", \"m_z\""<<
#endif
                      ", \"alpha\""<<", \"tag\""<<", \"hb_x\""<<", \"hb_y\""<<
#if AMREX_SPACEDIM==3
                      ", \"hb_z\""<<
#endif
                      ", \"ht_x\""<<", \"ht_y\""<<
#if AMREX_SPACEDIM==3
                      ", \"ht_z\""<<
#endif
                      ", \"kappa\""<<", \"rho\""<<", \"f_x\""<<", \"f_y\""<<
#if AMREX_SPACEDIM==3
                      ", \"f_z\""<<
#endif
                      "\n";

        for (int lev = 0; lev <= finest_level; ++lev) {
            bool m_use_cc_proj=v_incflo->m_use_cc_proj;
            auto& ld = *v_incflo->m_leveldata[lev];
            Geometry const& geom = v_incflo->Geom(lev);
            Box const& domain  = geom.Domain();
            auto const& dx     = geom.CellSizeArray();
            auto const& problo = geom.ProbLoArray();
            auto const& probhi = geom.ProbHiArray();
            auto const& ijk_min= domain.smallEnd();
            auto const& ijk_max= domain.bigEnd();
            const BoxArray& ba = ld.tracer.boxArray();  //cell-centered multifab
            int nb = ba.size();
            const DistributionMapping& dm = ld.tracer.DistributionMap();
            std::string IJK = "IJK";

          // Coordinate data
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

          //for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
             //Box const& bx = mfi.tilebox();
            for (MFIter mfi(ld.tracer); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.validbox();
                const auto lo = lbound(bx);
                const auto hi = ubound(bx);

                auto const& ijk_min= bx.smallEnd();
                auto const& ijk_max= bx.bigEnd();
                std::string zonetitle=("Level_"+std::to_string(lev)+"_Box_"+std::to_string(mfi.index())
                                        +"_Proc_"+std::to_string(myproc)+"_step_"+std::to_string(nstep));
                TecplotFile <<(std::string("ZONE T=")+zonetitle);
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
                  TecplotFile <<", "<<(IJK[dim]+std::string("="))<<(ijk_max[dim]-ijk_min[dim]+2);
                TecplotFile <<", DATAPACKING=BLOCK"<<", VARLOCATION=(["<<(m_use_cc_proj?AMREX_SPACEDIM+1:AMREX_SPACEDIM+2)<<"-"
                            <<(AMREX_SPACEDIM==3?24:19)<<"]=CELLCENTERED)"
                            <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//      AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";

                Array4<Real const> const& pa_nd = ld.p_nd.const_array(mfi);
                Array4<Real const> const& pa_cc = ld.p_cc.const_array(mfi);
                Array4<Real const> const& pa_mac = ld.mac_phi.const_array(mfi);
                Array4<Real const> const& tracer = ld.tracer.const_array(mfi);
                Array4<Real const> const& vel = ld.velocity.const_array(mfi);
                Array4<Real const> const& mv = normal[lev].const_array(mfi);
                Array4<Real const> const& al = alpha[lev].const_array(mfi);
                Array4<Real const> const& tag_arr = tag[lev].const_array(mfi);
                Array4<Real const> const& hb_arr = height[lev][0].const_array(mfi);
                Array4<Real const> const& ht_arr = height[lev][1].const_array(mfi);
                Array4<Real const> const& kappa_arr = kappa[lev].const_array(mfi);
                Array4<Real const> const& density_arr = ld.density.const_array(mfi);
                Array4<Real const> const& force_arr = force[lev].const_array(mfi);
                int nn=0, k=0;
                //write coordinate variables
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                  for (k = lo.z; k <= hi.z +1; ++k)
#endif
                  for (int j = lo.y; j <= hi.y +1; ++j) {
                  for (int i = lo.x; i <= hi.x +1; ++i) {
                      TecplotFile << (problo[dim]+dx[dim]*(dim<1?i:dim<2?j:k))<<" ";
                      ++nn;
                      if (nn > 100) {
                          TecplotFile <<"\n";
                          nn=0;
                      }
                  }
                  }

                }//
                //write pressure
                int nt=m_use_cc_proj?0:1;
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z+nt; ++k)
#endif
                for (int j = lo.y; j <= hi.y+nt; ++j) {
                for (int i = lo.x; i <= hi.x+nt; ++i) {
                    TecplotFile << (m_use_cc_proj?pa_cc(i,j,k):pa_nd(i,j,k))<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }

                //write VOF
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    TecplotFile << tracer(i,j,k,0)<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }

                //write velocity
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                    for (k = lo.z; k <= hi.z; ++k)
#endif
                    for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        TecplotFile << vel(i,j,k,dim)<<" ";
                        ++nn;
                        if (nn > 100) {
                            TecplotFile <<"\n";
                            nn=0;
                        }
                    }
                    }

                }//

                //write variables of the normal direction of the interface
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                    for (k = lo.z; k <= hi.z; ++k)
#endif
                    for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        TecplotFile << mv(i,j,k,dim)<<" ";
                        ++nn;
                        if (nn > 100) {
                            TecplotFile <<"\n";
                            nn=0;
                        }
                    }
                    }

                }//

                //write alpha of the interface
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    TecplotFile << al(i,j,k)<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }


                //write id of the droplets or bubbles
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    TecplotFile << tag_arr(i,j,k)<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }

                //write height function values
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                    for (k = lo.z; k <= hi.z; ++k)
#endif
                    for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        TecplotFile << hb_arr(i,j,k,dim)<<" ";
                        ++nn;
                        if (nn > 100) {
                            TecplotFile <<"\n";
                            nn=0;
                        }
                    }
                    }
                }//

                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                    for (k = lo.z; k <= hi.z; ++k)
#endif
                    for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        TecplotFile << ht_arr(i,j,k,dim)<<" ";
                        ++nn;
                        if (nn > 100) {
                            TecplotFile <<"\n";
                            nn=0;
                        }
                    }
                    }
                }//
                //write curvature
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    TecplotFile << kappa_arr(i,j,k)<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }

                //write density
#if AMREX_SPACEDIM==3
                for (k = lo.z; k <= hi.z; ++k)
#endif
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    TecplotFile << density_arr(i,j,k)<<" ";
                    ++nn;
                    if (nn > 100) {
                        TecplotFile <<"\n";
                        nn=0;
                    }
                }
                }


                //write force vector
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
#if AMREX_SPACEDIM==3
                    for (k = lo.z; k <= hi.z; ++k)
#endif
                    for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        TecplotFile << force_arr(i,j,k,dim)<<" ";
                        ++nn;
                        if (nn > 100) {
                            TecplotFile <<"\n";
                            nn=0;
                        }
                    }
                    }

                }//

                TecplotFile <<"\n";
            } // end MFIter

        } // end lev
    }


}

#include <queue>


#define ORTHOGONAL_COMPONENT(c)   (((c) + 1) % AMREX_SPACEDIM)

// @touch defines the touching connectivity. This function updates
// @touch with the info that region tagged with @tag1 touches the
// region tagged with @tag2  */
static void touching_regions (int tag1, int tag2, int * touch)
{
  if (tag2 < tag1) {
    int tmp = tag1;
    tag1 = tag2;
    tag2 = tmp;
  }
  else if (tag2 == tag1)
    return;
  int ntag = touch[tag2];
  if (ntag == tag1)
    return;
  if (ntag == 0)
    touch[tag2] = tag1;
  else {
    if (tag1 < ntag)
      touch[tag2] = tag1;
    touching_regions (tag1, ntag, touch);
  }
}

static void reduce_touching_regions (void * in, void * inout, int * len, MPI_Datatype * type)
{
  int * ltouch = (int *) in;
  int * gtouch = (int *) inout;
  int i;

  for (i = 1; i < *len; i++)
    if (ltouch[i] > 0)
      touching_regions (i, ltouch[i], gtouch);
}
////////////////////////////////////////////////////////////////////////////////////////
//                  domain_tag_droplets
// Simple flood fill algorithm is used to find the cells that belong to the same droplet.
//
// Fills the @tag variable of the cells of @domain with the (strictly
// positive) index of the droplet they belong to. The cells belonging
// to the background phase have an index of zero.
//
// Note that the volume fraction @c must be defined on all levels.
//
// Returns: the number of droplets.
////////////////////////////////////////////////////////////////////////////////

int VolumeOfFluid::domain_tag_droplets (int finest_level, Vector<BoxArray > const &grids, Vector<Geometry > const& geom,
                         Vector<MultiFab*> const& vof,Vector<MultiFab*> const& tag)
{

  //fix me: a temporary solution for one level mesh
  int ntag = 0;
  for (int lev = 0; lev <= finest_level; ++lev) {
    tag[lev]->setVal(0.);
    auto const& dx = geom[lev].CellSizeArray();
    bool touching = false;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vof[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       Box const& bx = mfi.tilebox();
//       const auto lo = lbound(bx);
//       const auto hi = ubound(bx);
       //note: we use the limit of the validbox to avoid include the ghost
       //cells in the FIFO queue
       Box const& bxv = mfi.validbox();
       auto const& ijk_min= bxv.smallEnd();
       auto const& ijk_max= bxv.bigEnd();
       Array4<Real const> const& vof_arr = vof[lev]->const_array(mfi);
       Array4<Real> const& tag_arr = tag[lev]->array(mfi);
       //fix me: not compatible with GPUs
       ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
       /* if(i==7&&j==11&&k==11){
            Print()<<"-----  "<<vof_arr(i,j,k)
            <<"  "<<"\n";
        }*/

        if(vof_arr(i,j,k)>EPS && tag_arr(i,j,k)==0.){
          std::queue <Array<int,3>> fifo;
          tag_arr(i,j,k)=++ntag;
          fifo.push({i,j,k});
          while (!fifo.empty()){
           Array<int,3> cell=fifo.front(),ncell=cell;
           for (int d = 0; d < AMREX_SPACEDIM; d++)
             for (int det =-1;det<=1;det+=2){
               ncell[d]=cell[d]+det;
               /*Print()<<vof_arr(ncell[0],ncell[1],ncell[2])<<"  "
                      <<(!CELL_IS_BOUNDARY(ncell,ijk_min,ijk_max))<<"  "
                      <<tag_arr(ncell[0],ncell[1],ncell[2])<<"\n";*/
               if(vof_arr(ncell[0],ncell[1],ncell[2])>EPS&&
                  !CELL_IS_BOUNDARY(ncell,ijk_min,ijk_max)&&
                  tag_arr(ncell[0],ncell[1],ncell[2])==0.){
                  tag_arr(ncell[0],ncell[1],ncell[2])=ntag;
                  fifo.push(ncell);
               }
             }
             fifo.pop();
          }// end while
        } // end if
       });

    }// end MFIter
 // the rest of the algorithm deals with periodic and parallel BCs
    if (ParallelDescriptor::NProcs() > 1){
      int myproc = ParallelDescriptor::MyProc();
      int nprocs = ParallelDescriptor::NProcs();
      int tags[nprocs];
      MPI_Allgather (&ntag, 1, MPI_INT, tags, 1, MPI_INT, MPI_COMM_WORLD);
    // tags[] now contains the 'ntag' value on each PE
      int i;
      ntag = 0;
      for (i = 0; i < nprocs; i++)
       ntag += tags[i];
    // shift tag values to get a single tag space across all PEs
      if (myproc > 0) {
        int tagshift = 0;
        for (i = 0; i < myproc; i++)
          tagshift += tags[i];
        //AllPrint()<<"tagshift---  "<<tagshift<<" processid "<<myproc<<"\n";

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vof[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
           Box const& bx = mfi.tilebox();
           Array4<Real> const& tag_arr = tag[lev]->array(mfi);
           ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
            if (tag_arr(i,j,k)>0)
              tag_arr(i,j,k) +=tagshift;
           });
        }
      }
    }//end if (algorithm to deal with the parallel process)
     //fix me: temporary solution for MPI boundary
    tag[lev]->FillBoundary(geom[lev].periodicity());
    int touch[ntag + 1]={};
//We search the cells in the box boundaries to determine if the tag value
//of the cell and tag value of its neighboring ghost cell are connected by
//the same droplet/bubble.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vof[lev]); mfi.isValid(); ++mfi)
    {
       Box const& bx = mfi.validbox();
       auto const& ijk_min= bx.smallEnd();
       auto const& ijk_max= bx.bigEnd();
       Array4<Real> const& tag_arr = tag[lev]->array(mfi);
       Array<IntVect,2> dim_limit={ijk_min,ijk_max};
       // we search the cells on the 6 boundaries of the 3D box.
       for (int d = 0; d < AMREX_SPACEDIM; d++){
        int ort1=ORTHOGONAL_COMPONENT(d)  // 1st transverse direction
           ,ort2=ORTHOGONAL_COMPONENT(ort1);//2nd transverse direction
        for (int n=0;n<2;n++){
         int k0=dim_limit[n][d],gd=k0+(n==0?-1:1);
         for(int i0=ijk_min[ort1];i0<=ijk_max[ort1];i0++){
#if AMREX_SPACEDIM==2 /*2D*/
            Real tag_cell=(d==0?tag_arr(k0,i0,0):tag_arr(i0,k0,0));
            if(tag_cell > 0){
              Real tag_gcell=(d==0?tag_arr(gd,i0,0):tag_arr(i0,gd,0));
              if(tag_gcell > 0)
                touching_regions (tag_cell, tag_gcell, touch);
            }
#else /*3D */
          for(int j0=ijk_min[ort2];j0<=ijk_max[ort2];j0++){
            Real tag_cell=(d==0?tag_arr(k0,i0,j0):
                           d==1?tag_arr(j0,k0,i0):
                                tag_arr(i0,j0,k0));
            if(tag_cell > 0){
              Real tag_gcell=(d==0?tag_arr(gd,i0,j0):
                              d==1?tag_arr(j0,gd,i0):
                                   tag_arr(i0,j0,gd));
              if(tag_gcell > 0)
                touching_regions (tag_cell, tag_gcell, touch);
            }
          }// end for-loop for searching cells in the boundaries.
#endif
         }
        }// end for-loop for low and high boundary
       }// end for-loop for AMREX_SPACEDIM
    }
    if (ParallelDescriptor::NProcs() > 1){
      //int gtouch[ntag + 1]={};
      MPI_Op op;
      MPI_Op_create (reduce_touching_regions, false, &op);
      //MPI_Allreduce (touch, gtouch, ntag + 1, MPI_INT, op, MPI_COMM_WORLD);
      ParallelDescriptor::detail::DoAllReduce<int>(touch,op,ntag+1);
      MPI_Op_free (&op);
      //std::memcpy(touch, gtouch, sizeof(gtouch));
    }
    /*Print()<<"-----touching-----  "<<"\n";
    for (int i=0;i<ntag+1;i++)
        Print()<<i<<" "<<touch[i]<<" "<<"\n";
    Print()<<"-----touching-----  "<<"\n";*/
// find region with smallest tag touching each region i.e. reduces
// the chain of touching tags
    int  maxtag = 0;
    for (int i = 1; i <= ntag; i++) {
      int itouch = touch[i];
      while (itouch > 0) {
       touch[i] = itouch;
       itouch = touch[itouch];
       touching = true;
      }
      if (touch[i] == 0 && i > maxtag)
       maxtag = i;
    }
//  fix touching regions
    if (touching) {
      int n_tag = 0; /* fresh tag index */
      int ntags[maxtag + 1];
      ntags[0] = 0;
      for (int i = 1; i <= maxtag; i++)
        if (touch[i] == 0) { /* this region is not touching any other */
         touch[i] = i;
         ntags[i] = ++n_tag;
        }
      maxtag = n_tag;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*vof[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& tag_arr = tag[lev]->array(mfi);
        //fix me: not compatible with GPUs
        ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           int ttag=tag_arr(i,j,k);
           tag_arr(i,j,k)=ntags[touch[ttag]];
        });

      }
     //fix me: temporary solution for MPI boundary
     tag[lev]->FillBoundary(geom[lev].periodicity());
     ntag=maxtag;
    }

  }//end lev

  return ntag;

}



//////////////////////////////////////////////////////////////////////////////
////
////              Computing sums for each droplet.
////
//////////////////////////////////////////////////////////////////////////////

void VolumeOfFluid::output_droplet (Real time, int nstep)
{
    static int first = 1;
    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();
    const std::string& filename = "droplet_his.dat";
    if(nstep%output_drop_frequence!=0)
        return;

    for (int lev = 0; lev <= finest_level; ++lev) {
      auto& ld = *v_incflo->m_leveldata[lev];
      Geometry const& geom = v_incflo->Geom(lev);
      auto const& dx = geom.CellSizeArray();
      auto const& problo = geom.ProbLoArray();
      auto const& probhi = geom.ProbHiArray();

      auto tvol = ld.tracer.sum()*AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
      //Print() <<" total vof--- "<<"   "<< tvol<<"\n";

      int n_tag=domain_tag_droplets (finest_level, v_incflo->grids,
                v_incflo->geom, v_incflo->get_tracer_new (),GetVecOfPtrs(tag));
      Print()<<"number of droplets  "<< n_tag<<"\n";
      // 'mcent' is mass center of each droplet
      Real vols[n_tag], vels[n_tag], mcent[AMREX_SPACEDIM][n_tag],surfA[n_tag];
      int ncell[n_tag];
      // find the max and min location of interface and curvature*/
      VofRange s[AMREX_SPACEDIM][n_tag], kappa_range[n_tag];
      // the range of location of interfacial cells */
      Real vtop[n_tag], range[AMREX_SPACEDIM][2][n_tag];
      for (int n = 0; n < n_tag; n++){
       ncell[n]=0; vols[n] = 0.; vels[n] = 0.; surfA[n]=0.; vtop[n] = 0.;
       range_init (kappa_range[n]);
       for(int d = 0; d < AMREX_SPACEDIM; d++) {
         mcent[d][n]=0.;
         range_init (s[d][n]);
         range[d][0][n] = 0.;
         range[d][1][n] = 0.;
       }
      }
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
       Box const& bx = mfi.tilebox();
       const auto lo = lbound(bx);
       const auto hi = ubound(bx);
       Array4<Real const> const& tag_arr = tag[lev].const_array(mfi);
       Array4<Real const> const& fv = ld.tracer.const_array(mfi);
       Array4<Real const> const& vel_arr = ld.velocity.const_array(mfi);
       Array4<Real const> const& mv = normal[lev].const_array(mfi);
       Array4<Real const> const& al =  alpha[lev].const_array(mfi);
       Array4<Real const> const& ka =  kappa[lev].const_array(mfi);
       //fix me: not compatible with GPUs
      // ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
     //  {
       for (int k = lo.z; k <= hi.z; ++k) {
       for (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) {
        int itag=tag_arr(i,j,k);
        if(itag > 0){
         /* count number */
         ncell[itag - 1]++;
         /* calculate the cell volume */
         Real dV = AMREX_D_TERM(dx[0],*dx[1],*dx[2])*fv(i,j,k);
         vols[itag - 1] += dV;
         /* calculate the momentum */
         Real vel = sqrt(AMREX_D_TERM(vel_arr(i,j,k,0)*vel_arr(i,j,k,0),
                                     +vel_arr(i,j,k,1)*vel_arr(i,j,k,1),
                                     +vel_arr(i,j,k,2)*vel_arr(i,j,k,2)));
         vels[itag - 1] += vel * dV;
         XDim3 center;
         for (int d = 0; d < AMREX_SPACEDIM; d++)
           (&center.x)[d] = problo[d] + dx[d]*((d<1?i:d<2?j:k)+Real(0.5));
         /* calculate the mass center */
         for(int d = 0; d < AMREX_SPACEDIM; d++)
          mcent[d][itag - 1] += (&center.x)[d] * dV;

         if (interface_cell (i,j,k, fv, fv(i,j,k))){
           Real   alpha;
           XDim3  m, p;
           for (int d = 0; d < AMREX_SPACEDIM; d++)
              (&m.x)[d]= mv(i,j,k,d);
           alpha= al(i,j,k,0);
           //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "<<" alpha "<<alpha<<" "<<fvol<<"\n";
           Real area=plane_area_center(m, alpha, p);
           surfA[itag-1] += area*AMREX_D_TERM(1,*dx[0],*dx[1]);
           /* convert the coord. of the center to the global sys*/
           for (int d = 0; d < AMREX_SPACEDIM; d++)
              (&p.x)[d] = problo[d] + dx[d]*((d<1?i:d<2?j:k)+(&p.x)[d]);
           for(int d = 0; d < AMREX_SPACEDIM; d++){
               //fixme: just for testing
               if (fabs(p.x)<dx[0])
              range_add_value (s[d][itag-1], (&p.x)[d]);
           }
           // do statistics of the curvature data
           range_add_value (kappa_range[itag-1], ka(i,j,k,0));
         }
        }
//       });
       }}}    //end of the ijk-loop

    }//end MFIter
    for (int n = 0; n < n_tag; n++)
      range_update (kappa_range[n]);
 // the rest of the algorithm deals with  parallel BCs
    if (ParallelDescriptor::NProcs()> 1){
      Real sum[n_tag];
      /*sum number of cells of each drop from different pid*/
      ParallelDescriptor::ReduceIntSum(ncell,n_tag);
      /*sum drop volume from different pid*/
      ParallelDescriptor::ReduceRealSum (vols, n_tag);
      //Print()<<"drople volume "<<vols[0]<<"\n";
      /*sum interfacial area from different pid*/
      ParallelDescriptor::ReduceRealSum (surfA, n_tag);
      /*sum drop momentum from different pid*/
      ParallelDescriptor::ReduceRealSum (vels, n_tag);

      /*sum mass center from different pid*/
      for (int d = 0; d < AMREX_SPACEDIM; d++) {
       ParallelDescriptor::ReduceRealSum (mcent[d], n_tag);
       for (int n = 0; n < n_tag; n++)
         domain_range_reduce(s[d][n]);
      }
      // sum curvature info.
      for (int n = 0; n < n_tag; n++){
         domain_range_reduce(kappa_range[n]);
         range_update (kappa_range[n]);
      }
    }
////////////////////////////////////////////////////////////////////
//////
//////   Estimate the phase error for cube case
//////
////////////////////////////////////////////////////////////////////

////Translation of a cube
Real error=0.;
if (0){
     Real lencube=5./16., cube_vol=lencube*lencube*lencube;
     Array<Real, AMREX_SPACEDIM> o0={0.1875},o,
                                 cube_min,cube_max;
     //theoretical centroid of regid body movement
     for (int d = 0; d < AMREX_SPACEDIM; d++){
        o[d] = o0[d]+1.0*time;
        cube_min[d]=o[d]-lencube*.5;
        cube_max[d]=o[d]+lencube*.5;
        int np=cube_min[d]/(probhi[d]-problo[d]);
        cube_min[d]-=np*(probhi[d]-problo[d]);
        np=cube_max[d]/(probhi[d]-problo[d]);
        cube_max[d]-=np*(probhi[d]-problo[d]);
     }
     Print()<<"cube center"<<o<<"\n";

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
       Box const& bx = mfi.tilebox();
       const auto lo = lbound(bx);
       const auto hi = ubound(bx);
       Array4<Real const> const& fv = ld.tracer.const_array(mfi);

       //fix me: not compatible with GPUs
       ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         if (fv(i,j,k)>0.){
           Real vertex[3][8];
           //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "<<" alpha "<<alpha<<" "<<fvol<<"\n";
           /* convert the coord. of the center to the global sys*/

           for (int d = 0; d < AMREX_SPACEDIM; d++){
            vertex[d][0]=problo[d] + dx[d]*((d<1?i:d<2?j:k));
            vertex[d][1]=problo[d] + dx[d]*((d<1?i+1:d<2?j:k));
            vertex[d][2]=problo[d] + dx[d]*((d<1?i+1:d<2?j+1:k));
            vertex[d][3]=problo[d] + dx[d]*((d<1?i:d<2?j+1:k));
            vertex[d][4]=problo[d] + dx[d]*((d<1?i:d<2?j:k+1));
            vertex[d][5]=problo[d] + dx[d]*((d<1?i+1:d<2?j:k+1));
            vertex[d][6]=problo[d] + dx[d]*((d<1?i+1:d<2?j+1:k+1));
            vertex[d][7]=problo[d] + dx[d]*((d<1?i:d<2?j+1:k+1));
           }
           Real vof=1.;
           for (int d = 0; d < AMREX_SPACEDIM; d++){
             int nn= 0;
             Real dd =0.;
             for (int n = 0; n < 8; n++){
               if (vertex[d][n]< cube_min[d]){
                 nn++;
                 dd+=(cube_min[d]-vertex[d][n])/dx[d];
               }
               else if(vertex[d][n]>cube_max[d]){
                 nn++;
                 dd+=(vertex[d][n]-cube_max[d])/dx[d];
               }
             }
             if (nn>0)
               vof*=(1.-dd/nn);
           }
           error+=fabs(fv(i,j,k)-vof)*AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
        }
       });

    }//end MFIter
   // if (time > .15){
   // Print()<<"test----"<<"\n";
   // }
    //AllPrint()<<"error ----"<< error<<"\n";
    ParallelDescriptor::ReduceRealSum (&error, 1);
    //Print()<<error<<"\n";
    error /=cube_vol;

}
/////Deformation of a sphere in sinusoidal velocity field
if (0){
   if (first)
    MultiFab::Copy(ld.tracer, ld.tracer, 0, 1, 1, 1);
//    MultiFab::Subtract (ld.tracer, ld.tracer, 0, 1, 1, 1);
    Real sphere_vol=4.*3.14159265357*.15*.15*.15/3;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
       Box const& bx = mfi.tilebox();
       Array4<Real const> const& fv = ld.tracer.const_array(mfi);

       //fix me: not compatible with GPUs
       ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         //if (fv(i,j,k)>0.){

       error+=fabs(fv(i,j,k,1)-fv(i,j,k,0))*AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
        //}
       });

      }//end MFIter
   ParallelDescriptor::ReduceRealSum (&error, 1);
   error /=sphere_vol;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* for parallel run, only master process outputs the data */
    if (ParallelDescriptor::IOProcessor()) {
    /* average values of each drop */
     for (int n = 0; n < n_tag; n++)
      if (vols[n] > 0.) {
        vels[n] /= vols[n];
        for (int d = 0; d < AMREX_SPACEDIM; d++)
          mcent[d][n] /= vols[n];
      }
    /* sort the drops */
     for (int n = 0; n < n_tag - 1; n++)
       for (int m = n + 1; m < n_tag; m++){
           Real xx,yy;
      //fix me: hard coded
         if (true){ /* sort according to the volume*/
           xx = vols[n], yy = vols[m];
         }
         else { /* sort according to the location*/
          //fix me: hard coded
       int c=0;
          if (true){ //positive direction
           xx = mcent[c][m],yy = mcent[c][n];
          }
          else {
           xx = mcent[c][n],yy = mcent[c][m];
          }
         }
         if (xx < yy) {
           int id = ncell[n];
           ncell[n] = ncell[m];
           ncell[m] = id;
           Real t = vols[n];
           vols[n] = vols[m]; vols[m] = t;
           t = vels[n];
           vels[n] = vels[m]; vels[m] = t;
           for (int d = 0; d < AMREX_SPACEDIM; d++) {
             t = mcent[d][n];
             mcent[d][n] = mcent[d][m];
             mcent[d][m] = t;
           }
           VofRange ss;
           for (int d = 0; d < AMREX_SPACEDIM; d++) {
             ss = s[d][n];
             s[d][n] = s[d][m];
             s[d][m] = ss;
           }
           t = surfA[n];
           surfA[n] = surfA[m]; surfA[m] = t;
         }
     } /*end of sort*/

     int ndrops = 0;
     Real sumtotal = 0.;
     for (int n = 0; n < n_tag; n++){
       ndrops++;
       sumtotal += vols[n];
     }
     for (int n = 0; n < n_tag; n++)
      for(int d = 0; d < AMREX_SPACEDIM; d++) {
        if (s[d][n].min<std::numeric_limits<Real>::max())
          range[d][0][n] = s[d][n].min;
        if (s[d][n].max>-std::numeric_limits<Real>::max())
          range[d][1][n] = s[d][n].max;
      }
     Real pos_limit[AMREX_SPACEDIM][2];
     for(int d = 0; d < AMREX_SPACEDIM; d++){
      pos_limit[d][0]=std::numeric_limits<Real>::max(),pos_limit[d][1]=std::numeric_limits<Real>::lowest();
      for (int n = 0; n < n_tag; n++){
        if (pos_limit[d][0]>range[d][0][n])
          pos_limit[d][0]=range[d][0][n];
        if (pos_limit[d][1]<range[d][1][n])
          pos_limit[d][1]=range[d][1][n];
      }
     }


     std::ofstream outputFile;
     if (first )
      outputFile.open(filename, std::ios_base::trunc);
     else
      outputFile.open(filename, std::ios_base::app);

// Check if file opened successfully
     if (!outputFile) {
       std::cerr << "Error opening file!" << std::endl;
       return;
     }
// If the file was just created, write the header
     if (outputFile.tellp() == 0) {
       outputFile << "The file contains the data of drops \n"
                 <<"Time, Num of Drops, Volume, Speed, Centroid, Range, Surface area, Curvature range (mean,min,max,stddev)\n";
     }

// Write data to file
     outputFile << std::fixed << std::setprecision(8) << time << ",   "<<ndrops<<",  "
                <<std::setprecision(8) << sumtotal<<",  "<<error;
     for (int n = 0; n < n_tag && n < 7; n++){
       outputFile <<", #"<<n+1<<",  "<<vols[n]<<",  "<<vels[n]<<",  "<<"L,  "
                  AMREX_D_TERM(<<mcent[0][n]<<",  ", <<mcent[1][n]<<",  ", <<mcent[2][n]<<",  ")
                  <<"R ";
       for(int d = 0; d < AMREX_SPACEDIM; d++)
         outputFile <<",  "<<range[d][0][n]<<",  "<<range[d][1][n];
       /*outputFile <<"s  "<<surfA[n]<<"  "
                  <<"k: "<<kappa_range[n].mean<<" "<<kappa_range[n].min<<" "<<kappa_range[n].max
                  <<" "<<kappa_range[n].stddev;    */
     }
     outputFile <<"\n";
     outputFile.close();
    }
    } // end lev
// Close the file
      first =0;
}

void VolumeOfFluid::apply_velocity_field (Real time, int nstep)
{


    for (int lev = 0; lev <= finest_level; ++lev) {
      auto& ld = *v_incflo->m_leveldata[lev];
      Geometry const& geom = v_incflo->Geom(lev);
      auto const& dx = geom.CellSizeArray();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();
        Array4<Real > const& vel = ld.velocity.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real x = Real(i+0.5)*dx[0];
          Real y = Real(j+0.5)*dx[1];
          Real z = Real(k+0.5)*dx[2];
          Real pi = 3.14159265357;
//          vel(i,j,k,0) = 2*sin(2.*pi*y)*sin(pi*x)*sin(pi*x)*sin(2*pi*z)*cos(pi*time/3.);
//          vel(i,j,k,1) = -sin(2.*pi*x)*sin(pi*y)*sin(pi*y)*sin(2*pi*z)*cos(pi*time/3.);
          vel(i,j,k,0) = sin(pi*x)*sin(pi*x)*sin(2*pi*y)*cos(pi*time/8.);
          vel(i,j,k,1) =-sin(pi*y)*sin(pi*y)*sin(2*pi*x)*cos(pi*time/8.);
#if (AMREX_SPACEDIM == 3)
          vel(i,j,k,2) = -sin(2.*pi*x)*sin(pi*z)*sin(pi*z)*sin(2*pi*y)*cos(pi*time/3.);
#endif
        });
     }//end MFIter
    } //end lev
}

