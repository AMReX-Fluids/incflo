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
#define vector_norm(v) (sqrt((v)->x*(v)->x + (v)->y*(v)->y + (v)->z*(v)->z))
#define CLAMP(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
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
    }

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
  int nnodes;               /* number of nodes (2, 3 or 4) */
#if AMREX_SPACEDIM==2
  XDim3 node[2];          /* node coordinates */
#else
  XDim3 node[4];
#endif
  XDim3 mv;
  Real alpha, vof, tag;
 // Constructor to initialize the Segment
  Segment(int n, Array<XDim3, 12> const& nodes, XDim3 m, Real a, Real f, Real t, int ns=0)
   : nnodes(n), mv (m), alpha (a), vof(f), tag(t) {
   for (int i = 0; i < n; ++i)
      node[i]= nodes[ns==0?i:(i + 3)%(n + 2)];
  }
};

static void add_segment (XDim3 const & center, GpuArray <Real,AMREX_SPACEDIM> const & dx,
                         Real alpha, XDim3 const & o, XDim3 const & m,
            Vector<Segment> & segments, int & nt, Real vof, Real tag)
{

   /*Print() <<" add_segment "<< *o<<"  "<<" vector "<<*m
                     <<"vof"<<"  "<<vof<<"  "<<"alpha " <<alpha<<"\n";  */

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
         Real length_diff1 = vector_norm (&diff1);
         if (length_diff1 < 1e-20) /*degenerated case*/
              return;
         Real max_sintheta = 0.;
     /*  cycle through all other nodes (jnode) where cut face intersects cell edges */
         for (inode2 = 1; inode2 < nnodecutface; inode2++) {
           int jnode = (inode + inode2)%nnodecutface;
           XDim3 diff2 = {nodecutface[jnode].x - node.x,
                            nodecutface[jnode].y - node.y,
                            nodecutface[jnode].z - node.z};
           Real length_diff2 = vector_norm (&diff2);
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
   segments.emplace_back(nnodecutface, nodecutface, m, alpha, vof, tag);
   //Print() << " normal direction " <<nnodecutface<<" "<<*m<<""<<vof<<"\n";
 }
 else { /* cut face must be divided into 2 quadrilateral/triangular faces */
   /* first face is quadrilateral */
   nt += 4;
   segments.emplace_back(4, nodecutface, m, alpha, vof, tag);
   /* second face is triangular if nnodecutface=5; otherwise quadrilateral */
//   WallFace * face2 = g_malloc (sizeof (WallFace));

   nnodecutface -= 2;

   nt += nnodecutface;
   segments.emplace_back(nnodecutface, nodecutface, m, alpha, vof, tag, 3);
 } /* cut face must be divided into 2 quadrilateral/triangular faces */

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


/**
 * plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #amrex::XDim3.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 *
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

#if AMREX_SPACEDIM == 2
#include "myc2D.h"
#else
#include "myc.h"
#endif

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
// we also need to
// calculate the normal and alpha of the cut plane.
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

void
VolumeOfFluid::tracer_vof_update(int lev, MultiFab & vof_mf)
{

    auto const& dx = v_incflo->geom[lev].CellSizeArray();
    for (MFIter mfi(vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real const> const& vof = vof_mf.const_array(mfi);
        Array4<Real> const& mv = normal[lev].array(mfi);
        Array4<Real> const& al = alpha[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          auto fvol = vof(i,j,k,0);
          THRESHOLD(fvol);
          if (!interface_cell (i,j,k, vof, fvol)) {
              AMREX_D_TERM(mv(i,j,k,0) = Real(0.);,
                           mv(i,j,k,1) = Real(0.);,
                           mv(i,j,k,2) = Real(0.););
              al(i,j,k) = fvol;
          }
          else {
             AMREX_D_PICK( ,Real f[3][3];, Real f[3][3][3];)
             XDim3 m;
             stencil (i,j,k, vof, f);
             mycs (f, &m.x);
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
       /*  Print() <<" normal direction "<< m.x<<" "<<m.y<<" "<<m.z<<"("<<i<<","<<j<<","<<k<<")"
                     <<"vof"<<"  "<<fvol<<"\n";*/
             //if (i==22&&j==20&&k==19)
        al(i,j,k)= plane_alpha (m, fvol);

        }
                  /*    if (i==12&&j==7&&k==7)
              Print() <<" normal direction "<<"("<<i<<","<<j<<","<<k<<") "
                     <<"vof"<<"  "<<fvol<<"  "<<"alpha " <<al(i,j,k)<<" "<<
                     interface_cell (i,j,k, vof, fvol)<<"\n";*/
        }); //  ParallelFor

    } // end MFIter

    //!!!!!!!!fix me: a temperary solution for the normal and alpha!!!!!!!!!!!
    // fill value of ghost cells (BCs, MPI infor)
    //vof_mf.FillBoundary();
    normal[lev].FillBoundary();
    alpha[lev].FillBoundary();
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


    // *************************************************************************************
    // Allocate space for the fluxes for vof advection
    // *************************************************************************************
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
                           dir < 2? v_mac[lev]:w_mac[lev];

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
        m_total_flux[lev].FillBoundary();
        vof_total_flux[lev].FillBoundary();


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
        tracer[lev]->FillBoundary();



        // update the normal and alpha of the plane in each interface cell
        // after each sweep
        tracer_vof_update(lev, *tracer[lev]);


     }// end i-,j-,k-sweep: calculation of vof advection


  }// end lev
  start = (start + 1) % AMREX_SPACEDIM;

    // determine the normal direction and alpha of the plane segment intersecting each interface cell.
   // tracer_vof_update(tracer);

}
////////////////////////////////////////////////////////////////////
///////
///////  Initialize the VOF value using the EB implicit surface function
///////
/////////////////////////////////////////////////////////////////////
void
VolumeOfFluid::tracer_vof_init_fraction(int lev, MultiFab& a_tracer)
{
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
            Array<Real,AMREX_SPACEDIM> center{AMREX_D_DECL((problo[0]+.35),
                                                           (problo[1]+.35),
                                                           (problo[2]+.35))};
            Real radius = .15; //5.0*dx[0];
            bool fluid_is_inside = true;

            EB2::SphereIF my_sphere(radius, center, fluid_is_inside);
          //  auto gshop = EB2::makeShop(my_sphere);



    // Initialise cylinder parameters
    int direction = 1;
    Real height = 14.5*dx[0];


    center[0]=0.5*(problo[0]+probhi[0]);
    center[1]=0.5*(problo[1]+probhi[1]);
    center[2]=0.75*(problo[1]+probhi[1]);
    // Build the Cylinder implficit function representing the curved walls
    EB2::CylinderIF my_cyl(radius, height, direction, center, false);
    radius = 8.0*dx[0];
    EB2::CylinderIF my_cyl_1(radius, height, direction, center, fluid_is_inside);

    //box
   /* Array<Real,AMREX_SPACEDIM> low{AMREX_D_DECL((problo[0]+10.3*dx[0]),
                                                (problo[1]+10.3*dx[1]),
                                                (problo[2]+10.3*dx[2]))};
    Array<Real,AMREX_SPACEDIM> high{AMREX_D_DECL((probhi[0]-11.2*dx[0]),
                                                 (probhi[1]-11.2*dx[1]),
                                                 (probhi[2]-11.2*dx[2]))};    */
   Array<Real,AMREX_SPACEDIM> low{AMREX_D_DECL( (problo[0]+.5/16.),
                                                (problo[1]+.5/16.),
                                                (problo[2]+.5/16.))};
    Array<Real,AMREX_SPACEDIM> high{AMREX_D_DECL((problo[0]+5.5/16.),
                                                 (problo[1]+5.5/16.),
                                                 (problo[2]+5.5/16.))};
    auto my_box= EB2::BoxIF( low,  high, fluid_is_inside);
    //auto my_box=  EB2::rotate(EB2::BoxIF( low,  high, fluid_is_inside), .3, 1);
    auto my_box1=  EB2::rotate(my_box, .3, 0);
    auto my_box2=  EB2::rotate(my_box1, .2, 2);
    //auto two =EB2::makeUnion(my_cyl_1, my_cyl);
   //auto two = EB2::makeComplement(EB2::makeUnion(my_cyl_1, my_cyl));

    // Generate GeometryShop
    //auto gshop = EB2::makeShop(two);
     auto gshop = EB2::makeShop(my_sphere);
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
    } else
#endif
    {
#ifdef AMRE_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a_tracer); mfi.isValid(); ++mfi)
        {
            Box const& vbx = mfi.validbox();
            auto const& tracer = a_tracer.array(mfi);
            amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
            });
        }
    }

    // Once vof tracer is initialized, we calculate the normal direction and alpha of the plane segment
    // intersecting each interface cell.
    v_incflo->p_volume_of_fluid->tracer_vof_update(lev, a_tracer);

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
                Vector<Segment> segments;
                int totalnodes = 0;
                for (int k = lo.z; k <= hi.z; ++k) {
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

                        /* Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "
                          " center "<< p <<" alpha "<<alpha<<" "<<fvol<<"\n";  */
                        //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" "<<fvol<<"\n";
                      //if (i==16 &&j==28&&k==3)
                        add_segment (center, dx, alpha, p, m, segments, totalnodes, fvol, tag_arr(i,j,k));
                    }
                }
                }
                }
                //           Print() << " Outside add_interface " << segments.size()<<" "<<totalnodes<<"\n";
                if (totalnodes > 0) {
                 // std::ofstream TecplotFile;
                //  TecplotFile.open(tecplotfilename, std::ios_base::trunc);
                   TecplotFile << "TITLE = \"incflow simulation from processer# " << myproc << "\" "<< "\n";
                   //spatial coordinates
                   TecplotFile << (AMREX_SPACEDIM== 2 ? "VARIABLES = \"X\", \"Y\"":"VARIABLES = \"X\", \"Y\", \"Z\"");
                   //output varibles
                   TecplotFile <<", \"F\""<<", \"m_x\""<<", \"m_y\""<<", \"m_z\""<<", \"alpha\""<<", \"tag\""<<"\n";
                   std::string zonetitle=("Level_"+std::to_string(lev)+
                                          "_Box_" +std::to_string(mfi.index())+
                                          "_Proc_"+std::to_string(myproc));
                   TecplotFile <<(std::string("ZONE T=")+zonetitle);
                   TecplotFile <<", DATAPACKING=POINT"<<", NODES="<<totalnodes<<", ELEMENTS="<<segments.size()
                               <<", ZONETYPE=FEQUADRILATERAL" <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//      AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";
                   int nn=0;
                   for (const auto& seg : segments) {
                       for (int in = 0; in < seg.nnodes; in++)
                           TecplotFile <<seg.node[in].x<<" "<<seg.node[in].y<<" "<<seg.node[in].z<<" "
                                       <<seg.vof<<" "<<seg.mv.x<<"  "<<seg.mv.y<<"  "<<seg.mv.z<<"  "
                                       <<seg.alpha<<"  "<<seg.tag<<"\n";

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

        TecplotFile << "TITLE = \"incflow simulation from processer# " << myproc << "\" "<< "\n";
        //spatial coordinates
        TecplotFile << (AMREX_SPACEDIM== 2 ? "VARIABLES = \"X\", \"Y\"":"VARIABLES = \"X\", \"Y\", \"Z\"");
        //output varibles
        TecplotFile <<", \"F\""<<", \"v_x\""<<", \"v_y\""<<", \"v_z\""<<
                  ", \"m_x\""<<", \"m_y\""<<", \"m_z\""<<", \"alpha\""<<", \"tag\""<<"\n";

        for (int lev = 0; lev <= finest_level; ++lev) {
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
                std::string zonetitle=("Level_"+std::to_string(lev)+"_Box_"+std::to_string(mfi.index()));
                TecplotFile <<(std::string("ZONE T=")+zonetitle);
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
                    TecplotFile <<", "<<(IJK[dim]+std::string("="))<<(ijk_max[dim]-ijk_min[dim]+2);
                TecplotFile <<", DATAPACKING=BLOCK"<<", VARLOCATION=(["<<(AMREX_SPACEDIM+1)<<"-"<<12<<"]=CELLCENTERED)"
                            <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//      AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";
                Array4<Real const> const& tracer = ld.tracer.const_array(mfi);
                Array4<Real const> const& vel = ld.velocity.const_array(mfi);
                Array4<Real const> const& mv = normal[lev].const_array(mfi);
                Array4<Real const> const& al = alpha[lev].const_array(mfi);
                Array4<Real const> const& tag_arr = tag[lev].const_array(mfi);
                int nn=0;
                //write coordinate variables
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    for (int k = lo.z; k <= hi.z +1; ++k) {
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
                    }
                }//

                for (int k = lo.z; k <= hi.z; ++k) {
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
                }
                //write velocity
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    for (int k = lo.z; k <= hi.z; ++k) {
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
                    }
                }//

                //write variables of the normal direction of the interface
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    for (int k = lo.z; k <= hi.z; ++k) {
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
                    }
                }//

                //write alpha of the interface

                    for (int k = lo.z; k <= hi.z; ++k) {
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
                    }

                //write alpha of the interface
                    for (int k = lo.z; k <= hi.z; ++k) {
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
                    }

                TecplotFile <<"\n";
            } // end MFIter

        } // end lev
    }


}

#include <queue>

#define CELL_IS_BOUNDARY(c,min,max) (c[0]< min[0]||c[1]< min[1]||c[2]< min[2]||c[0]>max[0]||c[1]> max[1]||c[2]> max[2])
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

int domain_tag_droplets (int finest_level, Vector<BoxArray > const &grids, Vector<Geometry > const& geom,
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
          std::queue <IntVect> fifo;
          tag_arr(i,j,k)=++ntag;
          fifo.push({i,j,k});
          while (!fifo.empty()){
           IntVect cell=fifo.front(),ncell=cell;
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
    tag[lev]->FillBoundary();
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
         for(int i0=ijk_min[ort1];i0<=ijk_max[ort1];i0++)
          for(int j0=ijk_min[ort2];j0<=ijk_max[ort2];j0++){
            Real tag_cell=(d==0?tag_arr(k0,i0,j0):
                           d==1?tag_arr(j0,k0,i0):
                                tag_arr(i0,j0,k0));
            if(tag_cell > 0){
              Real tag_gcell=(d==0?tag_arr(gd,i0,j0):
                              d==1?tag_arr(j0,gd,i0):
                                   tag_arr(i0,j0,gd));
              if(tag_gcell > 0){
                touching_regions (tag_cell, tag_gcell, touch);
              }
            }

          }// end for-loop for searching cells in the boundaries.
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
     tag[lev]->FillBoundary();
     ntag=maxtag;
    }

  }//end lev

  return ntag;

}

typedef struct {
  Real min, max, sum, sum2, mean, stddev;
  int n;
} VofRange;
//
// gts_range_init:
// @r: a #VofRange.
//
// Initializes a #VofRange.
//
void range_init (VofRange * r)
{
  r->max = - 10.e30;
  r->min = 10.e30;
  r->sum = r->sum2 = 0.0;
  r->n = 0;
}
/**
 * range_add_value:
 * @r: a #GtsRange.
 * @val: a value to add to @r.
 *
 * Adds @val to @r.
 */
void range_add_value (VofRange * r, Real val)
{
  if (val < r->min) r->min = val;
  if (val > r->max) r->max = val;
  r->sum += val;
  r->sum2 += val*val;
  r->n++;
}
static void range_reduce (void * i, void * o,
              int * len,
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

static void domain_range_reduce ( VofRange * s)
{

  double in[5];
  double out[5] = { 10.e30, - 10.e30, 0., 0., 0. };
  MPI_Op op;

  MPI_Op_create (range_reduce, true, &op);
  in[0] = s->min; in[1] = s->max; in[2] = s->sum; in[3] = s->sum2;
  in[4] = s->n;
  MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
  MPI_Op_free (&op);
  s->min = out[0]; s->max = out[1]; s->sum = out[2]; s->sum2 = out[3];
  s->n = out[4];

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
      // find the max and min location of interface */
      VofRange s[AMREX_SPACEDIM][n_tag];
      // the range of location of interfacial cells */
      Real vtop[n_tag], range[AMREX_SPACEDIM][2][n_tag];
      for (int n = 0; n < n_tag; n++){
       ncell[n]=0; vols[n] = 0.; vels[n] = 0.; surfA[n]=0.; vtop[n] = 0.;
       for(int d = 0; d < AMREX_SPACEDIM; d++) {
         mcent[d][n]=0.;
         range_init (&s[d][n]);
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
       //fix me: not compatable with GPUs
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
           for(int d = 0; d < AMREX_SPACEDIM; d++)
              range_add_value (&s[d][itag-1], (&p.x)[d]);
         }
        }
//       });
       }}}    //end of the ijk-loop

    }//end MFIter

 // the rest of the algorithm deals with  parallel BCs
    if (ParallelDescriptor::NProcs() > 1){
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
         domain_range_reduce(&s[d][n]);
      }
    }
////////////////////////////////////////////////////////////////////
//////
//////   Estimate the phase error for cube case
//////
////////////////////////////////////////////////////////////////////
Real error=0.;
if (false){
     Real lencube=5./16., cube_vol=lencube*lencube*lencube;
     Array<Real, AMREX_SPACEDIM> o0={0.1875,0.1875,0.1875},o,
                                 cube_min,cube_max;
     for (int d = 0; d < AMREX_SPACEDIM; d++){
        o[d] = o0[d]+1.0*time;
        cube_min[d]=o[d]-lencube*.5;
        cube_max[d]=o[d]+lencube*.5;
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

       //fix me: not compatable with GPUs
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
    if (time > .15){
   // Print()<<"test----"<<"\n";
    }
    //AllPrint()<<"error ----"<< error<<"\n";
    //ParallelDescriptor::ReduceRealSum (&error, 1);
    //Print()<<error<<"\n";
    error /=cube_vol;

}

if (true){
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

       //fix me: not compatable with GPUs
       ParallelFor(bx, [&] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
         //if (fv(i,j,k)>0.){

       error+=fabs(fv(i,j,k,1)-fv(i,j,k,0))*AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
        //}
       });

      }//end MFIter
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
        if (s[d][n].min<10.e30)
          range[d][0][n] = s[d][n].min;
        if (s[d][n].max>-10.e30)
          range[d][1][n] = s[d][n].max;
      }
     Real pos_limit[AMREX_SPACEDIM][2];
     for(int d = 0; d < AMREX_SPACEDIM; d++){
      pos_limit[d][0]=10.e30,pos_limit[d][1]=-10.e30;
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
                 <<"Time, Num of Drops, Volume, Speed, Centroid, Range, Surface area\n";
     }

// Write data to file
     outputFile << std::fixed << std::setprecision(8) << time << "   "<<ndrops<<"  "
                <<std::setprecision(8) << sumtotal<<"  "<<"error "<<error;
     for (int n = 0; n < n_tag && n < 7; n++){
       outputFile <<" #"<<n+1<<"  "<<vols[n]<<"  "<<vels[n]<<"  "<<"L  "
                  <<mcent[0][n]<<"  "<<mcent[1][n]<<"  "<<mcent[2][n]<<"  "<<"R ";
       for(int d = 0; d < AMREX_SPACEDIM; d++)
         outputFile <<range[d][0][n]<<"  "<<range[d][1][n]<<"  ";
       outputFile <<"s  "<<surfA[n]<<"  ";
     }
     outputFile <<"\n";
     outputFile.close();
      first =0;
    }
    } // end lev
// Close the file

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
          vel(i,j,k,0) = 2*sin(2.*pi*y)*sin(pi*x)*sin(pi*x)*sin(2*pi*z)*cos(pi*time/3.);
          vel(i,j,k,1) = -sin(2.*pi*x)*sin(pi*y)*sin(pi*y)*sin(2*pi*z)*cos(pi*time/3.);
#if (AMREX_SPACEDIM == 3)
          vel(i,j,k,2) = -sin(2.*pi*x)*sin(pi*z)*sin(pi*z)*sin(2*pi*y)*cos(pi*time/3.);
#endif
        });
     }//end MFIter
    } //end lev
}

