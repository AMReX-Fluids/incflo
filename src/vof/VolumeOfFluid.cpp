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

VolumeOfFluid::VolumeOfFluid (incflo* a_incflo) : v_incflo(a_incflo)
{
    finest_level = v_incflo->finestLevel();
    // *************************************************************************************
    // Allocate space for the information of the interface segments
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; ++lev){
        normal.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], AMREX_SPACEDIM, 1, MFInfo(), v_incflo->Factory(lev));
         alpha.emplace_back(v_incflo->grids[lev], v_incflo->dmap[lev], 1, 1, MFInfo(), v_incflo->Factory(lev));
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

static void cube_plane_intersection (XDim3 cell, GpuArray <Real,AMREX_SPACEDIM> dx,
				     XDim3 const * O, XDim3 const * n, XDim3 p[12],
				     int orient[12])
{
  XDim3 o;
  int i;
  for (i=0; i<AMREX_SPACEDIM; i++)
    (&o.x)[i] = (&cell.x)[i]-Real(0.5)*dx[i];
  for (i = 0; i < 12; i++) {
    XDim3 e, d;
	Real h= i<4? dx[0]:(i<8?dx[1]:dx[2]);
    d.x = o.x + h*edge[i][0].x; d.y = o.y + h*edge[i][0].y; d.z = o.z + h*edge[i][0].z;
    e.x = o.x + h*edge[i][1].x; e.y = o.y + h*edge[i][1].y; e.z = o.z + h*edge[i][1].z;
    Real den = n->x*(e.x - d.x) + n->y*(e.y - d.y) + n->z*(e.z - d.z);
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      Real t = (n->x*(O->x - d.x) + n->y*(O->y - d.y) + n->z*(O->z - d.z))/den;
      if (t >= 0. && t < 1.) {
	p[i].x = d.x + t*(e.x - d.x); p[i].y = d.y + t*(e.y - d.y); p[i].z = d.z + t*(e.z - d.z);
	orient[i] = (n->x*(e.x - O->x) + n->y*(e.y - O->y) + n->z*(e.z - O->z) > 0.);
      }
    }
  }
}
using NODE_CUT=Array<XDim3, AMREX_SPACEDIM*(AMREX_SPACEDIM - 1) + 1>;
/**
 * vof_cut_cube_vertices:
 * @cell: a #FttCell.
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
static int vof_cut_cube_vertices (XDim3 cell, GpuArray <Real,AMREX_SPACEDIM> dx,
                           XDim3 const * p,  XDim3 const *  n,
			   NODE_CUT &  v, int d[12])
{
  XDim3 a[12];
  int orient[12];
  int i;

  AMREX_ASSERT (p != NULL);



  cube_plane_intersection (cell, dx, p, n, a, orient);
  for (i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = connect[e][orient[e]];
      d[nv] = ne[3];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
	e = ne[m++];
    }
    if (nv > 2)
      return nv;
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
  Real alpha, vof;
 // Constructor to initialize the Segment
  Segment(int n, NODE_CUT const& nodes, XDim3 m, Real a, Real f, int ns=0)
   : nnodes(n), mv (m), alpha (a), vof(f) {
   for (int i = 0; i < n; ++i)
      node[i]= nodes[ns==0?i:(i + 3)%(n + 2)];
  }
};

static void add_segment (XDim3 const & cell, GpuArray <Real,AMREX_SPACEDIM> const & dx,
                         Real alpha, XDim3 const * o, XDim3 const * m,
			Vector<Segment> & segments, int & nt, Real vof)
{

   /*Print() <<" add_segment "<< *o<<"  "<<" vector "<<*m
		             <<"vof"<<"  "<<vof<<"  "<<"alpha " <<alpha<<"\n";  */
  
  int d[12];
  /* array of node coordinates for a cut face */
  NODE_CUT  nodecutface;
  int inode, inode2, jnode_max_sintheta = 0,
  nnodecutface = vof_cut_cube_vertices (cell, dx, o, m, nodecutface, d);
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
	     XDim3 diff1 = {o->x - node.x, o->y - node.y, o->z - node.z};
	     Real length_diff1 = vector_norm (&diff1);
	     Real max_sintheta = 0.;
	 /*  cycle through all other nodes (jnode) where cut face intersects cell edges */
	     for (inode2 = 1; inode2 < nnodecutface; inode2++) {
	       int jnode = (inode + inode2)%nnodecutface;
	       XDim3 diff2 = {nodecutface[jnode].x - node.x,
	 	 	              nodecutface[jnode].y - node.y,
	 	 	              nodecutface[jnode].z - node.z};
	       Real length_diff2 = vector_norm (&diff2);
	       if (length_diff2 < 1e-20) /*Hua Tan(11-1-2016)*/
	          return;
	       Real sintheta = ((diff1.y*diff2.z - diff1.z*diff2.y)*m->x +
	 	 	                (diff1.z*diff2.x - diff1.x*diff2.z)*m->y +
	 	 	                (diff1.x*diff2.y - diff1.y*diff2.x)*m->z)/
	                        (length_diff1*length_diff2);

	       if (sintheta > max_sintheta) {
	         max_sintheta = sintheta;
	         jnode_max_sintheta = jnode;
	       }
	      }
	     /* terminate if cannot find positive angle between cut face nodes */
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
   segments.emplace_back(nnodecutface, nodecutface, *m, alpha, vof);
   //Print() << " normal direction " <<nnodecutface<<" "<<*m<<""<<vof<<"\n";
 }
 else { /* cut face must be divided into 2 quadrilateral/triangular faces */
   /* first face is quadrilateral */
   nt += 4;
   segments.emplace_back(4, nodecutface, *m, alpha, vof);
   /* second face is triangular if nnodecutface=5; otherwise quadrilateral */
//   WallFace * face2 = g_malloc (sizeof (WallFace));

   nnodecutface -= 2;

   nt += nnodecutface;
   segments.emplace_back(nnodecutface, nodecutface, *m, alpha, vof, 3);
 } /* cut face must be divided into 2 quadrilateral/triangular faces */

}

/**
 * vof_line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #amrex::XDim3.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
Real vof_line_area_center (XDim3 const * m, Real alpha, XDim3 * p)
{
  XDim3 n;


  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  //made change from "alpha>=n.x+n.y" to consider 
  //the extreme cases
  if (alpha <= 0. || alpha > n.x + n.y) {
    p->x = p->y = 0.;
    return 0.;
  }

  if (n.x < EPS) {
    p->x = 0.5;
    p->y = m->y < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  if (n.y < EPS) {
    p->y = 0.5;
    p->x = m->x < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  p->x = p->y = 0.;

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  Real ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  p->x /= 2.;
  p->y /= 2.;

  THRESHOLD (p->x);
  THRESHOLD (p->y);

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;

  return sqrt (ax*ax + ay*ay);

}


/**
 * vof_plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #amrex::XDim3.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 *
 * Returns: the area of the facet.
 */
Real vof_plane_area_center (XDim3 const * m, Real alpha, XDim3 * p)
{

  if (fabs (m->x) < EPS) {
    XDim3 n, q;
    n.x = m->y;
    n.y = m->z;
    Real area = vof_line_area_center (&n, alpha, &q);
    p->x = 0.5;
    p->y = q.x;
    p->z = q.y;
    return area;
  }
  if (fabs (m->y) < EPS) {
    XDim3 n, q;
    n.x = m->z;
    n.y = m->x;
    Real area = vof_line_area_center (&n, alpha, &q);
    p->x = q.y;
    p->y = 0.5;
    p->z = q.x;
    return area;
  }
  if (fabs (m->z) < EPS) {
    Real area = vof_line_area_center (m, alpha, p);
    p->z = 0.5;
    return area;
  }

  XDim3 n = *m;
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
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  Real area = alpha*alpha;
  p->x = p->y = p->z = area*alpha;

  Real b = alpha - n.x;
  if (b > 0.) {
    area -= b*b;
    p->x -= b*b*(2.*n.x + alpha);
    p->y -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    area -= b*b;
    p->y -= b*b*(2.*n.y + alpha);
    p->x -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    area -= b*b;
    p->z -= b*b*(2.*n.z + alpha);
    p->x -= b*b*b;
    p->y -= b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    area += b*b;
    p->y += b*b*(2.*n.y + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.y);
    p->x += b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.x);
    p->y += b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.y);
    p->y += b*b*(2.*n.y + alpha - n.x);
    p->z += b*b*b;
  }

  area *= 3.;
  p->x /= area*n.x;
  p->y /= area*n.y;
  p->z /= area*n.z;

  THRESHOLD (p->x);
  THRESHOLD (p->y);
  THRESHOLD (p->z);

  if (m->x < 0.) p->x = 1. - p->x;
  if (m->y < 0.) p->y = 1. - p->y;
  if (m->z < 0.) p->z = 1. - p->z;

  return area*sqrt (1./(n.x*n.x*n.y*n.y) + 1./(n.x*n.x*n.z*n.z) + 1./(n.z*n.z*n.y*n.y))/6.;
}


/**
 * vof_plane_alpha:
 *
 * Returns: the value @alpha such that the volume of a cubic cell
 * lying under the plane defined by @m.@x = @alpha is equal to @c.
 */
Real vof_plane_alpha (XDim3 * m, Real c)
{

  AMREX_ASSERT(c >= 0. && c <= 1.);
  AMREX_ASSERT(m != NULL);
  // made change to avoid the numerical issue for
  // full cell that has empty cell in the face neighbor.
  if (c == 1. )
   c =1.-EPS;     
  
  Real alpha;
  XDim3 n;
//  m->x =1., m->y=0., m->z=0., c=1.;
//  Print()<<"vector"<<m->x<<" "<<m->y<<" "<<m->z<<"  "<<"vof "<<c<<"\n";
  n.x = fabs (m->x); n.y = fabs (m->y); n.z = fabs (m->z);

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
//	Print()<<"p q p12  "<<p<<" "<<q<<" "<<p12<<"  "<<"vof"<<c<<"\n";
    Real teta = acos(q/(p*p12))/3.;
    Real cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;
  if (m->z < 0.)
    alpha += m->z; 
//  Print()<<"alpha----  "<<alpha<<"\n";
  return alpha;
}

#if AMREX_SPACEDIM == 2
#include "myc2D.h"
#else
#include "myc.h"
#endif

void stencil (AMREX_D_DECL(int const i, int const j, int const k),
              Array4<Real const> const & v, 
              AMREX_D_PICK(      ,
                    Real fv[3][3],
                    Real fv[3][3][3]))                           
{
  int x, y, z = 0;
  AMREX_D_PICK(0,fv[1][1],fv[1][1][1]) = v (AMREX_D_DECL(i,j,k));
#if AMREX_SPACEDIM == 3
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
	if (x != 0 || y != 0 || z != 0)
	  AMREX_D_PICK(,fv[x + 1][y+1],fv[x + 1][y+1][z+1])=v(AMREX_D_DECL(i+x,j+y,k+z));

	  
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

bool interface_cell (AMREX_D_DECL(int const i, int const j, int const k),
                     Array4<Real const> const & v, Real fc)
{
  if (fc == 1.){
//when a full cell has an empty cell in its face neighbor (i.e., +-x,+-y,+-z)
// we also need to
// calculate the normal and alpha of the cut plane.  
  for (int dim=0; dim<AMREX_SPACEDIM; dim++)
    for (int step = -1; step <= 1; step+=2)
	if (dim==0?v(AMREX_D_DECL(i+step,j,k)) == 0.:
	    dim==1?v(AMREX_D_DECL(i,j+step,k)) == 0.:
	           v(AMREX_D_DECL(i,j,k+step)) == 0.){          
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
VolumeOfFluid::tracer_vof_update(Vector<MultiFab*> const& tracer)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& dx = v_incflo->geom[lev].CellSizeArray();
        auto const& problo = v_incflo->geom[lev].ProbLoArray();
        auto const& probhi = v_incflo->geom[lev].ProbHiArray();

        auto& vof_mf = tracer[lev];

        for (MFIter mfi(*vof_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vof = vof_mf->array(mfi);
            Array4<Real> const& mv = normal[lev].array(mfi);
            Array4<Real> const& al = alpha[lev].array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              auto fvol = vof(i,j,k,0);
              THRESHOLD(fvol);
              if (!interface_cell (AMREX_D_DECL(i,j,k), vof, fvol)) {
                 AMREX_D_TERM(mv(i,j,k,0) = Real(0.);,
                              mv(i,j,k,1) = Real(0.);,
                              mv(i,j,k,2) = Real(0.););
                 al(i,j,k) = fvol;
              }
              else {
                AMREX_D_PICK( ,Real f[3][3];, Real f[3][3][3];)
		XDim3 m;
                stencil (AMREX_D_DECL(i,j,k), vof, f);
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
		         //	Print() <<" normal direction "<< m.x<<" "<<m.y<<" "<<m.z<<"("<<i<<","<<j<<","<<k<<")"
		         //    <<"vof"<<"  "<<fvol<<"\n";		        
                //if (i==3&&j==6&&k==4)
		al(i,j,k)= vof_plane_alpha (&m, fvol);
		         
              }
              	/*	if (i==12&&j==7&&k==7)
		      Print() <<" normal direction "<<"("<<i<<","<<j<<","<<k<<") "
		             <<"vof"<<"  "<<fvol<<"  "<<"alpha " <<al(i,j,k)<<" "<<
		             interface_cell (AMREX_D_DECL(i,j,k), vof, fvol)<<"\n";*/
            }); //  ParallelFor


        } // end MFIter





    }// end lev
}



void
VolumeOfFluid::tracer_vof_advection(Vector<MultiFab*> const& tracer,
                                    AMREX_D_DECL(Vector<MultiFab const*> const& u_mac,
                                                 Vector<MultiFab const*> const& v_mac,
                                                 Vector<MultiFab const*> const& w_mac),
				      Real dt)
{

  amrex::Print() << " VOF Level#" << finest_level<<"\n";
  tracer_vof_update(tracer);
}
////////////////////////////////////////////////////////////////////
//////  Initialize the VOF value using the implicit surface function
/////////////////////////////////////////////////////////////////////
void
tracer_vof_init_fraction(int lev, MultiFab& a_tracer, incflo const* a_incflo)
{
    int vof_init_with_eb = 1;
    ParmParse pp("incflo");
    pp.query("vof_init_with_eb", vof_init_with_eb);

    Geometry const& geom = a_incflo->Geom(lev);
    auto const& dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
    auto const& probhi = geom.ProbHiArray();

#ifdef AMREX_USE_EB
    if (vof_init_with_eb) {
        if (lev == 0) {
            Array<Real,AMREX_SPACEDIM> center{AMREX_D_DECL(0.5*(problo[0]+probhi[0]),
                                                           0.5*(problo[1]+probhi[1]),
                                                           0.5*(problo[2]+probhi[2]))};
            Real radius = 5.0*dx[0];
            bool fluid_is_inside = true;

          //  EB2::SphereIF my_sphere(radius, center, fluid_is_inside);
          //  auto gshop = EB2::makeShop(my_sphere);

			

    // Initialise cylinder parameters
    int direction = 0;
    Real rotation  = 0, height = 11.*dx[0];
    int rotation_axe  = 0;
    rotation = (rotation/180.)*M_PI;


    // Build the Cylinder implficit function representing the curved walls
    EB2::CylinderIF my_cyl(radius, height, direction, center, fluid_is_inside);

   // auto my_cyl_rot = EB2::rotate(my_cyl, rotation, rotation_axe);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_cyl);			
			
            int max_level = a_incflo->maxLevel();
            EB2::Build(gshop, a_incflo->Geom(max_level), max_level, max_level);			
			
			
			
        }

        auto fact = amrex::makeEBFabFactory(geom, a_tracer.boxArray(), a_tracer.DistributionMap(),
                                            {1,1,0}, EBSupport::volume);
        auto const& volfrac = fact->getVolFrac();
        MultiFab::Copy(a_tracer, volfrac, 0, 0, 1, 1);

        if (lev == a_incflo->finestLevel()) {
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
    amrex::AllPrint() << " Output surface file at process#" << myproc<<"  " << nprocs << " at time " << time << std::endl;
    const std::string& tecplotfilename = amrex::Concatenate("tecplot_surface_", nstep)+"_";

    const int nfiles = 1;

    for (NFilesIter nfi(nfiles, tecplotfilename, false, true); nfi.ReadyToWrite(); ++nfi)
    {

//    Print()<<"surface_plot"<< nstep<<"  "<< tecplotfilename<<"\n";
        auto& TecplotFile = (std::ofstream&) nfi.Stream();

        TecplotFile << "TITLE = \"incflow simulation from processer# " << myproc << "\" "<< "\n";
        //spatial coordinates
        TecplotFile << (AMREX_SPACEDIM== 2 ? "VARIABLES = \"X\", \"Y\"":"VARIABLES = \"X\", \"Y\", \"Z\"");
        //output varibles
        TecplotFile <<", \"F\""<<", \"m_x\""<<", \"m_y\""<<", \"m_z\""<<", \"alpha\""<<"\n";

        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *v_incflo->m_leveldata[lev];
            Box const& domain  = v_incflo->geom[lev].Domain();
            auto const& dx     = v_incflo->geom[lev].CellSizeArray();
            auto const& problo = v_incflo->geom[lev].ProbLoArray();
            auto const& probhi = v_incflo->geom[lev].ProbHiArray();
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
                Array4<Real> const& mv = normal[lev].array(mfi);
                Array4<Real> const& al =  alpha[lev].array(mfi);
                Vector<Segment> segments;
                int totalnodes = 0;
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    auto fvol = vof(i,j,k,0);

                    if (interface_cell (AMREX_D_DECL(i,j,k), vof, fvol)){
                        Real   alpha;
                        XDim3	m, p, cell;
                        for (int d = 0; d < AMREX_SPACEDIM; d++) {
                            (&m.x)[d]= mv(i,j,k,d);
                        }
                        alpha= al(i,j,k,0);
                        //Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" vector "<<m<<" "<<" alpha "<<alpha<<" "<<fvol<<"\n";
                        vof_plane_area_center(&m, alpha, &p);
                        /* convert the coord. of the center to the global sys*/
                        for (int dim = 0; dim < AMREX_SPACEDIM; dim++){
                            (&p.x)[dim] = problo[dim] + dx[dim]*((dim<1?i:dim<2?j:k)+(&p.x)[dim]);
                            (&cell.x)[dim] = problo[dim] + dx[dim]*((dim<1?i:dim<2?j:k)+Real(0.5));
                        }
                       // Print() << " ijk index " <<"("<<i<<","<<j<<","<<k<<")"<<" "<<fvol<<"\n";
                      //if (i==12 &&j==7&&k==7)
                        add_segment (cell, dx, alpha, &p, &m, segments, totalnodes, fvol);
                    }
                }
                }
                }
                //           Print() << " Outside add_interface " << segments.size()<<" "<<totalnodes<<"\n";

                std::string zonetitle=("Level_"+std::to_string(lev)+
                                       "_Box_" +std::to_string(mfi.index())+
                                       "_Proc_"+std::to_string(myproc));
                TecplotFile <<(std::string("ZONE T=")+zonetitle);
                TecplotFile <<", DATAPACKING=POINT"<<", NODES="<<totalnodes<<", ELEMENTS="<<segments.size()
                            <<", ZONETYPE=FEQUADRILATERAL" <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//	  AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";
                int nn=0;
                for (const auto& seg : segments) {
                    for (int in = 0; in < seg.nnodes; in++)
                        TecplotFile <<seg.node[in].x<<" "<<seg.node[in].y<<" "<<seg.node[in].z<<" "
                                    <<seg.vof<<" "<<seg.mv.x<<"  "<<seg.mv.y<<"  "<<seg.mv.z<<"  "<<seg.alpha<<"\n";

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
            } // end MFIter
        } // end lev
    } // NFiles
}

void VolumeOfFluid::WriteTecPlotFile(Real time, int nstep)
{
    BL_PROFILE("incflo::WriteTecPlotFile()");
    std::string m_tecplot_file{"tecplot"};
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
        TecplotFile <<", \"F\""<<", \"m_x\""<<", \"m_y\""<<", \"m_z\""<<"\n";

        for (int lev = 0; lev <= finest_level; ++lev) {
            auto& ld = *v_incflo->m_leveldata[lev];
            Box const& domain  = v_incflo->geom[lev].Domain();
            auto const& dx     = v_incflo->geom[lev].CellSizeArray();
            auto const& problo = v_incflo->geom[lev].ProbLoArray();
            auto const& probhi = v_incflo->geom[lev].ProbHiArray();
            auto const& ijk_min= domain.smallEnd();
            auto const& ijk_max= domain.bigEnd();
            const BoxArray& ba = ld.tracer.boxArray();  //cell-centered multifab
            int nb = ba.size();
            const DistributionMapping& dm = ld.tracer.DistributionMap();
            std::string IJK = "IJK";
 //       amrex::Print() << " process#" << myproc<<"  " << ld.tracer.nGrow()<<" " << nb<<"\n";
//amrex::Print() << " process#" << myproc<<"  " << (IJK[0]+std::string("= "))<<"\n";
//     Print() << " process#" << myproc<<"  " << problo[0]<<" dx "<< dx[0]<<"   -------"<<"\n";
    //Output data for each box in boxarray according to Tecplot data format
//       for (int ibox = 0; ibox <nb; ++ibox) {
//		// only works on the boxes that are active on current processor
//		 if (myproc == dm[ibox]) {

//          Box const& bx = ba.get(ibox);
//          auto const& ijk_min= bx.smallEnd();
//          auto const& ijk_max= bx.bigEnd();
//          amrex::AllPrint() << " process#" << myproc<<"  " << ibox<<ijk_min<<ijk_max<<"\n";
//          //define Zone header for tecplot
//		  std::string zonetitle=("Level_"+std::to_string(lev)+"_Box_"+std::to_string(ibox));
//		  // Zone title
//		  TecplotFile <<(std::string("ZONE T=")+zonetitle);
//		  for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
//			 //note: tracer is cell-centered multifab, so the number of points in a dimension is one
//		     // more than the number of cells in that dimension
//             TecplotFile <<", "<<(IJK[dim]+std::string("="))<<(ijk_max[dim]-ijk_min[dim]+2);
//		  TecplotFile <<", DATAPACKING=BLOCK"<<", VARLOCATION=(["<<(AMREX_SPACEDIM+1)<<"]=CELLCENTERED)"
//		              <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";
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
                TecplotFile <<", DATAPACKING=BLOCK"<<", VARLOCATION=(["<<(AMREX_SPACEDIM+1)<<"-"<<7<<"]=CELLCENTERED)"
                            <<", SOLUTIONTIME="<<std::to_string(time)<<"\n";

//	  AllPrint() << " process#" << myproc<<"  " << lo << hi<<mfi.index()<<"\n";
                Array4<Real const> const& tracer = ld.tracer.const_array(mfi);
                Array4<Real const> const& mv = normal[lev].const_array(mfi);

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


                TecplotFile <<"\n";
            } // end MFIter

        } // end lev
    }
	std::rename("output.txt", "x.out");
}
