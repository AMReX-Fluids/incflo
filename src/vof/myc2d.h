#define NOT_ZERO 1.e-30

/*-----------------------------------------------------*
 *MYC - Mixed Youngs and Central Scheme (2D)           *
 *-----------------------------------------------------*/
static void mycs(double c[3][3],double mxy[2])
{
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;
  
  /* top, bottom, right and left sums of c values */
  c_t = c[0][2] + c[1][2] + c[2][2];
  c_b = c[0][0] + c[1][0] + c[2][0];
  c_r = c[2][0] + c[2][1] + c[2][2];
  c_l = c[0][0] + c[0][1] + c[0][2];

  /* consider two lines: sgn(my) Y =  mx0 X + alpha,
     and: sgn(mx) X =  my0 Y + alpha */ 
  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);

  /* minimum coefficient between mx0 and my0 wins */
  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }

  /* Youngs' normal to the interface */
  mm1 = c[0][0] + 2.0*c[0][1] + c[0][2];
  mm2 = c[2][0] + 2.0*c[2][1] + c[2][2];
  mx1 = mm1 - mm2;
  mm1 = c[0][0] + 2.0*c[1][0] + c[2][0];
  mm2 = c[0][2] + 2.0*c[1][2] + c[2][2];
  my1 = mm1 - mm2;

  /* choose between the best central and Youngs' scheme */ 
  if (ix) {
    mm1 = fabs(my1) + NOT_ZERO; 
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1) + NOT_ZERO; 
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
	
  /* normalize the set (mx0,my0): |mx0|+|my0|=1 and
     write the two components of the normal vector  */
  mm1 = fabs(mx0) + fabs(my0) + NOT_ZERO; 
  mxy[0] = mx0/mm1;
  mxy[1] = my0/mm1;
  
  return;
}
