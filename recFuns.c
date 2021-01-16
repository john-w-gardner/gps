#include <stdio.h>
#include <math.h>
#include "receiver.h"

int computeVehLocation(struct satInstance satarr[], struct vehInstance *veh, int nsat)
{
  int done = 0;
  double x[3] = {veh.x, veh.y, veh.z};
  int ndim = 3;
  double F[ndim];
  double J[ndim*ndim];

  while (done == 0)
    {
      computeF(satarr, nsat, *veh, &F);
      computeJ(satarr, *veh, &J);
      solve(J, 3, s, F);
    }
}

/* 
   compute gradient of least squares function as in (66)
   break this down into partial derivatives and  N_i, A_i
 */
void computeF(struct satInstance satarr[], int nsat, struct vehInstance veh, double (*F)[])
{
  // reset
  (*F)[0] = 0;
  (*F)[1] = 0;
  (*F)[2] = 0;
  double Ai, Xi, Yi, Zi;
  int i;
  
  while (i < nsat-2) // n-1 equations
    {
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Xi = computeXi(satarr[i], satarr[i+1], veh);
      Yi = computeYi(satarr[i], satarr[i+1], veh);
      Zi = computeZi(satarr[i], satarr[i+1], veh);
      (*F)[0] += 2*Ai*Xi;
      (*F)[1] += 2*Ai*Yi;
      (*F)[2] += 2*Ai*Zi;
    }
}


double computeAi(struct satInstance si, struct satInstance si1, struct vehInstance veh)
{
  double Ni = computeNi(si, veh);
  double Ni1 = computeNi(si1, veh);
  
  double result = Ni1 - Ni - c*(si.t - si1.t);
  return result;              
}

double computeNi(struct satInstance s, struct vehInstance v)
{
  // matops lib needs vectors
  double xs[] = {s.x, s.y, s.z};
  double xv[] = {v.x, v.y, v.z};

  double result = eucmetric(xs, xv, 3);
  return result;
}

double computeXi(struct satInstance si, struct satInstance si1, struct vehInstance v))
{
  double result;
  double Ni = computeNi(si, v);
  double Ni1 = computeNi(si1, v);

  result = -(si1.x - v.x)/Ni1 + (si.x - v.x)/Ni;
} 

double computeYi(struct satInstance si, struct satInstance si1, struct vehInstance v))
{
  double result;
  double Ni = computeNi(si, v);
  double Ni1 = computeNi(si1, v);

  result = -(si1.y - v.y)/Ni1 + (si.y - v.y)/Ni;
} 

double computeZi(struct satInstance si, struct satInstance si1, struct vehInstance v))
{
  double result;
  double Ni = computeNi(si, v);
  double Ni1 = computeNi(si1, v);

  result = -(si1.z - v.z)/Ni1 + (si.z - v.z)/Ni;
} 

// end computeF support functions


// compute Jacobian (matrix of second derivatives of f) from eqns (70), (71)
// J must have size 9..
void computeJ(double (*J)[], struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  (*J)[0] = d2fx2(satarr[], nsat, veh);
}

/* helpers:
   "fxx" = f_{xx} = \frac{\partial^2 f}{\partial x^2} ie second derivative of f wrt x
   "fxy" = f_{xy} = f_{yx} = \frac{\partial^2 f}{\partial x \partial y} ie derivative of (derivative of f wrt y) wrt x
   etc
   careful: these are all very similar but with important differences
*/

double fxx(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  double result = 0;
  double Xi, Ai, pXi
  int i;


  while (i < nsat-2)
    {
      Xi = computeXi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Xx = Xx(satarr[i], satarr[i+1], veh);
      result += 2*(pow(Xi,2) + Ai*Xx);
    }  
  return result;
}

// derivative of X_i wrt x
double Xx(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = (pow(Ni1, 2) - pow((si1.x - v.x), 2))/pow(Ni1, 3);
  b = (pow(Ni, 2) - pow((si.x - v.x), 2))/pow(Ni, 3);
  
  result = a-b;
  return result;
}

double fxy(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  int i;
  double result = 0;

  while (i < nsat-2)
    {
      Xi = computeXi(satarr[i], satarr[i+1], veh);
      Yi = computeYi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Xy = computeXy(satarr[i], satarr[i+1], veh);
      result += 2*(Xi*Yi + Ai*Xy);
    }
  
  return result;
}

double Xy(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.y - v.y)*(si1.x - v.x)/pow(Ni1,3);
  b = (si.y - v.y)*(si.x - v.x)/pow(Ni,3);

  result = a+b;
  return result;
}

double fxz(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  double Xi, Zi, Ai, Xz;
  int i;
  double result = 0;

  while (i < nsat-2)
    {
      Xi = computeXi(satarr[i], satarr[i+1], veh);
      Yi = computeZi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Xz = Xz(satarr[i], satarr[i+1], veh);
      result += 2*(Xi*Zi + Ai*Xz);
    }
  
  return result;
}

double Xz(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.z - v.z)*(si1.x - v.x)/pow(Ni1,3);
  b = (si.z - v.z)*(si.x - v.x)/pow(Ni,3);

  result = a+b;
  return result;
}

double fyy(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  double result = 0;
  double Yi, Ai, Yy
  int i;


  while (i < nsat-2)
    {
      Yi = computeYi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Yy = Yy(satarr[i], satarr[i+1], veh);
      result += 2*(pow(Yi,2) + Ai*Yy);
    }  
  return result;
}

double Yy(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = (pow(Ni1, 2) - pow((si1.y - v.y), 2))/pow(Ni1, 3);
  b = (pow(Ni, 2) - pow((si.y - v.y), 2))/pow(Ni, 3);
  
  result = a-b;
  return result;
}

double fyz(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  int i;
  double result = 0;

  while (i < nsat-2)
    {
      Zi = computeZi(satarr[i], satarr[i+1], veh);
      Yi = computeYi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Yz = Yz(satarr[i], satarr[i+1], veh);
      result += 2*(Zi*Yi + Ai*Zy);
    }
  
  return result;
}

double Yz(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.z - v.z)*(si1.y - v.y)/pow(Ni1,3);
  b = (si.z - v.z)*(si.y - v.y)/pow(Ni,3);

  result = a+b;
  return result;
}

double fzz(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  double result = 0;
  double Zi, Ai, Zz
  int i;


  while (i < nsat-2)
    {
      Zi = computeZi(satarr[i], satarr[i+1], veh);
      Ai = computeAi(satarr[i], satarr[i+1], veh);
      Zz = Zz(satarr[i], satarr[i+1], veh);
      result += 2*(pow(Zi,2) + Ai*Zz);
    }  
  return result;
}

double Zz(struct satInstance si, struct satInstance si1, struct vehInstance v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = (pow(Ni1, 2) - pow((si1.z - v.z), 2))/pow(Ni1, 3);
  b = (pow(Ni, 2) - pow((si.z - v.z), 2))/pow(Ni, 3);
  
  result = a-b;
  return result;
}

// end computeJ support functions



void readDataDatII(double *pi, double *c, double *R, double *s) 
{
  // note: strong dependency on organization of data.dat 
  FILE *dd;
  dd = fopen("data.dat", "r");

  // grab constants
  int NCHAR = 1000;
  char buf[NCHAR];
  int status;

  fgets(buf, NCHAR, dd);
  status = sscanf(buf, "%lf", pi); // ignore text
  fgets(buf, NCHAR, dd);
  status = sscanf(buf, "%lf", c); 
  fgets(buf, NCHAR, dd);
  status = sscanf(buf, "%lf", R); 
  fgets(buf, NCHAR, dd);
  status = sscanf(buf, "%lf", s); 

  fclose(dd);
}


