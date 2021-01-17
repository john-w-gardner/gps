#include <stdio.h>
#include <math.h>
#include "receiver.h"

/* newton iteration for nonlinear system, see (58)
   f is least squares function we want to minimize
   F is its gradient, J the jacobian of F
   until satisfied, solve Js = -F, set x = x + s
 */

int computeVehLocation(struct satInstance satarr[], struct vehInstance *veh, int nsat)
{
  int done = 0;
  double x[3] = {veh.x, veh.y, veh.z};
  double s[3];
  double tol = 0.01;
  int ndim = 3;
  double F[ndim];
  double J[ndim*ndim];

  while (stepsize < tol)
    {
      computeF(satarr, nsat, *veh, &F);
      computeJ(satarr, *veh, &J);
      solveAxb(J, 3, s, F);
      stepsize = eucnorm(s, 3);
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


/* compute Jacobian (matrix of second derivatives of f) from eqns (70), (71)
   J must have size 9.. 
   Reading into row-major form matrix, though irrelevant since since symmetric
*/
void computeJ(double (*J)[], struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  (*J)[0] = fxx(satarr[], nsat, veh);
  (*J)[1] = (*J)[3] = fxy(satarr[], nsat, veh);
  (*J)[2] = (*J)[6] = fxz(satarr[], nsat, veh);
  (*J)[4] = fyy(satarr[], nsat, veh);
  (*J)[5] = (*J)[7] = fyz(satarr[], nsat, veh);
  (*J)[8] = fzz(satarr[], nsat, veh);
}

/* 
   helpers:
   "fxx" = f_{xx} = \frac{\partial^2 f}{\partial x^2} ie second derivative of f wrt x
   "fxy" = f_{xy} = f_{yx} = \frac{\partial^2 f}{\partial x \partial y} ie derivative of (derivative of f wrt y) wrt x
   etc
   careful: these are all very similar but with important differences
*/

double fxx(struct satInstance satarr[], struct vehInstance veh, int nsat)
{
  double result = 0;
  double Xi, Ai, Xx
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
  double Xi, Yi, Ai, Xy;
  double result = 0;
  int i;
  
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
  int i;
  double Xi, Zi, Ai, Xz;
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
  double Zi, Yi, Ai, Yz
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


// compute t_s given vehicle location
// use all satellites and eqn (42)
double computeVehTime(struct satInstance satarr[], int n, struct vehInstance veh)
{
  int i;
  double result = 0;

  for (i=0; i<n; i++)
    {
      Ni = computeNi(satarr[i], veh);
      result += Ni + c*satarr[i].t;
    }

  result = result/(n*c);
  return result;
}


// going backwards: convert from R^3 to lat/long decomposed as in vehicle output (10)
// (testing purposes - only receiver actually uses this)
void convertCoords(struct vehInstance vehCart, struct latlongVeh *vehLL) 
{
  //printf("--constants pi: %lf s: %lf R: %lf\n", pi, s, R);
  //extern double pi, s, R;
  int EW, NS;
  int lmm, lmd;
  double lm;
  double lms;
  double h;
  double x = vehCart.x;
  double y = vehCart.y;
  double z = vehCart.z;
  double t = vehCart.t;

  unrotate(&vehCart);
  toLatitude(vehCart, vehLL);
  toLongitude(vehCart, vehLL);

  vehLL->alt = sqrt(pow(vehCart.x,2) + pow(vehCart.y,2) + pow(vehCart.z,2)) - R;
  vehLL->t = vehCart.t
}


// undo givens rotation (19): x = R(-\alpha)x
// note: edits vehicle inside convertCoords, NOT main
void unrotate(struct vehInstance *v)
{
  double tx, ty;
  double alpha = -2.0*pi*v.t/s; // sign flipped for undoing

  tx = v->x*cos(alpha) - v->y*sin(alpha);
  ty = v->x*sin(alpha) + v->y*cos(alpha);
  v->x = tx;
  v->y = ty;
}


void toLatitude(struct vehInstance v, struct latlong *ll)
{
  double ps;
  int psm, psd;
  //conversions on p.11
  if ( v.x*v.x+v.y*v.y != 0 ) {
    ps = atan(v.z/(sqrt(v.x*v.x+v.y*v.y)));
  } 
  else if (v.z < 0) 
    {
      ps = pi/2;
    } 
  else 
    ps = -pi/2;

 // hemisphere determinations
  ps = 360.0*ps/(2*pi);
  if (ps < 0) 
    {
      ll->NS = -1;
      ps = ps*(-1.0);
    }
  else 
    ll->NS = 1;

  // decompose into degrees, minutes, seconds
  // round down to get degrees, use remainder for minutes, seconds
  double temp;
  ll->psd = (int) ps; 
  temp = (ps - psd)*60.0; // minutes + remainder
  ll->psm = (int) temp;
  ll->pss = (temp - psm)*60; //seconds + remainder (as desired)

}

void toLongitude(struct vehInstance v, struct latlong ll)
{
  // FIX FOR RECEIVER: need to include x=0 cases
  if (x > 0.0) {
    if (y > 0.0) 
      lm = atan(y/x);
    else
      lm = 2.0*pi + atan(y/x);
  }  
  else if (x == 0.0) // longitude (lambda) could be anything
    {
      lm = 0.0; 
    }
  else
    lm = pi + atan(y/x);

  // hemisphere determination
  lm = 360.0*lm/(2.0*pi);
  if (lm > 180.0) {
    ll->EW = -1;
    lm = 360.0 - lm;
  }
  else 
    ll->EW = 1;

  //printf("ps: %lf lm: %lf\n", ps, lm);

  // decompose into degrees, minutes, seconds
  ll->lmd = (int) lm; 
  temp = (lm - lmd)*60.0; // minutes + remainder
  ll->lmm = (int) temp;
  ll->lms = (temp - lmm)*60; //seconds + remainder (as desired)
}


