#include <stdio.h>
#include <math.h>
#include "gps.h"
#include "mat.h"
/* 
helper functions for both satellite and receiver
*/

// read vehicle data from std_in 
// returns 1 if another line found, else 0
int getVeh(struct vehR3 *v){
  //extern double pi, R;
  struct latlong l;
  int scanCount; // 

  // grab data from stdin
  scanCount = scanf("%lf %d %d %lf %d %d %d %lf %d %lf", 
                    &l.t,
                    &l.psd, &l.psm, &l.pss, &l.NS,
                    &l.lmd, &l.lmm, &l.lms, &l.EW,
                    &l.alt);
  
  if (CHECKVEH == 1) {
    printf("vehicle data check:\n");
    printf("%lf %d %d %lf %d %d %d %lf %d %lf\n", 
                    l.t,
                    l.psd, l.psm, l.pss, l.NS,
                    l.lmd, l.lmm, l.lms, l.EW,
                    l.alt);    
  }
  
  // convert to R^3
  convertVeh(&l, v);
  
  // rotate to account for time
  rotateVeh(v);

  return scanCount;
}

// convert vehicle data from lat/long to R^3
void convertVeh(struct latlong *l, struct vehR3 *r) {

  //extern double pi, R;

  // gather lat/long values into one double each
  double ps = 2*pi* l->NS * (l->psd/360.0 + l->psm/(360.0*60.0) + l->pss/(360.0*60.0*60.0));
  double lm = 2*pi* l->EW * (l->lmd/360.0 + l->lmm/(360.0*60.0) + l->lms/(360.0*60.0*60.0));

  double mag = (R + l->alt);
  // convert to x,y,z
  r->x = mag*cos(ps)*cos(lm);
  r->y = mag*cos(ps)*sin(lm);
  r->z = mag*sin(ps);
  r->t = l->t; // time unchanged

  if (CHECKVEH == 1) {
    printf("vehicle in R3 (pre-rotate)\n");
    printf("R: %.15lf alt: %.15lf\n", R, l->alt);
    printf("psi: %.15lf lambda: %.15lf mag: %.15lf t: %.15lf\n", ps, lm, mag, r->t);
    printf("x: %.15lf y: %.15lf z: %.15lf\n", r->x, r->y, r->z);
  }
}

// rotate vehicle position in R^3 to account for time
// see matrix on handout p.11 
void rotateVeh(struct vehR3 *r) {
  //extern double pi, s;

  double alpha = (2*pi*r->t)/s;
  
  // need temp variables for matrix mult
  double x, y;

  // givens rotate
  x = r->x*cos(alpha) - r->y*sin(alpha); 
  y = r->x*sin(alpha) + r->y*cos(alpha);

  // update 
  r->x = x;
  r->y = y;

  if (CHECKVEH == 1) {
    printf("vehicle in R3 (post-rotate)\n");
    printf("pi: %.15lf s: %.15lf\n", pi, s);
    printf("alpha: %.15lf, t: %.15lf\n", alpha, r->t);
    printf("x: %.15lf y: %.15lf z: %.15lf\n", r->x, r->y, r->z);
  }

}


// read data.dat
void readDataDat(double *pi, double *c, double *R, double *s, struct satellite (*sArray)[] ) {
  
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

  // grab satellites
  int i;
  for (i = 0; i < NSAT; i++) {
    //printf("new satellite %d\n", i);

    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u1); // would be simpler with array of pointers..
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u2); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u3); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v1); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v2); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v3); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].per); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].alt); 
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].phase); 

    (*sArray)[i].index = i;

    if (CHECKSAT == 1) {
      double sR3[3];
      printf("satellite test in readDataDat: \n");
      x35((*sArray)[i], &sR3, 0.0);
      satelliteStatus((*sArray)[i], sR3, *R, *pi);
    }
  }
  fclose(dd);
}


/*
// satellite position in R3 at given time, eqn 35
double x35(struct satellite s, double t) {
  extern double pi, R;
*/
void x35(struct satellite s, double (*l)[], double t) {
  //extern double R, pi;
 
  double mag = R + s.alt; // satellite distance from origin
  double arg = (2.0*pi*t)/s.per + s.phase; // need to rotate u/v to account for time

  (*l)[0] = mag * (s.u1*cos(arg) + s.v1*sin(arg) );
  (*l)[1] = mag * (s.u2*cos(arg) + s.v2*sin(arg) );
  (*l)[2] = mag * (s.u3*cos(arg) + s.v3*sin(arg) );

  if (CHECKx35 == 1) {
    printf("inside x35:\n");
    satelliteStatus(s, *l, R, pi);
    printf("exit x35\n");
  }
}

// satellite visibility determination
// sat is above horizon if x_v \dot x_s > x_v \dot x_v
int isAbove(struct vehR3 v, struct satellite s) {
  // setup vehicle vector
  double xv[] = {v.x, v.y, v.z};
  double xs[3];
  x35(s, &xs, v.t);

  double vDOTv = dotprod(xv, xv, 3);
  double vDOTs = dotprod(xs, xv, 3);

  if (CHECKSAT == 1) {
    printf("sat coords: %.15lf %.15lf %.15lf\n", xs[0], xs[1], xs[2]);
    printf("veh coords: %.15lf %.15lf %.15lf\n", xv[0], xv[1], xv[2]);
    printf("vDOTs: %.15lf vDOTv: %.15lf\n", vDOTs, vDOTv);
  }

  if (vDOTs > vDOTv) 
    return 1;
  else 
    return -1;
}

/* how to compute t_s:
   starting at t0
   t_{k+1} = t_k - f(t_k)/f'(t_k) until convergence
   where
   f(t) = (x_s(t) - x_v)^T (x_s(t) - x_v) - c^2(t_v-t)^2 = 0  (eqn 37) 
   f'(t) = 4pi(R+h)/p (x_s(t) - x_v)^T (-u*sin(alpha) + v*cos(alpha)) + 2c^2(t_v-t)  (eqn 39)
   alpha = 2*pi*t/p + phase
   and u and v are vectors determining satellite orbit (normalized)

   initial guess t_0 = t_v - ||x_s(t_v) - x_v||/c
*/
 
double computet_s(struct vehR3 v, struct satellite s) {
  double f, fp, tnext, tlast;
  double err = 0.1;
  double tol = 0.01/c;
  int step = 0; // count iterations

  // x_s(t_v) - x_v
  tlast = computeInitialGuess(v, s);
  
  if (CHECKITER == 1) {
    f = f37(v, s, tlast); // numerator
    iterStatus(v, s, f, fp, tnext, tlast, err, tol, step);
  }

  // newton iteration
  while ((err > tol) && step++ < MAXSTEP) {
    f = f37(v, s, tlast); // numerator
    fp = fp39(v, s, tlast); // denominator
    tnext = tlast - f/fp; // iteration update
    err = fabs(tnext - tlast);
    tlast = tnext;
    if (CHECKITER == 1) {
      printf("iterStatus inside loop:\n");
      iterStatus(v, s, f, fp, tnext, tlast, err, tol, step);
    }
  } 

  // check and return
  if (err > tol) {
    if (CHECKITER == 1) {
      printf("iteration failed..\n f: %.15lf |t_k - t_{k-1}|: %.15lf\n", f, err);
      printf("final data: step: %d f: %.15lf  fp: %.15lf  err: %.15lf\n", step, f, fp, err);
      printf("final data: %d  -  tnext: %.15lf  tlast: %.15lf tol: %.15lf\n", step, tnext, tlast, tol);
    }
    return tnext;
  } else {
    if (CHECKITER == 1) {
      printf("newton iteration converged in %d steps.\n", step);
      printf("final data: step: %d f: %.15lf  fp: %.15lf  err: %.15lf\n", step, f, fp, err);
      printf("final data: %d  -  tnext: %.15lf  tlast: %.15lf tol: %.15lf\n", step, tnext, tlast, tol);
    }
    return tnext;
  }
}

// f(t) given on p.14 eqn 37
double f37(struct vehR3 v, struct satellite s, double t) {

  double result = 0;
  double xs[3]; // satellite in R^3
  double xv[] = {v.x, v.y, v.z}; // vehicle in R^3

  x35(s, &xs, t); //setup x_s(t) vector

  if (CHECKF37 == 1) {
    printf("pre-compute f37:\n");
    f37status(v, s, xs, R, pi, c, t, result);
  }
 
  //(x_s - x_v) \dot (x_s - x_v)
  scaleVec(&xv, -1.0, 3);
  addVec(&xs, xv, 3); // (x_s - x_v) stored in xs
  result = dotprod(xs, xs, 3);
  result = result - c*c*(v.t - t)*(v.t - t);

  if (CHECKF37 == 1) {
    printf("post-compute f37:\n");
    f37status(v, s, xs, R, pi, c, t, result);
  }

  return result;
}

// f'(t) given on p.14 eqn 39
double fp39(struct vehR3 v, struct satellite s, double t) {  
  //extern double pi, R, c;
  double result = 0;

  // need some vectors..
  double xs[3]; // satellite in R^3 i.e. x_s(t)
  x35(s, &xs, t); // compute x_s(t)
  double xv[] = {v.x, v.y, v.z}; // vehicle in R^3 i.e. x_v
  double U[3] = {s.u1, s.u2, s.u3}; // satellite direction vectors
  double V[3] = {s.v1, s.v2, s.v3};

  double lead = 4.0*pi*(R+s.alt)/s.per; // leading scalar product
  double arg = 2.0*pi*t/s.per + s.phase; // sin,cos argument

  scaleVec(&xv, -1.0, 3);
  addVec(&xs, xv, 3); // overwrite x_s with (x_s(t) - x_v)
  scaleVec(&xs, lead, 3); 

  // compute trig-scaled satellite direction vectors in []
  scaleVec(&U, -sin(arg), n);
  scaleVec(&V, cos(arg), n);
  addVec(&U, V, 3);

  // (x_s-x_v) \dot [-usin() + vcos()]
  result = dotprod(xs, U, 3);

  result += 2*c*c*(v.t - t);

  return result;
}

/* starting value for newton iteration */
double computeInitialGuess(struct vehR3 v, struct satellite s) {
  //extern double c, R, pi;

  // setup
  double xv[3] = {v.x, v.y, v.z};
  double xs[3];
  double diff[3]; //x_s - x_v
  x35(s, &xs, v.t);
  
  // subtract and ||.||_2
  scaleVec(&xv, -1.0, 3);
  addVec(&xs, xv, 3);
  double result = sqrt(dotprod(xs, xs, 3));
  
  result = v.t - result/c;
  return result;
}

/* some simple linear algebra helpers
   might use blas instead... */

/*
// dot product in R^3
// u and v MUST be dim-3 arrays 
double dotR3(double u[], double v[]) {
  double value = 0;

  int i;
  for (i=0; i < 3; i++) {
    value += u[i]*v[i];
  }
  return value;
}


// scale vector in R^3
void scaleR3(double (*u)[], double c) {
  int i;
  for (i=0; i < 3; i++) 
    (*u)[i] = c*(*u)[i];
}

// add vectors in R^3
void addR3(double u[], double v[], double (*res)[]) {
  int i;
  for (i=0; i < 3; i++) 
    (*res)[i] = u[i] + v[i];
}

*/

/* debugging routines */

// compare satellite struct to converted R^3
void satelliteStatus(struct satellite s, double r3[], double R, double pi) 
{
  printf("--- begin satelliteStatus ---\n");
  printf("constants: R %.15lf pi %.15lf\n", R, pi);
  printf("satellite %d data:\n", s.index);
  printf("U: %.15lf %.15lf %.15lf\n", s.u1, s.u2, s.u3);
  printf("V: %.15lf %.15lf %.15lf \n", s.v1, s.v2, s.v3);
  printf("per: %.15lf alt: %.15lf phase: %.15lf\n", s.per, s.alt, s.phase);
  printf("converted coordinates: %.15lf %.15lf %.15lf\n", r3[0], r3[1], r3[2]);
  printf("norm test: %.15lf\n", fabs(sqrt(dotprod(r3, r3, 3))) - (R+s.alt));
  printf("--- end satelliteStatus ---\n");
}

// check f37
void f37status(struct vehR3 v, struct satellite s, double xs[3], double R, double pi, double c, double t, double f)
{
    printf("--- debug f37 ---\n");
    printf("result: %.15lf \n", f);
    printf("consts: R %.15lf pi %.15lf c %.15lf t %.15lf\n", R, pi, c, t);
    printf("veh: v.x %.15lf v.y %.15lf v.z %.15lf v.t %.15lf \n", v.x, v.y, v.z, v.t);
    printf("satellite %d data: U: %.15lf %.15lf %.15lf\n", 
           s.index, s.u1, s.u2, s.u3);
    printf("V: %.15lf %.15lf %.15lf \n", s.v1, s.v2, s.v3);
    printf("per: %.15lf alt: %.15lf phase: %.15lf\n", s.per, s.alt, s.phase);
    printf("sat R^3: %.15lf %.15lf %.15lf\n", xs[0], xs[1], xs[2]);
}

// ping from inside newton iteration
void iterStatus(struct vehR3 v, struct satellite s, double f, double fp, double tnext, double tlast, double err, double tol, int step) 
{
  printf("--- ping at iteration step %d ---\n", step);
  printf("veh loc: %.15lf %.15lf %.15lf v.t: %.15lf\n", v.x, v.y, v.z, v.t);
  printf("sat: U %.15lf %.15lf %.15lf V %.15lf %.15lf %.15lf \n", 
         s.u1, s.u2, s.u3, s.v1, s.v2, s.v3);
  printf("sat: alt %.15lf per %.15lf phase %.15lf \n", s.alt, s.per, s.phase);
  printf("f: %.15lf fp: %.15lf tnext %.15lf tlast %.15lf\n", f, fp, tnext, tlast);
  printf("err %.15lf tol %.15lf\n", err, tol);
}



/* 
original cartesian to lat/long conversion for testing purposes in satellite
version in use found in receiver


// going backwards: convert from R^3 to lat/long decomposed as in vehicle output
// (testing purposes - only receiver actually uses this)
void unrotate(double x, double y, double z, double t) {
  //extern double pi, s, R;
  int EW, NS;
  int psm, psd, lmm, lmd;
  double ps, lm;
  double pss, lms;
  double h;

  printf("--constants pi: %lf s: %lf R: %lf\n", pi, s, R);
  // first undo givens rotation
  double alpha = -2.0*pi*t/s; // sign flipped for undoing

  double tx, ty;
  tx = x*cos(alpha) - y*sin(alpha);
  ty = x*sin(alpha) + y*cos(alpha);
  x = tx;
  y = ty;

  // altitude 
  h = sqrt(x*x + y*y + z*z) - R;

  //conversions on p.11
  if ( x*x+y*y != 0 ) {
    ps = atan(z/(sqrt(x*x+y*y)));
  } else if (z < 0) {
    ps = pi/2;
  } else 
    ps = -pi/2;

  // FIX FOR RECEIVER: need to include x=0 cases
  if (x > 0) {
    if (y > 0) 
      lm = atan(y/x);
    else
      lm = 2*pi + atan(y/x);
  }  else
    lm = pi + atan(y/x);

  // hemisphere determinations
  ps = 360.0*ps/(2*pi);
  if (ps < 0) {
    NS = -1;
    ps = ps*(-1.0);
  }
  else 
    NS = 1;

  lm = 360.0*lm/(2*pi);
  if (lm > 180.0) {
    EW = -1;
    lm = 360.0 - lm;
  }
  else 
    EW = 1;

  printf("ps: %lf lm: %lf\n", ps, lm);
  // decompose into degrees, minutes, seconds
  // round down to get degrees, use remainder for minutes, seconds
  double temp;
  psd = (int) ps; 
  temp = (ps - psd)*60.0; // minutes + remainder
  psm = (int) temp;
  pss = (temp - psm)*60; //seconds + remainder (as desired)

  // same for longitude
  lmd = (int) lm; 
  temp = (lm - lmd)*60.0; // minutes + remainder
  lmm = (int) temp;
  lms = (temp - lmm)*60; //seconds + remainder (as desired)


  printf("%lf %d %d %lf %d %d %d %lf %d %lf\n",
         t, psd, psm, pss, NS, lmd, lmm, lms, EW, h);

}
*/


void printSatellite(struct satR3 s)
{
  printf("satellite %d: \n", s.idx);
  printf("x %.15lf y %.15lf z %.15lf t %.15lf\n", s.x, s.y, s.z, s.t);
}

void printVehicle(struct latlong v)
{
  printf("%lf %d %d %lf %d %d %d %lf %d %lf \n", 
         v.t, v.psd, v.psm, v.pss, v.NS, v.lmd, v.lmm, v.lms, v.EW, v.alt);
}


/* computeVehLocation():
   newton iteration for nonlinear system, see (58)
   f is least squares function we want to minimize
   F is its gradient, J the jacobian of F
   until satisfied, solve Js = -F, set x = x + s
 */
int computeVehLocation(struct satR3 satarr[], struct vehR3 *veh, int nsat)
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
void computeF(struct satR3 satarr[], int nsat, struct vehR3 veh, double (*F)[])
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


double computeAi(struct satR3 si, struct satR3 si1, struct vehR3 veh)
{
  double Ni = computeNi(si, veh);
  double Ni1 = computeNi(si1, veh);
  
  double result = Ni1 - Ni - c*(si.t - si1.t);
  return result;              
}

double computeNi(struct satR3 s, struct vehR3 v)
{
  // matops lib needs vectors
  double xs[] = {s.x, s.y, s.z};
  double xv[] = {v.x, v.y, v.z};

  double result = eucmetric(xs, xv, 3);
  return result;
}

double computeXi(struct satR3 si, struct satR3 si1, struct vehR3 v)
{
  double result;
  double Ni = computeNi(si, v);
  double Ni1 = computeNi(si1, v);

  result = -(si1.x - v.x)/Ni1 + (si.x - v.x)/Ni;
} 

double computeYi(struct satR3 si, struct satR3 si1, struct vehR3 v)
{
  double result;
  double Ni = computeNi(si, v);
  double Ni1 = computeNi(si1, v);

  result = -(si1.y - v.y)/Ni1 + (si.y - v.y)/Ni;
} 

double computeZi(struct satR3 si, struct satR3 si1, struct vehR3 v)
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
void computeJ(double (*J)[], struct satR3 satarr[], struct vehR3 veh, int nsat)
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
   careful: these are all very similar 
*/
double fxx(struct satR3 satarr[], struct vehR3 veh, int nsat)
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
double Xx(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = (pow(Ni1, 2) - pow((si1.x - v.x), 2))/pow(Ni1, 3);
  b = (pow(Ni, 2) - pow((si.x - v.x), 2))/pow(Ni, 3);
  
  result = a-b;
  return result;
}

double fxy(struct satR3 satarr[], struct vehR3 veh, int nsat)
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

double Xy(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.y - v.y)*(si1.x - v.x)/pow(Ni1,3);
  b = (si.y - v.y)*(si.x - v.x)/pow(Ni,3);

  result = a+b;
  return result;
}

double fxz(struct satR3 satarr[], struct vehR3 veh, int nsat)
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

double Xz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
{
  double Ni, Ni1, a, b, result;

  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.z - v.z)*(si1.x - v.x)/pow(Ni1,3);
  b = (si.z - v.z)*(si.x - v.x)/pow(Ni,3);

  result = a+b;
  return result;
}

double fyy(struct satR3 satarr[], struct vehR3 veh, int nsat)
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

double Yy(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = (pow(Ni1, 2) - pow((si1.y - v.y), 2))/pow(Ni1, 3);
  b = (pow(Ni, 2) - pow((si.y - v.y), 2))/pow(Ni, 3);
  
  result = a-b;
  return result;
}

double fyz(struct satR3 satarr[], struct vehR3 veh, int nsat)
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

double Yz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
{
  double Ni, Ni1, a, b, result;
  Ni = computeNi(si, v);
  Ni1 = computeNi(si1, v);

  a = -(si1.z - v.z)*(si1.y - v.y)/pow(Ni1,3);
  b = (si.z - v.z)*(si.y - v.y)/pow(Ni,3);

  result = a+b;
  return result;
}

double fzz(struct satR3 satarr[], struct vehR3 veh, int nsat)
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

double Zz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat)
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
double computeVehTime(struct satR3 satarr[], int n, struct vehR3 veh)
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
void convertCoords(struct vehR3 vehCart, struct latlong *vehLL) 
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
void unrotate(struct vehR3 *v)
{
  double tx, ty;
  double alpha = -2.0*pi*v.t/s; // sign flipped for undoing

  tx = v->x*cos(alpha) - v->y*sin(alpha);
  ty = v->x*sin(alpha) + v->y*cos(alpha);
  v->x = tx;
  v->y = ty;
}


void toLatitude(struct vehR3 v, struct latlong *ll)
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

void toLongitude(struct vehR3 v, struct latlong ll)
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


void computeInitialGuess(struct vehR3 *v)
{
  // using b12 (11) as starting point 
  struct latlong ll = {0.0, 40, 45, 55.0, 1, 111, 50, 58.0, -1, 1372.0};

  convertVeh(&ll, v);
}
