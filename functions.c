#include <stdio.h>
#include <math.h>
#include "testSatFunctions.h"

/* satellite component of gps system 

1. read vehicle data and data.dat to set up satellites

2. determine which satellites are above the horizon
 
3. for each satellite above horizon, compute t_s and x_s(t_s)

4. write i_s t_s x_s to std_out and satellite.log


int main() {
  // 1. read data.dat, set up satellites

  // 2. read each vehicle line from stdin
  while(status = getVeh(&veh) > 0) {
  
  // 3. determine satellites above horizon 

    // 4. compute t_s, x_s for qualifying satellites
    for () {
      
      // 5. write to stdout, (satellite.log)

    }
  }
}

*/

//double R = 6.367444500000000000E+06;
//double S = 8.616408999999999651E+04;

// read vehicle data from std_in 
// returns 1 if another line found, else 0
int getVeh(struct locR3 *v){

  struct latlong l;
  int scanCount; // 

  // grab data from stdin
  scanCount = scanf("%lf %d %d %lf %d %d %d %lf %d %lf", 
                    &l.t,
                    &l.psd, &l.psm, &l.pss,
                    &l.NS,
                    &l.lmd, &l.lmm, &l.lms,
                    &l.EW,
                    &l.alt);
  
  // convert to R^3
  convertVeh(&l, v);
  
  // rotate to account for time
  rotateVeh(v);

  return scanCount;
}

// convert vehicle data from lat/long to R^3
void convertVeh(struct latlong *l, struct locR3 *r) {

  extern double pi, R;

  // gather lat/long values into one double each
  double ps = 2*pi* l->NS * (l->psd/360.0 + l->psm/(360.0*60.0) + l->pss/(360.0*60.0*60.0));
  double lm = 2*pi* l->EW * (l->lmd/360.0 + l->lmm/(360.0*60.0) + l->lms/(360.0*60.0*60.0));
  // printf("psi: %lf lambda: %lf\n", ps, lm);
  
  // convert to x,y,z
  r->x = (R + l->alt)*cos(ps)*cos(lm);
  r->y = (R + l->alt)*cos(ps)*sin(lm);
  r->z = (R + l->alt)*sin(ps);
  r->t = l->t; // time unchanged
  /*
  //test
  double mag = sqrt((r->x)*(r->x) + (r->y)*(r->y) + (r->z)*(r->z));
  printf("magnitude:%lf (R+h): %lf diff: %lf\n", mag, (R+l->alt), mag-(R+l->alt));
  */
}

// rotate vehicle position in R^3 to account for time
// see matrix on handout p.11 
void rotateVeh(struct locR3 *r) {
  extern double pi, s;

  double alpha = (2*pi*r->t)/s;
  
  // need temp variables for matrix mult
  double x, y;

  // givens rotate
  x = r->x*cos(alpha) - r->y*sin(alpha); 
  y = r->x*sin(alpha) + r->y*cos(alpha);

  // update 
  r->x = x;
  r->y = y;

}

// going backwards: convert from R^3 to lat/long decomposed as in vehicle output
// (testing purposes - only receiver uses this)
void unrotate(double x, double y, double z, double t) {
  extern double pi, s, R;
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
  printf("%lf %lf %lf %lf buf: %s status: %d\n", *pi, *c, *R, *s, buf, status);


  // grab satellites
  int i;
  for (i = 0; i < NSAT; i++) {
    //printf("new satellite %d\n", i);

    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u1); // would be simpler with array of pointers..
    //printf("value: %lf status: %d\n", (*sArray)[i].u1, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u2); 
    //printf("value: %lf status: %d\n", (*sArray)[i].u2, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].u3); 
    //printf("value: %lf status: %d\n", (*sArray)[i].u3, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v1); 
    //printf("value: %lf status: %d\n", (*sArray)[i].v1, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v2); 
    //printf("value: %lf status: %d\n", (*sArray)[i].v2, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].v3); 
    //printf("value: %lf status: %d\n", (*sArray)[i].v3, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].per); 
    //printf("value: %lf status: %d\n", (*sArray)[i].per, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].alt); 
    //printf("value: %lf status: %d\n", (*sArray)[i].alt, status);
    fgets(buf, NCHAR, dd);
    status = sscanf(buf, "%lf", &(*sArray)[i].phase); 
    //printf("value: %lf status: %d\n", (*sArray)[i].phase, status);

    (*sArray)[i].index = i;

  }
  fclose(dd);
}


// setellite to R3 conversion - see p14
void sattoR3(struct satellite s, struct locR3 *l, double t) {
  extern double R;
  extern double pi;

  double mag = R + s.alt; // satellite distance from origin
  double arg = (2.0*pi*t)/s.phase; // need to rotate u/v to account for time
  l->x = mag * (s.u1*cos(arg) + s.v1*sin(arg) );
  l->y = mag * (s.u2*cos(arg) + s.v2*sin(arg) );
  l->z = mag * (s.u3*cos(arg) + s.v3*sin(arg) );
}

// satellite visibility determination
// sat is above horizon if x_v \dot x_s > x_v \dot x_v
double dotProd(struct locR3 u, struct locR3 v) {
  double value;
  value = u.x*v.x + u.y*v.y + u.z*v.z;
  return value;
}

