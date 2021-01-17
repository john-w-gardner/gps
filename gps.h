/* header for satellite.c and receiver.c
includes functions from satFuns.c and recFuns.c
*/

// macros
#define NSAT 24 // total number of satellites
#define MAXSTEP 100 // max newton iterations
#define NARGS 10 // number of fields from vehicle 
#define LINELEN 200

// debug 
#define CHECKSAT 0
#define CHECKITER 0
#define CHECKF37 0
#define CHECKFP39 0
#define CHECKVEH 0
#define CHECKx35 0

/* global constants: (read from datad.dat)
   pi = radians in semicircle
   c = speed of light
   R = radius of earth
   s = sidereal day, seconds 
*/
double pi, c, R, s; 

// struct definitions
// position in R3 plus time
struct vehR3 {
  double x;
  double y;
  double z;
  double t;
};

// satellite position in R3 plus horizon indicator
struct satR3 {
  double x;
  double y;
  double z;
  double t;
  int above; 
};

// vehicle data before conversion 
struct latlong { 
  double t; // vehicle time
  int psd; // psi_d, indicating degree of latitude
  int psm; // latitude minutes
  double pss; // latitude seconds
  int NS; // north/south hemisphere
  int lmd; //lambda_d, degree of longitude
  int lmm; // minutes
  double lms; // seconds 
  int EW; // east/west hemisphere
  double alt; //altitude
};

// satellite read from data.dat
struct satellite {
  int index; 
  double u1; 
  double u2; 
  double u3; 
  double v1;
  double v2;
  double v3;
  double per; // periodicity of orbit - should be same for all but..
  double alt; // altitude - should be same for all but..
  double phase; // offsetting angle, in radians
};


// SATELLITE ONLY FUNCTIONS
int getVeh(struct vehR3 *v);
void rotateVeh(struct vehR3 *r) ;
void readDataDat(double *pi, double *c, double *R, double *s, struct satellite (*sArray)[]);
void x35(struct satellite s, double (*l)[], double t);
double computet_s(struct vehR3 v, struct satellite s);
double f37(struct vehR3 v, struct satellite s, double t);
double fp39(struct vehR3 v, struct satellite s, double t);
double computeInitialGuess(struct vehR3 v, struct satellite s);
//void findSatsAbove(struct vehR3 xv, struct vehR3 *satR3[]);
int isAbove(struct vehR3 v, struct satellite s);
//double dotR3(double u[], double v[]);
void scaleR3(double (*u)[], double c);
//void addR3(double u[], double v[], double (*res)[]);
void satelliteStatus(struct satellite s, double r3[], double R, double pi);
void f37status(struct vehR3 v, struct satellite s, double xs[3], double R, double pi, double c, double t, double f);
void iterStatus(struct vehR3 v, struct satellite s, double f, double fp, double tnext, double tlast, double err, double tol, int step);

// FUNCTIONS FOUND IN BOTH
void convertVeh(struct latlong *l, struct vehR3 *r) ;

// RECEIVER ONLY FUNCTIONS
void printSatellite(struct satR3 s);
void printVehicle(struct latlong v);
int computeVehLocation(struct satR3 satarr[], struct vehR3 *veh, int nsat);
void computeF(struct satR3 satarr[], int nsat, struct vehR3 veh, double (*F)[]);
double computeAi(struct satR3 si, struct satR3 si1, struct vehR3 veh);
double computeNi(struct satR3 s, struct vehR3 v);
double computeXi(struct satR3 si, struct satR3 si1, struct vehR3 v);
double computeYi(struct satR3 si, struct satR3 si1, struct vehR3 v);
double computeZi(struct satR3 si, struct satR3 si1, struct vehR3 v);
void computeJ(double (*J)[], struct satR3 satarr[], struct vehR3 veh, int nsat);
double fxx(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Xx(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
double fxy(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Xy(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
double fxz(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Xz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
double fyy(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Yy(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
double fyz(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Yz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
double fzz(struct satR3 satarr[], struct vehR3 veh, int nsat);
double Zz(struct satR3 si, struct satR3 si1, struct vehR3 v, int nsat);
void readDataDatII(double *pi, double *c, double *R, double *s);
double computeVehTime(struct satR3 satarr[], int n, struct vehR3 veh);
void convertCoords(struct vehR3 vehCart, struct latlong *vehLL);
void unrotate(struct vehR3 *v);
void toLatitude(struct vehR3 v, struct latlong *ll);
void toLongitude(struct vehR3 v, struct latlong ll);
void computeInitialGuess(struct vehR3 *v);


