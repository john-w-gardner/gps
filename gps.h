/* header for satellite.c and receiver.c
includes functions from satFuns.c and recFuns.c
*/

// macros
#define NSAT 24 // total number of satellites
#define MAXSTEP 100 // max newton iterations
#define NARGS 10 // number of fields from vehicle 
#define LINELEN 200


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
  int idx;
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

// gps.c
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
//void printSatellite(struct satR3 s);
void printVehicle(struct latlong v);
void printVehR3(struct vehR3 v);
void printSatR3(struct satR3 s);
void computeVehLocation(struct satR3 satarr[], struct vehR3 *veh, int nsat);
void computeF(struct satR3 satarr[], int nsat, struct vehR3 veh, double F[]);
double computeAi(struct satR3 si, struct satR3 si1, struct vehR3 veh);
double computeNi(struct satR3 s, struct vehR3 v);
double computeXi(struct satR3 si, struct satR3 si1, struct vehR3 v);
double computeYi(struct satR3 si, struct satR3 si1, struct vehR3 v);
double computeZi(struct satR3 si, struct satR3 si1, struct vehR3 v);
void computeJ(struct satR3 satarr[], int nsat, struct vehR3 veh, double J[]);
//void computeJ(struct satR3 satarr[], int nsat, struct vehR3 veh, double (*J)[]);
double fxx(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeXx(struct satR3 si, struct satR3 si1, struct vehR3 v);
double fxy(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeXy(struct satR3 si, struct satR3 si1, struct vehR3 v);
double fxz(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeXz(struct satR3 si, struct satR3 si1, struct vehR3 v);
double fyy(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeYy(struct satR3 si, struct satR3 si1, struct vehR3 v);
double fyz(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeYz(struct satR3 si, struct satR3 si1, struct vehR3 v);
double fzz(struct satR3 satarr[], int nsat, struct vehR3 veh);
double computeZz(struct satR3 si, struct satR3 si1, struct vehR3 v);
void updateVehLocation(struct vehR3 *v, double s[]);
void readDataDatII(double *pi, double *c, double *R, double *s);
double computeVehicleTime(struct satR3 satarr[], int n, struct vehR3 veh);
void convertCoords(struct vehR3 vehCart, struct latlong *vehLL);
void unrotate(struct vehR3 *v);
void toLatitude(struct vehR3 v, struct latlong *ll);
void toLongitude(struct vehR3 v, struct latlong *ll);
void computeInitialVehicle(struct vehR3 *v);
// end gps.c

// from getSatellite.c
int getSatelliteArray(struct satR3 *satarr, int *nsat); 
void updateSatellite(struct satR3 *s1, struct satR3 *s2);
int getSatellite(struct satR3 *s);
