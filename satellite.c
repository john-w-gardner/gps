#include <stdio.h>
#include <math.h>
#include "gps.h"
#include "mat.h"

/* satellite component of gps system 

1. read vehicle data and data.dat to set up satellites

2. determine which satellites are above the horizon
 
3. for each satellite above horizon, compute t_s and x_s(t_s)

4. write i_s t_s x_s to std_out and satellite.log
*/


// function prototypes

// go
int main() {
  
  // 1. read data.dat, set up satellites
  struct satellite satArr[NSAT]; 
  readDataDat(&pi, &c, &R, &s, &satArr);
  
  struct vehR3 veh;
  int status;
  double currentSat[3];
  double t_s;
  int i;

  // 2. read each vehicle line from stdin
  while((status = getVeh(&veh)) == NARGS) {
    
    // 3. find satellites above horizon
    for (i=0; i < NSAT; i++) {

      if (isAbove(veh, satArr[i]) == 1) {
        
        // 4. compute t_s, x_s 
        t_s = computet_s(veh, satArr[i]);
        x35(satArr[i], &currentSat, t_s); // update satellite location
        // 5. write to stdout (, satellite.log)
        printf("%d %.15lf %.15lf %.15lf %.15lf\n", i, t_s, currentSat[0], currentSat[1], currentSat[2]);
      }
    }
  }
}

