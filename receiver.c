#include <stdio.h>
#include <math.h>
#include "gps.h"
#include "mat.h"

/* receiver component of gps system
compute vehicle location and time using satellite data
1. collect satellite data from stdin, constants from data.dat
2. compute x_v and t_v
3. unrotate, print
*/



int main()
{  
  printf("testing..\n");

  // get constants
  readDataDatII(&pi, &c, &R, &s);

  // setup 
  int epoch, done, nsat, status;
  epoch = done = nsat = status = 0;
  int i;
  struct vehR3 vehicle;
  struct latlong cartVehicle;

  // initial position
  computeInitialVehicle(&vehicle);

  // set up satellite array
  struct satR3 satArray[NSAT]; 

  // while new epoch of satellites available
  while (done == 0)
    {
      status = getSatelliteArray(satArray, &nsat);
      epoch++; 
      printf("epoch: %d nsat: %d\n", epoch, nsat);

      // calculations
      computeVehLocation(satArray, &vehicle, nsat);
      vehicle.t = computeVehicleTime(satArray, nsat, vehicle); 

      // convert coords and print
      convertCoords(vehicle, &cartVehicle);
      printVehicle(cartVehicle);

      if (status == -1)
        done = 1;
    } 
}

