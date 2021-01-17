/* receiver component of gps system
compute vehicle location and time using satellite data
1. collect satellite data from stdin, constants from data.dat
2. compute x_v and t_v
3. unrotate, print
*/

#include <stdio.h>
#include <math.h>
#include "gps.h"
#include "mat.h"


int main()
{
  // get constants
  void readDataDatII(double *pi, double *c, double *R, double *s);

  // set up line buffer
  char line[LINELEN];
  int epoch, done, nsat, status;
  epoch = done = nsat = status = 0;
  int i;

  // initial position
  void computeInitialGuess(&vehicle);

  // set up satellite array
 struct satInstance satArray[NSAT]; 

  // while new epoch of satellites available
  while (done == 0)
    {
      status = getSatelliteArray(satArray, &nsat);
      epoch++; 
      printf("epoch: %d nsat: %d\n", epoch, nsat);

      // calculations
      computeVehLocation(satArray, &vehicle, nsat);
      vehicle.t = computeVehicleTime(satArray, vehicle); 

      // convert coords and print
      convertCoords(vehicle, &cartVehicle);
      printVehicle(cartVehicle);

      if (status == -1)
        done = 1;
    } 
}

