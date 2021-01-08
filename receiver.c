/* receiver component of gps system
compute vehicle location and time using satellite data
1. collect satellite data from stdin, constants from data.dat
2. compute x_v and t_v
3. unrotate, print

*/

#include <stdio.h>
#include <math.h>
#include "receiver.h"

#define LINELEN 200
#define NSAT 24/* length of array holding visible satellites
                   impossible for all 24 to be visible, but.. */

int main()
{
  // get constants
  //readDataDat();

  // set up line buffer
  char line[LINELEN];
  int epoch, done, nsat;
  epoch = done = nsat = 0;

  // set up satellite array
  struct satInstance satArray[NSAT]; 

  // while new epoch of satellites available
  while (getSatelliteArray(satArray, &nsat))
    {
      epoch++; 
      printSatellite(satArray[5]);
      /*
      // calculations
      computeVehLocation(satArray, veh);
      computet_v(satArray[0], veh); // test here: compare results using several different satellites 
      // convert coords and print
      unrotate(veh);
      printveh(veh);
     
      updateInitialGuess(); // for the next iteration/epoch
      */
    } 
}

void printSatellite(struct satInstance s)
{
  printf("satellite %d: \n", s.idx);
  printf("x %.15lf y %.15lf z %.15lf t %.15lf\n", s.x, s.y, s.z, s.t);
}
