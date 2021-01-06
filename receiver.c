/* receiver component of gps system
compute vehicle location and time using satellite data
1. collect satellite data from stdin, constants from data.dat
2. compute x_v and t_v
3. unrotate, print

*/

#include <stdio.h>
#include <math.h>

#define LINELEN 200
#define NSAT 24/* length of array holding visible satellites
                   impossible for all 24 to be visible, but.. */



struct satInstance {
  double x;
  double y;
  double z;
  double t;
  double idx;
};

struct vehInstance {
  double x;
  double y;
  double z;
  double t;
};

int main()
{
  // get constants
  readDataDat();

  // set up line buffer
  char line[LINELEN];
  int epoch, done, nsat;
  epoch = done = nsat = 0;

  // set up satellite array
  struct satInstance satArray[NSAT]; 

  // while new epoch of satellites available
  while (getSatelliteArray(&satArray, &nsat))
    {
      epoch++; 
      // calculations
      computeVehLocation(satArray, veh);
      computet_v(satArray[0], veh); /* test here: compare results using several different satellites */
      // convert coords and print
      unrotate(veh);
      printveh(veh);
     
      updateInitialGuess(); // for the next iteration/epoch
    } 
}


/* 
ungetSat: push satellite back onto input (buffer) 
 */
int ungetSat(struct satInstance *satBuf, struct satInstance *s)
{

}






// jam contents of s2 into s1
void updateSat(struct satInstance *s1, struct satInstance *s2)
{
  s1->idx = s2->idx;
  s1->x = s2->x;
  s1->y = s2->y;
  s1->z = s2->z;
  s1->t = s2->t;
}
