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
  int epoch = 0; // number of lines we end up returning
  char done = 0; // set to 1 when finished
  int nsat = 0;  // number of satellites visible in current epoch

  // set up satellite array
  struct satInstance satArray[NSAT]; 

  //  while (getSat(&satArray[0]) != NULL)
  while (getSat(&satArray[nsat]))
    {
      epoch++; 
      nsat = 0;
      getSat(&satArray[nsat]);

      while (getSat(&satArray[nsat+1]) > satArray[nsat].idx) { 
        // not incrementing nsat in while() bc it appears on both sides of >
        nsat++;
      }

      computeVehLocation(satArray, veh);
      computet_v(satArray[0], veh); /* test here: compare results using several different satellites */
      unrotate(veh);
      printveh(veh);
      //resetlocation();
    
      if (currentSat == -1) {
        //done = 1;
        break;
      }
      else 
        satArray[0] = satArray[nsat]; 
    }
   
}

// get next satellite from stdin
// return satellite index if available, -1 otherwise
int getSat(struct satInstance *satBuf, struct satInstance *s)
{
  char buf[LINELEN];  
  
  if ((*satBuf).idx > -1) 
    {
      // update s with satBuf, set satBuf.idx to -1
      updateSat(
    }
  else
    {
      // get new satellite from stdin, check if new
      if (fgets(buf, LINELEN, stdin) != NULL)    
        {
          sscanf(buf, "%d %lf %lf %lf %lf", &(s->idx), &(s->t), &(s->x), &(s->y), &(s->z));
          
        }
      else 
        

  
    
      
      if (status != NFIELD)
        printf()
          
  */
}

      // jam contents of s2 into s1
void updateSat(struct satInstance *satBuf, struct satInstance *s)
{
  
