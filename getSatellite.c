#include <stdio.h>
#include <stdlib.h>
#include "gps.h"

/* 
   populate satellite array (for receiver.c)
   return number of satellites found or 0 if EOF
   loosely modeled after k&r getch/ungetch

 */

#define MAXLEN 1000

static struct satR3 satbuf = {0,0,0,0,-1};

// read satellites for current epoch
int getSatelliteArray(struct satR3 *satarr, int *nsat)
{
  int i = 0;
  int status = 0;
  struct satR3 *iptr;
  iptr = satarr;

  if (satbuf.idx == -1) // buffer empty
    {
      status = getSatellite(satarr);
    }
  else // put buffer satellite into array
    updateSatellite(satarr, &satbuf);

  // populate array
  while ((status = getSatellite(&satbuf)) > (*satarr++).idx)
    {
      updateSatellite(satarr, &satbuf);
    }      
  
  *nsat = satarr - iptr;

  if (status == -1)
    return -1;
  else
    {
      updateSatellite(&satbuf, satarr);
      return 1; 
    }

}

// put contents of s2 into s1
void updateSatellite(struct satR3 *s1, struct satR3 *s2)
{
  s1->idx = s2->idx;
  s1->x = s2->x;
  s1->y = s2->y;
  s1->z = s2->z;
  s1->t = s2->t;
}

/* read satellite struct from stdin, 
   return index or -1 if EOF
*/
int getSatellite(struct satR3 *s)
{
  char buf[MAXLEN];

  if (fgets(buf, MAXLEN, stdin) != NULL)
    {
      sscanf(buf, "%d %lf %lf %lf %lf", &(s->idx), &(s->t), &(s->x), &(s->y), &(s->z));
      return s->idx;
    }
  else 
    return -1;
}        



/* 
   populate satellite array (for receiver.c)
   return number of satellites found or 0 if EOF
   loosely modeled after k&r getch/ungetch

   main calls getSatelliteArray, at first satbuf is empty.
   all intermediate calls will have nonzero buffer.
   final (nontrivial) call exhausts the buffer, 
   then main calls one more time, 
   getSatelliteArray() finds empty buffer, 
   getSatellite finds nothing, returns EOF, 
   which causes getSatelliteArray to return 0, 
   which main interprets as completion. 
   

   note: satbuf.idx = -1 is used in two ways:
   once at the beginning of main to let us know to go to stdin for first element,
   then at the end to let us know we're done.
 */

