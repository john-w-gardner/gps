#include <stdio.h>
#include <stdlib.h>

/* 
   populate satellite array (for receiver.c)
   return number of satellites found or 0 if EOF
   loosely modeled after k&r getch/ungetch
*/


static struct satInstance satelliteBuffer;
satelliteBuffer.idx = -1;

/* read satellites for current epoch
 */
int getSatelliteArray(struct satInstance *satbuf, struct satInstance *satarr[])
{
  int i = 0;
  struct satInstance tempsat; 

  if (satbuf->idx = -1) // buffer empty
    {
      getSatellite(satarr[i]);
    }
  else // put buffer satellite into array
    updateSatellite(satarr[i], satbuf);
  
  while (getSatellite(&satbuf) > (*satarr)[i++].idx)
    {
      updateSatellite(satarr[i], satbuf);
    }      
  



      if (tempsat.idx > (*satarr)[i].idx)
        {

        }
      else
        {
          
          
    }

}

// put contents of s2 into s1
int updateSatellite(struct satInstance *s1, struct satInstance *s2)
{

}

/* read satellite struct from stdin, 
   return index
*/
int getSatellite(struct satInstance *satBuf, struct satInstance *s)
{
  char buf[LINELEN];
  
  if ((*satBuf).idx > -1)
    {
      // put satBuf into s, reset satBuf
      updateSat(s, &satBuf);
      satBuf.idx = -1;
      return s.idx;
    }
  else
    {
      // get new satellite from stdin
      if (fgets(buf, LINELEN, stdin) != NULL)
        {
          sscanf(buf, "%d %lf %lf %lf %lf", &(s->idx), &(s->t), &(s->x), &(s->y), &(s->z));
          return s.idx
        }
      else 
        {
          return 0;
        }
    }
}
        
