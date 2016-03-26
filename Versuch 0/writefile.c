#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

int main ()
{
  double x;
  int i=0;
  FILE *fp;
  fp = fopen("OutputFile_c.txt", "w");
  for (x=0.0; x<6.29; x=x+0.1)
  {
    i=i+1;
    fprintf(fp,"%2d\t%11.5f\t%11.5e\n",i,x,sin(x));
  }
  fclose(fp);

  return 0;
}
