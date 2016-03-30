#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>
#include <limits>

/*
 The  <math.h>  header  shall  provide for the following constants.  The
   values are of type double and are accurate within the precision of  the
   double type.

   M_PI   Value of pi

   M_PI_2 Value of pi/2

   M_PI_4 Value of pi/4

   M_1_PI Value of 1/pi

   M_2_PI Value of 2/pi

   M_2_SQRTPI Value of 2/ sqrt pi
*/

using namespace std;
int main ()
{
  float pi_f = M_PI;
  double pi_d = M_PI;
  long double pi_l = M_PI;
  cout.precision(50);
  cout << "Const Pi:\t" << fixed << M_PI << endl;
  cout << "Float Pi:\t" << fixed << pi_f << endl;
  cout << "Double Pi:\t" << fixed << pi_d << endl;
  cout << "Long Pi:\t" << fixed << pi_l << endl;

  return 0;
}
