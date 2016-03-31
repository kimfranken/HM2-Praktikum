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

  fstream file;
  fstream file1;
  file.open("V_1.txt", ios::out);
  file1.open("V_1_scientific.txt", ios::out);

  long double V_1[41] = {M_2_PI, M_1_PI};

  file << setw( 2 ) << 0 << "\t" << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_1[0]  << endl;
  file << setw( 2 ) << 1 << "\t" << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_1[1]  << endl;

  file1 << setw( 2 ) << 0 << "\t" << scientific << setw( 50 ) << setprecision(45) << V_1[0]  << endl;
  file1 << setw( 2 ) << 1 << "\t" << scientific << setw( 50 ) << setprecision(45) << V_1[1]  << endl;

  for (int n = 2; n <= 40; n++) {
    V_1[n] = M_1_PI - V_1[n-2]*n*(n-1)/(M_PI*M_PI);
    file << setw( 2 ) << n << "\t"  << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_1[n]  << endl;
    file1 << setw( 2 ) << n << "\t"  << scientific << setw( 50 ) << setprecision(45) << V_1[n]  << endl;

  }

  file.close();
  /*
  fstream file_float;
  file_float.open("V_1_float.txt", ios::out);

  float V_1_float[41] = {M_2_PI, M_1_PI};

  file_float << setw( 2 ) << 0 << "\t" << fixed << setw( 40 ) << setprecision(std::numeric_limits<float>::digits10 + 1) << V_1_float[0]  << endl;
  file_float << setw( 2 ) << 1 << "\t" << fixed << setw( 40 ) << setprecision(std::numeric_limits<float>::digits10 + 1) << V_1_float[1]  << endl;

  for (int n = 2; n <= 40; n++) {
    V_1_float[n] = 1/M_PI - V_1[n-2]*n*(n-1)/(M_PI*M_PI);
    file_float << setw( 2 ) << n << "\t"  << fixed << setw( 40 ) << setprecision(std::numeric_limits<float>::digits10 + 1) << V_1_float[n]  << endl;
  }

  file_float.close();
  */

  fstream file_2;
  file_2.open("V_2.txt", ios::out);

  long double V_2[61];
  V_2[60] = 1./61.;
  V_2[59] = 1./60.;

  file_2 << setw( 2 ) << 60 << "\t" << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_2[60]  << endl;
  file_2 << setw( 2 ) << 59 << "\t" << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_2[59]  << endl;

  for (int n = 60; n >= 2; n--) {
    V_2[n-2] = M_PI*(1 - M_PI*V_2[n])/(n*(n-1));
    file_2 << setw( 2 ) << n-2 << "\t" << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_2[n-2]  << endl;
  }

  file_2.close();

  fstream file_3;
  file_3.open("V_2_V_1.txt", ios::out);

  long double V_3[41];

  for (int n = 0; n <= 40; n++) {
    V_3[n] = abs(V_2[n] - V_1[n]);
    file_3 << setw( 2 ) << n << "\t"  << setw( 40 ) << setprecision(std::numeric_limits<long double>::digits10 + 1) << V_3[n]  << endl;
  }

  file_3.close();

  return 0;
}
