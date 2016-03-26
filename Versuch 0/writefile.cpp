#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;
int main ()
{
  int i=0;
  fstream file;
  file.open("OutputFile_cpp.txt", ios::out);
  for (double x=0.0; x<6.29; x=x+0.1)
  {
    ++i;
    file << i << "\t" << fixed << setw( 11 )<< setprecision( 5 ) << x
    << "\t" << scientific << sin(x) << endl;
  }
  file.close();
  return 0;
}
