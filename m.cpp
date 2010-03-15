// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"

int main(int argc, char **argv)
{

 int n = 10;
 clMatrix A;
 A.Dim(n,n);
 A.Set(0,0,5);
 A.Set(1,1,5);

 for( int i=0;i<n;i++ )
  for( int j=0;j<n;j++ )
   A.Set(i,j,i*j);



 cout << "A = " << sizeof(A) << " bytes" << endl;


 return 0;
}


