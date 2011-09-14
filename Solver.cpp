// =================================================================== //
// this is file Solver.cpp, created at 10-Jun-2007               //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include "Solver.h"

Solver::Solver()
{
 n = 0;
 m = 0;      
 iteration = 0;  
 reason = "no reason";     
 solverName = "no solver"; 
 pcName = "no preconditioner";     
 initial = 0 ;    
 final = 0;      
 elapsed = 0;   
 residual = 0;   
}

Solver::~Solver(){}

void Solver::printAndSaveInfo()
{
 cout << " Preconditioner: " << color(none,yellow,black)
      << pcName << resetColor() << endl;
 
 if( strcmp( solverName.c_str(),"preonly") == 0 )
   cout << " Solver: " << color(none,yellow,black)
        << "mumps" << resetColor() << endl;
 else
   cout << " Solver: " << color(none,yellow,black)
        << solverName << resetColor() << endl;

 cout << " Size: " << color(none,yellow,black) 
      << n << resetColor() << endl;
 cout << " Iterations: " << color(none,yellow,black) 
      << iteration << resetColor() << endl;
 cout << " Reason: " << color(none,yellow,black) 
      << reason << resetColor() << endl;

 ostringstream sos;
 sos.precision(2);
 sos << fixed;
 sos << " Resolution Time: " << color(none,red,black) 
  << elapsed << resetColor() << "s" << endl;
 cout << sos.str();

 cout << " Residual Norm: " << color(none,yellow,black) 
      << residual << resetColor() << endl;

 if( residual > 10E-03 )
 {
  cerr << endl;
  cerr << endl;
  cerr << color(blink,red,black);
  cerr << "              *----------------------------------------*" << endl;
  cerr << "                The solver broke-up before convergence  " << endl;
  cerr << "              *----------------------------------------*" << endl;
  cerr << resetColor();
  cerr << endl;
  cerr << endl;
  exit(1);
 }


 /* saving file */
 string fileAux = "./dat/" + solverName + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
  testFile.close();
 else
 {
  file << "#Solution using " << solverName << " and " << pcName << endl;
  file << endl;
  file << "#size" << setw(13) 
       << "iterations" << setw(21)
	   << "reason" << setw(8)
	   << "time" << setw(10)
	   << "residual" << endl;
 }


 file << n << setw(12) 
      << iteration << setw(24) 
	  << reason << setw(7);

 ostringstream so;
 so.precision(2);
 so << fixed;
 so << elapsed;
 file << so.str() << setw(12);

 file << residual << endl;

 file.close();

}
