// =================================================================== //
// this is file Solver.h, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#ifndef Solver_H
#define Solver_H

#include <cstring>
#include <iomanip>
#include "clVector.h"
#include "clMatrix.h"
#include "Solver.h"
#include "colors.h"

class Solver
{
  public:
   Solver();
   virtual void solve(real eps, clMatrix & A, clVector & x, clVector & bb) = 0;
   virtual ~Solver();
   void printAndSaveInfo();

  protected:
   int n,m;            // matrix dimensions
   int iteration;      // number of iterations
   string reason;      // reason of solver stop
   string solverName;  // solver name
   string pcName;      // preconditioner name
   real initial;       // initial time 
   real final;         // final time
   real elapsed;       // difference between final and initial;
   real residual;      // residual norm
};


#endif


