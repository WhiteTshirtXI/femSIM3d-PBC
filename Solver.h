// =================================================================== //
// this is file Solver.h, created at 10-Jun-2007               //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#ifndef Solver_H
#define Solver_H

#include "clVector.h"
#include "clMatrix.h"
#include "Solver.h"

class Solver
{
  public:
   virtual void solve(real eps, clMatrix & A, clVector & x, clVector & bb) = 0;
   virtual ~Solver();
};


#endif


