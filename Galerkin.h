
#ifndef GALERKIN_H
#define GALERKIN_H

#include "Model3D.h"
#include "clVector.h"
#include "clMatrix.h"

class Galerkin
{
 public:
  Galerkin(Model3D &_m,clVector &_uSol,
	                   clVector &_vSol,
					   clVector &_wSol,
					   clVector &_cSol,
					   clMatrix &_gx,
					   clMatrix &_gy,
					   clMatrix &_gz);
  clVector compute(real dt);

  int numVerts,numNodes,numElems;
  real dt;
  clVector uc,vc,wc,pc;
  clVector X,Y,Z;
  clVector idbcu,idbcv,idbcw,idbcp;
  clVector *outflow;
  clMatrix IEN;
  clVector uSol,vSol,wSol,cSol;
  clMatrix gx,gy,gz;

};

#endif /* ifndef GALERKIN_H */

