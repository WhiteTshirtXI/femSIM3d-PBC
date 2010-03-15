
#ifndef SEMILAGRANGEAN_H
#define SEMILAGRANGEAN_H

#include "clVector.h"
#include "clMatrix.h"
#include "Model3D.h"
#include <iostream>
#include <math.h>
#include <vector>

class SemiLagrangean
{
 public:
  SemiLagrangean(Model3D &_m,clVector &_uSol,clVector &_vSol,
	             clVector &_wSol,clVector &_cSol);

  clVector compute(real dt);
  clVector computeFreeSurface(real dt);
  void getDepartElem(real dt);
  void getDepartElem2(real dt);
  void jumpToElem(int destElem,int iiVert,real R2X,real R2Y,real R2Z);
  bool testElement(int mele,int ii,real xP,real yP,real zP, 
	               real *l1,real *l2,real *l3,real *l4);
  void computeIntercept(int ii,real R2X,real R2Y,real R2Z,int ib1,int ib2,
	                    int ib3,real *l1,real *l2,real *l3);
  void setCentroid();
  void setBC();

 private:
  int numVerts,numNodes,numElems;
  real dt; 
  typedef list<int> listElem;
  vector<listElem> neighbourElem;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector uParticle,vParticle,wParticle,cParticle;
  clVector *X,*Y,*Z;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clMatrix *IEN;
  clVector uSol,vSol,wSol,cSol;
  clMatrix convLin;
  clMatrix *oFace;

};
#endif /* ifndef SEMILAGRANGEAN_H */

