
#ifndef SEMILAGRANGEAN_H
#define SEMILAGRANGEAN_H

#include "interpolations.h"
#include "clVector.h"
#include "clMatrix.h"
#include "Model3D.h"
#include "TElement.h"
#include <iostream>
#include <math.h>
#include <vector>

class SemiLagrangean
{
 public:
  SemiLagrangean(Model3D &_m);
  SemiLagrangean(Model3D &_m,clVector &_uSol,clVector &_vSol,
	             clVector &_wSol,
				 clVector &_velU,clVector &_velV,clVector &_velW,
				 clVector &_cSol);

  void compute(real dt);
  void computeFreeSurface(real dt);
  void getDepartElem(real dt);
  void getDepartElemQuad(real dt);
  void getDepartElem2(real dt);
  void jumpToElem(int destElem,int iiVert,real R2X,real R2Y,real R2Z);
  void jumpToElemQuad(int destElem,int iiVert,real R2X,real R2Y,real R2Z);
  void jumpToElem2(int destElem,int iiVert,real R2X,real R2Y,real R2Z);
  bool testElement(int mele,int ii,real xP,real yP,real zP, 
	               real *l1,real *l2,real *l3,real *l4);
  void computeIntercept(int ii,real R2X,real R2Y,real R2Z,int ib1,int ib2,
	                    int ib3,real *l1,real *l2,real *l3);
  void setCentroid();
  void setQuad();
  void setBC();
  clMatrix* getInterpLin();
  clVector* getUSL();
  clVector* getVSL();
  clVector* getWSL();
  clVector* getCSL();

 private:
  Model3D *m;
  int numVerts,numNodes,numElems;
  real dt; 
  vector< list<int> > *neighbourElem;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector uParticle,vParticle,wParticle,cParticle;
  clVector *X,*Y,*Z;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clMatrix *IEN;
  clVector uSol,vSol,wSol,cSol;
  clVector velU,velV,velW;
  clMatrix convLin,convQuad;
  clMatrix *oFace;
};
#endif /* ifndef SEMILAGRANGEAN_H */

