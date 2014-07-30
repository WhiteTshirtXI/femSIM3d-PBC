
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

  void compute(double dt);
  void computeFreeSurface(double dt);
  void getDepartElem(double dt);
  void getDepartElemQuad(double dt);
  void getDepartElem2(double dt);
  void jumpToElem(int destElem,int iiVert,double R2X,double R2Y,double R2Z);
  void jumpToElemQuad(int destElem,int iiVert,double R2X,double R2Y,double R2Z);
  void jumpToElem2(int destElem,int iiVert,double R2X,double R2Y,double R2Z);
  bool testElement(int mele,int ii,double xP,double yP,double zP, 
	               double *l1,double *l2,double *l3,double *l4);
  void computeIntercept(int ii,double R2X,double R2Y,double R2Z,int ib1,int ib2,
	                    int ib3,double *l1,double *l2,double *l3);
  void setCentroid();
  void setQuad();
  void setBC();
  clMatrix* getInterpLin();
 
  // PBC
  void computePBCFix(double dt);
  void getDepartElemPBCFix(double dt);
  void getDepartElemQuadPBCFix(double dt);
  clVector* getUSL();
  clVector* getVSL();
  clVector* getWSL();
  clVector* getCSL();

 private:
  Model3D *m;
  int numVerts,numNodes,numElems;
  double dt; 
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

