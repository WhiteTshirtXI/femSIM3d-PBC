// =================================================================== //
// this is file InterFace3D, created at 26-Mar-2009                    //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#ifndef INTERFACE3D
#define INTERFACE3D

#include <iostream>
#include <fstream>
#include <vector>
#include "Model3D.h"
#include "Simulator3D.h"
#include "Solver.h"
#include "PCGSolver.h"
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"

class Interface3D
 {
  public:
   Interface3D(Model3D &_m); // construtor padrao
   virtual ~Interface3D(); // destrutor padrao

  clVector curvature1();
  clVector curvature2();
  clVector smoothing(clMatrix &_AcTilde,clVector &_b1cTilde);
  clVector computeKappa1();
  clVector computeKappa2();
  void computeKappaTest();
  void computeKappa3();
  void plotKappa(clVector &_kappaAux);
  void setCloser();
  clDMatrix setKappaSurface(clVector &_kappaAux);
  clDMatrix setKappaSurface(clVector &_kappaNx,
	                        clVector &_kappaNy,
							clVector &_kappaNz);
  void setSolverSmooth(Solver *s);
  clVector getCloser();
  clVector crossProd(real x1,real y1,real z1,real x2,real y2,real z2);

  Solver *solverC;
 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  real xCenter,yCenter,zCenter,bubbleRadius;
  clVector *surface,xSurface,ySurface,zSurface;
  clVector closer,xCloser,yCloser,zCloser,closerViz;
  clDMatrix kappa;
  clVector *X,*Y,*Z;
  clMatrix *IEN;
  clVector *cc;
  clVector distance;

  vector< list<int> > *neighbourElem,*neighbourVert,*neighbourFace;
  vector< list<int> > *neighbourFaceVert,*surfaceViz,*faceIEN,*elemSurface;
};

#endif
