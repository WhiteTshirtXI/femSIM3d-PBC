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
   //virtual ~Interface2D(); // destrutor padrao

  clVector curvature1();
  clVector curvature2();
  clVector smoothing(clMatrix &_AcTilde,clVector &_b1cTilde);
  clVector computeKappa();
  clVector computeKappa2();
  clVector computeKappa3();
  clVector computeKappa4();
  void plotKappa(clVector &_kappaAux);
  clVector dsearchn(clVector _X,clVector _Y,clVector _Z,
	                clVector &_XI,clVector &_YI,clVector &_ZI);
  void setCloser();
  clDMatrix setKappaSurface(clVector &_kappaAux);
  void setSolverSmooth(Solver *s);
  void saveSurfaceVTK();
  clVector getCloser();


  Solver *solverC;
 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  real xCenter,yCenter,zCenter,bubbleRadius;
  clVector distance;
  clVector *surface,xSurface,ySurface,zSurface;
  clVector closer,xCloser,yCloser,zCloser,closerViz;
  clDMatrix kappa;
  clVector *X,*Y,*Z;
  clMatrix *IEN;
  clVector *cc;

  typedef list<int> listElem;
  vector<listElem> neighbourElem;
  vector<listElem> neighbourVert;
  vector<listElem> neighbourFace;
  vector<listElem> elemSurface;
  vector<listElem> neighbourFaceVert;
  vector<listElem> surfaceViz;
  vector<listElem> xSurfaceViz;
  vector<listElem> ySurfaceViz;
  vector<listElem> zSurfaceViz;
};

#endif
