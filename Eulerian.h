
#ifndef EULERIAN_H
#define EULERIAN_H

#include "clVector.h"
#include "clMatrix.h"
#include "Model3D.h"
#include <iostream>
#include <math.h>
#include <vector>

class Eulerian
{
 public:
  Eulerian(Model3D &_m);
  virtual ~Eulerian();

  clVector compute(real _dt);
  void stepSmooth();
  void stepSmoothTangent();
  void stepSmoothTangent2();
  void setBC();
  void setCentroid();
  void setSurface();
  int search(int node,real _XI, real _YI,real _ZI);

  clVector uSmooth,vSmooth,wSmooth;
  clVector uSmoothSurface,vSmoothSurface,wSmoothSurface;

 private:
  real dt;
  int numNodes,numVerts,numElems;
  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp;
  clMatrix *IEN;
  clVector *nonSurface,*surface,xSurface,ySurface,zSurface;
  clVector xAverage,yAverage,zAverage;

  typedef list<int> listElem;
  vector<listElem> neighbourVert;
  vector<listElem> surfaceViz;
};
#endif /* ifndef EULERIAN_H */

