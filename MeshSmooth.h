// =================================================================== // 
// this is file MeshSmooth.h, created at 20-Ago-2009                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //
//
#ifndef MESHSMOOTH_H
#define MESHSMOOTH_H

#include "clVector.h"
#include "clMatrix.h"
#include "Model3D.h"
#include <iostream>
#include <math.h>
#include <vector>

class MeshSmooth
{
 public:
  MeshSmooth(Model3D &_m);
  virtual ~MeshSmooth();

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
  Model3D *m;
  real dt;
  int numNodes,numVerts,numElems;
  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp;
  clMatrix *IEN;
  clVector *nonSurface,*surface,xSurface,ySurface,zSurface;
  clVector xAverage,yAverage,zAverage;

  vector< list<int> > *neighbourVert,*surfaceViz;
  list<int> *inVert;
};
#endif /* ifndef MESHSMOOTH_H */

