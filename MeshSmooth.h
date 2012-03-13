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
  MeshSmooth(Model3D &_m,real _dt);
  virtual ~MeshSmooth();

  clVector compute();
  void stepSmooth();
  void stepSmoothFujiwara();
  void stepSurfaceSmoothFujiwara();
  void stepSmooth(clVector &_uVel,clVector &_vVel,clVector &_wVel);
  void stepSmoothSurface();
  void stepSmoothSurface2();
  void stepSmoothSurface(clVector &_uVel,clVector &_vVel,clVector &_wVel);
  void setBC();
  void setCentroid();
  void setSurface();

  clVector* getUSmooth(); 
  clVector* getVSmooth(); 
  clVector* getWSmooth(); 
  clVector* getUSmoothSurface(); 
  clVector* getVSmoothSurface(); 
  clVector* getWSmoothSurface(); 
  clVector uSmooth,vSmooth,wSmooth;
  clVector uSmoothSurface,vSmoothSurface,wSmoothSurface;

 private:
  Model3D *m;
  real dt;
  int numNodes,numVerts,numElems;
  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp;
  clVector *heaviside;
  clMatrix *IEN;
  clVector *nonSurface,*surface,xSurface,ySurface,zSurface;
  SurfaceMesh *surfMesh;

  vector< list<int> > *neighbourVert,*neighbourPoint;
  list<int> *inVert;
};
#endif /* ifndef MESHSMOOTH_H */

