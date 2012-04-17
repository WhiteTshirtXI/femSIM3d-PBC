// =================================================================== //
// this is file Helmholtz3D.h, created at 21-Sep-2011                    //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#ifndef HELMHOLTZ3D_H
#define HELMHOLTZ3D_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "Model3D.h"
#include "FEMLinElement3D.h"
#include "Solver.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"
#include "colors.h"
#include "geometry.h"
#include "interpolations.h"
#include "searchInterp3D.h"
#include <algorithm>

class Helmholtz3D
{
 public:
  Helmholtz3D(); // construtor padrao
  Helmholtz3D( Model3D &_m ); // construtor 
  Helmholtz3D( Model3D &_m,Helmholtz3D &_d ); // construtor 
  //virtual ~Helmholtz3D(); // destrutor padrao

  void getModel3DAttrib(Model3D &_m);
  void allocateMemoryToAttrib();
  void init();
  void initMicro();
  void initSquareChannel();
  void initRisingBubble();
  void initSessile();
  void init2Bubbles();
  void assemble();
  void matMountC();

  void setCRHS();
  void setBC();
  void setUnCoupledCBC();
  void unCoupledC();
  void setModel3DEdgeSize();

  void setSolver(Solver *s);

  void saveVTK( const char* _dir,const char* _filename, int _iter );
  void saveChordalEdge( const char* _dir,const char* _filename, int _iter );

  void setk(real _k);
  real getk();
  void setCloserWall();
  bool isInsideLiquidFilm(int _node, real _thickness);

  clVector* getCSol();

 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  int numVertsOld,numElemsOld,numNodesOld;
  real k;
  vector<real> triEdge;
  clVector *X,*Y,*Z;
  clMatrix *IEN;
  clVector *heaviside,*interfaceDistance;
  clVector *edgeSize;
  clVector cc,idbcc;
  SurfaceMesh *surfMesh;
  Mesh3D *mesh3d;

  list<int> *boundaryVert;
  vector< list<int> > *neighbourPoint;

  clMatrix Kc,Mc;
  clMatrix matc,AcTilde;
  clDMatrix invC,invMcLumped,McLumped;
  clVector vcc,convC;
  clVector cTilde,b1c,ipc;
  clVector cSol,cSolOld;
  clVector closerWall;

  Solver *solver;
};

#endif

