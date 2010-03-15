/**=========================================================================
 *\file Id:Simulator3D.h, GESAR
 *\author Gustavo Rabello dos Anjos
 *\date   23-jun-2007
 *\comment Classe Simulator 3D para o modulo HYDRO
 *\references
 *\version 1.0
 *\updates
 	version    date         author             comment
	1.0        23/08/2007   Gustavo            incio de implementacao
 =========================================================================**/

#ifndef Simulator3D_H
#define Simulator3D_H

#include <iostream>
#include <fstream>
#include "Model3D.h"
#include "FEMMiniElement3D.h"
#include "FEMLinElement3D.h"
#include "Galerkin.h"
#include "SemiLagrangean.h"
#include "Eulerian.h"
#include "Solver.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"
#include "Interface3D.h"

class Simulator3D
{
 public:
  Simulator3D( Model3D &_m ); // construtor padrao
  virtual ~Simulator3D(); // destrutor padrao

  void init();
  void assembleK();
  void assembleM();
  void assemble();
  void assembleSlip();
  void assembleNuZ();
  void assembleNuCte();
  void matMount();
  void matMountC();

  void step();
  void stepSL();
  void stepLagrangian();
  void stepLagrangianZ();
  void stepALE();
  void stepSmooth();

  void setRHS();
  void setCRHS();
  void setGravity();
  void setGravityBoussinesq();
  void setInterface();
  void setInterfaceGeo();

  void coupled();
  void unCoupled();
  void unCoupledC();

  void setCoupledBC();
  void setUnCoupledBC();
  void setUnCoupledCBC();

  void setNuZ();
  void convergenceCriteria( real value );
  void setRe(real _Re);
  real getRe();
  void setSc(real _Sc);
  real getSc();
  void setFr(real _Fr);
  real getFr();
  void setWe(real _We);
  real getWe();
  void setAlpha(real _alpha);
  real getAlpha();
  void setBeta(real _beta);
  real getBeta();
  void setDt(real _dt);
  real getDt();
  real* getTime();
  void setCfl(real _cfl);
  void setCflDisk(real _cfl);
  void setCflBubble(real _cfl);
  real getCfl();
  void setUAnt(clVector &_uAnt);
  void setCSol(clVector &_cSol);

  void setSolverVelocity(Solver *s);
  void setSolverPressure(Solver *s);
  void setSolverConcentration(Solver *s);

  clVector getUSol()const;
  clVector* getPointerUSol();
  clVector getVSol()const;
  clVector* getPointerVSol();
  clVector getWSol()const;
  clVector* getPointerWSol();
  clVector getPSol()const;
  clVector* getPointerPSol();
  clVector getCSol()const;
  clVector* getPointerCSol();
  clVector getUAnt()const;
  clVector* getPointerUAnt();
  clVector* getPointerCAnt();
  clVector* getPointerDistance();
  clDMatrix* getPointerKappa();
  clMatrix getK()const;
  clMatrix* getPointerK();
  clMatrix getM()const;
  clMatrix* getPointerM();
  clMatrix getG()const;
  clMatrix* getPointerG();
  clMatrix getD()const;
  clMatrix* getPointerD();
  clMatrix getGx()const;
  clMatrix* getPointerGx();
  clMatrix getGy()const;
  clMatrix* getPointerGy();
  clMatrix getGz()const;
  clMatrix* getPointerGz();

  Solver *solverV,*solverP,*solverC;
  Model3D *m;
 private:
  int numVerts,numElems,numNodes;
  int numGLEU,numGLEP,numGLEC;
  real Re,Sc,Fr,We,alpha,beta,dt,cfl,time,sigma;
  real c1,c2,c3;

  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clVector *outflow;
  clMatrix *IEN;

  clMatrix K,Kc,M,Mc,G,D;
  clMatrix mat,matc;

  clDMatrix MLumped,McLumped;
  clVector velU,velV,velW,uSol,vSol,wSol,pSol,cSol;
  clVector uSL,vSL,wSL,cSL;
  clVector uALE,vALE,wALE;
  clVector uSmooth,vSmooth,wSmooth;

  clMatrix gx,gy,gz;
  clVector uAnt,cAnt;
  clVector nu;

  clVector va,vcc;
  clVector convUVW,convC;

  clMatrix A;
  clVector b;
  
  clMatrix ATilde,AcTilde,GTilde,DTilde,ETilde,E;
  clDMatrix invA,invC,invMLumped,invMcLumped;
  clVector uTilde,cTilde,pTilde,b1,b1c,b2,ip,ipc;
 
  clVector distance;
  clDMatrix kappa;
  clVector fint;
  clVector Fold;

};

#endif
