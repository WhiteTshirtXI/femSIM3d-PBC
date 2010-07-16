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
#include "MeshSmooth.h"
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
  Simulator3D( Model3D &_m, Simulator3D &_s );  // copia
  Simulator3D( Model3D &_mNew, Model3D &_mOld, Simulator3D &_s );  // copia
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
  void stepALE2();
  void stepALE3();
  void stepSmooth();
  void setInterfaceVel();
  void  stepMesh();

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

  void setHsmooth();
  void setNu(real nu0, real nu1);
  void setRho(real rho0, real rho1);
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
  void setSigma(real _sigma);
  real getSigma();
  void setDt(real _dt);
  void setTime(real _time);
  real getDt();
  real getTime2();
  real* getTime();
  void setCfl(real _cfl);
  void setCflDisk(real _cfl);
  void setCflBubble(real _cfl);
  real getCfl();
  void setUAnt(clVector &_uAnt);
  void setCSol(clVector &_cSol);
  void setUSol(clVector &_uSol);
  void setVSol(clVector &_vSol);
  void setWSol(clVector &_wSol);
  void updateIEN();
  void setCentroid();

  void setSolverVelocity(Solver *s);
  void setSolverPressure(Solver *s);
  void setSolverConcentration(Solver *s);

  clVector* getUSol();
  clVector* getVSol();
  clVector* getWSol();
  clVector* getPSol();
  clVector* getCSol();
  clVector* getUALE();
  clVector* getVALE();
  clVector* getWALE();
  clVector* getUAnt();
  clVector* getCAnt();
  clVector* getDistance();
  clVector* getFint();
  clDMatrix* getKappa();
  clMatrix* getK();
  clMatrix* getM();
  clMatrix* getG();
  clMatrix* getD();
  clMatrix* getGx();
  clMatrix* getGy();
  clMatrix* getGz();
  void operator=(Simulator3D &_s);
  int loadIteration();
  int loadIteration( const char* _dir,const char* _filename,int _iter );
  void loadSolution( const char* _dir,const char* _filename,int _iter );
  void applyLinearInterpolation(Model3D &_mOld);
  real getBubbleVelocity();
  void setALEVelBC();

  Solver *solverV,*solverP,*solverC;

  clVector cSol;
 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  int numGLEU,numGLEP,numGLEC;
  real Re,Sc,Fr,We,alpha,beta,dt,cfl,time,sigma;
  real c1,c2,c3;

  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clVector *outflow,*surface;
  clMatrix *IEN;

  clMatrix K,Kc,M,Mc,G,D;
  clMatrix mat,matc;

  clDMatrix MLumped,McLumped;
  clVector velU,velV,velW,uSol,vSol,wSol,pSol;
  //clVector velU,velV,velW,uSol,vSol,wSol,pSol,cSol;
  clVector uSL,vSL,wSL,cSL;
  clVector uALE,vALE,wALE;
  clVector uSmooth,vSmooth,wSmooth;
  clVector uALEOld,vALEOld,wALEOld;
  clVector uSolOld,vSolOld,wSolOld,pSolOld;

  clMatrix gx,gy,gz;
  clVector uAnt,cAnt;

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
  clVector Hsmooth,nu,rho;
  clVector Fold;

};

#endif
