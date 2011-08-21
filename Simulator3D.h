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
#include <iomanip>
#include "Model3D.h"
#include "FEMLinElement3D.h"
#include "FEMMiniElement3D.h"
#include "FEMQuadElement3D.h"
#include "Galerkin.h"
#include "SemiLagrangean.h"
#include "MeshSmooth.h"
#include "Solver.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "clVector.h"
#include "clMatrix.h"
#include "clDMatrix.h"
#include "colors.h"

class Simulator3D
{
 public:
  Simulator3D(); // construtor padrao
  Simulator3D( Model3D &_m ); // construtor 
  Simulator3D( Model3D &_m, Simulator3D &_s );  // copia
  Simulator3D( Model3D &_mNew, Model3D &_mOld, Simulator3D &_s );  // copia
  virtual ~Simulator3D(); // destrutor padrao

  void getModel3DAttrib(Model3D &_m);
  void allocateMemoryToAttrib();
  void init();
  void initFixedBubbleZ();
  void init2Bubbles();
  void assemble();
  void assembleK();
  void assembleM();
  void assembleC();
  void assembleSlip();
  void assembleNuCte();
  void assembleNuZ();
  void assembleNuC();
  void matMount();
  void matMountC();

  void step();
  void stepSL();
  void stepNoConvection();
  void stepLagrangian();
  void stepLagrangianZ();
  void stepALE();
  void stepALEVel();
  void stepSmooth();
  void setInterfaceVel();
  void setInterfaceVelNormal();

  void setRHS();
  void setCRHS();
  void setGravity(const char* _direction);
  void setGravityBoussinesq(const char* _direction);
  void setInterfaceGeo();
  void setInterfaceLevelSet();

  void coupled();
  void unCoupled();
  void unCoupledC();

  void setCoupledBC();
  void setUnCoupledBC();
  void setUnCoupledCBC();

  void setMuZ();
  void setMuC();
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
  void setDtDisk();
  void setDt();
  void setDt(real _dt);
  void setTime(real _time);
  real getDt();
  real getTime2();
  real* getTime();
  void setCfl(real _cfl);
  void setCflBubble(real _cfl);
  real getDtLagrangian();
  real getDtSemiLagrangian();
  real getDtSurfaceTension();
  real getDtGravity();
  real getCfl();
  void setIter(real _Iter);
  int getIter();
  void setC1(real _c1);
  void setC2(real _c2);
  void setC3(real _c3);
  void setC4(real _c4);
  real getC1();
  real getC2();
  real getC3();
  real getC4();
  void setMu(real _mu_in);
  void setMu(real _mu_in,real _mu_out);
  void setMuSmooth(real _mu_in,real _mu_out);
  real getMu_in();
  real getMu_out();
  void setRho(real _rho_in);
  void setRho(real _rho_in,real _rho_out);
  void setRrho_in(real _rho_in);
  void setRhoSmooth(real _rho_in,real _rho_out);
  void setHSmooth();
  real getRho_in();
  void setRho_out(real _rho_out);
  real getRho_out();
  void setUAnt(clVector &_uAnt);
  void setCSol(clVector &_cSol);
  void setUSol(clVector &_uSol);
  void setVSol(clVector &_vSol);
  void setWSol(clVector &_wSol);
  void updateIEN();
  clVector setCentroid(clVector &_vector);

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
  clVector* getFint();
  clVector* getGravity();
  real getGrav();
  clDMatrix* getKappa();
  clMatrix* getK();
  clMatrix* getM();
  clMatrix* getG();
  clMatrix* getD();
  clMatrix* getGx();
  clMatrix* getGy();
  clMatrix* getGz();
  clVector* getMu();
  clVector* getRho();
  clVector* getHSmooth();
  void operator=(Simulator3D &_s);
  void operator()(Model3D &_m);
  void operator()(Model3D &_m,Simulator3D &_s);
  int loadSolution( const char* _filename,int _iter );
  void applyLinearInterpolation(Model3D &_mOld);
  void getBubbleVelocity();
  void getBubbleVelocity(clVector _uVel,clVector _vVel,clVector _wVel);
  void setALEVelBC();
  void setLagrangianVelBC();

 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  int numVertsOld,numElemsOld,numNodesOld;
  real triEdge;
  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clVector *outflow,*surface;
  clMatrix *IEN;
  clDMatrix *curvature;
  SurfaceMesh *surfMesh;
  Mesh3D *mesh3d;
  clVector *heaviside,*interfaceDistance,*elemIdRegion;


  real Re,Sc,Fr,We,alpha,beta,dt,cfl,time;
  real c1,c2,c3,c4;
  real bubbleXVel,bubbleYVel,bubbleZVel;
  real sigma,g,rho_in,rho_out,mu_in,mu_out;
  real sigma_0,g_0,rho_0,mu_0;
  real sigmaAdimen,gAdimen,rho_inAdimen,rho_outAdimen,mu_inAdimen,mu_outAdimen;
  real bubbleVel;
  int iter;
  clMatrix K,Kc,Mrho,M,Mc,G,D,A;
  clMatrix mat,matc;
  clMatrix gx,gy,gz;
  clMatrix ATilde,AcTilde,GTilde,DTilde,ETilde,E;
  clDMatrix MrhoLumped,McLumped;
  clDMatrix invA,invC,invMrhoLumped,invMcLumped;
  clDMatrix MLumped,invMLumped;
  clDMatrix kappa;
  clVector va,vcc,b;
  clVector convUVW,convC;
  clVector uTilde,cTilde,pTilde,b1,b1c,b2,ip,ipc;
  clVector velU,velV,velW,uSol,vSol,wSol,pSol,cSol;
  clVector uSL,vSL,wSL,cSL,uALE,vALE,wALE;
  clVector uSmooth,vSmooth,wSmooth,uSmoothCoord,vSmoothCoord,wSmoothCoord;
  clVector uSmoothSurface,vSmoothSurface,wSmoothSurface;
  clVector fint,gravity,hSmooth,mu,rho;

  clDMatrix kappaOld;
  clVector uALEOld,vALEOld,wALEOld;
  clVector uSolOld,vSolOld,wSolOld,pSolOld,cSolOld;
  clVector fintOld,gravityOld,Fold,muOld,rhoOld,hSmoothOld;

  Solver *solverV,*solverP,*solverC;
};

#endif
