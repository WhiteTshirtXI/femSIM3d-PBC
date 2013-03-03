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
#include <algorithm>
#include "Model3D.h"
#include "TElement.h"
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
#include "geometry.h"
#include "searchInterp3D.h"

class Simulator3D
{
 public:
  Simulator3D(); // construtor padrao
  Simulator3D( Model3D &_m ); // construtor 
  Simulator3D( const Simulator3D &_sRight ); // construtor 
  Simulator3D( Model3D &_m, Simulator3D &_s );  // copia
  virtual ~Simulator3D(); // destrutor padrao

  void getModel3DAttrib(Model3D &_m);
  void allocateMemoryToAttrib();
  void init();
  void initDiskBaseState( const char* _dir,const char* _filename );
  void initAnnular();
  void initChannel();
  void initChannelSquare();
  void initHeatTransfer();
  void initFixedBubbleZ();
  void init2AxiBubbles();
  void init2Bubbles();
  void assemble();
  void assembleHeatTransfer();
  void assembleK();
  void assembleM();
  void assembleC();
  void assembleSlip();
  void assembleNuCte();
  void assembleNuZ(const char* _name);
  void assembleNuC();
  void matMount();
  void matMountC();

  void step();
  void stepSL();
  void stepNoConvection();
  void stepImposedPeriodicField(const char* _name, real T);
  void stepTimeHalf(const char* _name, real T);
  void copyALEtoSol();
  void stepLagrangian();
  void stepLagrangianZ();
  void stepALE();
  void movePoints();
  void movePoints(clVector *_uVel,
                  clVector *_vVel,
				  clVector *_wVel);
  void stepSmooth();
  void setInterfaceVelocity();

  void setRHS();
  void setCRHS();
  void setGravity(const char* _direction);
  void setGravityBoussinesq(const char* _direction);
  void setInterfaceGeo();
  void setInterfaceLevelSet();

  void coupled();
  void unCoupled();
  void unCoupledC();

  void saveOldData();

  void setCoupledBC();
  void setUnCoupledBC();
  void setUnCoupledCBC();

  void setNuZ(const char* _filename);
  void setNuC();
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
  void setDtLagrangian();
  void setDtLagrangianExtream();
  void setDtLagrangianNorberto();
  void setDtSemiLagrangian();
  void setDtGravity();
  void setDtSurfaceTension();
  void setDtEulerian();
  void setDtALESinglePhase();
  void setDtALETwoPhase();
  void setDt();
  void setDt(real _dt);
  void setTime(real _time);
  real getTime();
  void setCfl(real _cfl);
  void setCflBubble(real _cfl);
  real getDt();
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
  void setD1(real _d1);
  void setD2(real _d2);
  real getC1();
  real getC2();
  real getC3();
  real getD1();
  real getD2();
  void setURef(real _uRef);
  real getURef();
  void setVRef(real _vRef);
  real getVRef();
  void setWRef(real _wRef);
  real getWRef();
  void setXRef(real _XRef);
  real getXRef();
  void setYRef(real _YRef);
  real getYRef();
  void setZRef(real _ZRef);
  real getZRef();
  void setMu(real _mu_in);
  void setMu(real _mu_in,real _mu_out);
  void setMuSmooth(real _mu_in,real _mu_out);
  real getMu_in();
  real getMu_out();
  real getMu_inAdimen();
  real getMu_outAdimen();
  void setRho(real _rho_in);
  void setRho(real _rho_in,real _rho_out);
  void setRrho_in(real _rho_in);
  void setRhoSmooth(real _rho_in,real _rho_out);
  void setCp(real _cp_in);
  void setCp(real _cp_in,real _cp_out);
  void setCpSmooth(real _cp_in,real _cp_out);
  real getCp_in();
  real getCp_out();
  real getCp_inAdimen();
  real getCp_outAdimen();
  void setKt(real _kt_in);
  void setKt(real _kt_in,real _kt_out);
  void setKtSmooth(real _kt_in,real _kt_out);
  real getKt_in();
  real getKt_out();
  real getKt_inAdimen();
  real getKt_outAdimen();
  void setHSmooth();
  void setHeatFlux();
  real getRho_in();
  void setRho_out(real _rho_out);
  real getRho_out();
  real getRho_inAdimen();
  real getRho_outAdimen();
  void setUAnt(clVector &_uAnt);
  void setCSol(clVector &_cSol);
  void setUSol(clVector &_uSol);
  void setVSol(clVector &_vSol);
  void setWSol(clVector &_wSol);
  void updateIEN();
  clVector setCentroid(clVector &_vector);
  void setUSol(real _vel);
  void setVSol(real _vel);
  void setWSol(real _vel);

  void setSolverVelocity(Solver *s);
  void setSolverPressure(Solver *s);
  void setSolverConcentration(Solver *s);

  clVector* getUSol();
  clVector* getUSolOld();
  clVector* getVSol();
  clVector* getVSolOld();
  clVector* getWSol();
  clVector* getWSolOld();
  clVector* getPSol();
  clVector* getPSolOld();
  clVector* getCSol();
  clVector* getCSolOld();
  clVector* getUALE();
  clVector* getVALE();
  clVector* getUALEOld();
  clVector* getWALE();
  clVector* getVALEOld();
  clVector* getUAnt();
  clVector* getWALEOld();
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
  clVector* getCp();
  clVector* getKt();
  clVector* getHSmooth();
  clVector* getHeatFlux();
  void operator=(Simulator3D &_s);
  void operator()(Model3D &_m);
  int loadSolution( const char* _dir,const char* _filename, int _iter );
  void applyLinearInterpolation(Model3D &_mOld);
  void setCentroidVelPos();
  void setCentroidVelPosInterface();
  void setALEBC();
  void setAnnularALEBC();
  void setLagrangianVelBC();

  vector<real> getCentroidVelX();
  real getCentroidVelXAverage();
  vector<real> getCentroidVelY();
  real getCentroidVelYAverage();
  vector<real> getCentroidVelZ();
  real getCentroidVelZAverage();
  real getCentroidVelXMax();
  real getCentroidVelYMax();
  real getCentroidVelZMax();
  real getCentroidVelXMin();
  real getCentroidVelYMin();
  real getCentroidVelZMin();

  void setCentroidVelX(vector<real> _centroidVelX);
  void setCentroidVelY(vector<real> _centroidVelY);
  void setCentroidVelZ(vector<real> _centroidVelZ);

  vector<real> getCentroidPosX();
  real getCentroidPosXAverage();
  vector<real> getCentroidPosY();
  real getCentroidPosYAverage();
  vector<real> getCentroidPosZ();
  real getCentroidPosZAverage();
  void setCentroidPosX(vector<real> _centroidPosX);
  void setCentroidPosY(vector<real> _centroidPosY);
  void setCentroidPosZ(vector<real> _centroidPosZ);
  void setSurfaceTSat();
  void setMassTransfer();

 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  int numVertsOld,numElemsOld,numNodesOld;
  vector<real> triEdge;
  clVector *X,*Y,*Z;
  clVector *uc,*vc,*wc,*pc,*cc;
  clVector *idbcu,*idbcv,*idbcw,*idbcp,*idbcc;
  clVector *outflow,*surface;
  clMatrix *IEN;
  clDMatrix *curvature;
  SurfaceMesh *surfMesh;
  Mesh3D *mesh3d;
  clVector *heaviside,*interfaceDistance,*elemIdRegion;
  list<int> *boundaryVert;


  real Re,Sc,Fr,We,alpha,beta,cfl,time;
  real dt,dtLagrangian,dtSemiLagrangian,dtSurfaceTension,dtGravity;
  real c1,c2,c3,d1,d2;
  real bubbleXVel,bubbleYVel,bubbleZVel;
  real sigma,g,rho_in,rho_out,mu_in,mu_out;
  real sigma_0,g_0,rho_0,mu_0;
  real cp_in,cp_out,cp_0,cp_inAdimen,cp_outAdimen;
  real kt_in,kt_out,kt_0,kt_inAdimen,kt_outAdimen;
  real sigmaAdimen,gAdimen,rho_inAdimen,rho_outAdimen,mu_inAdimen,mu_outAdimen;
  real bubbleVel;
  int iter;
  vector<real> centroidVelX,centroidVelY,centroidVelZ;
  vector<real> centroidVelXOld,centroidVelYOld,centroidVelZOld;
  vector<real> centroidPosX,centroidPosY,centroidPosZ;
  vector<real> centroidPosXOld,centroidPosYOld,centroidPosZOld;

  // moving referential
  real uRef,vRef,wRef;
  real xRef,yRef,zRef;

  clMatrix K,Kc,Mrho,M,Mc,G,Gc,D,A;
  clMatrix mat,matc;
  clMatrix gx,gy,gz;
  clMatrix ATilde,AcTilde,GTilde,DTilde,ETilde,E;
  clMatrix interpLin;
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
  clVector fint,gravity,hSmooth,mu,rho,cp,kt;
  clVector heatFlux,heatFluxOld;

  clDMatrix kappaOld;
  clVector uALEOld,vALEOld,wALEOld;
  clVector uSolOld,vSolOld,wSolOld,pSolOld,cSolOld;
  clVector fintOld,gravityOld,Fold,muOld,rhoOld,hSmoothOld,cpOld,ktOld;

  Solver *solverV,*solverP,*solverC;
};

#endif
