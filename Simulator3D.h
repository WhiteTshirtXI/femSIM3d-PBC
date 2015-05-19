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
#include "Periodic3D.h"
#include "Linalg.h"

// Boost header file 
#include <boost/math/special_functions/erf.hpp>

class Simulator3D
{
 public:
  Simulator3D(); // construtor padrao
  Simulator3D( Model3D &_m ); // construtor 
  Simulator3D( Periodic3D &_pbc, Model3D &_m ); // PBC 
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

  // PBC ******* 
  void assemblePBC();
  void assemblePBCNew(); // idem assemblePBC with vector structure
  void assembleCPBC();
  void assembleCPBCNew(); // PBC
  void unCoupledPBC(); // modified solution velocity+pressure
  void unCoupledPBCNew(); 
  void unCoupledBeta();
  void unCoupledCPBC(); // modified solution scalar
  void unCoupledCPBCNew(); // modified solution scalar
  void getPeriodic3DToAttrib(Periodic3D &_pbc);
  void setPressureJump(double _pJump);
  void setBetaPressureLiquid(double _val);
  void setBetaFlowLiq(const char* _direction);
  void setRHSBeta();
  void inputVelocityPBC();
  void initTaylorVortex();
  void initTaylorGreenVortex();
  void initTanHJetProfile(); // hyperbolic tangent jet profile
  void initJetVelocity(double _vel); // inner phase with velocity
  void initPastCylinderFlow(double _U, double _V);
  void initDoubletFlow(double _force);
  void inputPurePressurePBC();
  void sumIndexPBCVel(clVector* _indexL, clVector* _indexR, clVector& _b);
  void sumIndexPBCVelNew(vector<int>* _indexL, vector<int>* _indexR, clVector& _b);
  void sumIndexPBCScalar(clVector* _indexL, clVector* _indexR, clVector& _b);
  void sumIndexPBCScalarNew(vector<int>* _indexL, vector<int>* _indexR, clVector& _b);
  void setCopyDirectionPBC(string _direction);
  void stepSLPBCFix();
  void setDirichletPressurePointPBC(string _method);
  double getTaylorVortexError();
  void initCTwoShearLayers(double _cLayerXBot, double _cLayerXTop); // PBC
  void initCGaussian(double _peak);
  
  // gets & sets
  double getPeriodicFaceVelXAverage();
  double getPeriodicFaceVelYAverage();
  double getPeriodicFaceVelZAverage();
  vector<double> getPeriodicFaceTimeAveragePressure(const char* _type); 
  void setBetaPressureLiquidTimeAverage(const char* _direction,const char* _type);
  double getBetaPressLiq();
  double getMeanPressureDomain(string _type);
  double calcTaylorVortexError();
  //******* 

  void step();
  void stepSL();
  void stepNoConvection();
  void stepImposedPeriodicField(const char* _name, double T,double _time);
  void stepImposedPeriodicField(const char* _name,double T,double _time,double _dt);
  void copyALEtoSol();
  void stepLagrangian();
  void stepLagrangianZ();
  void stepALE();
  void stepALEPBC(); // for PBC SL correction
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
  void timeStep();

  void setCoupledBC();
  void setUnCoupledBC();
  void setUnCoupledCBC();

  void setNuZ(const char* _filename);
  void setNuC();
  void setRe(double _Re);
  double getRe();
  void setSc(double _Sc);
  double getSc();
  void setFr(double _Fr);
  double getFr();
  void setWe(double _We);
  double getWe();
  void setAlpha(double _alpha);
  double getAlpha();
  void setBeta(double _beta);
  double getBeta();
  void setSigma(double _sigma);
  double getSigma();
  void setDtLagrangian();
  double setDtHeight(clVector &_uVel,clVector &_vVel,clVector &_wVel);
  void setDtLagrangianNorberto();
  void setDtSemiLagrangian();
  void setDtGravity();
  void setDtSurfaceTension();
  void setDtEulerian();
  void setDtALESinglePhase();
  void setDtALETwoPhase();
  void setDt();
  void setDt(double _dt);
  void setTime(double _time);
  double getTime();
  void setCfl(double _cfl);
  void setCflBubble(double _cfl);
  double getDt();
  double getDtLagrangian();
  double getDtSemiLagrangian();
  double getDtSurfaceTension();
  double getDtGravity();
  double getCfl();
  void setIter(double _Iter);
  int getIter();
  void setC1(double _c1);
  void setC2(double _c2);
  void setC3(double _c3);
  void setD1(double _d1);
  void setD2(double _d2);
  double getC1();
  double getC2();
  double getC3();
  double getD1();
  double getD2();
  void setURef(double _uRef);
  double getURef();
  void setVRef(double _vRef);
  double getVRef();
  void setWRef(double _wRef);
  double getWRef();
  void setXRef(double _XRef);
  double getXRef();
  void setYRef(double _YRef);
  double getYRef();
  void setZRef(double _ZRef);
  double getZRef();
  void setMu(double _mu_in);
  void setMu(double _mu_in,double _mu_out);
  void setMuSmooth(double _mu_in,double _mu_out);
  double getMu_in();
  double getMu_out();
  double getMu_inAdimen();
  double getMu_outAdimen();
  void setRho(double _rho_in);
  void setRho(double _rho_in,double _rho_out);
  void setRrho_in(double _rho_in);
  void setRhoSmooth(double _rho_in,double _rho_out);
  void setCp(double _cp_in);
  void setCp(double _cp_in,double _cp_out);
  void setCpSmooth(double _cp_in,double _cp_out);
  double getCp_in();
  double getCp_out();
  double getCp_inAdimen();
  double getCp_outAdimen();
  void setKt(double _kt_in);
  void setKt(double _kt_in,double _kt_out);
  void setKtSmooth(double _kt_in,double _kt_out);
  double getKt_in();
  double getKt_out();
  double getKt_inAdimen();
  double getKt_outAdimen();
  void setHSmooth();
  void setHeatFlux();
  double getRho_in();
  void setRho_out(double _rho_out);
  double getRho_out();
  double getRho_inAdimen();
  double getRho_outAdimen();
  void setUAnt(clVector &_uAnt);
  void setCSol(clVector &_cSol);
  void setUSol(clVector *_uSol);
  void setVSol(clVector *_vSol);
  void setWSol(clVector *_wSol);
  void setUALE(clVector *_uALE);
  void setVALE(clVector *_vALE);
  void setWALE(clVector *_wALE);
  void updateIEN();
  clVector setCentroid(clVector &_vector);
  void setUSol(double _vel);
  void setVSol(double _vel);
  void setWSol(double _vel);

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
  clVector* getBetaFlowLiq(); // PBC
  double getGrav();
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

  vector<double> getCentroidVelX();
  double getCentroidVelXAverage();
  vector<double> getCentroidVelY();
  double getCentroidVelYAverage();
  vector<double> getCentroidVelZ();
  double getCentroidVelZAverage();
  double getCentroidVelXMax();
  double getCentroidVelYMax();
  double getCentroidVelZMax();
  double getCentroidVelXMin();
  double getCentroidVelYMin();
  double getCentroidVelZMin();

  void setCentroidVelX(vector<double> _centroidVelX);
  void setCentroidVelY(vector<double> _centroidVelY);
  void setCentroidVelZ(vector<double> _centroidVelZ);

  vector<double> getCentroidPosX();
  double getCentroidPosXAverage();
  vector<double> getCentroidPosY();
  double getCentroidPosYAverage();
  vector<double> getCentroidPosZ();
  double getCentroidPosZAverage();
  void setCentroidPosX(vector<double> _centroidPosX);
  void setCentroidPosY(vector<double> _centroidPosY);
  void setCentroidPosZ(vector<double> _centroidPosZ);
  void setSurfaceTSat();
  void setMassTransfer();

 private:
  Model3D *m;
  int numVerts,numElems,numNodes;
  int numVertsOld,numElemsOld,numNodesOld;
  vector<double> triEdge;
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


  double Re,Sc,Fr,We,alpha,beta,cfl,time;
  double dt,dtLagrangian,dtSemiLagrangian,dtSurfaceTension,dtGravity;
  double c1,c2,c3,d1,d2;
  double bubbleXVel,bubbleYVel,bubbleZVel;
  double sigma,g,rho_in,rho_out,mu_in,mu_out;
  double sigma_0,g_0,rho_0,mu_0;
  double cp_in,cp_out,cp_0,cp_inAdimen,cp_outAdimen;
  double kt_in,kt_out,kt_0,kt_inAdimen,kt_outAdimen;
  double sigmaAdimen,gAdimen,rho_inAdimen,rho_outAdimen,mu_inAdimen,mu_outAdimen;
  double bubbleVel;
  int iter;
  vector<double> centroidVelX,centroidVelY,centroidVelZ;
  vector<double> centroidVelXOld,centroidVelYOld,centroidVelZOld;
  vector<double> centroidPosX,centroidPosY,centroidPosZ;
  vector<double> centroidPosXOld,centroidPosYOld,centroidPosZOld;

  // moving referential
  double uRef,vRef,wRef;
  double xRef,yRef,zRef;

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
  clVector fint,gravity,betaFlowLiq,hSmooth,mu,rho,cp,kt;
  clVector heatFlux,heatFluxOld;

  clDMatrix kappaOld;
  clVector uALEOld,vALEOld,wALEOld;
  clVector uSolOld,vSolOld,wSolOld,pSolOld,cSolOld;
  clVector fintOld,gravityOld,betaFlowLiqOld,Fold,muOld,rhoOld,hSmoothOld,cpOld,ktOld;

  Solver *solverV,*solverP,*solverC;

  // PBC
  Periodic3D *pbc;
  int nyPointsL;
  int IBR;
  double pJump;
  double betaPressLiq_0, betaPressGas_0;
  double betaPressLiq, betaPressGas;
  double betaPressLiqAdimen, betaPressGasAdimen;
  double betaPressTimeAverage;
  double pMaster, pSlave;
  string direction;
  clVector *VecXMin, *VecXMax,*VecXMid, *VecXMidVerts;
  clVector VecXMinGlob, VecXMaxGlob;
  vector<int> *MasterIndices;
  vector<int> *SlaveIndices;
  vector<int> *MasterElements;
  vector<int> *SlaveElements;
  vector<double> PeriodicFacePressures;

  double TaylorGreenError, TaylorGreenErrorPoint;
  double TaylorVortexError, TaylorVortexErrorPoint;
  double firstStokesProblemError, firstStokesProblemErrorPoint;
  double normAnalPoint, normNumPoint, uAnalPoint, vAnalPoint, wAnalPoint;
};

#endif
