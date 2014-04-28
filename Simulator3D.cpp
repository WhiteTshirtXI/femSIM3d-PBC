// =================================================================== //
// this is file Simulator3D.cpp, created at 23-Ago-2007                //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== // 
#include "Simulator3D.h"

Simulator3D::Simulator3D(){}

Simulator3D::Simulator3D( Periodic3D &_pbc, Model3D &_m )  
{
 getModel3DAttrib(_m);
 getPeriodic3DToAttrib(_pbc);
 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 Re    = 10;
 Sc    = 2;
 Fr    = 0.1;
 We    = 10;
 alpha = 1;
 beta  = 0;
 dt    = 0.01;
 dtSemiLagrangian = 0.01;
 dtLagrangian = 0.01;
 dtSurfaceTension = 0.01;
 dtGravity = 0.01;
 time  = 0.0;
 cfl   = 0.5;
 iter  = 0;

 c1    = 1.0;
 c2    = 0.0;
 c3    = 0.0;
 d1    = 1.0;
 d2    = 0.1;

 g     = 9.81;
 sigma = 1.0;
 mu_in  = 1.0;
 mu_out  = 1.0;
 rho_in = 1.0;
 rho_out = 1.0;
 cp_in = 1.0;
 cp_out = 1.0;
 kt_in = 1.0;
 kt_out = 1.0;
 uRef = 0.0;
 vRef = 0.0;
 wRef = 0.0;
 xRef = 0.0;
 yRef = 0.0;
 zRef = 0.0;

 betaPressLiq = 0.0;

 allocateMemoryToAttrib();

}


Simulator3D::Simulator3D( Model3D &_m )  
{
 getModel3DAttrib(_m);

 // Simulator3D pre-configures parameters
 // =================================================================== //
 // Re   <= 1   -> viscous flow
 // Re   >> 2000 -> turbulent flow
 // ------------------------------------------------------------------- //
 // alpha = 0   -> explicity
 // alpha = 0.5 -> crank-nicholson
 // alpha = 1   -> implicity
 // ------------------------------------------------------------------- //
 // dt          -> time between 2 iterations
 // ------------------------------------------------------------------- //
 // cfl         -> stability condiction (velocity x length)
 // ------------------------------------------------------------------- //
 // Solver      -> Pre-Conjugate Gradiente     :PCGSolver
 //             -> Conjugate Gradiente         :CGSolver
 //             -> GMRes                       :GMRes
 //             -> direct method (library: GSL):GSLSolver
 // =================================================================== //
 Re    = 10;
 Sc    = 2;
 Fr    = 0.1;
 We    = 10;
 alpha = 1;
 beta  = 0;
 dt    = 0.01;
 dtSemiLagrangian = 0.01;
 dtLagrangian = 0.01;
 dtSurfaceTension = 0.01;
 dtGravity = 0.01;
 time  = 0.0;
 cfl   = 0.5;
 iter  = 0;

 c1    = 1.0;
 c2    = 0.0;
 c3    = 0.0;
 d1    = 1.0;
 d2    = 0.1;

 g     = 9.81;
 sigma = 1.0;
 mu_in  = 1.0;
 mu_out  = 1.0;
 rho_in = 1.0;
 rho_out = 1.0;
 cp_in = 1.0;
 cp_out = 1.0;
 kt_in = 1.0;
 kt_out = 1.0;
 uRef = 0.0;
 vRef = 0.0;
 wRef = 0.0;
 xRef = 0.0;
 yRef = 0.0;
 zRef = 0.0;

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 allocateMemoryToAttrib();
}

/* the constructor will be initialize with &_sRight attributes and
 * &_sRight.m attributes, i.e. copying the pointer from _sRight to the
 * current object.
 * 
 * input: &_s
 * output: Simulator3D &_sRight and &_sRight.m
 *
 * */
Simulator3D::Simulator3D( const Simulator3D &_sRight )  
{
 getModel3DAttrib(*_sRight.m);

 Re = _sRight.Re;
 Sc = _sRight.Sc;
 Fr = _sRight.Fr;
 We = _sRight.We;
 alpha = _sRight.alpha;
 beta = _sRight.beta;
 dt = _sRight.dt;
 dtSemiLagrangian = _sRight.dtSemiLagrangian;
 dtLagrangian = _sRight.dtLagrangian;
 dtSurfaceTension = _sRight.dtSurfaceTension;
 dtGravity = _sRight.dtGravity;
 cfl = _sRight.cfl;
 time = _sRight.time;
 c1 = _sRight.c1;
 c2 = _sRight.c2;
 c3 = _sRight.c3;
 d1 = _sRight.d1;
 d2 = _sRight.d2;
 iter = _sRight.iter;

 g = _sRight.g;
 sigma = _sRight.sigma;
 rho_in = _sRight.rho_in;
 rho_out = _sRight.rho_out;
 mu_in = _sRight.mu_in;
 mu_out = _sRight.mu_out;
 cp_in = _sRight.cp_in;
 cp_out = _sRight.cp_out;
 kt_in = _sRight.kt_in;
 kt_out = _sRight.kt_out;
 
 uRef = _sRight.uRef;
 vRef = _sRight.vRef;
 wRef = _sRight.wRef;
 xRef = _sRight.xRef;
 yRef = _sRight.yRef;
 zRef = _sRight.zRef;

 g_0 = _sRight.g_0;
 sigma_0 = _sRight.sigma_0;
 rho_0 = _sRight.rho_0;
 mu_0 = _sRight.mu_0;
 cp_0 = _sRight.cp_0;
 kt_0 = _sRight.kt_0;

 gAdimen = _sRight.gAdimen;
 sigmaAdimen = _sRight.sigmaAdimen;
 rho_inAdimen = _sRight.rho_inAdimen;
 rho_outAdimen = _sRight.rho_outAdimen;
 mu_inAdimen = _sRight.mu_inAdimen;
 mu_outAdimen = _sRight.mu_outAdimen;
 cp_inAdimen = _sRight.cp_inAdimen;
 cp_outAdimen = _sRight.cp_outAdimen;
 kt_inAdimen = _sRight.kt_inAdimen;
 kt_outAdimen = _sRight.kt_outAdimen;

 K = _sRight.K;
 Kc = _sRight.Kc;
 Mrho = _sRight.Mrho;
 M = _sRight.M;
 Mc = _sRight.Mc;
 G = _sRight.G;
 Gc = _sRight.Gc;
 D = _sRight.D;
 mat = _sRight.mat;
 matc = _sRight.matc;
 MrhoLumped = _sRight.MrhoLumped;
 MLumped = _sRight.MLumped;
 McLumped = _sRight.McLumped;
 gx = _sRight.gx;
 gy = _sRight.gy;
 gz = _sRight.gz;
 A = _sRight.A;
 b = _sRight.b;
 ATilde = _sRight.ATilde;
 AcTilde = _sRight.AcTilde;
 GTilde = _sRight.GTilde;
 DTilde = _sRight.DTilde;
 ETilde = _sRight.ETilde;
 E = _sRight.E;
 invA = _sRight.invA;
 invC = _sRight.invC;
 invMrhoLumped = _sRight.invMrhoLumped;
 invMLumped = _sRight.invMLumped;
 invMcLumped = _sRight.invMcLumped;

 uSol = _sRight.uSol;
 vSol = _sRight.vSol;
 wSol = _sRight.wSol;
 pSol = _sRight.pSol;
 cSol = _sRight.cSol;
 velU = _sRight.velU;
 velV = _sRight.velV;
 velW = _sRight.velW;
 uALE = _sRight.uALE;
 vALE = _sRight.vALE;
 wALE = _sRight.wALE;
 uSL = _sRight.uSL;
 vSL = _sRight.vSL;
 wSL = _sRight.wSL;
 cSL = _sRight.cSL;
 uSmooth = _sRight.uSmooth;
 vSmooth = _sRight.vSmooth;
 wSmooth = _sRight.wSmooth;
 uSmoothCoord = _sRight.uSmoothCoord;
 vSmoothCoord = _sRight.vSmoothCoord;
 wSmoothCoord = _sRight.wSmoothCoord;
 uSmoothSurface = _sRight.uSmoothSurface;
 vSmoothSurface = _sRight.vSmoothSurface;
 wSmoothSurface = _sRight.wSmoothSurface;
 centroidVelX = _sRight.centroidVelX;
 centroidVelY = _sRight.centroidVelY;
 centroidVelZ = _sRight.centroidVelZ;
 centroidPosX = _sRight.centroidPosX;
 centroidPosY = _sRight.centroidPosY;
 centroidPosZ = _sRight.centroidPosZ;

 va = _sRight.va;
 vcc = _sRight.vcc;
 convUVW = _sRight.convUVW;
 convC = _sRight.convC;
 uTilde = _sRight.uTilde;
 cTilde = _sRight.cTilde;
 pTilde = _sRight.pTilde;
 b1 = _sRight.b1;
 b1c = _sRight.b1c;
 b2 = _sRight.b2;
 ip = _sRight.ip;
 ipc = _sRight.ipc;

 // two-phase vectors
 kappa   = _sRight.kappa;
 fint    = _sRight.fint;
 gravity = _sRight.gravity;
 betaFlowLiq = _sRight.betaFlowLiq;
 Fold    = _sRight.Fold;
 mu      = _sRight.mu;
 rho     = _sRight.rho;
 cp      = _sRight.cp;
 kt      = _sRight.kt;
 hSmooth = _sRight.hSmooth;
 heatFlux= _sRight.heatFlux;
 
 // old ints
 numVertsOld = _sRight.numVerts;
 numNodesOld = _sRight.numNodes;
 numElemsOld = _sRight.numElems;
 
 // oldSol vectors
 uSolOld    = _sRight.uSolOld;
 vSolOld    = _sRight.vSolOld;
 wSolOld    = _sRight.wSolOld;
 pSolOld    = _sRight.pSolOld;
 cSolOld    = _sRight.cSolOld;
 uALEOld    = _sRight.uALEOld;
 vALEOld    = _sRight.vALEOld;
 wALEOld    = _sRight.wALEOld;
 kappaOld   = _sRight.kappaOld;
 fintOld    = _sRight.fintOld;
 gravityOld = _sRight.gravityOld;
 betaFlowLiqOld = _sRight.betaFlowLiqOld;
 muOld      = _sRight.muOld;
 rhoOld     = _sRight.rhoOld;
 cpOld      = _sRight.cpOld;
 ktOld      = _sRight.ktOld;
 hSmoothOld = _sRight.hSmoothOld;
 heatFluxOld = _sRight.heatFluxOld;
 centroidVelXOld = _sRight.centroidVelXOld;
 centroidVelYOld = _sRight.centroidVelYOld;
 centroidVelZOld = _sRight.centroidVelZOld;
 centroidPosXOld = _sRight.centroidPosXOld;
 centroidPosYOld = _sRight.centroidPosYOld;
 centroidPosZOld = _sRight.centroidPosZOld;

 solverV = _sRight.solverV;
 solverP = _sRight.solverP;
 solverC = _sRight.solverC;
}

/* the constructor will be initialize with &_m attributes and &_sRight
 * attributes. numVertsOld,numNodesOld and numElentsOld will be
 * initialize with &_sRight.m attributes
 * 
 * input: &_m,&_s
 * output: Simulator3D with &_m and &_s
 *
 * */
Simulator3D::Simulator3D( Model3D &_m, Simulator3D &_sRight)  
{
 // mesh information vectors
 getModel3DAttrib(_m);

 Re    = _sRight.getRe();
 Sc    = _sRight.getSc();
 Fr    = _sRight.getFr();
 We    = _sRight.getWe();
 sigma = _sRight.getSigma();
 alpha = _sRight.getAlpha();
 beta  = _sRight.getBeta();
 dt    = _sRight.getDt();
 dtLagrangian = _sRight.getDtLagrangian();
 dtSemiLagrangian = _sRight.getDtSemiLagrangian();
 dtGravity = _sRight.getDtGravity();
 dtSurfaceTension = _sRight.getDtSurfaceTension();
 time  = _sRight.getTime();
 cfl   = _sRight.getCfl();
 g     = _sRight.getGrav();
 mu_in  = _sRight.getMu_in();
 mu_out  = _sRight.getMu_out();
 mu_inAdimen  = _sRight.getMu_inAdimen();
 mu_outAdimen  = _sRight.getMu_outAdimen();
 rho_in = _sRight.getRho_in();
 rho_out = _sRight.getRho_out();
 rho_inAdimen = _sRight.getRho_inAdimen();
 rho_outAdimen = _sRight.getRho_outAdimen();
 cp_in = _sRight.getCp_in();
 cp_out = _sRight.getCp_out();
 cp_inAdimen = _sRight.getCp_inAdimen();
 cp_outAdimen = _sRight.getCp_outAdimen();
 kt_in = _sRight.getKt_in();
 kt_out = _sRight.getKt_out();
 kt_inAdimen = _sRight.getKt_inAdimen();
 kt_outAdimen = _sRight.getKt_outAdimen();
 iter  = _sRight.getIter();
 c1    = _sRight.getC1();
 c2    = _sRight.getC2();
 c3    = _sRight.getC3();
 d1    = _sRight.getD1();
 d2    = _sRight.getD2();

 uRef = _sRight.getURef();
 vRef = _sRight.getVRef();
 wRef = _sRight.getWRef();
 xRef = _sRight.getXRef();
 yRef = _sRight.getYRef();
 zRef = _sRight.getZRef();

 numVertsOld = _sRight.m->getNumVerts();
 numNodesOld = _sRight.m->getNumNodes();
 numElemsOld = _sRight.m->getNumElems();

 allocateMemoryToAttrib();

 // recuperando campo de velocidade e pressao da malha antiga
 uSolOld    = *_sRight.getUSol();
 vSolOld    = *_sRight.getVSol();
 wSolOld    = *_sRight.getWSol();
 pSolOld    = *_sRight.getPSol();
 cSolOld    = *_sRight.getCSol();
 uALEOld    = *_sRight.getUALE();
 vALEOld    = *_sRight.getVALE();
 wALEOld    = *_sRight.getWALE();
 kappaOld   = *_sRight.getKappa();
 fintOld    = *_sRight.getFint();
 gravityOld = *_sRight.getGravity();
 betaFlowLiqOld = *_sRight.getBetaFlowLiq();
 muOld      = *_sRight.getMu();
 rhoOld     = *_sRight.getRho();
 cpOld     = *_sRight.getCp();
 ktOld     = *_sRight.getKt();
 hSmoothOld = *_sRight.getHSmooth();
 heatFluxOld = *_sRight.getHeatFlux();
 centroidVelXOld = _sRight.getCentroidVelX();
 centroidVelYOld = _sRight.getCentroidVelY();
 centroidVelZOld = _sRight.getCentroidVelZ();
 centroidPosXOld = _sRight.getCentroidPosX();
 centroidPosYOld = _sRight.getCentroidPosY();
 centroidPosZOld = _sRight.getCentroidPosZ();

 solverV = _sRight.solverV;
 solverP = _sRight.solverP;
 solverC = _sRight.solverC;
}

Simulator3D::~Simulator3D()
{
 //delete solverV;
 //delete solverP;
 //delete solverC;
}

void Simulator3D::init()
{
 uSolOld.CopyFrom( 0,*uc );
 vSolOld.CopyFrom( 0,*vc );
 wSolOld.CopyFrom( 0,*wc );
 pSolOld.CopyFrom( 0,*pc );
 cSolOld.CopyFrom( 0,*cc );

 uSol.CopyFrom( 0,*uc );
 vSol.CopyFrom( 0,*vc );
 wSol.CopyFrom( 0,*wc );
 pSol.CopyFrom( 0,*pc );
 cSol.CopyFrom( 0,*cc );
}

void Simulator3D::initDiskBaseState( const char* _dir,const char* _filename )
{
 init(); 

 double aux = 0.0;
 clMatrix solFile(2401,5); 

 string file = (string) _dir + (string) _filename;
 const char* filename = file.c_str();
 ifstream fileA( filename,ios::in );

 if( !fileA )
 {
  cerr << "Esta faltando o arquivo de perfis!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !fileA.eof() )
 {
  for( int i=0;i<solFile.DimI();i++ )
  {
   fileA >> aux;
   solFile.Set(i,0,aux);
   fileA >> aux;
   solFile.Set(i,1,aux);
   fileA >> aux;
   solFile.Set(i,2,aux);
   fileA >> aux;
   solFile.Set(i,3,aux);
   fileA >> aux;
   solFile.Set(i,4,aux);
  }
 }

 int j;
 double L,L1,L2;
 double omega = 1.0;
 double EPSlocal = 1e-04;
 for( int i=0;i<numNodes;i++ )
 {
  for( j=0;j<solFile.DimI()-1;j++ )
  {
   L = solFile(j+1,0)-solFile(j,0);
   L1 = ( Z->Get(i)-solFile(j,0) )/L;
   L2 = 1.0-L1;
   if( (L1>=0.0-EPSlocal) && (L1<=1.0+EPSlocal) &&
       (L2>=0.0-EPSlocal) && (L2<=1.0+EPSlocal) ) break;
  }
  // interpolant
  double interp = (Z->Get(i)-solFile(j,0))/(solFile(j+1,0)-solFile(j,0));
  double FA = solFile(j,1)+(solFile(j+1,1)-solFile(j,1))*interp;
  double GA = solFile(j,2)+(solFile(j+1,2)-solFile(j,2))*interp;
  double HA = solFile(j,3)+(solFile(j+1,3)-solFile(j,3))*interp;
  double CA = solFile(j,4)+(solFile(j+1,4)-solFile(j,4))*interp;

  aux = ( FA*X->Get(i)-GA*Y->Get(i) )*omega; // F
  uSol.Set(i,aux);
  uSolOld.Set(i,aux);
  aux = ( GA*X->Get(i)+FA*Y->Get(i) )*omega; // G
  vSol.Set(i,aux);
  vSolOld.Set(i,aux);
  aux = (-1)*HA; // H (positive on file)
  wSol.Set(i,aux);
  wSolOld.Set(i,aux);
  aux = CA; // C

  if( i < numVerts )
  {
   cSol.Set(i,aux);
   cSolOld.Set(i,aux);
  }
 }
}


void Simulator3D::initAnnular()
{
 init();

 //double p1y = Y->Min();
 //double p2y = Y->Max();
 //double p1z = Z->Min();
 //double p2z = Z->Max();

 for( int i=0;i<numNodes;i++ )
 {
  // gas
  if( heaviside->Get(i) < 1.0 )
  {
   double aux = 1.0;
   wSol.Set(i,aux);
   wSolOld.Set(i,aux);
  }
  else
  {
   double aux = 2.5;
   wSol.Set(i,aux);
   wSolOld.Set(i,aux);
  }
 }
}

void Simulator3D::initChannel()
{
 init();

 double p1y = Y->Min();
 double p2y = Y->Max();
 double p1z = Z->Min();
 double p2z = Z->Max();

 // calculating channel's diameter.
 double diameterYZ = ( dist(p1y,p2y) + 
                     dist(p1z,p2z) ) / 2.0;

 for( int i=0;i<numNodes;i++ )
 {
  double radius = sqrt( Y->Get(i)*Y->Get(i) + 
	                  Z->Get(i)*Z->Get(i) );

  // Parabolic profile
  double Umax = 1.0;
  double aux = 2*Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
	                                     (diameterYZ/2.0)) );
  aux = 1.0;
  uSol.Set(i,aux-1);
  uSolOld.Set(i,aux-1);
 }
}

void Simulator3D::initChannelSquare()
{
 init();

 double p1y = Y->Min();
 double p2y = Y->Max();
 double p1z = Z->Min();
 double p2z = Z->Max();

 // calculating channel's diameter.
 double diameterYZ = distance(p1y,p1z,p2y,p2z); 

 for( int i=0;i<numNodes;i++ )
 {
  double radius = sqrt( Y->Get(i)*Y->Get(i) + 
	                  Z->Get(i)*Z->Get(i) );

  // Parabolic profile
  double Umax = 1.0;
  double aux = 2*Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
	                                     (diameterYZ/2.0)) );

  aux = 1.0;
  uSol.Set(i,aux-1.0);
  uSolOld.Set(i,aux-1.0);
 }
}

void Simulator3D::initHeatTransfer()
{
//--------------------------------------------------
//  for( int i=0;i<numVerts;i++ )
//   if( Z->Get(i) > 0.4*Z->Max() )
//    cSolOld.Set(i,1.0);
//-------------------------------------------------- 

/* heat transfer */
 for( int i=0;i<numVerts;i++ )
  if( heaviside->Get(i) < 0.5 )
   cSolOld.Set(i,0.001);
}

void Simulator3D::initFixedBubbleZ()
{
 init();

 for( int i=0;i<numNodes;i++ )
 {
  double aux = 1.0;
  wSolOld.Set(i,aux);
 }
 for( int i=0;i<idbcw->Dim();i++ )
  wSolOld.Set( (int) idbcw->Get(i),0.0 ); 
}

void Simulator3D::init2AxiBubbles()
{
 init();

/* two bubbles */
 for( int i=0;i<numNodes;i++ )
 {
  double aux = X->Get(i);
  uSolOld.Set(i,aux);
  aux = -1.0*Y->Get(i);
  vSolOld.Set(i,aux);
  aux = 0.0;
  wSolOld.Set(i,aux);
 }
}

void Simulator3D::init2Bubbles()
{
 init();

/* two bubbles */
 for( int i=0;i<numNodes;i++ )
 {
  double aux = X->Get(i);
  uSolOld.Set(i,aux);
  aux = -1.0*Y->Get(i);
  vSolOld.Set(i,aux);
  aux = Z->Get(i);
  wSolOld.Set(i,aux);
 }
}

void Simulator3D::assemble()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyx( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzx( numNodes,numNodes );
 clMatrix Kzy( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Mx( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 //setMu( mu_in );
 //setRho( rho_in );
 setMu( mu_in,mu_out );
 setRho( rho_in,rho_out );

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  double muValue=0;
  double rhoValue=0;
  if( elemIdRegion->Get(mele) == 0.0 ) // out
  {
   muValue = mu_outAdimen;
   rhoValue = rho_outAdimen;
  }
  else
  {
   muValue = mu_inAdimen;
   rhoValue = rho_inAdimen;
  }

  miniElem.getM(*v);  // para problemas SEM deslizamento

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muValue*( 2*miniElem.kxx[i][j] + 
	                                   miniElem.kyy[i][j] + 
									   miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux); // matriz de massa

	aux = Mx.Get(ii,jj) + miniElem.massele[i][j];
	Mx.Set(ii,jj,aux); // matriz de massa sem rho

	// bloco 12
	aux = Kxy.Get(ii,jj) + muValue*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muValue*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                               2*miniElem.kyy[i][j] + 
							         miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muValue*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                                 miniElem.kyy[i][j] + 
								   2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);
   }
   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
 }
 
 /* GALERKIN */
 // gx = Gx;
 // gy = Gy;
 // gz = Gz;
 Mrho.CopyFrom(          0,          0,  Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes,  Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes,  Mx_rho );

 M.CopyFrom(          0,          0,     Mx );
 M.CopyFrom(   numNodes,   numNodes,     Mx );
 M.CopyFrom( 2*numNodes, 2*numNodes,     Mx );

 K.CopyFrom(          0,          0,     Kxx );
 K.CopyFrom(          0,   numNodes,     Kxy );
 K.CopyFrom(          0, 2*numNodes,     Kxz );
 K.CopyFrom(   numNodes,          0,     Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,     Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,     Kyz );
 K.CopyFrom( 2*numNodes,          0,     Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes,     Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,     Kzz );

 G.CopyFrom(          0,          0,     Gx );
 G.CopyFrom(   numNodes,          0,     Gy );
 G.CopyFrom( 2*numNodes,          0,     Gz );
 D.CopyFrom(          0,          0,     Gx.Transpose() );
 D.CopyFrom(          0,   numNodes,     Gy.Transpose() );
 D.CopyFrom(          0, 2*numNodes,     Gz.Transpose() );
} // fecha metodo ASSEMBLE

void Simulator3D::assembleHeatTransfer()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyx( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzx( numNodes,numNodes );
 clMatrix Kzy( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Mx( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );
 clMatrix GcMatx( numVerts,numVerts );
 clMatrix GcMaty( numVerts,numVerts );
 clMatrix GcMatz( numVerts,numVerts );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 FEMLinElement3D linElem(*X,*Y,*Z);

 //setMu( mu_in );
 //setRho( rho_in );
 //setCp( cp_in );
 //setKt( kt_in );
 setMu( mu_in,mu_out );
 setRho( rho_in,rho_out );
 setCp( cp_in,cp_out );
 setKt( kt_in,kt_out );

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  double muValue=0;
  double rhoValue=0;
  //double cpValue=0;
  double ktValue=0;
  if( elemIdRegion->Get(mele) == 0.0 ) // out
  {
   muValue = mu_outAdimen;
   rhoValue = rho_outAdimen;
   //cpValue = cp_outAdimen;
   ktValue = kt_outAdimen;
  }
  else
  {
   muValue = mu_inAdimen;
   rhoValue = rho_inAdimen;
   //cpValue = cp_inAdimen;
   ktValue = kt_inAdimen;
  }

  miniElem.getM(*v);  // para problemas SEM deslizamento
  linElem.getM(*v); 

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muValue*( 2*miniElem.kxx[i][j] + 
	                                   miniElem.kyy[i][j] + 
									   miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux); // matriz de massa

	aux = Mx.Get(ii,jj) + miniElem.massele[i][j];
	Mx.Set(ii,jj,aux); // matriz de massa sem rho

	// bloco 12
	aux = Kxy.Get(ii,jj) + muValue*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muValue*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                               2*miniElem.kyy[i][j] + 
							         miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muValue*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                                 miniElem.kyy[i][j] + 
								   2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);
   }
   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + ktValue*( linElem.kxxc[i][j] + 
	                                   linElem.kyyc[i][j] + 
							           linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);

	//aux = McMat.Get(ii,jj) + rhoValue*linElem.masselec[i][j];
	aux = McMat.Get(ii,jj) + linElem.masselec[i][j];
	McMat.Set(ii,jj,aux);

	aux = GcMatx.Get(ii,jj) + linElem.gxelec[i][j];
	GcMatx.Set(ii,jj,aux);

	aux = GcMaty.Get(ii,jj) + linElem.gyelec[i][j];
	GcMaty.Set(ii,jj,aux);

	aux = GcMatz.Get(ii,jj) + linElem.gzelec[i][j];
	GcMatz.Set(ii,jj,aux);
   }
  }
 }
 
 /* GALERKIN */
 // gx = Gx;
 // gy = Gy;
 // gz = Gz;
 Mrho.CopyFrom(          0,          0,  Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes,  Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes,  Mx_rho );

 M.CopyFrom(          0,          0,     Mx );
 M.CopyFrom(   numNodes,   numNodes,     Mx );
 M.CopyFrom( 2*numNodes, 2*numNodes,     Mx );

 K.CopyFrom(          0,          0,     Kxx );
 K.CopyFrom(          0,   numNodes,     Kxy );
 K.CopyFrom(          0, 2*numNodes,     Kxz );
 K.CopyFrom(   numNodes,          0,     Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,     Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,     Kyz );
 K.CopyFrom( 2*numNodes,          0,     Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes,     Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,     Kzz );

 G.CopyFrom(          0,          0,     Gx );
 G.CopyFrom(   numNodes,          0,     Gy );
 G.CopyFrom( 2*numNodes,          0,     Gz );
 D.CopyFrom(          0,          0,     Gx.Transpose() );
 D.CopyFrom(          0,   numNodes,     Gy.Transpose() );
 D.CopyFrom(          0, 2*numNodes,     Gz.Transpose() );

 Kc.CopyFrom(          0,         0,     KcMat );
 Mc.CopyFrom(          0,         0,     McMat );
 Gc.CopyFrom(          0,         0,     GcMatx );
 Gc.CopyFrom(   numVerts,         0,     GcMaty );
 Gc.CopyFrom( 2*numVerts,         0,     GcMatz );
 
} // fecha metodo ASSEMBLEHEATTRANSFER

void Simulator3D::assembleC()
{
 int i,j,ii,jj;
 int v[NUMGLEC];
 double aux;
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  double dif = 1.0;

  linElem.getM(*v); 

  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + dif*( linElem.kxxc[i][j] + 
	                               linElem.kyyc[i][j] + 
								   linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);

	aux = McMat.Get(ii,jj) + linElem.masselec[i][j];
	McMat.Set(ii,jj,aux);
   }
  }
 }
 
 Kc.CopyFrom( 0, 0, KcMat );
 Mc.CopyFrom( 0, 0, McMat );
 
} // fecha metodo ASSEMBLEC

void Simulator3D::assembleNuCte()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );
 clMatrix Dx( numVerts,numNodes );
 clMatrix Dy( numVerts,numNodes );
 clMatrix Dz( numVerts,numNodes );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 setMu(mu_in);
 setRho(rho_in);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  double muValue = mu_inAdimen;
  double rhoValue = rho_inAdimen;

  miniElem.getMSlip(*v);  // para problemas SEM deslizamento

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                                 miniElem.kyy[i][j] + 
							         miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux); // matriz de massa
   };

   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEP;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Dx.Get(ii,jj) + miniElem.dxele[i][j];
	Dx.Set(ii,jj,aux);
	// bloco 2
	aux = Dy.Get(ii,jj) + miniElem.dyele[i][j];
	Dy.Set(ii,jj,aux);
	// bloco 3
	aux = Dz.Get(ii,jj) + miniElem.dzele[i][j];
	Dz.Set(ii,jj,aux);
   }
  }
 }
 
 Mrho.CopyFrom(          0,          0, Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes, Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes, Mx_rho );

 // waste of memory!!!
 M = Mrho;

 K.CopyFrom(          0,          0,    Kxx );
 K.CopyFrom(   numNodes,   numNodes,    Kxx );
 K.CopyFrom( 2*numNodes, 2*numNodes,    Kxx );

 G.CopyFrom(          0,          0,    Gx );
 G.CopyFrom(   numNodes,          0,    Gy );
 G.CopyFrom( 2*numNodes,          0,    Gz );
 D.CopyFrom(          0,          0,    Dx );
 D.CopyFrom(          0,   numNodes,    Dy );
 D.CopyFrom(          0, 2*numNodes,    Dz );
//--------------------------------------------------
//  D.CopyFrom(          0,          0,  Gx.Transpose() );
//  D.CopyFrom(          0,   numNodes,  Gy.Transpose() );
//  D.CopyFrom(          0, 2*numNodes,  Gz.Transpose() );
//-------------------------------------------------- 

}; // fecha metodo ASSEMBLENuCte

void Simulator3D::assembleNuC()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );
 clMatrix Dx( numVerts,numNodes );
 clMatrix Dy( numVerts,numNodes );
 clMatrix Dz( numVerts,numNodes );
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 FEMLinElement3D linElem(*X,*Y,*Z);

 setRho(rho_in);

 for( int mele=0;mele<numElems;mele++ )
 {
  double c = 0;
  for( int n=0;n<NUMGLEU;n++ )
  {
   v[n] = (int) IEN->Get(mele,n);
   c += cSolOld.Get(v[n]);
  }
  c = c/NUMGLE;

  double eme = 0.81315;
  double muC = exp(eme*c);
  double dif = 1.0/muC;
  double rhoValue = 1.0;

  // updating mu
  for( int n=0;n<NUMGLE;n++ )
   mu.Set(v[n],muC);

  miniElem.getMSlip(*v);  // free-slip
  linElem.getM(*v); 

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muC*( 2*miniElem.kxx[i][j] + 
	                               miniElem.kyy[i][j] + 
								   miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux); // matriz de massa

	// bloco 12
	aux = Kxy.Get(ii,jj) + muC*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muC*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muC*( miniElem.kxx[i][j] + 
	                           2*miniElem.kyy[i][j] + 
							     miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muC*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muC*( miniElem.kxx[i][j] + 
	                             miniElem.kyy[i][j] + 
							   2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);

   };

   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEP;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Dx.Get(ii,jj) + miniElem.dxele[i][j];
	Dx.Set(ii,jj,aux);
	// bloco 2
	aux = Dy.Get(ii,jj) + miniElem.dyele[i][j];
	Dy.Set(ii,jj,aux);
	// bloco 3
	aux = Dz.Get(ii,jj) + miniElem.dzele[i][j];
	Dz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + dif*( linElem.kxxc[i][j] + 
	                               linElem.kyyc[i][j] + 
								   linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);

	aux = McMat.Get(ii,jj) + linElem.masselec[i][j];
	McMat.Set(ii,jj,aux);
   }
  }
 }
 
 /* GALERKIN */
 // gx = Gx;
 // gy = Gy;
 // gz = Gz;
 
 Mrho.CopyFrom(          0,          0,  Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes,  Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes,  Mx_rho );

 // waste of memory - find solution!!!
 M = Mrho;

 K.CopyFrom(          0,          0,     Kxx );
 K.CopyFrom(          0,   numNodes,     Kxy );
 K.CopyFrom(          0, 2*numNodes,     Kxz );
 K.CopyFrom(   numNodes,          0,     Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,     Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,     Kyz );
 K.CopyFrom( 2*numNodes,          0,     Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes,     Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,     Kzz );

 G.CopyFrom(          0,          0,     Gx );
 G.CopyFrom(   numNodes,          0,     Gy );
 G.CopyFrom( 2*numNodes,          0,     Gz );
 D.CopyFrom(          0,          0,     Dx );
 D.CopyFrom(          0,   numNodes,     Dy );
 D.CopyFrom(          0, 2*numNodes,     Dz );

 Kc.CopyFrom(         0,          0,     KcMat );
 Mc.CopyFrom(         0,          0,     McMat );
 
}; // fecha metodo ASSEMBLENUC

void Simulator3D::assembleSlip()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Mx( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );
 clMatrix Dx( numVerts,numNodes );
 clMatrix Dy( numVerts,numNodes );
 clMatrix Dz( numVerts,numNodes );
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  // muValue and rhoValue = mean value of element vertices
  double muValue = 0;
  double rhoValue = 0;
  for( int n=0;n<NUMGLE;n++ )
  {
   muValue += mu.Get(v[n]);
   rhoValue += rho.Get(v[n]);
  }
  muValue = muValue/NUMGLE;
  rhoValue = rhoValue/NUMGLE;

  //miniElem.getM(*v);  // no-slip
  miniElem.getMSlip(*v);  // free-slip
  linElem.getM(*v); 

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muValue*( 2*miniElem.kxx[i][j] + 
	                                   miniElem.kyy[i][j] + 
									   miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux); // matriz de massa

	aux = Mx.Get(ii,jj) + miniElem.massele[i][j];
	Mx.Set(ii,jj,aux); // matriz de massa sem rho

	// bloco 12
	aux = Kxy.Get(ii,jj) + muValue*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muValue*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muValue*(   miniElem.kxx[i][j] + 
	                                 2*miniElem.kyy[i][j] + 
							           miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muValue*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muValue*(   miniElem.kxx[i][j] + 
	                                   miniElem.kyy[i][j] + 
							         2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);

   };

   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEP;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Dx.Get(ii,jj) + miniElem.dxele[i][j];
	Dx.Set(ii,jj,aux);
	// bloco 2
	aux = Dy.Get(ii,jj) + miniElem.dyele[i][j];
	Dy.Set(ii,jj,aux);
	// bloco 3
	aux = Dz.Get(ii,jj) + miniElem.dzele[i][j];
	Dz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + ( linElem.kxxc[i][j] + 
	                           linElem.kyyc[i][j] + 
							   linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);

	aux = McMat.Get(ii,jj) + linElem.masselec[i][j];
	McMat.Set(ii,jj,aux);
   }
  }
 }
 
 /* GALERKIN */
 // gx = Gx;
 // gy = Gy;
 // gz = Gz;

 Mrho.CopyFrom(          0,          0, Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes, Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes, Mx_rho );

 M.CopyFrom(          0,          0,    Mx );
 M.CopyFrom(   numNodes,   numNodes,    Mx );
 M.CopyFrom( 2*numNodes, 2*numNodes,    Mx );

 K.CopyFrom(          0,          0,    Kxx );
 K.CopyFrom(          0,   numNodes,    Kxy );
 K.CopyFrom(          0, 2*numNodes,    Kxz );
 K.CopyFrom(   numNodes,          0,    Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,    Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,    Kyz );
 K.CopyFrom( 2*numNodes,          0,    Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes,    Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,    Kzz );

 G.CopyFrom(          0,          0,    Gx );
 G.CopyFrom(   numNodes,          0,    Gy );
 G.CopyFrom( 2*numNodes,          0,    Gz );

 D.CopyFrom(          0,          0,    Dx );
 D.CopyFrom(          0,   numNodes,    Dy );
 D.CopyFrom(          0, 2*numNodes,    Dz );

 Kc.CopyFrom(         0,          0,    KcMat );
 Mc.CopyFrom(         0,          0,    McMat );
 
}; // fecha metodo ASSEMBLESLIP

void Simulator3D::assembleNuZ(const char* _name)
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix Mx_rho( numNodes,numNodes );
 clMatrix Gx( numNodes,numVerts );
 clMatrix Gy( numNodes,numVerts );
 clMatrix Gz( numNodes,numVerts );
 clMatrix Dx( numVerts,numNodes );
 clMatrix Dy( numVerts,numNodes );
 clMatrix Dz( numVerts,numNodes );

 setNuZ(_name);      // carregando o arquivo de perfil nuZ
 setRho(rho_in);
 
#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);
  
  // muValue and rhoValue = mean value of element vertices
  double muValue = 0;
  double rhoValue = 0;
  for( int n=0;n<NUMGLE;n++ )
  {
   muValue += mu.Get(v[n]);
   rhoValue += rho.Get(v[n]);
  }
  muValue = muValue/NUMGLE;
  rhoValue = rhoValue/NUMGLE;

  //miniElem.getM(*v);      // no-slip
  miniElem.getMSlip(*v);  // free-slip

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muValue*( 2*miniElem.kxx[i][j] + 
	                                   miniElem.kyy[i][j] + 
								       miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	aux = Mx_rho.Get(ii,jj) + rhoValue*miniElem.massele[i][j];
	Mx_rho.Set(ii,jj,aux);

	// bloco 12
	aux = Kxy.Get(ii,jj) + muValue*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muValue*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                               2*miniElem.kyy[i][j] + 
								     miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muValue*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muValue*( miniElem.kxx[i][j] + 
	                                 miniElem.kyy[i][j] + 
								   2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);

   };
   for( j=0;j<NUMGLEP;j++ )
   {
	jj=v[j];
	
	// bloco 1
	aux = Gx.Get(ii,jj) + miniElem.gxele[i][j];
	Gx.Set(ii,jj,aux);
	// bloco 2
	aux = Gy.Get(ii,jj) + miniElem.gyele[i][j];
	Gy.Set(ii,jj,aux);
	// bloco 3
	aux = Gz.Get(ii,jj) + miniElem.gzele[i][j];
	Gz.Set(ii,jj,aux);
   }
  }
  for( i=0;i<NUMGLEP;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];
	// bloco 1
	aux = Dx.Get(ii,jj) + miniElem.dxele[i][j];
	Dx.Set(ii,jj,aux);
	// bloco 2
	aux = Dy.Get(ii,jj) + miniElem.dyele[i][j];
	Dy.Set(ii,jj,aux);
	// bloco 3
	aux = Dz.Get(ii,jj) + miniElem.dzele[i][j];
	Dz.Set(ii,jj,aux);
   }
  }
 };
 
 gx = Gx;
 gy = Gy;
 gz = Gz;

 Mrho.CopyFrom(          0,          0, Mx_rho );
 Mrho.CopyFrom(   numNodes,   numNodes, Mx_rho );
 Mrho.CopyFrom( 2*numNodes, 2*numNodes, Mx_rho );

 // waste of memory - find solution!!!
 M = Mrho;

 K.CopyFrom(          0,          0,    Kxx );
 K.CopyFrom(          0,   numNodes,    Kxy );
 K.CopyFrom(          0, 2*numNodes,    Kxz );
 K.CopyFrom(   numNodes,          0,    Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,    Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,    Kyz );
 K.CopyFrom( 2*numNodes,          0,    Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes,    Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,    Kzz );

 G.CopyFrom(          0,          0,    Gx );
 G.CopyFrom(   numNodes,          0,    Gy );
 G.CopyFrom( 2*numNodes,          0,    Gz );
 D.CopyFrom(          0,          0,    Dx );
 D.CopyFrom(          0,   numNodes,    Dy );
 D.CopyFrom(          0, 2*numNodes,    Dz );
 
}; // fecha metodo ASSEMBLENUZ

void Simulator3D::assembleK()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix Kxx( numNodes,numNodes );
 clMatrix Kxy( numNodes,numNodes );
 clMatrix Kxz( numNodes,numNodes );
 clMatrix Kyy( numNodes,numNodes );
 clMatrix Kyz( numNodes,numNodes );
 clMatrix Kzz( numNodes,numNodes );
 clMatrix KcMat( numVerts,numVerts );

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  // muValue and rhoValue = mean value of element vertices
  double c = 0;
  for( int n=0;n<NUMGLE;n++ )
   c += cSolOld.Get(v[n]);
  c = c/NUMGLE;

  double eme = 0.81315;
  double muC = exp(eme*c);
  double dif = 1.0/muC;

  miniElem.getK(*v);  
  linElem.getK(*v); 

  for( i=0;i<NUMGLEU;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEU;j++ )
   {
	jj=v[j];

	// bloco 11
	aux = Kxx.Get(ii,jj) + muC*( 2*miniElem.kxx[i][j] + 
	                               miniElem.kyy[i][j] + 
								   miniElem.kzz[i][j] ); 
	Kxx.Set(ii,jj,aux);

	// bloco 12
	aux = Kxy.Get(ii,jj) + muC*( miniElem.kxy[i][j] ); 
	Kxy.Set(ii,jj,aux);

	// bloco 13
	aux = Kxz.Get(ii,jj) + muC*( miniElem.kxz[i][j] ); 
	Kxz.Set(ii,jj,aux);

	// bloco 22
	aux = Kyy.Get(ii,jj) + muC*( miniElem.kxx[i][j] + 
	                           2*miniElem.kyy[i][j] + 
							     miniElem.kzz[i][j] ); 
	Kyy.Set(ii,jj,aux);

	// bloco 23
	aux = Kyz.Get(ii,jj) + muC*( miniElem.kyz[i][j] ); 
	Kyz.Set(ii,jj,aux);

	// bloco 33
	aux = Kzz.Get(ii,jj) + muC*( miniElem.kxx[i][j] + 
	                             miniElem.kyy[i][j] + 
							   2*miniElem.kzz[i][j] ); 
	Kzz.Set(ii,jj,aux);

   };
  }
  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + dif*( linElem.kxxc[i][j] + 
	                               linElem.kyyc[i][j] + 
								   linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);
   }
  }
 }
 
 K.CopyFrom(          0,          0,             Kxx );
 K.CopyFrom(          0,   numNodes,             Kxy );
 K.CopyFrom(          0, 2*numNodes,             Kxz );
 K.CopyFrom(   numNodes,          0, Kxy.Transpose() );
 K.CopyFrom(   numNodes,   numNodes,             Kyy );
 K.CopyFrom(   numNodes, 2*numNodes,             Kyz );
 K.CopyFrom( 2*numNodes,          0, Kxz.Transpose() );
 K.CopyFrom( 2*numNodes,   numNodes, Kyz.Transpose() );
 K.CopyFrom( 2*numNodes, 2*numNodes,             Kzz );

 Kc.CopyFrom(         0,          0,           KcMat );
 
}; // fecha metodo ASSEMBLEK

//este método cria os vetores uSL, vSL e wSL relacionados ao 
//termo convectivo da equação de Navier-Stokes. 
//Mais precisamente u \cdot \Nabla u.
//convUVW = uSL, vSL e wSL em um unico vetor.
//Como o metodo eh explicito para calculo do termo convectivo, 
//convUVW vai para setRHS.
void Simulator3D::stepSL()
{
 clVector velU = uSolOld-uALE;
 clVector velV = vSolOld-vALE;
 clVector velW = wSolOld-wALE;

//--------------------------------------------------
//  double test = centroidVelZ[1];
//  for ( int i=0;i<numVerts;i++ )
//  {
//   int node = i;
//   if( heaviside->Get(node) == 0.5 )
//   {
//    velU.Set(node,0.0);
//    velV.Set(node,0.0);
//    velW.Set(node,0.0);
//    //velW.Set(node,velW.Get(node) - test);
//    //wSolOld.Set(node,wSolOld.Get(node) + test);
//   }
//  }
//-------------------------------------------------- 

 SemiLagrangean sl(*m,uSolOld,vSolOld,wSolOld,velU,velV,velW,cSolOld);

 sl.compute(dt);
 uSL = *sl.getUSL();
 vSL = *sl.getVSL();
 wSL = *sl.getWSL();

 convUVW.CopyFrom(0,uSL);
 convUVW.CopyFrom(numNodes,vSL);
 convUVW.CopyFrom(2*numNodes,wSL);
 convC = *sl.getCSL();

//--------------------------------------------------
//  clVector convSurface = sl.computeSurface(dt);
//  convSurface.CopyTo(0,uSLSurface);
//  convSurface.CopyTo(numNodes,vSLSurface);
//  convSurface.CopyTo(2*numNodes,wSLSurface);
//-------------------------------------------------- 
 
} // fecha metodo stepSL 

void Simulator3D::stepSLPBCFix()
{
   clVector velU = uSolOld-uALE;
   clVector velV = vSolOld-vALE;
   clVector velW = wSolOld-wALE;
 
   SemiLagrangean slp (*m,uSolOld,vSolOld,wSolOld,velU,velV,velW,cSolOld);

   slp.computePBCFix(dt);
   uSL = *slp.getUSL();
   vSL = *slp.getVSL();
   wSL = *slp.getWSL();

   convUVW.CopyFrom(0,uSL);
   convUVW.CopyFrom(numNodes,vSL);
   convUVW.CopyFrom(2*numNodes,wSL);
   convC = *slp.getCSL();
 
} // fecha metodo stepSLPBCFix 

void Simulator3D::stepNoConvection()
{
 convUVW.CopyFrom(0,uSolOld);
 convUVW.CopyFrom(numNodes,vSolOld);
 convUVW.CopyFrom(2*numNodes,wSolOld);
 convC = cSolOld;
} // fecha metodo stepNoConvection

/* 
 * This method should be used only to move the mesh, thus Navier-Stokes
 * equations will not be solved. To do so, mainVortex.cpp is an example
 * of how-to.

 * input: X, Y and Z coordinates
 * output: Sol velocity field
 * */
void Simulator3D::stepImposedPeriodicField(const char* _name,
                                           double _T,
										   double _time)
{
 double aux;
 double pi = 3.14159265358;
 /*
  * Reference:
  *
  * */
 if( strcmp( _name,"2d") == 0 || 
     strcmp( _name,"2D") == 0 )
 {
  for( int i=0;i<numVerts;i++ )
  {
   aux = (-1.0)*sin(pi*X->Get(i))*
	            sin(pi*X->Get(i))*
				sin(2*pi*Y->Get(i))*
				cos(_time/_T);
   uSol.Set(i,aux);
   aux = sin(pi*Y->Get(i))*
	     sin(pi*Y->Get(i))*
		 sin(2*pi*X->Get(i))*
		 cos(_time/_T);
   vSol.Set(i,aux);
   aux = 0.0;
   wSol.Set(i,aux);
  }
 }
 /* Front tracking with moving-least-squares surfaces
  * Joao Paulo Gois, Anderson Nakano, Luis Gustavo Nonato, Gustavo C.
  * Buscaglia
  * */
 else if( strcmp( _name,"3d") == 0 || 
          strcmp( _name,"3D") == 0 ) 
 {
  for( int i=0;i<numVerts;i++ )
  {
   aux = (2.0)*sin(pi*X->Get(i))*
	           sin(pi*X->Get(i))*
			   sin(2*pi*Y->Get(i))*
			   sin(2*pi*Z->Get(i))*
			   cos(pi*_time/_T);
   uSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*X->Get(i))*
	            sin(pi*Y->Get(i))*
				sin(pi*Y->Get(i))*
				sin(2*pi*Z->Get(i))*
				cos(pi*_time/_T);
   vSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*X->Get(i))*
	            sin(2*pi*Y->Get(i))*
				sin(pi*Z->Get(i))*
				sin(pi*Z->Get(i))*
				cos(pi*_time/_T);
   wSol.Set(i,aux);
  }
 }
 /* A simple package for front tracking
  * Jian Du, Brian Fix, James Glimm, Xicheng Jia, Xiaolin Li, Yuanhua
  * Li, Lingling Wu
  * */
 else if( strcmp( _name,"shear3d") == 0 || 
          strcmp( _name,"shear3D") == 0 ) 
 {
  double R = 0.5;
  double x0 = 0.5;
  double y0 = 0.5;
  for( int i=0;i<numVerts;i++ )
  {
   aux = sin(pi*X->Get(i))*
		 sin(pi*X->Get(i))*
		 sin(2.0*pi*Y->Get(i))*
		 cos(pi*_time/_T);
   uSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*X->Get(i))*
	            sin(pi*Y->Get(i))*
				sin(pi*Y->Get(i))*
				cos(pi*_time/_T);
   vSol.Set(i,aux);
   double r0 = sqrt( (X->Get(i)-x0)*(X->Get(i)-x0)+
	               (Y->Get(i)-y0)*(Y->Get(i)-y0) );
   aux = ( 1.0-r0/R )*( 1.0-r0/R )*cos(pi*_time/_T);
   wSol.Set(i,aux);
  }
 }
 else if( strcmp( _name,"one") == 0 || 
          strcmp( _name,"One") == 0 ) 
 {
  for( int i=0;i<numVerts;i++ )
  {
   aux = 1.0*cos(pi*_time/_T);

   uSol.Set(i,aux);
   vSol.Set(i,aux);
   wSol.Set(i,aux);
  }
 }
 else if( strcmp( _name,"rotating") == 0 || 
          strcmp( _name,"Rotating") == 0 ) 
 {
  for( int i=0;i<numNodes;i++ )
  {
   double omega=1.0;

   aux = (-1.0)*Y->Get(i)*omega;
   uSol.Set(i,aux);
   aux = X->Get(i)*omega;
   vSol.Set(i,aux);
   aux = 0.0;
   wSol.Set(i,aux);
  }
 }
 else
 {
  cerr << "Periodic field not defined!" << endl;
  exit(1);
 }
} // fecha metodo stepImposedPeriodicField

/* 1/2 time step is computed from ALE velocity field
 * 
 * input:  SolOld velocity
 * output: SolOld velocity
 * */
void Simulator3D::stepImposedPeriodicField(const char* _name,
                                           double _T,
										   double _time,
										   double _dt)
{
 double aux;
 double pi = 3.14159265358;
 /*
  * Reference:
  *
  * */
 if( strcmp( _name,"2d") == 0 || 
     strcmp( _name,"2D") == 0 )
 {
  for( int i=0;i<numVerts;i++ )
  {
   double Xp = X->Get(i)+(uSolOld.Get(i)*_dt);
   double Yp = Y->Get(i)+(vSolOld.Get(i)*_dt);
   aux = (-1.0)*sin(pi*Xp)*
	            sin(pi*Xp)*
				sin(2*pi*Yp)*
				cos(_time/_T);
   uSol.Set(i,aux);
   aux = sin(pi*Yp)*
	     sin(pi*Yp)*
		 sin(2*pi*Xp)*
		 cos(_time/_T);
   vSol.Set(i,aux);
   aux = 0.0;
   wSol.Set(i,aux);
  }
 }
 /* Front tracking with moving-least-squares surfaces
  * Joao Paulo Gois, Anderson Nakano, Luis Gustavo Nonato, Gustavo C.
  * Buscaglia
  * */
 else if( strcmp( _name,"3d") == 0 || 
          strcmp( _name,"3D") == 0 ) 
 {
  for( int i=0;i<numVerts;i++ )
  {
   double Xp = X->Get(i)+(uSolOld.Get(i)*_dt);
   double Yp = Y->Get(i)+(vSolOld.Get(i)*_dt);
   double Zp = Z->Get(i)+(wSolOld.Get(i)*_dt);
   aux = (2.0)*sin(pi*Xp)*
	           sin(pi*Xp)*
			   sin(2*pi*Yp)*
			   sin(2*pi*Zp)*
			   cos(pi*_time/_T);
   uSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*Xp)*
	            sin(pi*Yp)*
				sin(pi*Yp)*
				sin(2*pi*Zp)*
				cos(pi*_time/_T);
   vSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*Xp)*
	            sin(2*pi*Yp)*
				sin(pi*Zp)*
				sin(pi*Zp)*
				cos(pi*_time/_T);
   wSol.Set(i,aux);
  }
 }
 /* A simple package for front tracking
  * Jian Du, Brian Fix, James Glimm, Xicheng Jia, Xiaolin Li, Yuanhua
  * Li, Lingling Wu
  * */
 else if( strcmp( _name,"shear3d") == 0 || 
          strcmp( _name,"shear3D") == 0 ) 
 {
  double R = 0.5;
  double x0 = 0.5;
  double y0 = 0.5;
  for( int i=0;i<numVerts;i++ )
  {
   double Xp = X->Get(i)+(uSolOld.Get(i)*dt/2.0);
   double Yp = Y->Get(i)+(vSolOld.Get(i)*dt/2.0);
   //double Zp = Z->Get(i)+(wSolOld.Get(i)*dt/2.0);
   aux = sin(pi*Xp)*
		 sin(pi*Xp)*
		 sin(2.0*pi*Yp)*
		 cos(pi*_time/_T);
   uSol.Set(i,aux);
   aux = (-1.0)*sin(2*pi*Xp)*
	            sin(pi*Yp)*
				sin(pi*Yp)*
				cos(pi*_time/_T);
   vSol.Set(i,aux);
   double r0 = sqrt( (Xp-x0)*(Xp-x0)+
	               (Yp-y0)*(Yp-y0) );
   aux = ( 1.0-r0/R )*( 1.0-r0/R )*cos(pi*_time/_T);
   wSol.Set(i,aux);
  }
 }
 else if( strcmp( _name,"one") == 0 || 
          strcmp( _name,"One") == 0 ) 
 {
  for( int i=0;i<numVerts;i++ )
  {
   //double Xp = X->Get(i)+(uSolOld.Get(i)*dt/2.0);
   //double Yp = Y->Get(i)+(vSolOld.Get(i)*dt/2.0);
   //double Zp = Z->Get(i)+(wSolOld.Get(i)*dt/2.0);
   aux = 1.0*cos(pi*_time/_T);

   uSol.Set(i,aux);
   vSol.Set(i,aux);
   wSol.Set(i,aux);
  }
 }
 else if( strcmp( _name,"rotating") == 0 || 
          strcmp( _name,"Rotating") == 0 ) 
 {
  for( int i=0;i<numNodes;i++ )
  {
   double Xp = X->Get(i)+(uSolOld.Get(i)*dt/2.0);
   double Yp = Y->Get(i)+(vSolOld.Get(i)*dt/2.0);
   //double Zp = Z->Get(i)+(wSolOld.Get(i)*dt/2.0);
   double omega=1.0;

   aux = (-1.0)*Yp*omega;
   uSol.Set(i,aux);
   aux = Xp*omega;
   vSol.Set(i,aux);
   aux = 0.0;
   wSol.Set(i,aux);
  }
 }
 else
 {
  cerr << "Periodic field not defined!" << endl;
  exit(1);
 }
} // fecha metodo stepImposedPeriodicField

void Simulator3D::copyALEtoSol()
{
 uSolOld = uSolOld-uALE;
 vSolOld = vSolOld-vALE;
 wSolOld = wSolOld-wALE;
}

void Simulator3D::step()
{
 Galerkin galerkin(*m,uSolOld,vSolOld,wSolOld,cSolOld,gx,gy,gz);

 galerkin.compute(dt);
 clVector uGalerkin = *galerkin.getConvU();
 clVector vGalerkin = *galerkin.getConvV();
 clVector wGalerkin = *galerkin.getConvW();
 convUVW.CopyFrom(0,uGalerkin);
 convUVW.CopyFrom(numNodes,vGalerkin);
 convUVW.CopyFrom(2*numNodes,wGalerkin);
 convC = *galerkin.getConvC();

 clVector uvwSol(0);
 uvwSol.Append(uSolOld);
 uvwSol.Append(vSolOld);
 uvwSol.Append(wSolOld);
} // fecha metodo step 

// metodo para movimentacao dos pontos da malha atraves da velocidade do
// escoamento (uSol, vSol e wSol), caracterizando a movimentacao
// puramente lagrangiana
void Simulator3D::stepLagrangian()
{
 // impoe velocidade SolOld = 0 no contorno
 setLagrangianVelBC();

 convUVW.CopyFrom(0,uSolOld);
 convUVW.CopyFrom(numNodes,vSolOld);
 convUVW.CopyFrom(2*numNodes,wSolOld);
 convC = cSolOld;

 m->moveXPoints(uSolOld,dt);
 m->moveYPoints(vSolOld,dt);
 m->moveZPoints(wSolOld,dt);

} // fecha metodo stepLagrangian

// metodo para movimentacao dos pontos da malha na direcao Z atraves da 
// velocidade do escoamento (wSol), caracterizando a movimentacao
// puramente lagrangiana em Z. Para os pontos nas direcoes X e Y eh
// utilizado o metodo explicito semi lagrangiano
void Simulator3D::stepLagrangianZ()
{
 m->moveZPoints(wSolOld,dt);
 m->centroidPositionCorrection();

 SemiLagrangean sl(*m,uSolOld,vSolOld,wSolOld,velU,velV,velW,cSolOld);

 sl.computeFreeSurface(dt);
 uSL = *sl.getUSL();
 vSL = *sl.getVSL();
 //wSL = *sl.getWSL();
 convC = *sl.getCSL();

 convUVW.CopyFrom(0,uSL);
 convUVW.CopyFrom(numNodes,vSL);
 convUVW.CopyFrom(2*numNodes,wSol);
} // fecha metodo stepLagragianZ

/* compute ALE velocity according to mesh parameters c1,c2,c3,d1 and d3
 * 
 * input:  SolOld velocity
 * output: ALE velocity
 * */
void Simulator3D::stepALE()
{
 // vertice velocity (uVert,vVert)
 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 uSolOld.CopyTo(0,uVert);
 vSolOld.CopyTo(0,vVert);
 wSolOld.CopyTo(0,wVert);

 setInterfaceVelocity();
 //setMassTransfer();

 if( c2 > 0 )
 {
  // smoothing - velocidade
  MeshSmooth e1(*m,dt); // criando objeto MeshSmooth
  e1.stepSmooth(uALE,vALE,wALE);
  uSmooth = *e1.getUSmooth();
  vSmooth = *e1.getVSmooth();
  wSmooth = *e1.getWSmooth();

  uSmooth=uALE;
  vSmooth=vALE;
  wSmooth=wALE;
  MeshSmooth e3(*m,dt); // criando objeto MeshSmooth
  for( int i=0;i<20;i++ )
  {
   // smoothing - velocidade
   e3.stepSmoothLonger(uSmooth,vSmooth,wSmooth);
   uSmooth = *e3.getUSmooth();
   vSmooth = *e3.getVSmooth();
   wSmooth = *e3.getWSmooth();

   for( int i=0;i<surfMesh->numVerts;i++ )
   {
	uSmooth.Set(i,uALE.Get(i));
	vSmooth.Set(i,vALE.Get(i));
	wSmooth.Set(i,wALE.Get(i));
   }
  }
 }

 if( c3 > 0.0 )
 {
  // smoothing coords
  MeshSmooth e2(*m,1.0); // criando objeto MeshSmooth
  e2.stepSmoothFujiwara();
  uSmoothCoord = *e2.getUSmooth();
  vSmoothCoord = *e2.getVSmooth();
  wSmoothCoord = *e2.getWSmooth();
 }

 // compute ALE
 uALE = c1*uVert+c2*uSmooth+c3*uSmoothCoord;
 vALE = c1*vVert+c2*vSmooth+c3*vSmoothCoord;
 wALE = c1*wVert+c2*wSmooth+c3*wSmoothCoord;
 
 // impoe velocidade do fluido na interface
 setInterfaceVelocity();
 //setMassTransfer();

 clVector zeroAux(numNodes-numVerts);
 uALE.Append(zeroAux);
 vALE.Append(zeroAux);
 wALE.Append(zeroAux);
#if NUMGLEU == 5
 uALE = setTetCentroid(*IEN,uALE);
 vALE = setTetCentroid(*IEN,vALE);
 wALE = setTetCentroid(*IEN,wALE);
#else
 uALE = setTetQuad(*IEN,uALE);
 vALE = setTetQuad(*IEN,vALE);
 wALE = setTetQuad(*IEN,wALE);
#endif

 // impoe velocidade ALE = 0 no contorno
 setALEBC();
 //setAnnularALEBC();

 // calcula velocidade do fluido atraves do metodo semi-lagrangeano
 // comment if using mainVortex.cpp
 //stepSL();
 stepSLPBCFix(); // for PBC
} // fecha metodo stepALE

/* move nodes according to ALE velocity 
 * 
 * input:  ALE velocity
 * output: ----
 * */
void Simulator3D::movePoints()
{
 movePoints(&uALE,&vALE,&wALE);
}

/* move nodes according to _?vel velocity, which are passed by the user
 * 
 * input:  _?Vel velocity
 * output: ----
 * */
void Simulator3D::movePoints(clVector *_uVel,
                             clVector *_vVel,
							 clVector *_wVel)
{
 // movimentando os vertices pontos da malha com velocidade ALE
 m->moveXPoints(*_uVel,dt);
 m->moveYPoints(*_vVel,dt);
 m->moveZPoints(*_wVel,dt);
 m->centroidPositionCorrection();

 // correcao do volume da bolha
 m->applyBubbleVolumeCorrection();
}

void Simulator3D::setInterfaceVelocity()
{
 if( d2 > 0 )
 {
  // smoothing - coordenadas
  MeshSmooth e1(*m,dt); // criando objeto MeshSmooth
  e1.stepSurfaceSmoothFujiwara();
  uSmoothSurface = *e1.getUSmoothSurface();
  vSmoothSurface = *e1.getVSmoothSurface();
  wSmoothSurface = *e1.getWSmoothSurface();
 }

 for( int i=0;i<surface->Dim();i++ )
 {
  int surfaceNode = surface->Get(i);

  // unitario do vetor normal (ponderado com a area) resultante
  double xNormalUnit = surfMesh->xNormal.Get(surfaceNode);
  double yNormalUnit = surfMesh->yNormal.Get(surfaceNode);
  double zNormalUnit = surfMesh->zNormal.Get(surfaceNode);

  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  double prod = (uSolOld.Get(surfaceNode)+1.3*uRef)*xNormalUnit+ 
              (vSolOld.Get(surfaceNode)+1.3*vRef)*yNormalUnit + 
			  (wSolOld.Get(surfaceNode)+1.3*wRef)*zNormalUnit;
  double uSolNormal = xNormalUnit*prod;
  double vSolNormal = yNormalUnit*prod;
  double wSolNormal = zNormalUnit*prod;

  // 1.3 is a pragmatic number which fits the velocity for the bhaga5
  // and the moving frame technique. Still don't know why!
  double uSolTangent = uSolOld.Get(surfaceNode) + 1.3*uRef - uSolNormal;
  double vSolTangent = vSolOld.Get(surfaceNode) + 1.3*vRef - vSolNormal;
  double wSolTangent = wSolOld.Get(surfaceNode) + 1.3*wRef - wSolNormal;

  // tratamento da superficie
  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  double prod2 = uSmoothSurface.Get(surfaceNode)*xNormalUnit + 
               vSmoothSurface.Get(surfaceNode)*yNormalUnit + 
			   wSmoothSurface.Get(surfaceNode)*zNormalUnit;
  double uSmoothNormal = xNormalUnit*prod2;
  double vSmoothNormal = yNormalUnit*prod2;
  double wSmoothNormal = zNormalUnit*prod2;

  double uSmoothTangent = uSmoothSurface.Get(surfaceNode) - uSmoothNormal;
  double vSmoothTangent = vSmoothSurface.Get(surfaceNode) - vSmoothNormal;
  double wSmoothTangent = wSmoothSurface.Get(surfaceNode) - wSmoothNormal;

  double uALESurface =   uSolOld.Get(surfaceNode) 
                     - d1*uSolTangent 
					 + d2*uSmoothTangent;
  double vALESurface =   vSolOld.Get(surfaceNode) 
                     - d1*vSolTangent 
					 + d2*vSmoothTangent;
  double wALESurface =   wSolOld.Get(surfaceNode) 
                     - d1*wSolTangent 
					 + d2*wSmoothTangent;

  uALE.Set(surfaceNode,uALESurface);
  vALE.Set(surfaceNode,vALESurface);
  wALE.Set(surfaceNode,wALESurface);
 }
} // fecha metodo setInterfaceVelocity

//setRHS eh o metodo que cria o vetor do lado direito (condicao de
//contorno + termo convectivo)
void Simulator3D::setRHS()
{
 // sem correcao na pressao
 va = ( (1.0/dt) * Mrho + (1-alpha) * -(1.0/Re) * K ) * convUVW;

 // com correcao na pressao
 //va = ( (1.0/dt) * Mrho - (1-alpha) * (1.0/Re) * K ) * convUVW - (G*pSol);
}

void Simulator3D::setCRHS()
{
 //vcc = ( (1.0/dt) * Mc - (1-alpha) * (1.0/(Re*Sc)) * Kc ) * convC;
 //vcc = ( (1.0/dt) * McLumped - (1-alpha) * (1.0/(Re*Sc)) * Kc ) * convC;
 vcc = ( (1.0/dt) * McLumped ) * convC;
}

void Simulator3D::setGravity(const char* _direction)
{
 g = 9.81;
 g_0 = g;
 gAdimen = g/g_0;

 clVector gx(numNodes);
 clVector gy(numNodes);
 clVector gz(numNodes);

 if( strcmp( _direction,"x") == 0 || 
     strcmp( _direction,"+x") == 0 || 
     strcmp( _direction,"X") == 0 || 
     strcmp( _direction,"+X") == 0 ) 
 {
  gx.SetAll(gAdimen);
  gy.SetAll(0.0);
  gz.SetAll(0.0);
 }
 else if( strcmp( _direction,"y") == 0 || 
          strcmp( _direction,"+y") == 0 || 
		  strcmp( _direction,"Y") == 0 || 
		  strcmp( _direction,"+Y") == 0 ) 
 {
  gx.SetAll(0.0);
  gy.SetAll(gAdimen);
  gz.SetAll(0.0);
 }
 else if( strcmp( _direction,"z") == 0 || 
          strcmp( _direction,"+z") == 0 || 
		  strcmp( _direction,"Z") == 0 ||
		  strcmp( _direction,"+Z") == 0 ) 
 {
  gx.SetAll(0.0);
  gy.SetAll(0.0);
  gz.SetAll(gAdimen);
 }
 else if( strcmp( _direction,"-x") == 0 || 
          strcmp( _direction,"-X") == 0 ) 
 {
  gx.SetAll(-1.0*gAdimen);
  gy.SetAll(0.0);
  gz.SetAll(0.0);
 }
 else if( strcmp( _direction,"-y") == 0 || 
          strcmp( _direction,"-Y") == 0 ) 
 {
  gx.SetAll(0.0);
  gy.SetAll(-1.0*gAdimen);
  gz.SetAll(0.0);
 }
 else if( strcmp( _direction,"-z") == 0 || 
          strcmp( _direction,"-Z") == 0 ) 
 {
  gx.SetAll(0.0);
  gy.SetAll(0.0);
  gz.SetAll(-1.0*gAdimen);
 }
 else 
 {
  gx.SetAll(0.0);
  gy.SetAll(0.0);
  gz.SetAll(0.0);
 }

 gravity.Dim(0);
 gravity.Append(gx);
 gravity.Append(gy);
 gravity.Append(gz);

 //gravity = -( 1.0/(Fr*Fr) )*( Mrho*gUnit );
 //gravity = -( 1.0/(Fr*Fr) )*( (Mrho - M)*gUnit );

 // update RHS
 //va = va + gravity;
}


/* Sets dimensionless pressure gradient for PBC applications, which
 * acts on the flow in the streamwise direction only. Its construction
 * is similar to the vector "gravity", except for the value set as
 * input not always unitary.
 *
 * \param[in]: streamwise direction which the periodicity applies to.  
 * \param[out]: pressure gradient 
 * */
void Simulator3D::setBetaFlowLiq(const char* _direction)
{
 betaPressLiqAdimen = betaPressLiq;

 clVector px(numNodes);
 clVector py(numNodes);
 clVector pz(numNodes);

 if( strcmp( _direction,"x") == 0 || 
     strcmp( _direction,"+x") == 0 || 
     strcmp( _direction,"X") == 0 || 
     strcmp( _direction,"+X") == 0 ) 
 {
  px.SetAll(betaPressLiqAdimen);
  py.SetAll(0.0);
  pz.SetAll(0.0);
 }
 else if( strcmp( _direction,"y") == 0 || 
     strcmp( _direction,"+y") == 0 || 
     strcmp( _direction,"Y") == 0 || 
     strcmp( _direction,"+Y") == 0 ) 
 {
  px.SetAll(0.0);
  py.SetAll(betaPressLiqAdimen);
  pz.SetAll(0.0);
 }
 else if( strcmp( _direction,"z") == 0 || 
     strcmp( _direction,"+z") == 0 || 
     strcmp( _direction,"Z") == 0 || 
     strcmp( _direction,"+Z") == 0 ) 
 {
  px.SetAll(0.0);
  py.SetAll(0.0);
  pz.SetAll(betaPressLiqAdimen);
 }
 else if( strcmp( _direction,"-x") == 0 || 
     strcmp( _direction,"-X") == 0 ) 
 {
  px.SetAll(-1.0*betaPressLiqAdimen);
  py.SetAll(0.0);
  pz.SetAll(0.0);
 }
 else if( strcmp( _direction,"-y") == 0 || 
     strcmp( _direction,"-Y") == 0 ) 
 {
  px.SetAll(0.0);
  py.SetAll(-1.0*betaPressLiqAdimen);
  pz.SetAll(0.0);  
 }
 else if( strcmp( _direction,"-z") == 0 || 
     strcmp( _direction,"-Z") == 0 ) 
 {
  px.SetAll(0.0);
  py.SetAll(0.0);
  pz.SetAll(-1.0*betaPressLiqAdimen);
 }
 else 
 {
  px.SetAll(0.0);
  py.SetAll(0.0);
  pz.SetAll(0.0);
 }

 betaFlowLiq.Dim(0);
 betaFlowLiq.Append(px);
 betaFlowLiq.Append(py);
 betaFlowLiq.Append(pz);

}


void Simulator3D::setGravityBoussinesq(const char* _direction)
{
 clVector zeroAux(numNodes-numVerts);
 clVector forceC(numNodes);

 if( strcmp( _direction,"x") == 0 || strcmp( _direction,"X") == 0 ) 
 {
  forceC.Append(convC); 
  forceC.Append(zeroAux);
  forceC.Append(zeroAux);
 }
 if( strcmp( _direction,"y") == 0 || strcmp( _direction,"Y") == 0 ) 
 {
  forceC.Append(zeroAux);
  forceC.Append(convC); 
  forceC.Append(zeroAux);
 }
 if( strcmp( _direction,"z") == 0 || strcmp( _direction,"Z") == 0 ) 
 {
  forceC.Append(zeroAux);
  forceC.Append(zeroAux);
  forceC.Append(convC); 
 }

 gravity =  beta*( Mrho*forceC );

 // update RHS
 va = va + gravity;
}

/* Set kappa vector according to Model3D's curvature vector:
 *   kappa = curvature
 * Compyute Surface Tension Force fint as: 
 * fint = kappa*GH
 * 
 * input: X, Y and Z coordinates
 * output: Sol velocity field
 * */
void Simulator3D::setInterfaceGeo()
{
 m->setNormalAndKappa();
 m->setKappaSurface();
 kappa = *m->getCurvature();

 fint = kappa*(G*(*heaviside));
 //fint = (1.0/We) * ( kappa*(G*(*heaviside)) );
 //fint = (1.0/We) * ( kappa*(GTilde*(*heaviside)) );
 
 //va = va + fint;
} // fecha metodo setInterface 

void Simulator3D::setInterfaceLevelSet()
{
 clVector half(numVerts);half.SetAll(0.5);
 clVector zeroLevel = ((*heaviside)-half)*2;
 clVector levelSet(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  double aux = zeroLevel.Get(i)*interfaceDistance->Get(i);
  levelSet.Set(i,aux);
 }

 // ----------- smoothing step ------------ //
 matMountC();
 setUnCoupledCBC();
 vcc = ( (1.0/dt) * McLumped ) * levelSet;
 unCoupledC();
 //clVector kappaAux = invMcLumped*(Kc * vcc);
 //clVector kappaAux = invMcLumped*(Kc * levelSet);
 clVector kappaAux = invMcLumped*(Kc * cTilde);
 // --------------------------------------- //

 m->setKappaSurface(kappaAux);
 kappa = *m->getCurvature();

 fint = (1.0/We) * ( kappa*(G*(*heaviside)) );
 //fint = (1.0/We) * ( kappa*(GTilde*(*heaviside)) );

 //va = va + invA*fint;
} // fecha metodo setInterface 

void Simulator3D::matMount()
{
 mat = ( (1.0/dt) * Mrho) + (alpha * (1.0/Re) * K );

 for( int i = 0; i < 3*numNodes; i++ )
 {
  double sumMat = mat.SumLine(i);
  invA.Set( i,1.0/sumMat );

  double sumM = M.SumLine(i);
  //MLumped.Set( i,sumM );
  invMLumped.Set( i,1.0/sumM );

  double sumMrho = Mrho.SumLine(i);
  //MrhoLumped.Set( i,sumMrho );
  invMrhoLumped.Set( i,1.0/sumMrho );
 }
}

void Simulator3D::matMountC()
{
 // set para Matriz Lumped da concentracao
 for( int i=0;i<numVerts;i++ )
 {
  double sumMc = Mc.SumLine(i);
  McLumped.Set( i,sumMc );
  invMcLumped.Set( i,1.0/sumMc );
 }

 //matc = ((1.0/dt) * Mc) + (alpha * (1.0/(Sc*Re)) * Kc);
 matc = ((1.0/dt) * McLumped) + (alpha * (1.0/(Sc*Re)) * Kc);

 for( int i=0;i<numVerts;i++ )
 {
  double sumMatc = matc.SumLine(i);
  invC.Set( i,1.0/sumMatc );
 }
}

void Simulator3D::coupled()
{
 // A=[K G]=[matU  0 Gx][ ]=[us]  [us]=[ ]
 //   [D 0] [ 0 matV Gy][u]=[vs]  [vs]=[b]
 //         [Dx  Dy   0][ ]=[ps]  [ps]=[ ]

 A.CopyFrom(         0,         0,mat);
 A.CopyFrom(         0,3*numNodes,G  );
 A.CopyFrom(3*numNodes,         0,D  );

 clVector zeros(numVerts);

 b = va + (1.0/We)*fint + (1.0/(Fr*Fr)*gravity);
 b.Append(zeros);
 setCoupledBC();

 clVector u(3*numNodes+numVerts);

 cout << " --> solving velocity and pressure -- " << endl;
 solverV->solve(1E-15,A,u,b);
 cout << " ------------------------------------ " << endl;
 
 u.CopyTo(0,uSol);
 u.CopyTo(numNodes,vSol);
 u.CopyTo(numNodes*2,wSol);
 u.CopyTo(numNodes*3,pSol);
} // fecha metodo coupled 

void Simulator3D::unCoupled()
{
 clVector uvw(3*numNodes);
 clVector vaIp(3*numNodes);
 clVector b1Tilde;
 clVector b2Tilde;

 // applying boundary conditions to b1Tilde
 vaIp = va.MultVec(ip); // operacao vetor * vetor (elemento a elemento)

 b1Tilde = b1 + vaIp;

 // resolve sitema ATilde uTilde = b
 cout << " --------> solving velocity --------- " << endl;
 solverV->solve(1E-15,ATilde,uTilde,b1Tilde);
 cout << " ------------------------------------ " << endl;

 // uvw = uTilde + surface tension + gravity
 if( rho_in <= rho_out ) // BUBBLE
 {
  cout << setw(70) << "BUBBLE SIMULATION" << endl;
  uvw = uTilde + 
        dt*invMrhoLumped*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		dt*invMLumped*( ((1.0/We)*(fint) ).MultVec(ip) );

 }
 else // DROPLET
 {
  cout << setw(70) << "DROPLET SIMULATION" << endl;
  uvw = uTilde + 
        invA*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		invA*( ((1.0/We)*(fint) ).MultVec(ip) );
 }
 
 //uvw = uTilde;


 /* 
  * Mass Transfer
  * */
 //clVector massTransfer = dt*invMcLumped*( heatFlux );
 //clVector massTransfer = ( invMcLumped*heatFlux );
 clVector massTransfer = invC*
                         (1.0/rho_inAdimen - 1.0/rho_outAdimen)
						 *heatFlux;

 b2Tilde = (-1.0)*( b2 - (DTilde * uvw) ); 
 //b2Tilde = (-1.0)*( b2 - (DTilde * uvw) + (massTransfer) );


 // resolve sistema E pTilde = b2
 cout << " --------> solving pressure --------- " << endl;
 solverP->solve(1E-15,ETilde,pTilde,b2Tilde);
 cout << " ------------------------------------ " << endl;
 
 //uvw = uTilde - dt*(invMrhoLumped * GTilde * pTilde);
 uvw = uvw - (invA * GTilde * pTilde);

 uvw.CopyTo(         0,uSol);
 uvw.CopyTo(  numNodes,vSol);
 uvw.CopyTo(numNodes*2,wSol);

 pSol = pTilde;       // sem correcao na pressao
 //pSol = pSol + pTilde;  // com correcao na pressao

 // compute bubble's centroid velocity
 if( surfMesh->numInterfaces > 0 )
  setCentroidVelPos();
 //setCentroidVelPosInterface();
 // NOT WORKING!
 //setCentroidVelPosInterface();
//--------------------------------------------------
//  for ( int i=0;i<numVerts;i++ )
//  {
//   int node = i;
//   uSol.Set(node,uSol.Get(node) - getCentroidVelXMax());
//  }
//-------------------------------------------------- 
} // fecha metodo unCoupled 

void Simulator3D::unCoupledC()
{
 clVector vcIp(numVerts);
 clVector b1cTilde;

 vcIp = vcc.MultVec(ipc); // operacao vetor * vetor (elemento a elemento)

 b1cTilde = b1c + vcIp;

 b1cTilde = b1cTilde + 
        //dt*invMcLumped*( ((1.0/(Re*Sc))*( heatFlux )).MultVec(ipc) );  
        invC*( ((1.0/(Re*Sc))*( heatFlux )).MultVec(ipc) );  

 // resolve sitema ATilde uTilde = b
 cout << endl;
 cout << " --------> solving scalar ----------- " << endl;
 solverC->solve(1E-15,AcTilde,cTilde,b1cTilde);
 cout << " ------------------------------------ " << endl;
 cout << endl;

 // comentar em caso de utilizacao de setInterface()
 // pois cSol nao pode ser atualizado
 cSol = cTilde;
}

void Simulator3D::unCoupledCPBC()
{
 clVector vcIp(numVerts);
 clVector b1cTilde;
 Periodic3D pbc;

 vcIp = vcc.MultVec(ipc); // operacao vetor * vetor (elemento a elemento)
 b1cTilde = b1c + vcIp; //<<< Por que duplicado??

 b1cTilde = b1cTilde + 
        //dt*invMcLumped*( ((1.0/(Re*Sc))*( heatFlux )).MultVec(ipc) );  
        invC*( ((1.0/(Re*Sc))*( heatFlux )).MultVec(ipc) );  

 //*** modifying r.h.s. for scalar
 sumIndexPBCScalar(VecXMin,VecXMax,b1cTilde);

 //*** Reassemble (scalar field)
 assembleCPBC();

 // solving scalar modified system AcTilde cTilde = b1cTilde
 cout << endl;
 cout << " --------> solving scalar ----------- " << endl;
 solverC->solve(1E-15,AcTilde,cTilde,b1cTilde);
 cout << " ------------------------------------ " << endl;
 cout << endl;

 //*** copying scalar 
 pbc.SetPureScalarPBC(cTilde,VecXMinGlob,VecXMaxGlob,nyPointsL,"RL");

 // comentar em caso de utilizacao de setInterface()
 // pois cSol nao pode ser atualizado
 cSol = cTilde;
}



void Simulator3D::saveOldData()
{
 uSolOld     = uSol;
 vSolOld     = vSol;
 wSolOld     = wSol;
 pSolOld     = pSol;
 cSolOld     = cSol;
 uALEOld     = uALE;
 vALEOld     = vALE;
 wALEOld     = wALE;
 fintOld     = fint;
 gravityOld  = gravity;
 betaFlowLiqOld  = betaFlowLiq;
 kappaOld    = kappa;
 muOld       = mu;
 rhoOld      = rho;
 cpOld       = cp;
 ktOld       = kt;
 hSmoothOld  = hSmooth;
 heatFluxOld = heatFluxOld;

 numVertsOld = numVerts;
 numNodesOld = numNodes;
 numElemsOld = numElems;
}

void Simulator3D::timeStep()
{
 time = time + dt;
 iter++;
}

/**
 * @brief ESTRATEGIA PARA APLICACAO DA CONDICAO DE CONTORNO PARA SISTEMA
 * RESOLVIDO DIRETAMENTE
 * 1-localizar o vertice onde ha condicao de contorno
 * 2-somar a coluna onde o vertice se encontra no vetor do lado direto 
 * 3-zerar a coluna e a linha colocando 1 no vertice da matriz
 * 4-repetir a operacao para outros vertices condicao de contorno 
 *
 * @return VOID (VAZIO)
 **/
void Simulator3D::setCoupledBC()
{
 int nbc,i,j;
 clVector vetColSizeUVWP(3*numNodes+numVerts);
 clVector UVWPC = *uc;
 UVWPC.Append(*vc);
 UVWPC.Append(*wc);
 UVWPC.Append(*pc);
 
 nbc = idbcu->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcu->Get(i);
  b.CopyMult(j,A,UVWPC);
  A.Set(j,j,1);
  b.Set(j,uc->Get(j));
 }
 cout << endl;
 cout << " boundary condition ";
 cout << color(none,red,black) << "U ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcv->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcv->Get(i);
  b.CopyMult(j+numNodes,A,UVWPC);
  A.Set( j+numNodes, j+numNodes, 1 );
  b.Set( j+numNodes, vc->Get(j) );
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "V ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcw->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcw->Get(i);
  b.CopyMult(j+numNodes*2,A,UVWPC);
  A.Set( j+numNodes*2,j+numNodes*2, 1);
  b.Set( j+numNodes*2,wc->Get(j) );
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "W ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcp->Dim();
 for( i=0;i<nbc;i++ )
 {
  j = (int) idbcp->Get(i);
  b.CopyMult(j+numNodes*3,A,UVWPC);
  A.Set( j+3*numNodes,j+3*numNodes, -1 );
  b.Set( j+3*numNodes, -pc->Get(j) ); // sem correcao na pressao
  //b.Set( j+3*numNodes,-pc->Get(j)*0 ); // com correcao na pressao
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "P ";
 cout << resetColor() << "--> SET " << endl;
 cout << endl;
 
} // fecha metodo setCoupledBC 

void Simulator3D::setUnCoupledBC()
{
 int nbc,i,j;
 // ----------- nova versao ------------ //
 clVector UVWC = *uc;
 UVWC.Append(*vc);
 UVWC.Append(*wc);
 // ------------------------------------ //
 
 ATilde = mat;
 DTilde = D;
 GTilde = G;

 b1.Dim(3*numNodes,0);  // zerando os vetores b1 e b2
 b2.Dim(numVerts,0);
 ip.Dim(3*numNodes,1); // inicializando vetor ip com 1

 nbc = idbcu->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcu->Get(i);
  b1.CopyMult(j,ATilde,UVWC);
  b2.CopyMult(j,DTilde,GTilde,UVWC);
  ATilde.Set(j,j,1);
  b1.Set(j,uc->Get(j));
  ip.Set(j,0);
 }
 cout << endl;
 cout << " boundary condition ";
 cout << color(none,red,black) << "U ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcv->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcv->Get(i);
  b1.CopyMult(j+numNodes,ATilde,UVWC);
  b2.CopyMult(j+numNodes,DTilde,GTilde,UVWC);
  ATilde.Set(j+numNodes,j+numNodes,1);
  b1.Set(j+numNodes,vc->Get(j));
  ip.Set(j+numNodes,0);
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "V ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcw->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcw->Get(i);
  b1.CopyMult(j+numNodes*2,ATilde,UVWC);
  b2.CopyMult(j+numNodes*2,DTilde,GTilde,UVWC);
  ATilde.Set(j+numNodes*2,j+numNodes*2,1);
  b1.Set(j+numNodes*2,wc->Get(j));
  ip.Set(j+numNodes*2,0);
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "W ";
 cout << resetColor() << "--> SET " << endl;

 E.Dim(numVerts,numVerts);
 nbc = idbcp->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcp->Get(i);
  b1.CopyMult(j,GTilde,DTilde,*pc);
  E.Set(j,j,-1);
  b2.Set(j,-pc->Get(j));  // sem correcao na pressao
  //b2.Set(j,-pc->Get(j)*0);  // com correcao na pressao
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "P ";
 cout << resetColor() << "--> SET " << endl;
 cout << endl;

 ETilde = E - ((DTilde * invA) * GTilde); 
 //ETilde = E - dt*((DTilde * invMrhoLumped) * GTilde); 
} // fecha metodo setUnCoupledBC 


/* Idem explicacao para setUnCoupledBC(), exceto por chamada final
 * que impoe ETilde = E. */
void Simulator3D::setUnCoupledPBC()
{
 int nbc,i,j;
 // ----------- nova versao ------------ //
 clVector UVWC = *uc;
 UVWC.Append(*vc);
 UVWC.Append(*wc);
 // ------------------------------------ //
 
 ATilde = mat;
 DTilde = D;
 GTilde = G;

 b1.Dim(3*numNodes,0);  // zerando os vetores b1 e b2
 b2.Dim(numVerts,0);
 ip.Dim(3*numNodes,1); // inicializando vetor ip com 1

 nbc = idbcu->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcu->Get(i);
  b1.CopyMult(j,ATilde,UVWC);
  b2.CopyMult(j,DTilde,GTilde,UVWC);
  ATilde.Set(j,j,1);
  b1.Set(j,uc->Get(j));
  ip.Set(j,0);
 }
 cout << endl;
 cout << " boundary condition ";
 cout << color(none,red,black) << "U ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcv->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcv->Get(i);
  b1.CopyMult(j+numNodes,ATilde,UVWC);
  b2.CopyMult(j+numNodes,DTilde,GTilde,UVWC);
  ATilde.Set(j+numNodes,j+numNodes,1);
  b1.Set(j+numNodes,vc->Get(j));
  ip.Set(j+numNodes,0);
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "V ";
 cout << resetColor() << "--> SET " << endl;

 nbc = idbcw->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcw->Get(i);
  b1.CopyMult(j+numNodes*2,ATilde,UVWC);
  b2.CopyMult(j+numNodes*2,DTilde,GTilde,UVWC);
  ATilde.Set(j+numNodes*2,j+numNodes*2,1);
  b1.Set(j+numNodes*2,wc->Get(j));
  ip.Set(j+numNodes*2,0);
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "W ";
 cout << resetColor() << "--> SET " << endl;

 E.Dim(numVerts,numVerts);
 nbc = idbcp->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcp->Get(i);
  b1.CopyMult(j,GTilde,DTilde,*pc);
  E.Set(j,j,-1);
  b2.Set(j,-pc->Get(j));  // sem correcao na pressao
  //b2.Set(j,-pc->Get(j)*0);  // com correcao na pressao
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "P ";
 cout << resetColor() << "--> SET " << endl;
 cout << endl;

 ETilde = E; 
} // fecha metodo setUnCoupledPBC 


void Simulator3D::setUnCoupledCBC()
{
 int nbc,i,j;
 b1c.Dim(numVerts,0);  // zerando o vetore b1c
 ipc.Dim(numVerts,1); // inicializando vetor ipc com 1
 
 AcTilde = matc;

 nbc = idbcc->Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcc->Get(i);
  b1c.CopyMult(j,AcTilde,*cc);
  AcTilde.Set(j,j,1);
  b1c.Set(j,cc->Get(j));
  ipc.Set(j,0);
 }
 cout << " boundary condition ";
 cout << color(none,red,black) << "C ";
 cout << resetColor() << "--> SET " << endl;
 cout << endl;

} // fecha metodo setUnCoupledCBC 

/* cfl based on the booklet:
 * Direct Numerical Simulation of Gas-Liquid Multiphase Flows
 * G. Tryggvason, R. Scardovelli and S. Zaleski
 *
 * - cfl of advection term:
 *
 *  Umax*dt
 * --------- <= 1
 *  triEdge
 *  
 *
 * - smallest resolved wave:
 *         pi
 * k = ---------
 *      triEdge
 *
 * - phase velocity of a capillary wave:
 * 
 *            sigma*kappa
 * c = sqrt( ------------- )
 *             rho1+rho2
 *
 *          (rho1+rho2)*triEdge*triEdge*triEdge
 * dt <= ( ------------------------------------- )
 *                       pi*sigma
 *
 * cfl based on the paper:
 * A Continuum Method for Modeling Surface Tension
 * J.U. Brackbill, D.B. Kothe and C. Zemach
 *
 * - cfl of advection term:
 *
 *  Umax*dt      1
 * --------- <= ---
 *  minEdge      2
 *
 *          0.5*(rho1+rho2)*minEdgeTri*minEdgeTri*minEdgeTri
 * dt <= ( -------------------------------------------------- )
 *                            2*pi*sigma
 * */
void Simulator3D::setDtSurfaceTension()
{
 if( rho_in >= rho_out )
  rho_0 = rho_in; 
 else
  rho_0 = rho_out; 

 rho_inAdimen = rho_in/rho_0; 
 rho_outAdimen = rho_out/rho_0;

 vector<double> minLength = m->getMinLength();
 double minEdge = *min_element(minLength.begin()+1,minLength.end());

 dtSurfaceTension = sqrt( ( We*0.5*(rho_inAdimen+rho_outAdimen)*
                    minEdge*minEdge*minEdge)/
                    (2*3.141592) );
}

/*
 *
 *                                 v1              
 *                                  o              ---
 *                                 / \                |
 *                                /   \               | y1
 *                               /     \              |
 *                   -------    /   x   \     --------
 *               y2 |          /         \            |
 *                  |         /           \           | y3
 *                   -----   o ----------- o       ---
 *                          v2             v3
 *                
 *                           |      |
 *                            ------
 *                              x2
 *                                  |      |
 *                                   ------
 *                                     x3
 *                
 *                
 *                   o vertices
 *                   x centroid
 *                
 * */
void Simulator3D::setDtLagrangianExtream()
{
 int idMinVolume = *min_element(m->getIdMinVolume().begin(),m->getIdMinVolume().end());
 double minEdge = m->getMinEdge();

 //cout << idMinVolume << " " << minEdge << endl;

 int v1 = IEN->Get(idMinVolume,0);
 int v2 = IEN->Get(idMinVolume,1);
 int v3 = IEN->Get(idMinVolume,2);
 int v4 = IEN->Get(idMinVolume,3);

 /* X-component */
 double p1x = X->Get(v1);
 double p2x = X->Get(v2);
 double p3x = X->Get(v3);
 double p4x = X->Get(v4);

 double xCentroid = ( p1x+p2x+p3x+p4x )*0.25;

 double d1x = dist(p1x,xCentroid);
 double d2x = dist(p2x,xCentroid);
 double d3x = dist(p3x,xCentroid);
 double d4x = dist(p4x,xCentroid);
 
 double minXdist = min(d1x,d2x);
 minXdist = min(minXdist,d3x);
 minXdist = min(minXdist,d4x);

 double xVelMax = max( 1.0,uALE.Abs().Max() );
 double minDtx = minXdist/xVelMax;

 /* Y-component */
 double p1y = Y->Get(v1);
 double p2y = Y->Get(v2);
 double p3y = Y->Get(v3);
 double p4y = Y->Get(v4);

 double yCentroid = ( p1y+p2y+p3y+p4y )*0.25;

 double d1y = dist(p1y,yCentroid);
 double d2y = dist(p2y,yCentroid);
 double d3y = dist(p3y,yCentroid);
 double d4y = dist(p4y,yCentroid);

 double minYdist = min(d1y,d2y);
 minYdist = min(minYdist,d3y);
 minYdist = min(minYdist,d4y);

 double yVelMax = max( 1.0,vALE.Abs().Max() );
 double minDty = minYdist/yVelMax;

 /* Z-component */
 double p1z = Z->Get(v1);
 double p2z = Z->Get(v2);
 double p3z = Z->Get(v3);
 double p4z = Z->Get(v4);

 double zCentroid = ( p1z+p2z+p3z+p4z )*0.25;

 double d1z = dist(p1z,zCentroid);
 double d2z = dist(p2z,zCentroid);
 double d3z = dist(p3z,zCentroid);
 double d4z = dist(p4z,zCentroid);

 double minZdist = min(d1z,d2z);
 minZdist = min(minZdist,d3z);
 minZdist = min(minZdist,d4z);

 double zVelMax = max( 1.0,wALE.Abs().Max() );
 double minDtz = minZdist/zVelMax;

 double minDt1 = min(minDtx,minDty);
 minDt1 = min(minDt1,minDtz);

 double velMax = max(xVelMax,yVelMax);
 velMax = max(velMax,zVelMax);
 double minDt2 = 0.5*minEdge/velMax;

 dtLagrangian = min(minDt1,minDt2);

 //cout << "minDt: " <<  minDt1 << " " << minDt2 << endl;
 //cout << minXdist << " " << minYdist << " " << minZdist << endl;
 //cout << xVelMax << " " << yVelMax << " " << zVelMax << endl;
 //cout << minDtx << " " << minDty << " " << minDtz << endl;
}

/*     
 *                l*3^(1/3)
 *   height = h = --------- = l*0.86602
 *                    2
 * */
void Simulator3D::setDtLagrangian()
{
 double minEdge = m->getMinEdge();

 double length = minEdge*0.86602;

 double velMax = max( 1.0,uALE.Abs().Max() );
 velMax = max( velMax,vALE.Abs().Max() );
 velMax = max( velMax,wALE.Abs().Max() );

//--------------------------------------------------
//  cout << "length: " << length << endl; 
//  cout << "uMax: " << uALE.Abs().Max() << endl; 
//  cout << "vMax: " << vALE.Abs().Max() << endl; 
//  cout << "wMax: " << wALE.Abs().Max() << endl; 
//  cout << "velMax: " << velMax << endl;
//-------------------------------------------------- 
 
 dtLagrangian = 0.5*length/velMax;
}

/*
 * Set lagrangian dt.
 * This method checks all the tetrahedron's edge length and respective
 * velocities at both extremities, then the gradient module |u| is used to
 * compute dt. 
 *
 * 1D example:
 *
 *
 *                  u1=+1                u2=-1               
 * case1:           x--->            <---x       |u|=|u1-u2|=2
 * (contract)       o--------------------o       dt=l/|u| = 0.5 (ok!) 
 *                 v1                    v2      
 *
 *              
 *                  u1=-1                u2=+1
 * case2:       <---x                    x--->   |u1-u2|=2
 * (stretch)        o--------------------o       dt=l/|u| = 0.5 (not ok!)
 *                 v1                    v2      (should be) dt = inf
 *              
 *
 *                  u1=+1                u2=+1
 * case3:           x--->                x--->   |u1-u2|=0
 * (displace)       o--------------------o       dt=l/|u| = inf (ok!)
 *                 v1                    v2
 *              
 *                  |                    |
 *                  |<------ l=1 ------->|
 *
 * OBS.: for case2, the dt could be anything, since the vertices are
 * stretching and they would never touch themselves. But is difficult to
 * predict such a case and for this reason we treat such a case as we
 * treat case1.
 *
 * */
void Simulator3D::setDtLagrangianNorberto()
{
 clMatrix* mapEdge = m->getMapEdge();

 dtLagrangian = 0.1;
 for( int edge=0;edge<mapEdge->DimI();edge++ )
 {
  // v1
  int v1 = mapEdge->Get(edge,4);
  double p1x=X->Get(v1);
  double p1y=Y->Get(v1);
  double p1z=Z->Get(v1);

  // v2
  int v2 = mapEdge->Get(edge,5);
  double p2x=X->Get(v2);
  double p2y=Y->Get(v2);
  double p2z=Z->Get(v2);

  double length = distance(p1x,p1y,p1z,p2x,p2y,p2z);

  // bubble.py - 146 iterations
  double xVel = fabs(uALE.Get(v1)) - fabs(uALE.Get(v2));
  double yVel = fabs(vALE.Get(v1)) - fabs(vALE.Get(v2));
  double zVel = fabs(wALE.Get(v1)) - fabs(wALE.Get(v2));

  double vel = vectorLength(xVel,yVel,zVel);
  //double vel = distance(fabs(xVel),fabs(yVel),fabs(zVel),
  //                    fabs(x),fabs(y),fabs(z));

  double a = 0.2; // security parameter

  double minDt = a*length/vel;

  if( minDt < dtLagrangian && minDt > 0 )
   dtLagrangian = minDt;
 }
}

/* 
 * Method to compute the Semi-Lagrangian dt using the biggest velocity
 * and the smallest tetrahedron edge length as:
 *
 * dt = l_min/v_max
 *
 * OBS.: This will limit the dt too much! Should be the same as the
 * setDtLagrangianNorberto!
 *
 * */
void Simulator3D::setDtSemiLagrangian()
{
 clMatrix* mapEdge = m->getMapEdge();
 clVector velU = uSolOld-uALE;
 clVector velV = vSolOld-vALE;
 clVector velW = wSolOld-wALE;

 dtSemiLagrangian = 1.0;
 for( int edge=0;edge<mapEdge->DimI();edge++ )
 {
  // v1
  int v1 = mapEdge->Get(edge,4);
  double p1x=X->Get(v1);
  double p1y=Y->Get(v1);
  double p1z=Z->Get(v1);

  // v2
  int v2 = mapEdge->Get(edge,5);
  double p2x=X->Get(v2);
  double p2y=Y->Get(v2);
  double p2z=Z->Get(v2);

  double length = distance(p1x,p1y,p1z,p2x,p2y,p2z);

  // bubble.py - 146 iterations
  double xVel = fabs(velU.Get(v1)) - fabs(velU.Get(v2));
  double yVel = fabs(velV.Get(v1)) - fabs(velV.Get(v2));
  double zVel = fabs(velW.Get(v1)) - fabs(velW.Get(v2));

  double vel = vectorLength(xVel,yVel,zVel);
  //double vel = distance(fabs(xVel),fabs(yVel),fabs(zVel),
  //                    fabs(x),fabs(y),fabs(z));

  double a = 0.2; // security parameter

  double minDt = a*length/vel;

  if( minDt < dtSemiLagrangian && minDt > 0 )
   dtSemiLagrangian = minDt;
 }
}

void Simulator3D::setDtGravity()
{
 double minEdge = m->getMinEdge();
 double velMax = max( 1.0,gravity.Abs().Max() );

 dtGravity = sqrt(velMax/minEdge);
}

/*
 * Set Dt of the current simulation.
 * Explicity terms:
 *  - Semi-Lagrangian;
 *  - gravity
 *  */
void Simulator3D::setDtEulerian()
{
 // setting required dt
 setDtSemiLagrangian();
 setDtGravity();
 dtSurfaceTension=0.0;

 double dtEulerian = min(1.0,getDtSemiLagrangian());
 dtEulerian = min(dtEulerian,getDtGravity());

 dt = cfl*dtEulerian;
}

/*
 * Set Dt of the current simulation.
 * Explicity terms:
 *  - Semi-Lagrangian;
 *  - Surface Tension;
 *  - Gravity.
 *  */
void Simulator3D::setDtALETwoPhase()
{
 // setting required dt
 setDtLagrangianNorberto();
 setDtGravity();
 dtSemiLagrangian = dtLagrangian;
 setDtSurfaceTension();

 double dtALETwoPhase = min(getDtLagrangian(),getDtSurfaceTension());
 dtALETwoPhase = min(dtALETwoPhase,getDtGravity());

 dt = cfl*dtALETwoPhase;
}

/*
 * Set Dt of the current simulation.
 * Explicity terms:
 *  - Semi-Lagrangian;
 *  - Surface Tension;
 *  - Gravity.
 *  */
void Simulator3D::setDtALESinglePhase()
{
 // setting required dt
 setDtLagrangianNorberto();
 setDtGravity();
 dtSemiLagrangian=dtLagrangian;
 dtSurfaceTension = 0.0;

 double dtALETwoPhase = min(getDtLagrangian(),getDtGravity());

 dt = cfl*dtALETwoPhase;
}

void Simulator3D::setDt()
{
 setDtALETwoPhase();
}

void Simulator3D::setCSol(clVector &_cSol)
{
 cSol = _cSol;
}

void Simulator3D::setNuZ(const char* _filename)
{
 // -- Leitura do perfil de nu variavel em Z para os nos da malha -- //

 double aux;
 double dist1,dist2;
 clMatrix muFile(1001,2); // number of points in the file

 ifstream file( _filename,ios::in );

 if( !file )
 {
  cerr << "Esta faltando o arquivo de perfis!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !file.eof() )
 {
  for( int i=0;i<muFile.DimI();i++ )
  {
   file >> aux;
   muFile.Set(i,0,aux);
   file >> aux;
   muFile.Set(i,1,aux);
  }
 }

 int j;
 for( int i=0;i<numVerts;i++ )
 {
  for( j=0;j<muFile.DimI()-1;j++ )
  {
   dist1 = fabs( Z->Get(i) - muFile(j,0) );
   dist2 = fabs( Z->Get(i) - muFile(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = muFile(j,1);
  mu.Set(i,aux);
 }
}

void Simulator3D::setNuC()
{
 for( int mele=0;mele<numElems;mele++ )
 {
  int v1 = (int) IEN->Get(mele,0);
  int v2 = (int) IEN->Get(mele,1);
  int v3 = (int) IEN->Get(mele,2);
  int v4 = (int) IEN->Get(mele,3);
  double c = ( cSol.Get(v1)+
	         cSol.Get(v2)+
	         cSol.Get(v3)+
	         cSol.Get(v4) )/4.0;

  double eme = 0.81315;
  double muC = exp(eme*c);

  // updating mu
  mu.Set(v1,muC);
  mu.Set(v2,muC);
  mu.Set(v3,muC);
  mu.Set(v4,muC);
 }
}

void Simulator3D::setSigma(double _sigma)
{
 sigma = _sigma;
 sigma_0 = sigma;
 sigmaAdimen = sigma/sigma_0;
}

void Simulator3D::setSolverVelocity(Solver *s){solverV = s;}
void Simulator3D::setSolverPressure(Solver *s){solverP = s;}
void Simulator3D::setSolverConcentration(Solver *s){solverC = s;}
void Simulator3D::setRe(double _Re){Re = _Re;}
double Simulator3D::getRe(){return Re;}
void Simulator3D::setSc(double _Sc){Sc = _Sc;}
double Simulator3D::getSc(){return Sc;}
void Simulator3D::setFr(double _Fr){Fr = _Fr;}
double Simulator3D::getFr(){return Fr;}
void Simulator3D::setWe(double _We){We = _We;}
double Simulator3D::getWe(){return We;}
void Simulator3D::setAlpha(double _alpha){alpha = _alpha;}
double Simulator3D::getAlpha(){return alpha;}
void Simulator3D::setBeta(double _beta){beta = _beta;}
double Simulator3D::getBeta(){return beta;}
double Simulator3D::getSigma(){return sigma;}
void Simulator3D::setDt(double _dt){dt = _dt;}
void Simulator3D::setTime(double _time){time = _time;}
double Simulator3D::getDt(){return dt;}
double Simulator3D::getDtLagrangian(){return dtLagrangian;}
double Simulator3D::getDtSemiLagrangian(){return dtSemiLagrangian;}
double Simulator3D::getDtGravity(){return dtGravity;}
double Simulator3D::getDtSurfaceTension(){return dtSurfaceTension;}
void Simulator3D::setIter(double _iter){iter = _iter;}
int  Simulator3D::getIter(){return iter;}
double Simulator3D::getCfl(){return cfl;}
void Simulator3D::setC1(double _c1){c1 = _c1;}
void Simulator3D::setC2(double _c2){c2 = _c2;}
void Simulator3D::setC3(double _c3){c3 = _c3;}
void Simulator3D::setD1(double _d1){d1 = _d1;}
void Simulator3D::setD2(double _d2){d2 = _d2;}
double Simulator3D::getC1(){return c1;}
double Simulator3D::getC2(){return c2;}
double Simulator3D::getC3(){return c3;}
double Simulator3D::getD1(){return d1;}
double Simulator3D::getD2(){return d2;}

// reference frame velocity
void Simulator3D::setURef(double _uRef)
{
 uRef = _uRef;
 xRef = uRef*time;
}
void Simulator3D::setVRef(double _vRef)
{
 vRef = _vRef;
 yRef = vRef*time;
}
void Simulator3D::setWRef(double _wRef)
{
 wRef = _wRef;
 zRef = wRef*time;
}
double Simulator3D::getURef(){return uRef;}
double Simulator3D::getVRef(){return vRef;}
double Simulator3D::getWRef(){return wRef;}
void Simulator3D::setXRef(double _xRef){xRef = _xRef;}
double Simulator3D::getXRef(){return xRef;}
void Simulator3D::setYRef(double _yRef){yRef = _yRef;}
double Simulator3D::getYRef(){return yRef;}
void Simulator3D::setZRef(double _zRef){zRef = _zRef;}
double Simulator3D::getZRef(){return zRef;}

void Simulator3D::setMu(double _mu_in)
{ 
 mu_in = _mu_in;
 mu_0 = _mu_in;
 mu_inAdimen = mu_in/mu_0;

 mu.SetAll(mu_inAdimen); 
}

void Simulator3D::setRho(double _rho_in)
{ 
 rho_in = _rho_in;
 rho_0 = rho_in; 
 rho_inAdimen = rho_in/rho_0; 

 rho.SetAll(rho_inAdimen); 
}

void Simulator3D::setCp(double _cp_in)
{ 
 cp_in = _cp_in;
 cp_0 = cp_in; 
 cp_inAdimen = cp_in/cp_0; 

 cp.SetAll(cp_inAdimen); 
}

void Simulator3D::setKt(double _kt_in)
{ 
 kt_in = _kt_in;
 kt_0 = kt_in; 
 kt_inAdimen = kt_in/kt_0; 

 kt.SetAll(kt_inAdimen); 
}

/* Since we are working here with non-dimensional equations and also
 * because the referential Reynolds used in two-phase flow is the one
 * from the liquid phase, we set the liquid viscosity (mu_fluid) as
 * equal to 1.0.
 *
 *          rho_0 * v_0 * D     
 * Re_ref = ---------------     
 *               mu_0            
 *
 *               rho_in                   rho_out
 * rho_inAdimen = -----     rho_outAdimen = -----
 *               rho_0                   rho_0
 *
 *               mu_in                   mu_out
 * muAdimen_l =  ----      muAdimen_g = ----
 *               mu_0                   mu_0
 *
 *
 *        rho_in * v_l * D          rho_out * v_g * D 
 * Re_l = ---------------   Re_g = --------------- 
 *             mu_in                     mu_out       
 *
 * Navier-Stokes Viscous Term (liquid):
 *
 *    mu_0                 [            (                       ) ]
 * --------- * \nabla \dot [ muAdimen * ( \nabla u + \nabla u^T ) ]
 * rho_0*v*D               [            (                       ) ]
 *
 *
 * Navier-Stokes Viscous Term (gas phase): [INSIDE BUBBLE]
 *
 *    mu_0                 [ mu_out   (                       ) ]
 * --------- * \nabla \dot [ ---- * ( \nabla u + \nabla u^T ) ]
 * rho_0*v*D               [ mu_0   (                       ) ]
 *
 *
 * Navier-Stokes Viscous Term (liquid phase): [OUTSIDE BUBBLE]
 *
 *    mu_0                 [ mu_in   (                       ) ]
 * --------- * \nabla \dot [ ---- * ( \nabla u + \nabla u^T ) ]
 * rho_0*v*D               [ mu_0   (                       ) ]
 *
 * */
void Simulator3D::setMu(double _mu_in,double _mu_out)
{ 
 mu_in = _mu_in;
 mu_out = _mu_out;

 if( mu_in >= mu_out )
  mu_0 = mu_in; 
 else
  mu_0 = mu_out;

 mu_inAdimen = mu_in/mu_0;
 mu_outAdimen = mu_out/mu_0;

 clVector one(numVerts);one.SetAll(1.0);
 mu = mu_inAdimen*(*heaviside) + mu_outAdimen*(one-(*heaviside));

 //double rMax = 1.0;
//--------------------------------------------------
//  for (list<int>::iterator it=boundaryVert->begin(); 
//                           it!=boundaryVert->end(); 
// 						  ++it)
//   //if (X->Get(*it)*X->Get(*it)+Y->Get(*it)*Y->Get(*it) > rMax*rMax - 0.001) 
//    mu.Set(*it,0.0); // zero viscosity at wall
//-------------------------------------------------- 


//--------------------------------------------------
//  mu.Dim(numNodes);
//  mu.SetAll(mu_inAdimen);
//  list<int> *inElem;
//  inElem = m->getInElem();
//  for (list<int>::iterator it=inElem->begin(); it!=inElem->end(); ++it)
//   mu.Set(*it,mu_outAdimen);
//-------------------------------------------------- 
}

void Simulator3D::setRho(double _rho_in,double _rho_out)
{ 
 rho_in = _rho_in;
 rho_out = _rho_out;

 if( rho_in >= rho_out )
  rho_0 = rho_in; 
 else
  rho_0 = rho_out; 

 rho_inAdimen = rho_in/rho_0; 
 rho_outAdimen = rho_out/rho_0;

 clVector one(numVerts);one.SetAll(1.0);
 rho = rho_inAdimen*(*heaviside) + rho_outAdimen*(one-(*heaviside));

 //double rMax = 1.0;
//--------------------------------------------------
//  for (list<int>::iterator it=boundaryVert->begin(); 
//                           it!=boundaryVert->end(); 
// 						  ++it)
//   //if (X->Get(*it)*X->Get(*it)+Y->Get(*it)*Y->Get(*it) > rMax*rMax - 0.001) 
//    rho.Set(*it,0.0); // zero density at wall
//-------------------------------------------------- 

//--------------------------------------------------
//  rho.Dim(numNodes);
//  rho.SetAll(rho_inAdimen);
//  list<int> *inElem;
//  inElem = m->getInElem();
//  for (list<int>::iterator it=inElem->begin(); it!=inElem->end(); ++it)
//   rho.Set(*it,rho_outAdimen);
//-------------------------------------------------- 
}

void Simulator3D::setCp(double _cp_in,double _cp_out)
{ 
 cp_in = _cp_in;
 cp_out = _cp_out;

 if( cp_in >= cp_out )
  cp_0 = cp_in; 
 else
  cp_0 = cp_out;

 cp_inAdimen = cp_in/cp_0; 
 cp_outAdimen = cp_out/cp_0;

 clVector one(numVerts);one.SetAll(1.0);
 cp = cp_inAdimen*(*heaviside) + cp_outAdimen*(one-(*heaviside));

//--------------------------------------------------
//  for (list<int>::iterator it=boundaryVert->begin(); 
//                           it!=boundaryVert->end(); 
// 						  ++it)
//   cp.Set(*it,0.0); // set cp = 0 in the Wall
//-------------------------------------------------- 
}

void Simulator3D::setKt(double _kt_in,double _kt_out)
{ 
 kt_in = _kt_in;
 kt_out = _kt_out;

 if( kt_in >= kt_out )
  kt_0 = kt_in; 
 else
  kt_0 = kt_out;

 kt_inAdimen = kt_in/kt_0; 
 kt_outAdimen = kt_out/kt_0;

 clVector one(numVerts);one.SetAll(1.0);
 kt = kt_inAdimen*(*heaviside) + kt_outAdimen*(one-(*heaviside));

//--------------------------------------------------
//  for (list<int>::iterator it=boundaryVert->begin(); 
//                           it!=boundaryVert->end(); 
// 						  ++it)
//   kt.Set(*it,0.0); // set kt = 0 in the Wall
//-------------------------------------------------- 
}

void Simulator3D::setMuSmooth(double _mu_in,double _mu_out)
{ 
 mu_in = _mu_in;
 mu_out = _mu_out;

 if( mu_in >= mu_out )
  mu_0 = mu_in; 
 else
  mu_0 = mu_out;

 mu_inAdimen = mu_in/mu_0; 
 mu_outAdimen = mu_out/mu_0;

 clVector one(numVerts);one.SetAll(1.0);
 mu = mu_inAdimen*hSmooth + mu_outAdimen*(one-hSmooth);
}

void Simulator3D::setRhoSmooth(double _rho_in,double _rho_out)
{ 
 rho_in = _rho_in;
 rho_out = _rho_out;

 if( rho_in >= rho_out )
  rho_0 = rho_in; 
 else
  rho_0 = rho_out; 

 rho_inAdimen = rho_in/rho_0; 
 rho_outAdimen = rho_out/rho_0;

 clVector one(numVerts);one.SetAll(1.0);
 rho = rho_inAdimen*hSmooth + rho_outAdimen*(one-hSmooth);
}

void Simulator3D::setHSmooth()
{
 hSmooth.Dim(numVerts);
 clVector half(numVerts);half.SetAll(0.5);
 clVector zeroLevel = ((*heaviside)-half)*2;
 triEdge = m->getTriEdge();
 double triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));
 for( int i=0;i<numVerts;i++ )
 {
  double len = 1.3*triEdgeMin;
  double d = interfaceDistance->Get(i);
  double aux = zeroLevel.Get(i)*d;

  if( aux < -len )
   hSmooth.Set( i,0.0 );
  else if( aux > len )
   hSmooth.Set( i,1.0 );
  else
  {
   //double func = 0.5;
   double func = 0.5*( 1.0+(aux/len)+(1.0/3.1415)*sin(3.1415*(aux/len)) );
   hSmooth.Set( i,func );
  }
 }
}

//--------------------------------------------------
// double Simulator3D::computeReynolds()
// {
//  double D = 1.0;
//  double L = D;
//  double U = sqrt(g*D);
//  Re = rho_in*L*U/mu_in;
// }
// 
// double Simulator3D::computeFroud()
// {
//  double D = 1.0;
//  double L = D;
//  double U = sqrt(g*D);
//  Fr = U/sqrt(g*L);
// }
// 
// double Simulator3D::computeWebber()
// {
//  double D = 1.0;
//  double L = D;
//  double U = sqrt(g*D);
//  We = rho_in*L*U*U/sigma;
// }
// 
// double Simulator3D::computeEotvos()
// {
//  double D = 1.0;
//  Eo = rho_in*g*D*D/sigma;
//  We = Eo;
// }
// 
// double Simulator3D::computeGalileo()
// {
//  double D = 1.0;
//  double L = D;
//  double U = sqrt(g*D);
//  N = rho_in*U*L/mu_in;
//  Re = N;
// }
// 
// double Simulator3D::computeMorton()
// {
//  double D = 1.0;
//  double L = D;
//  double U = sqrt(g*D);
//  Mo = (rho_in-rho_out)*mu_in*mu_in*mu_in*mu_in*g/(rho_in*rho_in*sigma*sigma*sigma); 
// }
//-------------------------------------------------- 

double Simulator3D::getMu_in(){return mu_in;}
double Simulator3D::getMu_out(){return mu_out;}
double Simulator3D::getRho_in(){return rho_in;}
double Simulator3D::getRho_out(){return rho_out;}
double Simulator3D::getCp_in(){return cp_in;}
double Simulator3D::getCp_out(){return cp_out;}
double Simulator3D::getKt_in(){return kt_in;}
double Simulator3D::getKt_out(){return kt_out;}
double Simulator3D::getMu_inAdimen(){return mu_inAdimen;}
double Simulator3D::getMu_outAdimen(){return mu_outAdimen;}
double Simulator3D::getRho_inAdimen(){return rho_inAdimen;}
double Simulator3D::getRho_outAdimen(){return rho_outAdimen;}
double Simulator3D::getCp_inAdimen(){return cp_inAdimen;}
double Simulator3D::getCp_outAdimen(){return cp_outAdimen;}
double Simulator3D::getKt_inAdimen(){return kt_inAdimen;}
double Simulator3D::getKt_outAdimen(){return kt_outAdimen;}
double Simulator3D::getTime(){return time;}
clVector* Simulator3D::getUSol(){return &uSol;} 
clVector* Simulator3D::getUSolOld(){return &uSolOld;} 
void Simulator3D::setUSol(clVector *_uSol){uSol = *_uSol;}
void Simulator3D::setUALE(clVector *_uALE){uALE = *_uALE;}
clVector* Simulator3D::getVSol(){return &vSol;} 
clVector* Simulator3D::getVSolOld(){return &vSolOld;} 
void Simulator3D::setVSol(clVector *_vSol){vSol = *_vSol;}
void Simulator3D::setVALE(clVector *_vALE){vALE = *_vALE;}
clVector* Simulator3D::getWSol(){return &wSol;}
clVector* Simulator3D::getWSolOld(){return &wSolOld;}
void Simulator3D::setWSol(clVector *_wSol){wSol = *_wSol;}
void Simulator3D::setWALE(clVector *_wALE){wALE = *_wALE;}
clVector* Simulator3D::getPSol(){return &pSol;}
clVector* Simulator3D::getPSolOld(){return &pSolOld;}
clVector* Simulator3D::getCSol(){return &cSol;}
clVector* Simulator3D::getCSolOld(){return &cSolOld;}
clVector* Simulator3D::getUALE(){return &uALE;} 
clVector* Simulator3D::getUALEOld(){return &uALEOld;} 
clVector* Simulator3D::getVALE(){return &vALE;} 
clVector* Simulator3D::getVALEOld(){return &vALEOld;} 
clVector* Simulator3D::getWALE(){return &wALE;} 
clVector* Simulator3D::getWALEOld(){return &wALEOld;} 
clVector* Simulator3D::getFint(){return &fint;}
clVector* Simulator3D::getGravity(){return &gravity;}
clVector* Simulator3D::getBetaFlowLiq(){return &betaFlowLiq;}
double Simulator3D::getGrav(){return g;}
clDMatrix* Simulator3D::getKappa(){return &kappa;}
clMatrix* Simulator3D::getK(){return &K;}
clMatrix* Simulator3D::getM(){return &Mrho;}
clMatrix* Simulator3D::getGx(){return &gx;}
clMatrix* Simulator3D::getGy(){return &gy;}
clMatrix* Simulator3D::getGz(){return &gz;}
clMatrix* Simulator3D::getG(){return &G;}
clMatrix* Simulator3D::getD(){return &D;}
clVector* Simulator3D::getMu(){return &mu;}
clVector* Simulator3D::getRho(){return &rho;}
clVector* Simulator3D::getCp(){return &cp;}
clVector* Simulator3D::getKt(){return &kt;}
clVector* Simulator3D::getHSmooth(){return &hSmooth;}
clVector* Simulator3D::getHeatFlux(){return &heatFlux;}
void Simulator3D::updateIEN(){IEN = m->getIEN();}
void Simulator3D::setCfl(double _cfl){cfl = _cfl;}
vector<double> Simulator3D::getCentroidVelX(){return centroidVelX;}
vector<double> Simulator3D::getCentroidVelY(){return centroidVelY;}
vector<double> Simulator3D::getCentroidVelZ(){return centroidVelZ;}
void Simulator3D::setCentroidVelX(vector<double> _centroidVelX)
{centroidVelX = _centroidVelX;}
void Simulator3D::setCentroidVelY(vector<double> _centroidVelY)
{centroidVelY = _centroidVelY;}
void Simulator3D::setCentroidVelZ(vector<double> _centroidVelZ)
{centroidVelZ = _centroidVelZ;}

vector<double> Simulator3D::getCentroidPosX(){return centroidPosX;}
vector<double> Simulator3D::getCentroidPosY(){return centroidPosY;}
vector<double> Simulator3D::getCentroidPosZ(){return centroidPosZ;}
void Simulator3D::setCentroidPosX(vector<double> _centroidPosX)
{centroidPosX = _centroidPosX;}
void Simulator3D::setCentroidPosY(vector<double> _centroidPosY)
{centroidPosY = _centroidPosY;}
void Simulator3D::setCentroidPosZ(vector<double> _centroidPosZ)
{centroidPosZ = _centroidPosZ;}


// set do centroide para o elemento mini apos a interpolacao linear
// (applyLinearInterpolation)
clVector Simulator3D::setCentroid(clVector &_vector)
{
 double aux;
 clVector zerosConv(numNodes-numVerts);
 _vector.Append(zerosConv);

 for( int mele=0;mele<numElems;mele++ )
 {
  int v1 = (int) IEN->Get(mele,0);
  int v2 = (int) IEN->Get(mele,1);
  int v3 = (int) IEN->Get(mele,2);
  int v4 = (int) IEN->Get(mele,3);
  int v5 = (int) IEN->Get(mele,4);
  aux = ( _vector.Get(v1)+_vector.Get(v2)+_vector.Get(v3)+_vector.Get(v4) );
  aux = aux/4.0;
  _vector.Set(v5,aux);
 }
  return _vector;
}// fim do metodo compute -> setCentroid

// Atribui o Simulator3D do argumento no corrente
void Simulator3D::operator=(Simulator3D &_sRight) 
{
 getModel3DAttrib(*_sRight.m);

 Re = _sRight.Re;
 Sc = _sRight.Sc;
 Fr = _sRight.Fr;
 We = _sRight.We;
 alpha = _sRight.alpha;
 beta = _sRight.beta;
 dt = _sRight.dt;
 dtSemiLagrangian = _sRight.dtSemiLagrangian;
 dtLagrangian = _sRight.dtLagrangian;
 dtSurfaceTension = _sRight.dtSurfaceTension;
 dtGravity = _sRight.dtGravity;
 cfl = _sRight.cfl;
 time = _sRight.time;
 c1 = _sRight.c1;
 c2 = _sRight.c2;
 c3 = _sRight.c3;
 d1 = _sRight.d1;
 d2 = _sRight.d2;
 iter = _sRight.iter;
 uRef = _sRight.uRef;
 vRef = _sRight.vRef;
 wRef = _sRight.wRef;
 xRef = _sRight.xRef;
 yRef = _sRight.yRef;
 zRef = _sRight.zRef;

 g = _sRight.g;
 sigma = _sRight.sigma;
 rho_in = _sRight.rho_in;
 rho_out = _sRight.rho_out;
 mu_in = _sRight.mu_in;
 mu_out = _sRight.mu_out;
 cp_in = _sRight.cp_in;
 cp_out = _sRight.cp_out;
 kt_in = _sRight.kt_in;
 kt_out = _sRight.kt_out;

 g_0 = _sRight.g_0;
 sigma_0 = _sRight.sigma_0;
 rho_0 = _sRight.rho_0;
 mu_0 = _sRight.mu_0;
 cp_0 = _sRight.cp_0;
 kt_0 = _sRight.kt_0;

 gAdimen = _sRight.gAdimen;
 sigmaAdimen = _sRight.sigmaAdimen;
 rho_inAdimen = _sRight.rho_inAdimen;
 rho_outAdimen = _sRight.rho_outAdimen;
 mu_inAdimen = _sRight.mu_inAdimen;
 mu_outAdimen = _sRight.mu_outAdimen;
 cp_inAdimen = _sRight.cp_inAdimen;
 cp_outAdimen = _sRight.cp_outAdimen;
 kt_inAdimen = _sRight.kt_inAdimen;
 kt_outAdimen = _sRight.kt_outAdimen;

 K = _sRight.K;
 Kc = _sRight.Kc;
 Mrho = _sRight.Mrho;
 M = _sRight.M;
 Mc = _sRight.Mc;
 G = _sRight.G;
 Gc = _sRight.Gc;
 D = _sRight.D;
 mat = _sRight.mat;
 matc = _sRight.matc;
 MrhoLumped = _sRight.MrhoLumped;
 MLumped = _sRight.MLumped;
 McLumped = _sRight.McLumped;
 gx = _sRight.gx;
 gy = _sRight.gy;
 gz = _sRight.gz;
 A = _sRight.A;
 b = _sRight.b;
 ATilde = _sRight.ATilde;
 AcTilde = _sRight.AcTilde;
 GTilde = _sRight.GTilde;
 DTilde = _sRight.DTilde;
 ETilde = _sRight.ETilde;
 E = _sRight.E;
 invA = _sRight.invA;
 invC = _sRight.invC;
 invMrhoLumped = _sRight.invMrhoLumped;
 invMLumped = _sRight.invMLumped;
 invMcLumped = _sRight.invMcLumped;

 uSol = _sRight.uSol;
 vSol = _sRight.vSol;
 wSol = _sRight.wSol;
 pSol = _sRight.pSol;
 cSol = _sRight.cSol;
 velU = _sRight.velU;
 velV = _sRight.velV;
 velW = _sRight.velW;
 uALE = _sRight.uALE;
 vALE = _sRight.vALE;
 wALE = _sRight.wALE;
 uSL = _sRight.uSL;
 vSL = _sRight.vSL;
 wSL = _sRight.wSL;
 cSL = _sRight.cSL;
 uSmooth = _sRight.uSmooth;
 vSmooth = _sRight.vSmooth;
 wSmooth = _sRight.wSmooth;
 uSmoothCoord = _sRight.uSmoothCoord;
 vSmoothCoord = _sRight.vSmoothCoord;
 wSmoothCoord = _sRight.wSmoothCoord;
 uSmoothSurface = _sRight.uSmoothSurface;
 vSmoothSurface = _sRight.vSmoothSurface;
 wSmoothSurface = _sRight.wSmoothSurface;
 centroidVelX = _sRight.centroidVelX;
 centroidVelY = _sRight.centroidVelY;
 centroidVelZ = _sRight.centroidVelZ;
 centroidPosX = _sRight.centroidPosX;
 centroidPosY = _sRight.centroidPosY;
 centroidPosZ = _sRight.centroidPosZ;

 va = _sRight.va;
 vcc = _sRight.vcc;
 convUVW = _sRight.convUVW;
 convC = _sRight.convC;
 uTilde = _sRight.uTilde;
 cTilde = _sRight.cTilde;
 pTilde = _sRight.pTilde;
 b1 = _sRight.b1;
 b1c = _sRight.b1c;
 b2 = _sRight.b2;
 ip = _sRight.ip;
 ipc = _sRight.ipc;

 // two-phase vectors
 kappa   = _sRight.kappa;
 fint    = _sRight.fint;
 gravity = _sRight.gravity;
 betaFlowLiq = _sRight.betaFlowLiq;
 Fold    = _sRight.Fold;
 mu      = _sRight.mu;
 rho     = _sRight.rho;
 cp      = _sRight.cp;
 kt      = _sRight.kt;
 hSmooth = _sRight.hSmooth;
 heatFlux= _sRight.heatFlux;
 
 // old ints
 numVertsOld = _sRight.numVerts;
 numNodesOld = _sRight.numNodes;
 numElemsOld = _sRight.numElems;
 
 // oldSol vectors
 uSolOld    = _sRight.uSolOld;
 vSolOld    = _sRight.vSolOld;
 wSolOld    = _sRight.wSolOld;
 pSolOld    = _sRight.pSolOld;
 cSolOld    = _sRight.cSolOld;
 uALEOld    = _sRight.uALEOld;
 vALEOld    = _sRight.vALEOld;
 wALEOld    = _sRight.wALEOld;
 kappaOld   = _sRight.kappaOld;
 fintOld    = _sRight.fintOld;
 gravityOld = _sRight.gravityOld;
 betaFlowLiqOld = _sRight.betaFlowLiqOld;
 muOld      = _sRight.muOld;
 rhoOld     = _sRight.rhoOld;
 cpOld      = _sRight.cpOld;
 ktOld      = _sRight.ktOld;
 hSmoothOld = _sRight.hSmoothOld;
 heatFluxOld= _sRight.heatFluxOld;
 centroidVelXOld = _sRight.centroidVelXOld;
 centroidVelYOld = _sRight.centroidVelYOld;
 centroidVelZOld = _sRight.centroidVelZOld;
 centroidPosXOld = _sRight.centroidPosXOld;
 centroidPosYOld = _sRight.centroidPosYOld;
 centroidPosZOld = _sRight.centroidPosZOld;

 solverV = _sRight.solverV;
 solverP = _sRight.solverP;
 solverC = _sRight.solverC;
}

void Simulator3D::operator()(Model3D &_m) 
{
 // mesh information vectors
 getModel3DAttrib(_m);

 Re    = 10;
 Sc    = 2;
 Fr    = 0.1;
 We    = 10;
 sigma = 1;
 alpha = 1;
 beta  = 0;
 dt    = 0.01;
 dtSemiLagrangian = 0.01;
 dtLagrangian = 0.01;
 dtSurfaceTension = 0.01;
 dtGravity = 0.01;
 time  = 0.0;
 cfl   = 0.5;
 iter  = 0;
 c1    = 1.0;
 c2    = 0.0;
 c3    = 0.0;
 d1    = 1.0;
 d2    = 0.1;
 uRef  = 0.0;
 vRef  = 0.0;
 wRef  = 0.0;
 xRef  = 0.0;
 yRef  = 0.0;
 zRef  = 0.0;
 g     = 9.81;
 mu_in  = 1.0;
 mu_out  = 1.0;
 rho_in = 1.0;
 rho_out = 1.0;
 cp_in = 1.0;
 cp_out = 1.0;
 kt_in = 1.0;
 kt_out = 1.0;

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 allocateMemoryToAttrib();
}

int Simulator3D::loadSolution( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // read numVerts and numNodes and all the simulation parameters
 // iteration
 string file = (string) _dir + "vtk/" + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ifstream fileP( filename,ios::in ); 
 if( !fileP )
 {
  cerr << "VTK file is missing for reading!" << endl;
  exit(1);
 }

 char auxstr[255];

 while( ( !fileP.eof())&&(strcmp(auxstr,"TIME") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> dt;
 fileP >> cfl;
 fileP >> time;

 while( ( !fileP.eof())&&(strcmp(auxstr,"ITERATION") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> iter;

 while( ( !fileP.eof())&&(strcmp(auxstr,"NODES") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> numVertsOld;
 fileP >> numNodesOld;
 fileP >> numElemsOld;

 while( ( !fileP.eof())&&(strcmp(auxstr,"PARAMETERS") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> Re;
 fileP >> Sc;
 fileP >> Fr;
 fileP >> We;

 while( ( !fileP.eof())&&(strcmp(auxstr,"PROPERTIES") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> mu_in;
 fileP >> mu_out;
 fileP >> rho_in;
 fileP >> rho_out;
 fileP >> sigma;

 while( ( !fileP.eof())&&(strcmp(auxstr,"REFERENCEVELOCITY") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> uRef;
 fileP >> vRef;
 fileP >> wRef;
 fileP >> xRef;
 fileP >> yRef;
 fileP >> zRef;

 while( ( !fileP.eof())&&(strcmp(auxstr,"COEFFICIENTS") != 0) )
  fileP >> auxstr;

 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> auxstr;
 fileP >> c1;
 fileP >> c2;
 fileP >> c3;
 fileP >> d1;
 fileP >> d2;
 fileP >> alpha;
 fileP >> beta;

 while( ( !fileP.eof())&&(strcmp(auxstr,"CHARACTERISTICLENGTH") != 0) )
  fileP >> auxstr;

 int nRegions = 0;
 fileP >> auxstr; // 1
 fileP >> nRegions;
 fileP >> auxstr; // float
 for( int nb=0;nb<nRegions;nb++ )
  fileP >> triEdge[nb];

//--------------------------------------------------
//  cout << dt << " " << cfl << " " << time << endl;
//  cout << iter << endl;
//  cout << numVertsOld << " " << numNodesOld << " " << numElemsOld << endl;
//  cout << Re << " " << Sc << " " << Fr << " " << We << " " << endl;
//  cout << mu_in << " " << mu_out << " " << rho_in << " " << rho_out << " " << endl;
//  cout << c1 << " " << c2 << " " << c3 << " " << d1 << " " 
//       << d2 << " "<< alpha << " " << beta << endl;
//  cout << lineEdge << endl;
//-------------------------------------------------- 

 // read solution of U, V, W and P
 string fileUVWPC = (string) _dir + "bin/" + (string) _filename + "-" + str + ".bin";
 const char* filenameUVPC = fileUVWPC.c_str();

 clVector aux2(6*numNodesOld+2*numVertsOld); // vetor tambem carrega a concentracao

 ifstream UVPC_file( filenameUVPC,ios::in | ios::binary ); 

 if( !UVPC_file)
 {
  cerr << "Solution file is missing for reading!" << endl;
  exit(1);
 }

 UVPC_file.read( (char*) aux2.GetVec(),aux2.Dim()*sizeof(double) );

 UVPC_file.close();

 uSol.Dim(numNodesOld);
 vSol.Dim(numNodesOld);
 wSol.Dim(numNodesOld);
 pSol.Dim(numVertsOld);
 cSol.Dim(numVertsOld);
 uALE.Dim(numNodesOld);
 vALE.Dim(numNodesOld);
 wALE.Dim(numNodesOld);

 aux2.CopyTo(0,uSol);
 aux2.CopyTo(numNodesOld,vSol);
 aux2.CopyTo(2*numNodesOld,wSol);
 aux2.CopyTo(3*numNodesOld,pSol);
 aux2.CopyTo(3*numNodesOld+numVertsOld,cSol);
 aux2.CopyTo(3*numNodesOld+2*numVertsOld,uALE);
 aux2.CopyTo(4*numNodesOld+2*numVertsOld,vALE);
 aux2.CopyTo(5*numNodesOld+2*numVertsOld,wALE);

 // copy to SolOld
 uSolOld = uSol;
 vSolOld = vSol;
 wSolOld = wSol;
 pSolOld = pSol;
 cSolOld = cSol;
 uALEOld = uALE;
 vALEOld = vALE;
 wALEOld = wALE;

 cout << "Solution No. " << iter << " read" << endl;

 return iter+1;
} // fecha metodo loadSol 

// interpolacao linear dos vetores velocidade e pressao calculados na
// malha antiga (Model3D mold) para malha nova (this). Note que a
// interpolacao eh feita somente nos vertices da malha deixando para o
// metodo setCentroid a configuracao dos valores de cada centroide
// pertencente aos elementos
void Simulator3D::applyLinearInterpolation(Model3D &_mOld)
{
 // xVert da malha nova
 clVector xVert(numVerts);
 clVector yVert(numVerts);
 clVector zVert(numVerts);
 X->CopyTo(0,xVert);
 Y->CopyTo(0,yVert);
 Z->CopyTo(0,zVert);

 /*
  * LINEAR interpolation of mesh - numVerts - (lib/interpolation.h)
  * interpLin is a interpolation matrix, considering the previous Model3D 
  * passed as argument and the new mesh coordinates (xVert,yVert and zVert).
  * Due to its linear kink, the iterpolation is only considering the
  * vertices and it is NOT including the others nodes (centroid, mid
  * edge point, etc.)
  * This matrix is ready to be applied to some vector.
  * */ 
 clMatrix interpLin = meshInterp(_mOld,xVert,yVert,zVert);

 uSol.Dim( numVerts );
 vSol.Dim( numVerts );
 wSol.Dim( numVerts );
 pSol.Dim( numVerts );
 cSol.Dim( numVerts );
 uALE.Dim( numVerts );
 vALE.Dim( numVerts );
 wALE.Dim( numVerts );
 mu.Dim( numVerts );
 rho.Dim( numVerts );
 cp.Dim( numVerts );
 kt.Dim( numVerts );
 hSmooth.Dim( numVerts );
 heatFlux.Dim( numVerts );

 // the edgeSize vector is not part of the Simulator3D, but because we
 // are not saving the matrix interpLin, we should interpolate edgeSize
 // and send it back to Model3D. 
 // It is better to find another solution!!!
 clVector edgeSize( numVerts );

 // only 1 component because the 2 others are exactly the same
 clVector xKappaOld(_mOld.getNumVerts());
 for( int i=0;i<_mOld.getNumVerts();i++ )
 {
  double aux = kappaOld.Get(i);
  xKappaOld.Set(i,aux);
 }
 clVector xFintOld(_mOld.getNumVerts());
 clVector yFintOld(_mOld.getNumVerts());
 clVector zFintOld(_mOld.getNumVerts());
 fintOld.CopyTo(_mOld.getNumNodes()*0,xFintOld);
 fintOld.CopyTo(_mOld.getNumNodes()*1,yFintOld);
 fintOld.CopyTo(_mOld.getNumNodes()*2,zFintOld);

 clVector xGravityOld(_mOld.getNumVerts());
 clVector yGravityOld(_mOld.getNumVerts());
 clVector zGravityOld(_mOld.getNumVerts());
 gravityOld.CopyTo(_mOld.getNumNodes()*0,xGravityOld);
 gravityOld.CopyTo(_mOld.getNumNodes()*1,yGravityOld);
 gravityOld.CopyTo(_mOld.getNumNodes()*2,zGravityOld);
 
 // setting dimension of numVertsOld for all vectors.
 int numVertsOld =  _mOld.getNumVerts();
 clVector uSolOldVert(numVertsOld);
 clVector vSolOldVert(numVertsOld);
 clVector wSolOldVert(numVertsOld);
 clVector pSolOldVert(numVertsOld);
 clVector cSolOldVert(numVertsOld);
 clVector uALEOldVert(numVertsOld);
 clVector vALEOldVert(numVertsOld);
 clVector wALEOldVert(numVertsOld);
 clDMatrix kappaOldVert(numVertsOld);
 clVector fintOldVert(numVertsOld);
 clVector gravityOldVert(numVertsOld);
 clVector muOldVert(numVertsOld);
 clVector rhoOldVert(numVertsOld);
 clVector cpOldVert(numVertsOld);
 clVector ktOldVert(numVertsOld);
 clVector hSmoothOldVert(numVertsOld);
 clVector heatFluxOldVert(numVertsOld);
 
 // shrinking the solution vectors with size numNodes to numVerts. 
 uSolOld.CopyTo(0,uSolOldVert);
 vSolOld.CopyTo(0,vSolOldVert);
 wSolOld.CopyTo(0,wSolOldVert);
 pSolOld.CopyTo(0,pSolOldVert);
 cSolOld.CopyTo(0,cSolOldVert);
 uALEOld.CopyTo(0,uALEOldVert);
 vALEOld.CopyTo(0,vALEOldVert);
 wALEOld.CopyTo(0,wALEOldVert);
 fintOld.CopyTo(0,fintOldVert);
 gravityOld.CopyTo(0,gravityOldVert);
 muOld.CopyTo(0,muOldVert);
 rhoOld.CopyTo(0,rhoOldVert);
 cpOld.CopyTo(0,cpOldVert);
 ktOld.CopyTo(0,ktOldVert);
 hSmoothOld.CopyTo(0,hSmoothOldVert);
 heatFluxOld.CopyTo(0,heatFluxOldVert);
 clVector edgeSizeOld = *_mOld.getEdgeSize();

 // interpolation process between the old mesh with numVertsOld to the new
 // numVert.
 pSol = interpLin*(pSolOld);
 cSol = interpLin*(cSolOld);
 mu = interpLin*(muOld);
 rho = interpLin*(rhoOld);
 cp  = interpLin*(cpOld);
 kt  = interpLin*(ktOld);
 hSmooth = interpLin*(hSmoothOld);
 heatFlux = interpLin*(heatFluxOld);
 edgeSize = interpLin*(edgeSizeOld);

 // For velocities, kappa, fint and gravity it is mandatory to
 // reallocate these vectors using numNodes. For this we should first
 // append zeros (numNodes-numVerts) and then setCentroi or setQuad.
 clVector zeros(numNodes-numVerts);
 uSol = interpLin*(uSolOld);uSol.Append(zeros);
 vSol = interpLin*(vSolOld);vSol.Append(zeros);
 wSol = interpLin*(wSolOld);wSol.Append(zeros);
 uALE = interpLin*(uALEOld);uALE.Append(zeros);
 vALE = interpLin*(vALEOld);vALE.Append(zeros);
 wALE = interpLin*(wALEOld);wALE.Append(zeros);
 clVector xKappa = interpLin*(xKappaOld);xKappa.Append(zeros);
 clVector xFint = interpLin*(xFintOld);xFint.Append(zeros);
 clVector yFint = interpLin*(yFintOld);yFint.Append(zeros);
 clVector zFint = interpLin*(zFintOld);zFint.Append(zeros);
 clVector xGravity = interpLin*(xGravityOld);xGravity.Append(zeros);
 clVector yGravity = interpLin*(yGravityOld);yGravity.Append(zeros);
 clVector zGravity = interpLin*(zGravityOld);zGravity.Append(zeros);

#if NUMGLEU == 5 // for the MINI element
  // set do centroid
  uSol = setTetCentroid(*IEN,uSol);
  vSol = setTetCentroid(*IEN,vSol);
  wSol = setTetCentroid(*IEN,wSol);
  uALE = setTetCentroid(*IEN,uALE);
  vALE = setTetCentroid(*IEN,vALE);
  wALE = setTetCentroid(*IEN,wALE);
  xKappa = setTetCentroid(*IEN,xKappa);
  xFint = setTetCentroid(*IEN,xFint);
  yFint = setTetCentroid(*IEN,yFint);
  zFint = setTetCentroid(*IEN,zFint);
  xGravity = setTetCentroid(*IEN,xGravity);
  yGravity = setTetCentroid(*IEN,yGravity);
  zGravity = setTetCentroid(*IEN,zGravity);
#else // for the QUAD element
  // set do quad
  uSol = setTetQuad(*IEN,uSol);
  vSol = setTetQuad(*IEN,vSol);
  wSol = setTetQuad(*IEN,wSol);
  uALE = setTetQuad(*IEN,uALE);
  vALE = setTetQuad(*IEN,vALE);
  wALE = setTetQuad(*IEN,wALE);
  xKappa = setTetQuad(*IEN,xKappa);
  xFint = setTetQuad(*IEN,xFint);
  yFint = setTetQuad(*IEN,yFint);
  zFint = setTetQuad(*IEN,zFint);
  xGravity = setTetQuad(*IEN,xGravity);
  yGravity = setTetQuad(*IEN,yGravity);
  zGravity = setTetQuad(*IEN,zGravity);
#endif

 /* fint and gravity (as well as kappa) are calculated every time step,
  * thus their interpolation to the new mesh is not mandatory, but since
  * we are saving up all the simulation after the re-meshing process, it
  * is necessary to interpolate then so they can be shown on the VTK
  * file. 
  * */
 // kappaX,kappaY and kappaZ have the same values, so we need to set up 
 // these values on the vector kappa(3*numNodes)
 kappa.Dim(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  double aux = xKappa.Get(i);
  kappa.Set(i,aux);
  kappa.Set(i+numNodes,aux);
  kappa.Set(i+numNodes*2,aux);
 }
 
 // setting up fint and gravity.
 fint.Dim(3*numNodes);
 fint.CopyFrom(numNodes*0,xFint);
 fint.CopyFrom(numNodes*1,yFint);
 fint.CopyFrom(numNodes*2,zFint);
 gravity.Dim(3*numNodes);
 gravity.CopyFrom(numNodes*0,xGravity);
 gravity.CopyFrom(numNodes*1,yGravity);
 gravity.CopyFrom(numNodes*2,zGravity);
 
 // updating setEdgeSize on Model3D.
 m->setEdgeSize(edgeSize);

 // updating centroidVelPos
 if( surfMesh->numInterfaces > 0 )
  setCentroidVelPos();
 //setCentroidVelPosInterface();

 // updating old data vectors with the new mesh.
 uSolOld     = uSol;
 vSolOld     = vSol;
 wSolOld     = wSol;
 pSolOld     = pSol;
 cSolOld     = cSol;
 uALEOld     = uALE;
 vALEOld     = vALE;
 wALEOld     = wALE;
 fintOld     = fint;
 gravityOld  = gravity;
 kappaOld    = kappa;
 muOld       = mu;
 rhoOld      = rho;
 cpOld       = cp;
 ktOld       = kt;
 hSmoothOld  = hSmooth;
 heatFluxOld = heatFluxOld;

 numVertsOld = numVerts;
 numNodesOld = numNodes;
 numElemsOld = numElems;
} // fecha metodo applyLinearInterpolation


// impoe velocidade Lagrangian = 0 no contorno
void Simulator3D::setLagrangianVelBC()
{
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
 {
  uSolOld.Set(*it,0.0);
  vSolOld.Set(*it,0.0);
  wSolOld.Set(*it,0.0);
 }
}

// impoe velocidade ALE = 0 no contorno
//-------------------------------------------------- 
void Simulator3D::setALEBC()
{
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
 {
  uALE.Set(*it,0.0);
  vALE.Set(*it,0.0);
  wALE.Set(*it,0.0);
 }
}

void Simulator3D::setAnnularALEBC()
{
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
 {
  if( surfMesh->phyBounds.at(*it) == "\"wallNormalU\"" )
   uALE.Set(*it,0.0);
  else if( surfMesh->phyBounds.at(*it) == "\"wallNormalV\"" )
   vALE.Set(*it,0.0);
  else if( surfMesh->phyBounds.at(*it) == "\"wallNormalW\"" )
   wALE.Set(*it,0.0);
  else
  {
   uALE.Set(*it,0.0);
   vALE.Set(*it,0.0);
   wALE.Set(*it,0.0);
  }
 }
}

/* method to copy the pointer of each attribute of _m to the Simulator3D
 * 
 * input: &_m
 * output: numVerts,Elems,Nodes,X,Y,Z,uc,vc,wc,pc,cc,heaviside etc.
 *
 * */
void Simulator3D::getModel3DAttrib(Model3D &_m)
{
 m = &_m;
 // mesh information vectors
 numVerts = m->getNumVerts();
 numElems = m->getNumElems();
 numNodes = m->getNumNodes();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 heaviside = m->getHeaviside();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 outflow = m->getOutflow();
 IEN = m->getIEN();
 curvature = m->getCurvature();
 surface = m->getSurface();
 surfMesh = m->getSurfMesh();
 mesh3d = m->getMesh3d();
 interfaceDistance = m->getInterfaceDistance();
 elemIdRegion = m->getElemIdRegion();
 boundaryVert = m->getBoundaryVert();
 triEdge = m->getTriEdge();
}


void Simulator3D::allocateMemoryToAttrib()
{
 // assembly matrix
 K.Dim( 3*numNodes,3*numNodes );
 Kc.Dim( numVerts,numVerts );
 Mrho.Dim( 3*numNodes,3*numNodes );
 M.Dim( 3*numNodes,3*numNodes );
 Mc.Dim( numVerts,numVerts );
 MrhoLumped.Dim( 3*numNodes );
 MLumped.Dim( 3*numNodes );
 McLumped.Dim( numVerts );
 G.Dim( 3*numNodes,numVerts );
 Gc.Dim( 3*numVerts,numVerts );
 D.Dim( numVerts,3*numNodes );
 gx.Dim( numNodes,numVerts );
 gy.Dim( numNodes,numVerts );
 gz.Dim( numNodes,numVerts );

 // COUPLED method matrix and vector
 A.Dim( 3*numNodes+numVerts,3*numNodes+numVerts );
 b.Dim( 3*numNodes+numVerts );

 // right hand side vectors
 va.Dim( 3*numNodes );
 vcc.Dim( numVerts );
 b1.Dim( 3*numNodes );
 b1c.Dim( numVerts );
 b2.Dim( numVerts );

 // boundary condition configured matrix
 ATilde.Dim( 3*numNodes,3*numNodes );
 AcTilde.Dim( numVerts,numVerts );
 GTilde.Dim( 3*numNodes,numVerts );
 DTilde.Dim( numVerts,3*numNodes );
 ETilde.Dim( numVerts,numVerts );
 E.Dim( numVerts, numVerts );

 // K + M matrix set
 mat.Dim( 3*numNodes,3*numNodes );
 matc.Dim( numVerts,numVerts );
 invA.Dim( 3*numNodes );
 invMrhoLumped.Dim( 3*numNodes );
 invMLumped.Dim( 3*numNodes );
 invC.Dim( numVerts );
 invMcLumped.Dim( numVerts );

 // solution vectors 
 // vetores solucao
 uTilde.Dim( 3*numNodes );
 pTilde.Dim( numVerts );
 cTilde.Dim( numVerts );
 uSol.Dim( numNodes );
 vSol.Dim( numNodes );
 wSol.Dim( numNodes );
 pSol.Dim( numVerts );
 cSol.Dim( numVerts );

 // auxiliar vectors
 Fold.Dim( 3*numNodes+numVerts );
 ip.Dim( 3*numNodes,1 );
 ipc.Dim( numVerts,1 );

 // convective term vectors
 convUVW.Dim( 3*numNodes );
 convC.Dim( numVerts );

 // convective term vectors (semi-lagrangian)
 uSL.Dim( numNodes );
 vSL.Dim( numNodes );
 wSL.Dim( numNodes );
 cSL.Dim( numVerts );

 // convective term vectors (ALE)
 uALE.Dim( numNodes );
 vALE.Dim( numNodes );
 wALE.Dim( numNodes );
 uSmooth.Dim( numVerts );
 vSmooth.Dim( numVerts );
 wSmooth.Dim( numVerts );
 uSmoothCoord.Dim( numVerts );
 vSmoothCoord.Dim( numVerts );
 wSmoothCoord.Dim( numVerts );
 uSmoothSurface.Dim( numVerts );
 vSmoothSurface.Dim( numVerts );
 wSmoothSurface.Dim( numVerts );
 
 // oldSol vectors
 uSolOld.Dim( numNodes );
 vSolOld.Dim( numNodes );
 wSolOld.Dim( numNodes );
 pSolOld.Dim( numVerts );
 cSolOld.Dim( numVerts );
 uALEOld.Dim( numNodes );
 vALEOld.Dim( numNodes );
 wALEOld.Dim( numNodes );
 kappaOld.Dim( numNodes );
 fintOld.Dim( numNodes );
 gravityOld.Dim( numVerts );
 betaFlowLiqOld.Dim( numVerts );
 muOld.Dim( numVerts );
 rhoOld.Dim( numVerts );
 cpOld.Dim( numVerts );
 ktOld.Dim( numVerts );
 hSmoothOld.Dim( numVerts );
 heatFluxOld.Dim( numVerts );

 // interface vectors (two-phase)
 fint.Dim ( 3*numNodes );
 gravity.Dim( 3*numNodes );
 betaFlowLiq.Dim( 3*numNodes );
 mu.Dim( numVerts );
 rho.Dim( numVerts );
 cp.Dim( numVerts );
 kt.Dim( numVerts );
 kappa.Dim( numNodes );
 hSmooth.Dim( numVerts );
 heatFlux.Dim( numVerts );

 int numSurface = surfMesh->numInterfaces+1; 
 centroidVelX.resize(numSurface);
 centroidVelY.resize(numSurface);
 centroidVelZ.resize(numSurface);
 centroidVelXOld.resize(numSurface);
 centroidVelYOld.resize(numSurface);
 centroidVelZOld.resize(numSurface);
 centroidPosX.resize(numSurface);
 centroidPosY.resize(numSurface);
 centroidPosZ.resize(numSurface);
 centroidPosXOld.resize(numSurface);
 centroidPosYOld.resize(numSurface);
 centroidPosZOld.resize(numSurface);
}

void Simulator3D::setCentroidVelPos()
{
 clVector _uVel = uSol;
 clVector _vVel = vSol;
 clVector _wVel = wSol;

 int v = elemIdRegion->Max()+1;
 vector<double> velX(v,0);
 vector<double> velY(v,0);
 vector<double> velZ(v,0);
 vector<double> posX(v,0);
 vector<double> posY(v,0);
 vector<double> posZ(v,0);
 vector<double> volume(v,0);
 vector<double> sumVolume(v,0);
 vector<double> sumXVelVolume(v,0);
 vector<double> sumYVelVolume(v,0);
 vector<double> sumZVelVolume(v,0);
 vector<double> sumXPosVolume(v,0);
 vector<double> sumYPosVolume(v,0);
 vector<double> sumZPosVolume(v,0);

 list<int> *inElem;
 inElem = m->getInElem();
 //for( int elem=0;elem<surfMesh->numElems;elem ) 
 for (list<int>::iterator it=inElem->begin(); it!=inElem->end(); ++it)
 {
  int v1 = IEN->Get(*it,0);
  int v2 = IEN->Get(*it,1);
  int v3 = IEN->Get(*it,2);
  int v4 = IEN->Get(*it,3);

  int elemID = elemIdRegion->Get(*it);

  velX[elemID] = ( _uVel.Get(v1)+
	               _uVel.Get(v2)+
				   _uVel.Get(v3)+
				   _uVel.Get(v4) )/4.0;

  velY[elemID] = ( _vVel.Get(v1)+
                   _vVel.Get(v2)+
				   _vVel.Get(v3)+
				   _vVel.Get(v4) )/4.0;

  velZ[elemID] = ( _wVel.Get(v1)+
                   _wVel.Get(v2)+
				   _wVel.Get(v3)+
				   _wVel.Get(v4) )/4.0;

  posX[elemID] = ( X->Get(v1)+
                   X->Get(v2)+
				   X->Get(v3)+
				   X->Get(v4) )/4.0;

  posY[elemID] = ( Y->Get(v1)+
                   Y->Get(v2)+
				   Y->Get(v3)+
				   Y->Get(v4) )/4.0;

  posZ[elemID] = ( Z->Get(v1)+
                   Z->Get(v2)+
				   Z->Get(v3)+
				   Z->Get(v4) )/4.0;

  volume[elemID] = m->getVolume(*it);

  sumXVelVolume[elemID] += velX[elemID] * volume[elemID];
  sumYVelVolume[elemID] += velY[elemID] * volume[elemID];
  sumZVelVolume[elemID] += velZ[elemID] * volume[elemID];

  sumXPosVolume[elemID] += posX[elemID] * volume[elemID];
  sumYPosVolume[elemID] += posY[elemID] * volume[elemID];
  sumZPosVolume[elemID] += posZ[elemID] * volume[elemID];
 }

 centroidVelX.clear();centroidVelY.clear();centroidVelZ.clear();
 centroidPosX.clear();centroidPosY.clear();centroidPosZ.clear();
 vector<double> surfaceVolume = m->getSurfaceVolume();
 for( int nb=0;nb<v;nb++ )
 {
  centroidVelX.push_back(sumXVelVolume[nb]/surfaceVolume[nb]);
  centroidVelY.push_back(sumYVelVolume[nb]/surfaceVolume[nb]);
  centroidVelZ.push_back(sumZVelVolume[nb]/surfaceVolume[nb]);

  centroidPosX.push_back(sumXPosVolume[nb]/surfaceVolume[nb]);
  centroidPosY.push_back(sumYPosVolume[nb]/surfaceVolume[nb]);
  centroidPosZ.push_back(sumZPosVolume[nb]/surfaceVolume[nb]);
 }
}

double Simulator3D::getCentroidVelXAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelX[nb];
 return sum/v;
}

double Simulator3D::getCentroidVelYAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelY[nb];
 return sum/v;
}

double Simulator3D::getCentroidVelZAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelZ[nb];
 return sum/v;
}

double Simulator3D::getCentroidPosXAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosX[nb];
 return sum/v;
}

double Simulator3D::getCentroidPosYAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosY[nb];
 return sum/v;
}

double Simulator3D::getCentroidPosZAverage()
{
 double sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosZ[nb];
 return sum/v;
}

void Simulator3D::setCentroidVelPosInterface()
{
 clVector _uVel = uALEOld;
 clVector _vVel = vALEOld;
 clVector _wVel = wALEOld;

 int v = surfMesh->elemIdRegion.Max()+1;
 vector<double> velX(v,0);
 vector<double> velY(v,0);
 vector<double> velZ(v,0);
 vector<double> posX(v,0);
 vector<double> posY(v,0);
 vector<double> posZ(v,0);
 vector<double> area(v,0);
 vector<double> sumArea(v,0);
 vector<double> sumXVelArea(v,0);
 vector<double> sumYVelArea(v,0);
 vector<double> sumZVelArea(v,0);
 vector<double> sumXPosArea(v,0);
 vector<double> sumYPosArea(v,0);
 vector<double> sumZPosArea(v,0);

 for( int elem=0;elem<surfMesh->numElems;elem++ ) 
 {
  // P1
  int v1 = surfMesh->IEN.Get(elem,0);
  double p1x = surfMesh->X.Get(v1);
  double p1y = surfMesh->Y.Get(v1);
  double p1z = surfMesh->Z.Get(v1);

  // P2
  int v2 = surfMesh->IEN.Get(elem,1);
  double p2x = surfMesh->X.Get(v2);
  double p2y = surfMesh->Y.Get(v2);
  double p2z = surfMesh->Z.Get(v2);

  // P3
  int v3 = surfMesh->IEN.Get(elem,2);
  double p3x = surfMesh->X.Get(v3);
  double p3y = surfMesh->Y.Get(v3);
  double p3z = surfMesh->Z.Get(v3);

  int elemID = surfMesh->elemIdRegion.Get(elem);

  velX[elemID] = ( _uVel.Get(v1)+
	               _uVel.Get(v2)+
				   _uVel.Get(v3) )/3.0;

  velY[elemID] = ( _vVel.Get(v1)+
                   _vVel.Get(v2)+
				   _vVel.Get(v3) )/3.0;

  velZ[elemID] = ( _wVel.Get(v1)+
                   _wVel.Get(v2)+
				   _wVel.Get(v3) )/3.0;

  posX[elemID] = ( surfMesh->X.Get(v1)+
                   surfMesh->X.Get(v2)+
				   surfMesh->X.Get(v3) )/3.0;

  posY[elemID] = ( surfMesh->Y.Get(v1)+
                   surfMesh->Y.Get(v2)+
				   surfMesh->Y.Get(v3) )/3.0;

  posZ[elemID] = ( surfMesh->Z.Get(v1)+
                   surfMesh->Z.Get(v2)+
				   surfMesh->Z.Get(v3) )/3.0;

  area[elemID] = getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z);

  sumXVelArea[elemID] += velX[elemID] * area[elemID];
  sumYVelArea[elemID] += velY[elemID] * area[elemID];
  sumZVelArea[elemID] += velZ[elemID] * area[elemID];

  sumXPosArea[elemID] += posX[elemID] * area[elemID];
  sumYPosArea[elemID] += posY[elemID] * area[elemID];
  sumZPosArea[elemID] += posZ[elemID] * area[elemID];
 }

 centroidVelX.clear();centroidVelY.clear();centroidVelZ.clear();
 centroidPosX.clear();centroidPosY.clear();centroidPosZ.clear();
 vector<double> surfaceArea = m->getSurfaceArea();
 for( int nb=0;nb<v;nb++ )
 {
  centroidVelX.push_back(sumXVelArea[nb]/surfaceArea[nb]);
  centroidVelY.push_back(sumYVelArea[nb]/surfaceArea[nb]);
  centroidVelZ.push_back(sumZVelArea[nb]/surfaceArea[nb]);

  centroidPosX.push_back(sumXPosArea[nb]/surfaceArea[nb]);
  centroidPosY.push_back(sumYPosArea[nb]/surfaceArea[nb]);
  centroidPosZ.push_back(sumZPosArea[nb]/surfaceArea[nb]);
 }
}

double Simulator3D::getCentroidVelXMax()
{
 return *max_element(centroidVelX.begin()+1,centroidVelX.end());
}
double Simulator3D::getCentroidVelYMax()
{
 return *max_element(centroidVelY.begin()+1,centroidVelY.end());
}
double Simulator3D::getCentroidVelZMax()
{
 return *max_element(centroidVelZ.begin()+1,centroidVelZ.end());
}
double Simulator3D::getCentroidVelXMin()
{
 return *min_element(centroidVelX.begin()+1,centroidVelX.end());
}
double Simulator3D::getCentroidVelYMin()
{
 return *min_element(centroidVelY.begin()+1,centroidVelY.end());
}
double Simulator3D::getCentroidVelZMin()
{
 return *min_element(centroidVelZ.begin()+1,centroidVelZ.end());
}

void Simulator3D::setUSol(double _vel)
{
 for( int i=0;i<numNodes;i++ )
 {
  uSol.Set(i,uSol.Get(i)-_vel);
  uSolOld.Set(i,uSolOld.Get(i)-_vel);
 }
}
void Simulator3D::setVSol(double _vel)
{
 for( int i=0;i<numNodes;i++ )
 {
  vSol.Set(i,vSol.Get(i)-_vel);
  vSolOld.Set(i,vSolOld.Get(i)-_vel);
 }
}
void Simulator3D::setWSol(double _vel)
{
 for( int i=0;i<numNodes;i++ )
 {
  wSol.Set(i,wSol.Get(i)-_vel);
  wSolOld.Set(i,wSolOld.Get(i)-_vel);
 }
}

void Simulator3D::setSurfaceTSat()
{
 double Tsat = 0.0;
 for( int i=0;i<surface->Dim();i++ )
 {
  int surfaceNode = surface->Get(i);
  cSol.Set(surfaceNode,Tsat);
  cSolOld.Set(surfaceNode,Tsat);
 }
}

//--------------------------------------------------
// void Simulator3D::setMassTransferVelocity()
// {
//  for( int i=0;i<surface->Dim();i++ )
//  {
//   int surfaceNode = surface->Get(i);
// 
//   // unitario do vetor normal (ponderado com a area) resultante
//   double xNormalUnit = surfMesh->xNormal.Get(surfaceNode);
//   double yNormalUnit = surfMesh->yNormal.Get(surfaceNode);
//   double zNormalUnit = surfMesh->zNormal.Get(surfaceNode);
// 
//   double uSurface = xNormalUnit;
//   double vSurface = yNormalUnit;
//   double wSurface = zNormalUnit;
// 
//   uMassTransfer.Set(surfaceNode,uSurface);
//   vMassTransfer.Set(surfaceNode,vSurface);
//   wMassTransfer.Set(surfaceNode,wSurface);
//  }
// } // fecha metodo setInterfaceVelocity
//-------------------------------------------------- 

void Simulator3D::setMassTransfer()
{
 setSurfaceTSat();

 // compute q
 //clVector q(numVerts);q.SetAll(-0.1);
 //clVector q = (-1.0)*invMcLumped*Kc*cSolOld;
 clVector q = (1.0)*Kc*cSolOld;
//--------------------------------------------------
//  for( int i=0;i<surface->Dim();i++ )
//  {
//   int surfaceNode = surface->Get(i);
//    cout << cSolOld.Get(surfaceNode) << " " << q.Get(surfaceNode) << endl;
//  }
//-------------------------------------------------- 

 clDMatrix distrib(numVerts);
 clVector* closer = m->getCloser();
 for( int i=0;i<numVerts;i++ )
 {
  int aux = closer->Get(i);
  distrib.Set( i,q.Get(aux) );
 }

 clVector GH = Gc*(*heaviside);
 //clVector GT = Gc*(cSolOld);
 //clVector *vertIdRegion = m->getVertIdRegion();
 for( int i=0;i<numVerts;i++ )
 {
  double aux = sqrt( GH.Get(i)*GH.Get(i) + 
                   GH.Get(i+numVerts)*GH.Get(i+numVerts) +
                   GH.Get(i+2*numVerts)*GH.Get(i+2*numVerts) );
  heatFlux.Set(i,distrib.Get(i)*aux);
  //heatFlux.Set(i,1.0*aux);


//--------------------------------------------------
// // Uni-dimensional problem in Y
// //-------------------------------------------------- 
//   //if( Y->Get(i) < Y->Max() )
//   if( vertIdRegion->Get(i) == 1 && i < lineMesh->numVerts )
//   {
//    double aux = fabs( GH.Get(i+numVerts) );
//    heatFlux.Set(i,0.0*aux);
//   }
//-------------------------------------------------- 
 }


 /* 
  * Interface velocity + Mesh Smooth
  *
  * */

 // smoothing - coordenadas
 MeshSmooth e1(*m,dt); // criando objeto MeshSmooth
 e1.stepSurfaceSmoothFujiwara();
 uSmoothSurface = *e1.getUSmoothSurface();
 vSmoothSurface = *e1.getVSmoothSurface();
 wSmoothSurface = *e1.getWSmoothSurface();

 for( int i=0;i<surface->Dim();i++ )
 {
  int surfaceNode = surface->Get(i);

  // unitario do vetor normal (ponderado com a area) resultante
  double xNormalUnit = surfMesh->xNormal.Get(surfaceNode);
  double yNormalUnit = surfMesh->yNormal.Get(surfaceNode);
  double zNormalUnit = surfMesh->zNormal.Get(surfaceNode);

  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  double prod = (uSolOld.Get(surfaceNode)+1.3*uRef)*xNormalUnit+ 
              (vSolOld.Get(surfaceNode)+1.3*vRef)*yNormalUnit + 
			  (wSolOld.Get(surfaceNode)+1.3*wRef)*zNormalUnit;
  double uSolNormal = xNormalUnit*prod;
  double vSolNormal = yNormalUnit*prod;
  double wSolNormal = zNormalUnit*prod;

  // 1.3 is a pragmatic number which fits the velocity for the bhaga5
  // and the moving frame technique. Still don't know why!
  double uSolTangent = uSolOld.Get(surfaceNode) + 1.3*uRef - uSolNormal;
  double vSolTangent = vSolOld.Get(surfaceNode) + 1.3*vRef - vSolNormal;
  double wSolTangent = wSolOld.Get(surfaceNode) + 1.3*wRef - wSolNormal;

  // tratamento da superficie
  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  double prod2 = uSmoothSurface.Get(surfaceNode)*xNormalUnit + 
               vSmoothSurface.Get(surfaceNode)*yNormalUnit + 
			   wSmoothSurface.Get(surfaceNode)*zNormalUnit;
  double uSmoothNormal = xNormalUnit*prod2;
  double vSmoothNormal = yNormalUnit*prod2;
  double wSmoothNormal = zNormalUnit*prod2;

  double uSmoothTangent = uSmoothSurface.Get(surfaceNode) - uSmoothNormal;
  double vSmoothTangent = vSmoothSurface.Get(surfaceNode) - vSmoothNormal;
  double wSmoothTangent = wSmoothSurface.Get(surfaceNode) - wSmoothNormal;

  double uALESurface =   uSolOld.Get(surfaceNode) 
                     - d1*uSolTangent 
					 + d2*uSmoothTangent;
  double vALESurface =   vSolOld.Get(surfaceNode) 
                     - d1*vSolTangent 
					 + d2*vSmoothTangent;
  double wALESurface =   wSolOld.Get(surfaceNode) 
                     - d1*wSolTangent 
					 + d2*wSmoothTangent;


  // MASS TRASNFER at SURFACE NODES
  double rho1 = ( (1.0/(rho_inAdimen*rho_inAdimen))+
                (1.0/(rho_outAdimen*rho_outAdimen)) )/
                ( (1.0/rho_inAdimen)+(1.0/rho_outAdimen) ) ;

  double rho2 = 0.5*( 1.0/rho_inAdimen + 1.0/rho_outAdimen );

  cout << surfaceNode << " " 
       << q.Get(surfaceNode) << " "
       << rho1 << " "
	   << rho2 << endl;

  double uMassTransfer = uALESurface
                       - (q.Get(surfaceNode)*xNormalUnit)*rho2;

  double vMassTransfer = vALESurface
                       - (q.Get(surfaceNode)*yNormalUnit)*rho2;

  double wMassTransfer = wALESurface
                       - (q.Get(surfaceNode)*zNormalUnit)*rho2;

  uALE.Set(surfaceNode,uMassTransfer);
  vALE.Set(surfaceNode,vMassTransfer);
  wALE.Set(surfaceNode,wMassTransfer);
 }
} // fecha metodo setMassTransfer


// PBC
void Simulator3D::getPeriodic3DToAttrib(Periodic3D &_pbc)
{
	pbc = &_pbc;
	nyPointsL = pbc->GetNyPointsL();
	nyPointsM = pbc->GetNyPointsM();
	NumVertsMid = pbc->GetNumVertsMid();
	VecXMin = pbc->GetVecXMin();
	VecXMax = pbc->GetVecXMax();
	VecXMid = pbc->GetVecXMid();
	VecXMidVerts = pbc->GetVecXMidVerts();
	VecXMinGlob.Dim(0);
	VecXMaxGlob.Dim(0);
	betaFlowLiq.Dim(3*numNodes);

	MasterIndices = pbc->GetMasterIndices();
	SlaveIndices = pbc->GetSlaveIndices();

} // fecha metodo


/* VelPressMatrixModifierPBC() */
void Simulator3D::VelPressMatrixModifierPBC()
{
   size_t i,j;
   int ibL,ibR;

   // indices to which PBC will be applied
   vector<int> indxMasterToModify(0);
   vector<int> indxSlaveToModify(0);
   
   // finding PBC nodes set in Model
   for ( list<int>::iterator it = boundaryVert->begin(); it != boundaryVert->end(); ++it )
   {   
	if ( surfMesh->phyBounds.at(*it) == "\"wallLeft\"" )
	{
	  indxMasterToModify.push_back(*it);	  
	}

    if ( surfMesh->phyBounds.at(*it) == "\"wallRight\"" )
	{
	  indxSlaveToModify.push_back(*it);	  
	}
  }

   // indxMasterSlave must have the same size. Hence, ny suffices.
   size_t ny = indxMasterToModify.size();

  /* Copying rows and columns from ibL into ibR:
   * i) Remove the contributions from left and overloads in right;
   * ii) Now, the positions that stayed opened receive the same 
   * values at right. */
  if ( direction == "RL" )
  {
	  // ATilde
	  for ( i = 0; i < ny; i++ ) // loop paired points
	  {
	     // left and right nodes
	     ibL = indxMasterToModify.at(i);
	     ibR = indxSlaveToModify.at(i);

		  for ( j = 0; j < 3*numNodes; j++ ) // loop rows
		  {
			 // x-direction
			 double ATildeRow = ATilde.Get(ibR,j);
			 ATildeRow += ATilde.Get(ibL,j);
			 ATilde.Set(ibR,j,ATildeRow);
			  
			 // y-direction
			 double ATildeRowN = ATilde.Get(ibR + numNodes,j);
			 ATildeRowN += ATilde.Get(ibL + numNodes,j);
			 ATilde.Set(ibR + numNodes,j,ATildeRowN);
			 
			 // z-direction
			 ATildeRowN = ATilde.Get(ibR + 2*numNodes,j);
			 ATildeRowN += ATilde.Get(ibL + 2*numNodes,j);
			 ATilde.Set(ibR + 2*numNodes,j,ATildeRowN);
	      }

		  for ( j = 0; j < 3*numNodes; j++ ) // loop columns
		  {
		   	 // x-direction
			 double ATildeColumn = ATilde.Get(j,ibR);
			 ATildeColumn += ATilde.Get(j,ibL);
			 ATilde.Set(j,ibR,ATildeColumn);

		   	 // y-direction
			 double ATildeColumnN = ATilde.Get(j,ibR + numNodes);
			 ATildeColumnN += ATilde.Get(j,ibL + numNodes);
			 ATilde.Set(j,ibR + numNodes,ATildeColumnN);

		   	 // z-direction
			 ATildeColumnN = ATilde.Get(j,ibR + 2*numNodes);
			 ATildeColumnN += ATilde.Get(j,ibL + 2*numNodes);
			 ATilde.Set(j,ibR + 2*numNodes,ATildeColumnN);
		  }
	  }
	
	  for ( i = 0; i < ny; i++ )
	  {

		 // eliminating rows and columns
		 ibL = indxMasterToModify.at(i);
		 ATilde.SetRowZero(ibL);
		 ATilde.SetColumnZero(ibL);
		 ATilde.SetRowZero(ibL + numNodes);
		 ATilde.SetColumnZero(ibL + numNodes);
		 ATilde.SetRowZero(ibL + 2*numNodes);
		 ATilde.SetColumnZero(ibL + 2*numNodes);

		 // changed diagonal's values
	     ATilde.Set(ibL,ibL,1.0);
	     ATilde.Set(ibL + numNodes,ibL + numNodes,1.0);
	     ATilde.Set(ibL + 2*numNodes,ibL + 2*numNodes,1.0);

	  }

	  // DTilde
	  for ( i = 0; i < ny; i++ )
	  {
	     ibL = indxMasterToModify.at(i);
		 ibR = indxSlaveToModify.at(i);
		
		   for ( j = 0; j < 3*numNodes; j++ ) // loop rows
		   {
				double DTildeRow = DTilde.Get(ibR,j);
				DTildeRow += DTilde.Get(ibL,j);
				DTilde.Set(ibR,j,DTildeRow);
		   }
	  
           for ( j = 0; j < numVerts; j++ ) // loop columns
	       {
				double DTildeColumn = DTilde.Get(j,ibR);
				DTildeColumn += DTilde.Get(j,ibL);
				DTilde.Set(j,ibR,DTildeColumn);

				DTildeColumn = DTilde.Get(j, ibR + numNodes);
				DTildeColumn += DTilde.Get(j, ibL + numNodes);
				DTilde.Set(j, ibR + numNodes, DTildeColumn);

				DTildeColumn = DTilde.Get(j, ibR + 2*numNodes);
				DTildeColumn += DTilde.Get(j, ibL + 2*numNodes);
				DTilde.Set(j, ibR + 2*numNodes, DTildeColumn);
	       }

	  }	
	  
	  for ( i = 0; i < ny; i++ )
	  {

		  ibL = indxMasterToModify.at(i);
		  DTilde.SetRowZero(ibL);
		  DTilde.SetColumnZero(ibL);
		  DTilde.SetColumnZero(ibL + numNodes);
		  DTilde.SetColumnZero(ibL + 2*numNodes);

	  }

	  // GTilde
	  for ( i = 0; i < ny; i++ )
	  {
	      ibL = indxMasterToModify.at(i);
		  ibR = indxSlaveToModify.at(i);

		     for ( j = 0; j < numVerts; j++ ) // loop rows
			 {
				double GTildeRow = GTilde.Get(ibR,j);
				GTildeRow += GTilde.Get(ibL,j);
				GTilde.Set(ibR,j,GTildeRow);

			 }
	     
			 for ( j = 0; j < numNodes; j++ ) // loop columns
			 {
			    double GTildeColumn = GTilde.Get(j,ibR);
				GTildeColumn += GTilde.Get(j,ibL);
				GTilde.Set(j,ibR,GTildeColumn);

				GTildeColumn = GTilde.Get(j + numNodes,ibR);
				GTildeColumn += GTilde.Get(j + numNodes,ibL);
				GTilde.Set(j + numNodes,ibR,GTildeColumn);
				
				GTildeColumn = GTilde.Get(j + 2*numNodes,ibR);
				GTildeColumn += GTilde.Get(j + 2*numNodes,ibL);
				GTilde.Set(j + 2*numNodes,ibR,GTildeColumn);

			 }

	  }

	  for ( i = 0; i < ny; i++ )
	  {

		 ibL = indxMasterToModify.at(i);
		 GTilde.SetRowZero(ibL);
		 GTilde.SetRowZero(ibL + numNodes);
		 GTilde.SetRowZero(ibL + 2*numNodes);
		 GTilde.SetColumnZero(ibL);
		 GTilde.SetColumnZero(ibL + numNodes);
		 GTilde.SetColumnZero(ibL + 2*numNodes);


	  }
	  
	  //*** ETilde call
	  ETilde = E - (( DTilde * invA) * GTilde );

	  for ( i = 0; i < ny; i++ )
	  {
		  ibL = indxMasterToModify.at(i);
		  ETilde.SetRowZero(ibL);
		  ETilde.SetColumnZero(ibL);
		  ETilde.Set(ibL,ibL,1.0);		  
	  }

	}

	/* Copying rows and columns from ibR to ibL. Idem, with inversed
	 * indices. */
	else
	{
		/* copy block! */
	}

} // fecha metodo


void Simulator3D::assemblePBC()
{
 	register int i,j;
	int ibL,ibR;

	clVector VecXMinAux, VecXMaxAux;

	VecXMinAux.Dim(nyPointsL);
	VecXMaxAux.Dim(nyPointsL);

	VecXMin->CopyTo(0,VecXMinAux);
	VecXMax->CopyTo(0,VecXMaxAux);

	/* Copying rows and columns from ibL into ibR:
	 * i) Remove the contributions from left and overloads in right;
	 * ii) Now, the positions that stayed opened receive the same 
	 * values at right. */

	if ( direction == "RL" )
	{
	// ATilde
	  for ( i = 0; i < nyPointsL; i++ ) // loop paired points
	  {
	     // left and right nodes
	     ibL = VecXMinAux.Get(i);
	     ibR = VecXMaxAux.Get(i);

		  for ( j = 0; j < 3*numNodes; j++ ) // loop rows
		  {
			 // x-direction
			 double ATildeRow = ATilde.Get(ibR,j);
			 ATildeRow += ATilde.Get(ibL,j);
			 ATilde.Set(ibR,j,ATildeRow);
			  
			 // y-direction
			 double ATildeRowN = ATilde.Get(ibR + numNodes,j);
			 ATildeRowN += ATilde.Get(ibL + numNodes,j);
			 ATilde.Set(ibR + numNodes,j,ATildeRowN);
			 
			 // z-direction
			 ATildeRowN = ATilde.Get(ibR + 2*numNodes,j);
			 ATildeRowN += ATilde.Get(ibL + 2*numNodes,j);
			 ATilde.Set(ibR + 2*numNodes,j,ATildeRowN);
	      }

		  for ( j = 0; j < 3*numNodes; j++ ) // loop columns
		  {
		   	 // x-direction
			 double ATildeColumn = ATilde.Get(j,ibR);
			 ATildeColumn += ATilde.Get(j,ibL);
			 ATilde.Set(j,ibR,ATildeColumn);

		   	 // y-direction
			 double ATildeColumnN = ATilde.Get(j,ibR + numNodes);
			 ATildeColumnN += ATilde.Get(j,ibL + numNodes);
			 ATilde.Set(j,ibR + numNodes,ATildeColumnN);

		   	 // z-direction
			 ATildeColumnN = ATilde.Get(j,ibR + 2*numNodes);
			 ATildeColumnN += ATilde.Get(j,ibL + 2*numNodes);
			 ATilde.Set(j,ibR + 2*numNodes,ATildeColumnN);
		  }
	  }
	
	  for ( i = 0; i < nyPointsL; i++ )
	  {

		 // eliminating rows and columns
		 ibL = VecXMinAux.Get(i);
		 ATilde.SetRowZero(ibL);
		 ATilde.SetColumnZero(ibL);
		 ATilde.SetRowZero(ibL + numNodes);
		 ATilde.SetColumnZero(ibL + numNodes);
		 ATilde.SetRowZero(ibL + 2*numNodes);
		 ATilde.SetColumnZero(ibL + 2*numNodes);

		 // changed diagonal's values
	     ATilde.Set(ibL,ibL,1.0);
	     ATilde.Set(ibL + numNodes,ibL + numNodes,1.0);
	     ATilde.Set(ibL + 2*numNodes,ibL + 2*numNodes,1.0);

	  }

	  // DTilde
	  for ( i = 0; i < nyPointsL; i++ )
	  {
	     ibL = VecXMinAux.Get(i);
		 ibR = VecXMaxAux.Get(i);
		
		   for ( j = 0; j < 3*numNodes; j++ ) // loop rows
		   {
				double DTildeRow = DTilde.Get(ibR,j);
				DTildeRow += DTilde.Get(ibL,j);
				DTilde.Set(ibR,j,DTildeRow);
		   }
	  
           for ( j = 0; j < numVerts; j++ ) // loop columns
	       {
				double DTildeColumn = DTilde.Get(j,ibR);
				DTildeColumn += DTilde.Get(j,ibL);
				DTilde.Set(j,ibR,DTildeColumn);

				DTildeColumn = DTilde.Get(j, ibR + numNodes);
				DTildeColumn += DTilde.Get(j, ibL + numNodes);
				DTilde.Set(j, ibR + numNodes, DTildeColumn);

				DTildeColumn = DTilde.Get(j, ibR + 2*numNodes);
				DTildeColumn += DTilde.Get(j, ibL + 2*numNodes);
				DTilde.Set(j, ibR + 2*numNodes, DTildeColumn);
	       }

	  }	
	  
	  for ( i = 0; i < nyPointsL; i++ )
	  {

		  ibL = VecXMinAux.Get(i);
		  DTilde.SetRowZero(ibL);
		  DTilde.SetColumnZero(ibL);
		  DTilde.SetColumnZero(ibL + numNodes);
		  DTilde.SetColumnZero(ibL + 2*numNodes);

	  }

	  // GTilde
	  for ( i = 0; i < nyPointsL; i++ )
	  {
	      ibL = VecXMinAux.Get(i);
		  ibR = VecXMaxAux.Get(i);

		     for ( j = 0; j < numVerts; j++ ) // loop rows
			 {
				double GTildeRow = GTilde.Get(ibR,j);
				GTildeRow += GTilde.Get(ibL,j);
				GTilde.Set(ibR,j,GTildeRow);

			 }
	     
			 for ( j = 0; j < numNodes; j++ ) // loop columns
			 {
			    double GTildeColumn = GTilde.Get(j,ibR);
				GTildeColumn += GTilde.Get(j,ibL);
				GTilde.Set(j,ibR,GTildeColumn);

				GTildeColumn = GTilde.Get(j + numNodes,ibR);
				GTildeColumn += GTilde.Get(j + numNodes,ibL);
				GTilde.Set(j + numNodes,ibR,GTildeColumn);
				
				GTildeColumn = GTilde.Get(j + 2*numNodes,ibR);
				GTildeColumn += GTilde.Get(j + 2*numNodes,ibL);
				GTilde.Set(j + 2*numNodes,ibR,GTildeColumn);

			 }

	  }

	  for ( i = 0; i < nyPointsL; i++ )
	  {

		 ibL = VecXMinAux.Get(i);
		 GTilde.SetRowZero(ibL);
		 GTilde.SetRowZero(ibL + numNodes);
		 GTilde.SetRowZero(ibL + 2*numNodes);
		 GTilde.SetColumnZero(ibL);
		 GTilde.SetColumnZero(ibL + numNodes);
		 GTilde.SetColumnZero(ibL + 2*numNodes);


	  }
	  
	  //*** ETilde call
	  ETilde = E - (( DTilde * invA) * GTilde );

	  for ( i = 0; i < nyPointsL; i++ )
	  {
		  ibL = VecXMinAux.Get(i);
		  ETilde.SetRowZero(ibL);
		  ETilde.SetColumnZero(ibL);
		  ETilde.Set(ibL,ibL,1.0);		  
	  }

	}

	/* Copying rows and columns from ibR to ibL. Idem, with inversed
	 * indices. */
	else
	{
		/* copy block! */
	}

} // fecha metodo


/* Method should be used after SetCopyDirectionPBC() */
void Simulator3D::assembleCPBC()
{
    register int i,j;
	int ibL,ibR;
  	
	clVector VecXMinAux, VecXMaxAux;

  	VecXMinAux.Dim(nyPointsL);
  	VecXMaxAux.Dim(nyPointsL);
    
	VecXMin->CopyTo(0,VecXMinAux);
  	VecXMax->CopyTo(0,VecXMaxAux);

    /* Copying rows and columns from ibL into ibR: 
	 * i) Remove the contributions from left and overloads in right;
	 * ii) Now, the positions that stayed opened receive the same
	 * values at right. */
	  
	if ( direction == "RL" )
	{	 
	// AcTilde
      for ( i = 0; i < nyPointsL; i++ ) // loop paired points
      {
	    // left and right nodes
        ibL = VecXMinAux.Get(i);
        ibR = VecXMaxAux.Get(i);
        
          for( j = 0; j < numVerts; j++ ) // loop rows
          {
              double AcTildeRow = AcTilde.Get(ibR,j);
              AcTildeRow += AcTilde.Get(ibL,j);
              AcTilde.Set(ibR,j,AcTildeRow);
            
		  }
        
          for( j = 0; j < numVerts; j++ ) // loop columns
          {
              double AcTildeColumn = AcTilde.Get(j,ibR);
              AcTildeColumn += AcTilde.Get(j,ibL);
              AcTilde.Set(j,ibR,AcTildeColumn);
          }
      }
    
	  for ( i = 0; i < nyPointsL; i++  )	 
	  {

		// eliminating rows and columns
        ibL = VecXMinAux.Get(i);
        AcTilde.SetRowZero(ibL);
        AcTilde.SetColumnZero(ibL);
		
		// changed diagonal's values
		AcTilde.Set(ibL,ibL,1.0);

	  }
    }

	/* Copying rows and columns from ibR to ibL. Idem, with inversed
	 * indices. */
	else
	{
	// AcTilde
      for ( i = 0; i < nyPointsL; i++ ) // loop paired points
      {
	    // left and right nodes
        ibL = VecXMinAux.Get(i);
        ibR = VecXMaxAux.Get(i);
        
          for( j = 0; j < numVerts; j++ ) // loop rows
          {
              double AcTildeRow = AcTilde.Get(ibL,j);
              AcTildeRow += AcTilde.Get(ibR,j);
              AcTilde.Set(ibL,j,AcTildeRow);
            
          }
        
          for( j = 0; j < numVerts; j++ ) // loop columns
          {
              double AcTildeColumn = AcTilde.Get(j,ibL);
              AcTildeColumn += AcTilde.Get(j,ibR);
              AcTilde.Set(j,ibL,AcTildeColumn);
            
          }
      }
    
	  for ( i = 0; i < nyPointsL; i++  )	 
	  {

		// eliminating rows and columns
        ibR = VecXMaxAux.Get(i);
        AcTilde.SetRowZero(ibR);
        AcTilde.SetColumnZero(ibR);
		
		// changed diagonal's values
		AcTilde.Set(ibR,ibR,1.0);
	  }
    
     }

} // close method assembleCPBC


/* unCoupledPBC with <vector> structure */
void Simulator3D::unCoupledPBCVector()
{
 clVector uvw(3*numNodes);
 clVector vaIp(3*numNodes);
 clVector b1Tilde;
 clVector b2Tilde;
 clVector uTildeU, uTildeV, uTildeW;
 Periodic3D pbc;

 vaIp = va.MultVec(ip); // operacao vetor * vetor (elemento a elemento)

 // indices to which PBC will be applied
 vector<int> indxMasterToModify(0);
 vector<int> indxSlaveToModify(0);
   
 // finding PBC nodes set in Model
 for ( list<int>::iterator it = boundaryVert->begin(); it != boundaryVert->end(); ++it )
 {   
  if ( surfMesh->phyBounds.at(*it) == "\"wallLeft\"" )
  {
	  indxMasterToModify.push_back(*it);	  
  } 

  if ( surfMesh->phyBounds.at(*it) == "\"wallRight\"" )
  {
      indxSlaveToModify.push_back(*it);	  
  }
 }

 // reallocating indices PBC 
 MasterIndices = indxMasterToModify;
 SlaveIndices = indxSlaveToModify;
 nyPointsL = indxMasterToModify.size();

 b1Tilde = b1 + vaIp;

 //*** sums ibL to ibR on rhs vector - velocity
 sumIndexPBCVelVector(MasterIndices,SlaveIndices,b1Tilde);

 //*** Reassemble
 VelPressMatrixModifierPBC();

 // resolve sistema ATilde uTilde = b1Tilde
 cout << " --------> solving velocity --------- " << endl;
 solverV->solve(1E-15,ATilde,uTilde,b1Tilde);
 cout << " ------------------------------------ " << endl;

 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 //*** copy uTilde to set it periodic
 uTildeU = uTilde.Copy(0,numNodes - 1);
 uTildeV = uTilde.Copy(numNodes,2*numNodes - 1);
 uTildeW = uTilde.Copy(2*numNodes,3*numNodes - 1);
 pbc.SetVelocityPBCVector(uTildeU,uTildeV,uTildeW,MasterIndices,SlaveIndices,nyPointsL,"RL");
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 //*** updated periodic velocity
 uTilde = uTildeU;
 uTilde.Append(uTildeV);
 uTilde.Append(uTildeW);
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 if( rho_in <= rho_out ) // BUBBLE
 {
  cout << setw(70) << "BUBBLE SIMULATION" << endl;
  uvw = uTilde + 
        dt*invMrhoLumped*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		dt*invMLumped*( ((1.0/We)*(fint) ).MultVec(ip) );

 }
 else // DROPLET
 {
  cout << setw(70) << "DROPLET SIMULATION" << endl;
  uvw = uTilde + 
        invA*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		invA*( ((1.0/We)*(fint) ).MultVec(ip) );
 }
 
 //uvw = uTilde;

 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 /* 
  * Mass Transfer
  * */
 //clVector massTransfer = dt*invMcLumped*( heatFlux );
 //clVector massTransfer = ( invMcLumped*heatFlux );
 clVector massTransfer = invC*
                         (1.0/rho_inAdimen - 1.0/rho_outAdimen)
						 *heatFlux;

 ///*** setting b2 periodic, because DTilde*uvw already is.
 sumIndexPBCPressVector(MasterIndices,SlaveIndices,b2);
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 b2Tilde = (-1.0)*( b2 - (DTilde * uvw) ); 
 //b2Tilde = (-1.0)*( b2 - (DTilde * uvw) + (massTransfer) );
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 // resolve sistema E pTilde = b2Tilde
 cout << " --------> solving pressure --------- " << endl;
 solverP->solve(1E-15,ETilde,pTilde,b2Tilde);
 cout << " ------------------------------------ " << endl;
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;
 
 //*** copying pressure
 pbc.SetPurePressurePBCVector(pTilde,MasterIndices,SlaveIndices,nyPointsL,"RL");
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 uvw = uvw - (invA * GTilde * pTilde);
 cout << uTilde.Get(164) << endl;
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;
 uTildeU = uvw.Copy(0,numNodes - 1);
 uTildeV = uvw.Copy(numNodes,2*numNodes - 1);
 uTildeW = uvw.Copy(2*numNodes,3*numNodes - 1);

 cout << uTilde.Get(164) << endl;
 pbc.SetVelocityPBCVector(uTildeU,uTildeV,uTildeW,MasterIndices,SlaveIndices,nyPointsL,"RL");
 cout << uvw.Get(164) << endl;
 cout << getUSol()->Get(164) << endl;

 //*** updated periodic velocity
 uSol = uTildeU;
 vSol = uTildeV;
 wSol = uTildeW;

 cout << getUSol()->Get(164) << endl;

 pSol = pTilde;       // sem correcao na pressao
 //pSol = pSol + pTilde;  // com correcao na pressao

 // compute bubble's centroid velocity
 if( surfMesh->numInterfaces > 0 )
  setCentroidVelPos();
} // fecha metodo unCoupledPBC 



void Simulator3D::unCoupledPBC()
{
 clVector uvw(3*numNodes);
 clVector vaIp(3*numNodes);
 clVector b1Tilde;
 clVector b2Tilde;
 clVector uTildeU, uTildeV, uTildeW;
 Periodic3D pbc;

 vaIp = va.MultVec(ip); // operacao vetor * vetor (elemento a elemento)

 VecXMaxGlob.Dim(nyPointsL);
 VecXMax->CopyTo(0,VecXMaxGlob);
 VecXMinGlob.Dim(nyPointsL);
 VecXMin->CopyTo(0,VecXMinGlob);
 
 b1Tilde = b1 + vaIp;

 //*** sums ibL to ibR on rhs vector - velocity
 sumIndexPBCVel(VecXMin,VecXMax,b1Tilde);

 //*** Reassemble
 assemblePBC();

 // resolve sistema ATilde uTilde = b1Tilde
 cout << " --------> solving velocity --------- " << endl;
 solverV->solve(1E-15,ATilde,uTilde,b1Tilde);
 cout << " ------------------------------------ " << endl;

 //*** copy uTilde to set it periodic
 uTildeU = uTilde.Copy(0,numNodes - 1);
 uTildeV = uTilde.Copy(numNodes,2*numNodes - 1);
 uTildeW = uTilde.Copy(2*numNodes,3*numNodes - 1);
 pbc.SetVelocityPBC(uTildeU,uTildeV,uTildeW,VecXMinGlob,VecXMaxGlob,nyPointsL,"RL");

 //*** updated periodic velocity
 uTilde = uTildeU;
 uTilde.Append(uTildeV);
 uTilde.Append(uTildeW);

 if( rho_in <= rho_out ) // BUBBLE
 {
  cout << setw(70) << "BUBBLE SIMULATION" << endl;
  uvw = uTilde + 
        dt*invMrhoLumped*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		dt*invMLumped*( ((1.0/We)*(fint) ).MultVec(ip) );

 }
 else // DROPLET
 {
  cout << setw(70) << "DROPLET SIMULATION" << endl;
  uvw = uTilde + 
        invA*( ((1.0/(Fr*Fr))*( Mrho*gravity )).MultVec(ip) ) + 
		invA*( ((1.0/We)*(fint) ).MultVec(ip) );
 }
 
 //uvw = uTilde;

 /* 
  * Mass Transfer
  * */
 //clVector massTransfer = dt*invMcLumped*( heatFlux );
 //clVector massTransfer = ( invMcLumped*heatFlux );
 clVector massTransfer = invC*
                         (1.0/rho_inAdimen - 1.0/rho_outAdimen)
						 *heatFlux;

 ///*** setting b2 periodic, because DTilde*uvw already is.
 sumIndexPBCPress(VecXMin,VecXMax,b2);

 b2Tilde = (-1.0)*( b2 - (DTilde * uvw) ); 
 //b2Tilde = (-1.0)*( b2 - (DTilde * uvw) + (massTransfer) );

 // resolve sistema E pTilde = b2Tilde
 cout << " --------> solving pressure --------- " << endl;
 solverP->solve(1E-15,ETilde,pTilde,b2Tilde);
 cout << " ------------------------------------ " << endl;
 
 //*** copying pressure
 pbc.SetPurePressurePBC(pTilde,VecXMinGlob,VecXMaxGlob,nyPointsL,"RL");

 uvw = uvw - (invA * GTilde * pTilde);
 uTildeU = uvw.Copy(0,numNodes - 1);
 uTildeV = uvw.Copy(numNodes,2*numNodes - 1);
 uTildeW = uvw.Copy(2*numNodes,3*numNodes - 1);

 pbc.SetVelocityPBC(uTildeU,uTildeV,uTildeW,VecXMinGlob,VecXMaxGlob,nyPointsL,"RL");

 //*** updated periodic velocity
 uSol = uTildeU;
 vSol = uTildeV;
 wSol = uTildeW;

 pSol = pTilde;       // sem correcao na pressao
 //pSol = pSol + pTilde;  // com correcao na pressao

 // compute bubble's centroid velocity
 if( surfMesh->numInterfaces > 0 )
  setCentroidVelPos();
} // fecha metodo unCoupledPBC 


/* Sets the dimensionless pressure gradient for PBC application.
 * 
 * \param[in]: pressure drop, channel length, density of liquid,
 * velocity (Quantities of reference).
 * 
 * 
 * \remark Standard values for code verification against the Poiseuille 
 * flow and recovery of the flow rate are Re = 1.0; betaPressLiq = 12.0.
 *
 */
void Simulator3D::setBetaPressureLiquid()
{
     betaPressLiq = 32.0/Re;
	
}


void Simulator3D::inputVelocityPBC()
{
	Periodic3D pbc;

	clVector VecXMinAux, VecXMaxAux;

	VecXMinAux.Dim(nyPointsL);
	VecXMaxAux.Dim(nyPointsL);

	VecXMin->CopyTo(0,VecXMinAux);
	VecXMax->CopyTo(0,VecXMaxAux);

	pbc.SetVelocityPBC(uSol, vSol, wSol, VecXMinAux, VecXMaxAux, nyPointsL, direction);

	*uc = uSol;
	*vc = vSol;
	*wc = wSol;

} // fecha metodo

void Simulator3D::inputPurePressurePBC()
{
	Periodic3D pbc;

	clVector VecXMinAux, VecXMaxAux;

	VecXMinAux.Dim(nyPointsL);
	VecXMaxAux.Dim(nyPointsL);

	VecXMin->CopyTo(0,VecXMinAux);
	VecXMax->CopyTo(0,VecXMaxAux);

	pbc.SetPurePressurePBC(pSol, VecXMinAux, VecXMaxAux, nyPointsL, direction);

	*pc = uSol;

} // fecha metodo


void Simulator3D::setRHS_PBC()
{
	va = ( (1.0/dt) * Mrho + (1-alpha) * -(1.0/Re) * K ) * convUVW 
	 	 + M*betaFlowLiq;

} // fecha metodo

void Simulator3D::sumIndexPBCVel(clVector* _indexL, clVector* _indexR, clVector& _b)
{
    for (int i = 0; i < _indexL->Dim(); i++)
    {
       int ibL = _indexL->Get(i);
       int ibR = _indexR->Get(i);
 
       double uL = _b.Get(ibL);
       double uR = _b.Get(ibR);
       _b.Set(ibR, uL + uR);
       _b.Set(ibL,0);
 
       uL = _b.Get(ibL + numNodes);
       uR = _b.Get(ibR + numNodes);
       _b.Set(ibR + numNodes, uL + uR);
       _b.Set(ibL + numNodes,0);
 
       uL = _b.Get(ibL + 2*numNodes);
       uR = _b.Get(ibR + 2*numNodes);
       _b.Set(ibR + 2*numNodes, uL + uR);
       _b.Set(ibL + 2*numNodes,0);

     }
} // fecha metodo


void Simulator3D::sumIndexPBCVelVector(vector<int> _indexL, vector<int> _indexR, clVector& _b)
{
    for (int i = 0; i < _indexL.size(); i++)
    {
       int ibL = _indexL.at(i);
       int ibR = _indexR.at(i);
 
       double uL = _b.Get(ibL);
       double uR = _b.Get(ibR);
       _b.Set(ibR, uL + uR);
       _b.Set(ibL,0);
 
       uL = _b.Get(ibL + numNodes);
       uR = _b.Get(ibR + numNodes);
       _b.Set(ibR + numNodes, uL + uR);
       _b.Set(ibL + numNodes,0);
 
       uL = _b.Get(ibL + 2*numNodes);
       uR = _b.Get(ibR + 2*numNodes);
       _b.Set(ibR + 2*numNodes, uL + uR);
       _b.Set(ibL + 2*numNodes,0);

     }
} // fecha metodo




void Simulator3D::sumIndexPBCPress(clVector* _indexL, clVector* _indexR, clVector& _p)
{
   for (int i = 0; i < _indexL->Dim(); i++)
   {
      int ibL = _indexL->Get(i);
      int ibR = _indexR->Get(i);

      double pL = _p.Get(ibL);
      double pR = _p.Get(ibR);
      _p.Set(ibR,pL + pR);
      _p.Set(ibL,0);

    }
} // fecha metodo


void Simulator3D::sumIndexPBCPressVector(vector<int> _indexL, vector<int> _indexR, clVector& _p)
{
   for (int i = 0; i < _indexL.size(); i++)
   {
      int ibL = _indexL.at(i);
      int ibR = _indexR.at(i);

      double pL = _p.Get(ibL);
      double pR = _p.Get(ibR);
      _p.Set(ibR,pL + pR);
      _p.Set(ibL,0);

    }
} // fecha metodo


void Simulator3D::sumIndexPBCScalar(clVector* _indexL, clVector* _indexR, clVector& _s)
{
   for (int i = 0; i < _indexL->Dim(); i++)
   {
      int ibL = _indexL->Get(i);
      int ibR = _indexR->Get(i);

      double qL = _s.Get(ibL);
      double qR = _s.Get(ibR);
      _s.Set(ibR,qL + qR);
      _s.Set(ibL,0);

    }

} // fecha metodo


void Simulator3D::sumIndexPBCScalarVector(vector<int> _indexL, vector<int> _indexR, clVector& _s)
{
   for (int i = 0; i < _indexL.size(); i++)
   {
      int ibL = _indexL.at(i);
      int ibR = _indexR.at(i);

      double qL = _s.Get(ibL);
      double qR = _s.Get(ibR);
      _s.Set(ibR,qL + qR);
      _s.Set(ibL,0);

    }

} // fecha metodo



void Simulator3D::setCopyDirectionPBC(string _direction)
{
 	direction = _direction;

} // fecha metodo


/* \brief Initializes a Taylor vortex in the flow. */
void Simulator3D::initTaylorVortex()
{
 	init();
	#if NUMGLEU == 5
 	double numBCPoints = numVerts;
	#else
	double numBCPoints = numNodes;
	#endif

	double xM = 0.25*( 3.0*X->Max() + X->Min() );
	double yM = 0.5*( Y->Max() + Y->Min() );

	for ( int i = 0; i < numBCPoints; i++ )
	{
		double x = X->Get(i) - xM;
		double y = Y->Get(i) - yM;
	
		double r = sqrt( x*x + y*y );
		double theta = atan2(y,x);

		double r2 = ( Y->Max() - Y->Min() )/40;

		double c1 = 0.5;

		double vtheta = c1*r*( exp ( -r*r/r2 ) );

		double U = 0.0;
		double V = 0.0;
		double W = 0.0;

		uSol.Set(i, U - vtheta*sin(theta));
		vSol.Set(i, V + vtheta*cos(theta));
		wSol.Set(i, W);

		uSolOld.Set(i, U - vtheta*sin(theta));
		vSolOld.Set(i, V + vtheta*cos(theta));
		wSolOld.Set(i, W);
	}

}

void Simulator3D::initTanHJetProfile()
{
    init();
	#if NUMGLEU == 5
 		double numBCPoints = numVerts;
	#else
		double numBCPoints = numNodes;
	#endif

	for ( int i = 0; i < numBCPoints; i++ )
	{
		double y = Y->Get(i);
		double z = Z->Get(i);
	
		double r = sqrt( y*y + z*z );
		double thetaZero = 1.0;

		double U = 0.5 - 0.5*tanh( (0.25/thetaZero)*( r - 1.0/r ) );
		double V = 0.0;
		double W = 0.0;

		uSol.Set(i, U);
		vSol.Set(i, V);
		wSol.Set(i, W);

		uSolOld.Set(i, U);
		vSolOld.Set(i, V);
		wSolOld.Set(i, W);
	}


}

/* \brief Intializes a Taylor-Green  vortex in the flow. */ 
void Simulator3D::initTaylorGreenVortex()
{
 	init();
	#if NUMGLEU == 4
 	double numBCPoints = numVerts;
	#else 
 	double numBCPoints = numNodes;
	#endif

	for ( int i = 0; i < numBCPoints; i++ )
	{
	   double x = X->Get(i);
	   double y = Y->Get(i);

	   double u =   sin(PI_CONSTANT*x)*cos(PI_CONSTANT*y);
	   double v = - cos(PI_CONSTANT*x)*sin(PI_CONSTANT*y);
	   double w = 0.0;

	   double U = 0.0;
	   double V = 0.0;
	   double W = 0.0;
	  
	   uSol.Set(i, U + u);
	   vSol.Set(i, V + v);
	   wSol.Set(i, W + w);

	   uSolOld.Set(i, U + u);
	   vSolOld.Set(i, V + v);
	   wSolOld.Set(i, W + w);
	}
}

