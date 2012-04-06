// =================================================================== //
// this is file Simulator3D.cpp, created at 23-Ago-2007                //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include "Simulator3D.h"

Simulator3D::Simulator3D(){}

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

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 allocateMemoryToAttrib();
}

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

 g_0 = _sRight.g_0;
 sigma_0 = _sRight.sigma_0;
 rho_0 = _sRight.rho_0;
 mu_0 = _sRight.mu_0;

 gAdimen = _sRight.gAdimen;
 sigmaAdimen = _sRight.sigmaAdimen;
 rho_inAdimen = _sRight.rho_inAdimen;
 rho_outAdimen = _sRight.rho_outAdimen;
 mu_inAdimen = _sRight.mu_inAdimen;
 mu_outAdimen = _sRight.mu_outAdimen;

 K = _sRight.K;
 Kc = _sRight.Kc;
 Mrho = _sRight.Mrho;
 M = _sRight.M;
 Mc = _sRight.Mc;
 G = _sRight.G;
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
 Fold    = _sRight.Fold;
 mu      = _sRight.mu;
 rho     = _sRight.rho;
 hSmooth = _sRight.hSmooth;
 
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
 muOld      = _sRight.muOld;
 rhoOld     = _sRight.rhoOld;
 hSmoothOld = _sRight.hSmoothOld;

 solverV = _sRight.solverV;
 solverP = _sRight.solverP;
 solverC = _sRight.solverC;
}

Simulator3D::Simulator3D( Model3D &_m, Simulator3D &_s)  
{
 // mesh information vectors
 getModel3DAttrib(_m);

 Re    = _s.getRe();
 Sc    = _s.getSc();
 Fr    = _s.getFr();
 We    = _s.getWe();
 sigma = _s.getSigma();
 alpha = _s.getAlpha();
 beta  = _s.getBeta();
 dt    = _s.getDt();
 dtLagrangian = _s.getDtLagrangian();
 dtSemiLagrangian = _s.getDtSemiLagrangian();
 dtGravity = _s.getDtGravity();
 dtSurfaceTension = _s.getDtSurfaceTension();
 time  = _s.getTime();
 cfl   = _s.getCfl();
 g     = _s.getGrav();
 mu_in  = _s.getMu_in();
 mu_out  = _s.getMu_out();
 mu_inAdimen  = _s.getMu_inAdimen();
 mu_outAdimen  = _s.getMu_outAdimen();
 rho_in = _s.getRho_in();
 rho_out = _s.getRho_out();
 rho_inAdimen = _s.getRho_inAdimen();
 rho_outAdimen = _s.getRho_outAdimen();
 iter  = _s.getIter();
 c1    = _s.getC1();
 c2    = _s.getC2();
 c3    = _s.getC3();
 d1    = _s.getD1();
 d2    = _s.getD2();

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 numVertsOld = _s.m->getNumVerts();
 numNodesOld = _s.m->getNumNodes();
 numElemsOld = _s.m->getNumElems();

 allocateMemoryToAttrib();

 // recuperando campo de velocidade e pressao da malha antiga
 uSolOld    = *_s.getUSol();
 vSolOld    = *_s.getVSol();
 wSolOld    = *_s.getWSol();
 pSolOld    = *_s.getPSol();
 cSolOld    = *_s.getCSol();
 uALEOld    = *_s.getUALE();
 vALEOld    = *_s.getVALE();
 wALEOld    = *_s.getWALE();
 kappaOld   = *_s.getKappa();
 fintOld    = *_s.getFint();
 gravityOld = *_s.getGravity();
 muOld      = *_s.getMu();
 rhoOld     = *_s.getRho();
 hSmoothOld = *_s.getHSmooth();
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
}

void Simulator3D::initHeatTransfer()
{
 init();

 for( int i=0;i<numVerts;i++ )
  if( Z->Get(i) > 0.4*Z->Max() )
   cSolOld.Set(i,1.0);
}

void Simulator3D::initFixedBubbleZ()
{
 init();

 for( int i=0;i<numNodes;i++ )
 {
  real aux = 1.0;
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
  real aux = X->Get(i);
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
  real aux = X->Get(i);
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
 real aux;
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

#if NUMGLEU == 5
 FEMMiniElement3D miniElem(*X,*Y,*Z);
#else
 FEMQuadElement3D miniElem(*X,*Y,*Z);
#endif

 FEMLinElement3D linElem(*X,*Y,*Z);

 //setMu( mu_in );
 //setRho( rho_in );
 setMu( mu_in,mu_out );
 setRho( rho_in,rho_out );

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  real muValue=0;
  real rhoValue=0;
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

 Kc.CopyFrom(         0,          0,     KcMat );
 Mc.CopyFrom(         0,          0,     McMat );
 
} // fecha metodo ASSEMBLE

void Simulator3D::assembleC()
{
 int i,j,ii,jj;
 int v[NUMGLEC];
 real aux;
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  real dif = 1.0;

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
 real aux;
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

  real muValue = mu_inAdimen;
  real rhoValue = rho_inAdimen;

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
 real aux;
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
  real c = 0;
  for( int n=0;n<NUMGLEU;n++ )
  {
   v[n] = (int) IEN->Get(mele,n);
   c += cSolOld.Get(v[n]);
  }
  c = c/NUMGLE;

  real eme = 0.81315;
  real muC = exp(eme*c);
  real dif = 1.0/muC;
  real rhoValue = 1.0;

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
 real aux;
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
  real muValue = 0;
  real rhoValue = 0;
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

void Simulator3D::assembleNuZ()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 real aux;
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

 setMuZ();      // carregando o arquivo de perfil nuZ
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
  real muValue = 0;
  real rhoValue = 0;
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
 real aux;
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
  real c = 0;
  for( int n=0;n<NUMGLE;n++ )
   c += cSolOld.Get(v[n]);
  c = c/NUMGLE;

  real eme = 0.81315;
  real muC = exp(eme*c);
  real dif = 1.0/muC;

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
//  for ( int i=0;i<numVerts;i++ )
//  {
//   int node = i;
//   //if( heaviside->Get(node) != 0.5 )
//   {
//    velW.Set(node,velW.Get(node) - centroidVelZ);
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

void Simulator3D::stepNoConvection()
{
 convUVW.CopyFrom(0,uSolOld);
 convUVW.CopyFrom(numNodes,vSolOld);
 convUVW.CopyFrom(2*numNodes,wSolOld);
 convC = cSolOld;
} // fecha metodo stepSL 

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

} // fecha metodo stepLagrangian

// metodo para movimentacao dos pontos da malha na direcao Z atraves da 
// velocidade do escoamento (wSol), caracterizando a movimentacao
// puramente lagrangiana em Z. Para os pontos nas direcoes X e Y eh
// utilizado o metodo explicito semi lagrangiano
void Simulator3D::stepLagrangianZ()
{
 m->moveZPoints(wSol,dt);
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

void Simulator3D::stepALE()
{
 setInterfaceVelocity();

 // smoothing - coordenadas
 MeshSmooth e1(*m,1.0); // criando objeto MeshSmooth
 e1.stepSmoothFujiwara();

#if NUMGLEU == 5
  uSmooth = setTetCentroid(*IEN,*e1.getUSmooth());
  vSmooth = setTetCentroid(*IEN,*e1.getVSmooth());
  wSmooth = setTetCentroid(*IEN,*e1.getWSmooth());
#else
  uSmooth = setTetQuad(*IEN,*e1.getUSmooth());
  vSmooth = setTetQuad(*IEN,*e1.getVSmooth());
  wSmooth = setTetQuad(*IEN,*e1.getWSmooth());
#endif

 uALE = c1*uSolOld+c3*uSmoothCoord;
 vALE = c1*vSolOld+c3*vSmoothCoord;
 wALE = c1*wSolOld+c3*wSmoothCoord;

 // impoe velocidade (componente normal) do fluido na interface
 setInterfaceVelocity();

 // impoe velocidade ALE = 0 no contorno
 setALEVelBC();
 //setAnnularALEVelBC();

 // calcula velocidade do fluido atraves do metodo semi-lagrangeano
 stepSL();
} // fecha metodo stepALE

void Simulator3D::stepALEVel()
{
 // vertice velocity (uVert,vVert)
 clVector uVert(numVerts);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 uSol.CopyTo(0,uVert);
 vSol.CopyTo(0,vVert);
 wSol.CopyTo(0,wVert);

 setInterfaceVelocity();

//--------------------------------------------------
//  // smoothing - velocidade
//  MeshSmooth e1(*m,dt); // criando objeto MeshSmooth
//  e1.stepSmooth(uALE,vALE,wALE);
//  uSmooth = *e1.getUSmooth();
//  vSmooth = *e1.getVSmooth();
//  wSmooth = *e1.getWSmooth();
//-------------------------------------------------- 

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

 // smoothing coords
 MeshSmooth e2(*m,1.0); // criando objeto MeshSmooth
 e2.stepSmoothFujiwara();
 uSmoothCoord = *e2.getUSmooth();
 vSmoothCoord = *e2.getVSmooth();
 wSmoothCoord = *e2.getWSmooth();

//--------------------------------------------------
//  // smoothing coords
//  MeshSmooth e4(*m,1.0); // criando objeto MeshSmooth
//  e4.stepSmoothFujiwaraByHeight();
//  clVector uSmoothHeight = *e4.getUSmooth();
//  clVector vSmoothHeight = *e4.getVSmooth();
//  clVector wSmoothHeight = *e4.getWSmooth();
//-------------------------------------------------- 

 // compute ALE
 uALE = c1*uVert+c2*uSmooth+c3*uSmoothCoord;
 vALE = c1*vVert+c2*vSmooth+c3*vSmoothCoord;
 wALE = c1*wVert+c2*wSmooth+c3*wSmoothCoord;
 
 // impoe velocidade do fluido na interface
 setInterfaceVelocity();

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
 setALEVelBC();
 //setAnnularALEVelBC();

 // calcula velocidade do fluido atraves do metodo semi-lagrangeano
 stepSL();
} // fecha metodo stepALEVel

void Simulator3D::movePoints()
{
 // movimentando os vertices pontos da malha com velocidade ALE
 m->moveXPoints(uALE,dt);
 m->moveYPoints(vALE,dt);
 m->moveZPoints(wALE,dt);
 m->centroidPositionCorrection();

 // correcao do volume da bolha
 m->applyBubbleVolumeCorrection();
}

void Simulator3D::setInterfaceVelocity()
{
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
  real xNormalUnit = surfMesh->xNormal.Get(surfaceNode);
  real yNormalUnit = surfMesh->yNormal.Get(surfaceNode);
  real zNormalUnit = surfMesh->zNormal.Get(surfaceNode);

  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  real prod = uSolOld.Get(surfaceNode)*xNormalUnit+ 
              vSolOld.Get(surfaceNode)*yNormalUnit + 
			  wSolOld.Get(surfaceNode)*zNormalUnit;
  real uSolNormal = xNormalUnit*prod;
  real vSolNormal = yNormalUnit*prod;
  real wSolNormal = zNormalUnit*prod;

  real uSolTangent = uSolOld.Get(surfaceNode) - uSolNormal;
  real vSolTangent = vSolOld.Get(surfaceNode) - vSolNormal;
  real wSolTangent = wSolOld.Get(surfaceNode) - wSolNormal;

  // tratamento da superficie
  // produto escalar --> projecao do vetor normalUnit no segmento de reta
  // | Unit.RetaUnit | . RetaUnit
  // resultado = vetor normal a reta situado na superficie
  real prod2 = uSmoothSurface.Get(surfaceNode)*xNormalUnit + 
               vSmoothSurface.Get(surfaceNode)*yNormalUnit + 
			   wSmoothSurface.Get(surfaceNode)*zNormalUnit;
  real uSmoothNormal = xNormalUnit*prod2;
  real vSmoothNormal = yNormalUnit*prod2;
  real wSmoothNormal = zNormalUnit*prod2;

  real uSmoothTangent = uSmoothSurface.Get(surfaceNode) - uSmoothNormal;
  real vSmoothTangent = vSmoothSurface.Get(surfaceNode) - vSmoothNormal;
  real wSmoothTangent = wSmoothSurface.Get(surfaceNode) - wSmoothNormal;

  real uALESurface =   uSolOld.Get(surfaceNode) 
                     - d1*uSolTangent 
					 + d2*uSmoothTangent;
  real vALESurface =   vSolOld.Get(surfaceNode) 
                     - d1*vSolTangent 
					 + d2*vSmoothTangent;
  real wALESurface =   wSolOld.Get(surfaceNode) 
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
  real aux = zeroLevel.Get(i)*interfaceDistance->Get(i);
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
  real sumMat = mat.SumLine(i);
  invA.Set( i,1.0/sumMat );

  real sumM = M.SumLine(i);
  MLumped.Set( i,sumM );
  invMLumped.Set( i,1.0/sumM );

  real sumMrho = Mrho.SumLine(i);
  MrhoLumped.Set( i,sumMrho );
  invMrhoLumped.Set( i,1.0/sumMrho );
 }
}

void Simulator3D::matMountC()
{
 // set para Matriz Lumped da concentracao
 for( int i=0;i<numVerts;i++ )
 {
  real sumMc = Mc.SumLine(i);
  McLumped.Set( i,sumMc );
  invMcLumped.Set( i,1.0/sumMc );
 }

 //matc = ((1.0/dt) * Mc) + (alpha * (1.0/(Sc*Re)) * Kc);
 matc = ((1.0/dt) * McLumped) + (alpha * (1.0/(Sc*Re)) * Kc);

 for( int i=0;i<numVerts;i++ )
 {
  real sumMatc = matc.SumLine(i);
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

 uSolOld = uSol;
 vSolOld = vSol;
 wSolOld = wSol;
 pSolOld = pSol;
 
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

 b2Tilde = (-1.0)*( b2 - (DTilde * uvw) ); 

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

 time = time + dt;
 iter++;
} // fecha metodo unCoupled 

void Simulator3D::unCoupledC()
{
 clVector vcIp(numVerts);
 clVector b1cTilde;

 vcIp = vcc.MultVec(ipc); // operacao vetor * vetor (elemento a elemento)

 b1cTilde = b1c + vcIp;

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

void Simulator3D::saveOldData()
{
 uSolOld = uSol;
 vSolOld = vSol;
 wSolOld = wSol;
 pSolOld = pSol;
 cSolOld = cSol;
 uALEOld    = uALE;
 vALEOld    = vALE;
 wALEOld    = wALE;
 fintOld    = fint;
 gravityOld = gravity;
 kappaOld   = kappa;
 muOld      = mu;
 rhoOld     = rho;
 hSmoothOld = hSmooth;
 centroidVelXOld = centroidVelX;
 centroidVelYOld = centroidVelY;
 centroidVelZOld = centroidVelZ;
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

 vector<real> minLength = m->getMinLength();
 real minEdge = *min_element(minLength.begin()+1,minLength.end());

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
 real minEdge = m->getMinEdge();

 //cout << idMinVolume << " " << minEdge << endl;

 int v1 = IEN->Get(idMinVolume,0);
 int v2 = IEN->Get(idMinVolume,1);
 int v3 = IEN->Get(idMinVolume,2);
 int v4 = IEN->Get(idMinVolume,3);

 /* X-component */
 real p1x = X->Get(v1);
 real p2x = X->Get(v2);
 real p3x = X->Get(v3);
 real p4x = X->Get(v4);

 real xCentroid = ( p1x+p2x+p3x+p4x )*0.25;

 real d1x = dist(p1x,xCentroid);
 real d2x = dist(p2x,xCentroid);
 real d3x = dist(p3x,xCentroid);
 real d4x = dist(p4x,xCentroid);
 
 real minXdist = min(d1x,d2x);
 minXdist = min(minXdist,d3x);
 minXdist = min(minXdist,d4x);

 real xVelMax = max( 1.0,uALE.Abs().Max() );
 real minDtx = minXdist/xVelMax;

 /* Y-component */
 real p1y = Y->Get(v1);
 real p2y = Y->Get(v2);
 real p3y = Y->Get(v3);
 real p4y = Y->Get(v4);

 real yCentroid = ( p1y+p2y+p3y+p4y )*0.25;

 real d1y = dist(p1y,yCentroid);
 real d2y = dist(p2y,yCentroid);
 real d3y = dist(p3y,yCentroid);
 real d4y = dist(p4y,yCentroid);

 real minYdist = min(d1y,d2y);
 minYdist = min(minYdist,d3y);
 minYdist = min(minYdist,d4y);

 real yVelMax = max( 1.0,vALE.Abs().Max() );
 real minDty = minYdist/yVelMax;

 /* Z-component */
 real p1z = Z->Get(v1);
 real p2z = Z->Get(v2);
 real p3z = Z->Get(v3);
 real p4z = Z->Get(v4);

 real zCentroid = ( p1z+p2z+p3z+p4z )*0.25;

 real d1z = dist(p1z,zCentroid);
 real d2z = dist(p2z,zCentroid);
 real d3z = dist(p3z,zCentroid);
 real d4z = dist(p4z,zCentroid);

 real minZdist = min(d1z,d2z);
 minZdist = min(minZdist,d3z);
 minZdist = min(minZdist,d4z);

 real zVelMax = max( 1.0,wALE.Abs().Max() );
 real minDtz = minZdist/zVelMax;

 real minDt1 = min(minDtx,minDty);
 minDt1 = min(minDt1,minDtz);

 real velMax = max(xVelMax,yVelMax);
 velMax = max(velMax,zVelMax);
 real minDt2 = 0.5*minEdge/velMax;

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
 real minEdge = m->getMinEdge();

 real length = minEdge*0.86602;

 real velMax = max( 1.0,uALE.Abs().Max() );
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
  real p1x=X->Get(v1);
  real p1y=Y->Get(v1);
  real p1z=Z->Get(v1);

  // v2
  int v2 = mapEdge->Get(edge,5);
  real p2x=X->Get(v2);
  real p2y=Y->Get(v2);
  real p2z=Z->Get(v2);

  real length = distance(p1x,p1y,p1z,p2x,p2y,p2z);

  // bubble.py - 146 iterations
  real xVel = fabs(uALE.Get(v1)) - fabs(uALE.Get(v2));
  real yVel = fabs(vALE.Get(v1)) - fabs(vALE.Get(v2));
  real zVel = fabs(wALE.Get(v1)) - fabs(wALE.Get(v2));

  real vel = vectorLength(xVel,yVel,zVel);
  //real vel = distance(fabs(xVel),fabs(yVel),fabs(zVel),
  //                    fabs(x),fabs(y),fabs(z));

  real a = 0.2; // security parameter

  real minDt = a*length/vel;

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
 setDtLagrangianNorberto();
 dtSemiLagrangian=dtLagrangian;
}

void Simulator3D::setDtGravity()
{
 real minEdge = m->getMinEdge();
 real velMax = max( 1.0,gravity.Abs().Max() );

 dtGravity = sqrt(velMax/minEdge);
}

/*
 * Set Dt of the current simulation.
 * Explicity terms:
 *  - Semi-Lagrangian;
 *  */
void Simulator3D::setDtEulerian()
{
 // setting required dt
 setDtSemiLagrangian();
 setDtGravity();
 dtSurfaceTension=0.0;

 real dtEulerian = min(1.0,getDtSemiLagrangian());
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

 real dtALETwoPhase = min(getDtLagrangian(),getDtSurfaceTension());
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

 real dtALETwoPhase = min(getDtLagrangian(),getDtGravity());

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

void Simulator3D::setMuZ()
{
 // -- Leitura do perfil de nu variavel em Z para os nos da malha -- //

 real aux;
 real dist1,dist2;
 clMatrix muFile(1002,2); // vetor super dimensionado!

 const char* _filename = "../../db/baseState/nuZ/nuZ.dat";
 ifstream file( _filename,ios::in );

 if( !file )
 {
  cerr << "Esta faltando o arquivo de perfis!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !file.eof() )
 {
  for( int i=0;i<1001;i++ )
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
  for( j=0;j<1001;j++ )
  {
   dist1 = fabs( Z->Get(i) - muFile(j,0) );
   dist2 = fabs( Z->Get(i) - muFile(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = muFile(j,1);
  mu.Set(i,aux);
 }
}

void Simulator3D::setMuC()
{
 // -- Leitura do perfil de nu variavel em Z para os nos da malha -- //
 for( int mele=0;mele<numElems;mele++ )
 {
  int v1 = (int) IEN->Get(mele,0);
  int v2 = (int) IEN->Get(mele,1);
  int v3 = (int) IEN->Get(mele,2);
  int v4 = (int) IEN->Get(mele,3);
  real c = ( cSol.Get(v1)+
	         cSol.Get(v2)+
	         cSol.Get(v3)+
	         cSol.Get(v4) )/4.0;

  real eme = 0.81315;
  real muC = exp(eme*c);

  // updating mu
  mu.Set(v1,muC);
  mu.Set(v2,muC);
  mu.Set(v3,muC);
  mu.Set(v4,muC);
 }
}

void Simulator3D::setSigma(real _sigma)
{
 sigma = _sigma;
 sigma_0 = sigma;
 sigmaAdimen = sigma/sigma_0;
}

void Simulator3D::setSolverVelocity(Solver *s){solverV = s;}
void Simulator3D::setSolverPressure(Solver *s){solverP = s;}
void Simulator3D::setSolverConcentration(Solver *s){solverC = s;}
void Simulator3D::setRe(real _Re){Re = _Re;}
real Simulator3D::getRe(){return Re;}
void Simulator3D::setSc(real _Sc){Sc = _Sc;}
real Simulator3D::getSc(){return Sc;}
void Simulator3D::setFr(real _Fr){Fr = _Fr;}
real Simulator3D::getFr(){return Fr;}
void Simulator3D::setWe(real _We){We = _We;}
real Simulator3D::getWe(){return We;}
void Simulator3D::setAlpha(real _alpha){alpha = _alpha;}
real Simulator3D::getAlpha(){return alpha;}
void Simulator3D::setBeta(real _beta){beta = _beta;}
real Simulator3D::getBeta(){return beta;}
real Simulator3D::getSigma(){return sigma;}
void Simulator3D::setDt(real _dt){dt = _dt;}
void Simulator3D::setTime(real _time){time = _time;}
real Simulator3D::getDt(){return dt;}
real Simulator3D::getDtLagrangian(){return dtLagrangian;}
real Simulator3D::getDtSemiLagrangian(){return dtSemiLagrangian;}
real Simulator3D::getDtGravity(){return dtGravity;}
real Simulator3D::getDtSurfaceTension(){return dtSurfaceTension;}
void Simulator3D::setIter(real _iter){iter = _iter;}
int  Simulator3D::getIter(){return iter;}
real Simulator3D::getCfl(){return cfl;}
void Simulator3D::setC1(real _c1){c1 = _c1;}
void Simulator3D::setC2(real _c2){c2 = _c2;}
void Simulator3D::setC3(real _c3){c3 = _c3;}
void Simulator3D::setD1(real _d1){d1 = _d1;}
void Simulator3D::setD2(real _d2){d2 = _d2;}
real Simulator3D::getC1(){return c1;}
real Simulator3D::getC2(){return c2;}
real Simulator3D::getC3(){return c3;}
real Simulator3D::getD1(){return d1;}
real Simulator3D::getD2(){return d2;}

void Simulator3D::setMu(real _mu_in)
{ 
 mu_in = _mu_in;
 mu_0 = _mu_in;
 mu_inAdimen = mu_in/mu_0;

 mu.SetAll(mu_inAdimen); 
}

void Simulator3D::setRho(real _rho_in)
{ 
 rho_in = _rho_in;
 rho_0 = rho_in; 
 rho_inAdimen = rho_in/rho_0; 

 rho.SetAll(rho_inAdimen); 
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
void Simulator3D::setMu(real _mu_in,real _mu_out)
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

 //real rMax = 1.0;
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
  //if (X->Get(*it)*X->Get(*it)+Y->Get(*it)*Y->Get(*it) > rMax*rMax - 0.001) 
   mu.Set(*it,0.0); // zero viscosity at wall


//--------------------------------------------------
//  mu.Dim(numNodes);
//  mu.SetAll(mu_inAdimen);
//  list<int> *inElem;
//  inElem = m->getInElem();
//  for (list<int>::iterator it=inElem->begin(); it!=inElem->end(); ++it)
//   mu.Set(*it,mu_outAdimen);
//-------------------------------------------------- 
}

void Simulator3D::setRho(real _rho_in,real _rho_out)
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

 //real rMax = 1.0;
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
  //if (X->Get(*it)*X->Get(*it)+Y->Get(*it)*Y->Get(*it) > rMax*rMax - 0.001) 
   rho.Set(*it,0.0); // zero density at wall

//--------------------------------------------------
//  rho.Dim(numNodes);
//  rho.SetAll(rho_inAdimen);
//  list<int> *inElem;
//  inElem = m->getInElem();
//  for (list<int>::iterator it=inElem->begin(); it!=inElem->end(); ++it)
//   rho.Set(*it,rho_outAdimen);
//-------------------------------------------------- 
}

void Simulator3D::setMuSmooth(real _mu_in,real _mu_out)
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

void Simulator3D::setRhoSmooth(real _rho_in,real _rho_out)
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
 real triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));
 for( int i=0;i<numVerts;i++ )
 {
  real len = 1.3*triEdgeMin;
  real d = interfaceDistance->Get(i);
  real aux = zeroLevel.Get(i)*d;

  if( aux < -len )
   hSmooth.Set( i,0.0 );
  else if( aux > len )
   hSmooth.Set( i,1.0 );
  else
  {
   //real func = 0.5;
   real func = 0.5*( 1.0+(aux/len)+(1.0/3.1415)*sin(3.1415*(aux/len)) );
   hSmooth.Set( i,func );
  }
 }
}

//--------------------------------------------------
// real Simulator3D::computeReynolds()
// {
//  real D = 1.0;
//  real L = D;
//  real U = sqrt(g*D);
//  Re = rho_in*L*U/mu_in;
// }
// 
// real Simulator3D::computeFroud()
// {
//  real D = 1.0;
//  real L = D;
//  real U = sqrt(g*D);
//  Fr = U/sqrt(g*L);
// }
// 
// real Simulator3D::computeWebber()
// {
//  real D = 1.0;
//  real L = D;
//  real U = sqrt(g*D);
//  We = rho_in*L*U*U/sigma;
// }
// 
// real Simulator3D::computeEotvos()
// {
//  real D = 1.0;
//  Eo = rho_in*g*D*D/sigma;
//  We = Eo;
// }
// 
// real Simulator3D::computeGalileo()
// {
//  real D = 1.0;
//  real L = D;
//  real U = sqrt(g*D);
//  N = rho_in*U*L/mu_in;
//  Re = N;
// }
// 
// real Simulator3D::computeMorton()
// {
//  real D = 1.0;
//  real L = D;
//  real U = sqrt(g*D);
//  Mo = (rho_in-rho_out)*mu_in*mu_in*mu_in*mu_in*g/(rho_in*rho_in*sigma*sigma*sigma); 
// }
//-------------------------------------------------- 

real Simulator3D::getMu_in(){return mu_in;}
real Simulator3D::getMu_out(){return mu_out;}
real Simulator3D::getRho_in(){return rho_in;}
real Simulator3D::getRho_out(){return rho_out;}
real Simulator3D::getMu_inAdimen(){return mu_inAdimen;}
real Simulator3D::getMu_outAdimen(){return mu_outAdimen;}
real Simulator3D::getRho_inAdimen(){return rho_inAdimen;}
real Simulator3D::getRho_outAdimen(){return rho_outAdimen;}
real Simulator3D::getTime(){return time;}
clVector* Simulator3D::getUSol(){return &uSol;} 
clVector* Simulator3D::getUSolOld(){return &uSolOld;} 
void Simulator3D::setUSol(clVector &_uSol){uSol = _uSol;}
clVector* Simulator3D::getVSol(){return &vSol;} 
clVector* Simulator3D::getVSolOld(){return &vSolOld;} 
void Simulator3D::setVSol(clVector &_vSol){vSol = _vSol;}
clVector* Simulator3D::getWSol(){return &wSol;}
clVector* Simulator3D::getWSolOld(){return &wSolOld;}
void Simulator3D::setWSol(clVector &_wSol){wSol = _wSol;}
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
real Simulator3D::getGrav(){return g;}
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
clVector* Simulator3D::getHSmooth(){return &hSmooth;}
void Simulator3D::updateIEN(){IEN = m->getIEN();}
void Simulator3D::setCfl(real _cfl){cfl = _cfl;}
vector<real> Simulator3D::getCentroidVelX(){return centroidVelX;}
vector<real> Simulator3D::getCentroidVelY(){return centroidVelY;}
vector<real> Simulator3D::getCentroidVelZ(){return centroidVelZ;}
void Simulator3D::setCentroidVelX(vector<real> _centroidVelX)
{centroidVelX = _centroidVelX;}
void Simulator3D::setCentroidVelY(vector<real> _centroidVelY)
{centroidVelY = _centroidVelY;}
void Simulator3D::setCentroidVelZ(vector<real> _centroidVelZ)
{centroidVelZ = _centroidVelZ;}

vector<real> Simulator3D::getCentroidPosX(){return centroidPosX;}
vector<real> Simulator3D::getCentroidPosY(){return centroidPosY;}
vector<real> Simulator3D::getCentroidPosZ(){return centroidPosZ;}
void Simulator3D::setCentroidPosX(vector<real> _centroidPosX)
{centroidPosX = _centroidPosX;}
void Simulator3D::setCentroidPosY(vector<real> _centroidPosY)
{centroidPosY = _centroidPosY;}
void Simulator3D::setCentroidPosZ(vector<real> _centroidPosZ)
{centroidPosZ = _centroidPosZ;}


// set do centroide para o elemento mini apos a interpolacao linear
// (applyLinearInterpolation)
clVector Simulator3D::setCentroid(clVector &_vector)
{
 real aux;
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
 m = _sRight.m;
 numVerts = _sRight.numVerts;
 numNodes = _sRight.numNodes;
 numElems = _sRight.numElems;
 X = _sRight.X;
 Y = _sRight.Y;
 Z = _sRight.Z;
 uc = _sRight.uc;
 vc = _sRight.vc;
 wc = _sRight.wc;
 pc = _sRight.pc;
 cc = _sRight.cc;
 heaviside = _sRight.heaviside;
 idbcu = _sRight.idbcu;
 idbcv = _sRight.idbcv;
 idbcw = _sRight.idbcw;
 idbcp = _sRight.idbcp;
 idbcc = _sRight.idbcc;
 outflow = _sRight.outflow;
 surface = _sRight.surface;
 IEN = _sRight.IEN;
 interfaceDistance = _sRight.interfaceDistance;
 elemIdRegion = _sRight.elemIdRegion;
 triEdge = _sRight.triEdge;

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

 g_0 = _sRight.g_0;
 sigma_0 = _sRight.sigma_0;
 rho_0 = _sRight.rho_0;
 mu_0 = _sRight.mu_0;

 gAdimen = _sRight.gAdimen;
 sigmaAdimen = _sRight.sigmaAdimen;
 rho_inAdimen = _sRight.rho_inAdimen;
 rho_outAdimen = _sRight.rho_outAdimen;
 mu_inAdimen = _sRight.mu_inAdimen;
 mu_outAdimen = _sRight.mu_outAdimen;

 K = _sRight.K;
 Kc = _sRight.Kc;
 Mrho = _sRight.Mrho;
 M = _sRight.M;
 Mc = _sRight.Mc;
 G = _sRight.G;
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
 Fold    = _sRight.Fold;
 mu      = _sRight.mu;
 rho     = _sRight.rho;
 hSmooth = _sRight.hSmooth;
 
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
 muOld      = _sRight.muOld;
 rhoOld     = _sRight.rhoOld;
 hSmoothOld = _sRight.hSmoothOld;
 centroidVelXOld = _sRight.centroidVelXOld;
 centroidVelYOld = _sRight.centroidVelYOld;
 centroidVelZOld = _sRight.centroidVelZOld;

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
 g     = 9.81;
 mu_in  = 1.0;
 mu_out  = 1.0;
 rho_in = 1.0;
 rho_out = 1.0;

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 allocateMemoryToAttrib();
}

void Simulator3D::operator()(Model3D &_m,Simulator3D &_s) 
{
 // mesh information vectors
 getModel3DAttrib(_m);

 Re    = _s.getRe();
 Sc    = _s.getSc();
 Fr    = _s.getFr();
 We    = _s.getWe();
 sigma = _s.getSigma();
 alpha = _s.getAlpha();
 beta  = _s.getBeta();
 dt    = _s.getDt();
 dtLagrangian = _s.getDtLagrangian();
 dtSemiLagrangian = _s.getDtSemiLagrangian();
 dtGravity = _s.getDtGravity();
 dtSurfaceTension = _s.getDtSurfaceTension();
 time  = _s.getTime();
 cfl   = _s.getCfl();
 iter  = _s.getIter();
 c1    = _s.getC1();
 c2    = _s.getC2();
 c3    = _s.getC3();
 d1    = _s.getD1();
 d2    = _s.getD2();
 g     = _s.getGrav();
 mu_in  = _s.getMu_in();
 mu_out  = _s.getMu_out();
 rho_in = _s.getRho_in();
 rho_out = _s.getRho_out();
 mu_inAdimen  = _s.getMu_inAdimen();
 mu_outAdimen  = _s.getMu_outAdimen();
 rho_inAdimen = _s.getRho_inAdimen();
 rho_outAdimen = _s.getRho_outAdimen();

 setSolverVelocity( new PCGSolver() );
 setSolverPressure( new PCGSolver() );
 setSolverConcentration( new PCGSolver() );

 allocateMemoryToAttrib();

 centroidVelX = _s.getCentroidVelX();
 centroidVelY = _s.getCentroidVelY();
 centroidVelZ = _s.getCentroidVelZ();

 numVertsOld = _s.m->getNumVerts();
 numNodesOld = _s.m->getNumNodes();
 numElemsOld = _s.m->getNumElems();

 // recuperando campo de velocidade e pressao da malha antiga
 uSolOld    = *_s.getUSol();
 vSolOld    = *_s.getVSol();
 wSolOld    = *_s.getWSol();
 pSolOld    = *_s.getPSol();
 cSolOld    = *_s.getCSol();
 uALEOld    = *_s.getUALE();
 vALEOld    = *_s.getVALE();
 wALEOld    = *_s.getWALE();
 muOld      = *_s.getMu();
 rhoOld     = *_s.getRho();
 hSmoothOld = *_s.getHSmooth();
 kappaOld   = *_s.getKappa();
 fintOld    = *_s.getFint();
 gravityOld = *_s.getGravity();
 centroidVelXOld = _s.getCentroidVelX();
 centroidVelYOld = _s.getCentroidVelY();
 centroidVelZOld = _s.getCentroidVelZ();
}

int Simulator3D::loadSolution( const char* _filename,int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // read numVerts and numNodes and all the simulation parameters
 // iteration
 string file = "./vtk/" + (string) _filename + "-" + str + ".vtk";
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
 string fileUVWPC = "./bin/" + (string) _filename + "-" + str + ".bin";
 const char* filenameUVPC = fileUVWPC.c_str();

 clVector aux2(6*numNodesOld+2*numVertsOld); // vetor tambem carrega a concentracao

 ifstream UVPC_file( filenameUVPC,ios::in | ios::binary ); 

 if( !UVPC_file)
 {
  cerr << "Solution file is missing for reading!" << endl;
  exit(1);
 }

 UVPC_file.read( (char*) aux2.GetVec(),aux2.Dim()*sizeof(real) );

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
 hSmooth.Dim( numVerts );

 // the edgeSize vector is not part of the Simulator3D, but because we
 // are not saving the matrix interpLin, we should interpolate edgeSize
 // and send it back to Model3D. 
 // It is better to find another solution!!!
 clVector edgeSize( numVerts );

 // only 1 component because the 2 others are exactly the same
 clVector xKappaOld(_mOld.getNumVerts());
 for( int i=0;i<_mOld.getNumVerts();i++ )
 {
  real aux = kappaOld.Get(i);
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
 clVector hSmoothOldVert(numVertsOld);
 
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
 hSmoothOld.CopyTo(0,hSmoothOldVert);
 clVector edgeSizeOld = *_mOld.getEdgeSize();

 // interpolation process between the old mesh with numVertsOld to the new
 // numVert.
 pSol = interpLin*(pSolOld);
 cSol = interpLin*(cSolOld);
 mu = interpLin*(muOld);
 rho = interpLin*(rhoOld);
 hSmooth = interpLin*(hSmoothOld);
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
  real aux = xKappa.Get(i);
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

 // updating old data vectors with the new mesh.
 saveOldData();
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
void Simulator3D::setALEVelBC()
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

void Simulator3D::setAnnularALEVelBC()
{
 real rMax = 1.0;
 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); 
						  ++it)
 {
  if( (Z->Get(*it) == Z->Max() || Z->Get(*it) == Z->Min() ) &&
      (X->Get(*it)*X->Get(*it)+Y->Get(*it)*Y->Get(*it) < rMax*rMax - 0.001) )
   wALE.Set(*it,0.0);
  else
  {
   uALE.Set(*it,0.0);
   vALE.Set(*it,0.0);
   wALE.Set(*it,0.0);
  }
 }
}

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
 muOld.Dim( numVerts );
 rhoOld.Dim( numVerts );
 hSmoothOld.Dim( numVerts );

 // interface vectors (two-phase)
 fint.Dim ( 3*numNodes );
 gravity.Dim( 3*numNodes );
 mu.Dim( numVerts );
 rho.Dim( numVerts );
 hSmooth.Dim( numVerts );
}

void Simulator3D::setCentroidVelPos()
{
 clVector _uVel = uSolOld;
 clVector _vVel = vSolOld;
 clVector _wVel = wSolOld;

 int v = elemIdRegion->Max()+1;
 vector<real> velX(v,0);
 vector<real> velY(v,0);
 vector<real> velZ(v,0);
 vector<real> posX(v,0);
 vector<real> posY(v,0);
 vector<real> posZ(v,0);
 vector<real> volume(v,0);
 vector<real> sumVolume(v,0);
 vector<real> sumXVelVolume(v,0);
 vector<real> sumYVelVolume(v,0);
 vector<real> sumZVelVolume(v,0);
 vector<real> sumXPosVolume(v,0);
 vector<real> sumYPosVolume(v,0);
 vector<real> sumZPosVolume(v,0);

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
 vector<real> surfaceVolume = m->getSurfaceVolume();
 for( int nb=0;nb<=v;nb++ )
 {
  centroidVelX.push_back(sumXVelVolume[nb]/surfaceVolume[nb]);
  centroidVelY.push_back(sumYVelVolume[nb]/surfaceVolume[nb]);
  centroidVelZ.push_back(sumZVelVolume[nb]/surfaceVolume[nb]);

  centroidPosX.push_back(sumXPosVolume[nb]/surfaceVolume[nb]);
  centroidPosY.push_back(sumYPosVolume[nb]/surfaceVolume[nb]);
  centroidPosZ.push_back(sumZPosVolume[nb]/surfaceVolume[nb]);
 }
}

real Simulator3D::getCentroidVelXAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelX[nb];
 return sum/v;
}

real Simulator3D::getCentroidVelYAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelY[nb];
 return sum/v;
}

real Simulator3D::getCentroidVelZAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidVelZ[nb];
 return sum/v;
}

real Simulator3D::getCentroidPosXAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosX[nb];
 return sum/v;
}

real Simulator3D::getCentroidPosYAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosY[nb];
 return sum/v;
}

real Simulator3D::getCentroidPosZAverage()
{
 real sum=0;
 int v = elemIdRegion->Max();
 for( int nb=1;nb<=v;nb++ )
  sum+=centroidPosZ[nb];
 return sum/v;
}
