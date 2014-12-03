/* \file mainMovingFrameCrossflowPBC.cpp
 * \author Gustavo Peixoto
 * \email gustavo.oliveira@uerj.br 
 * \date November 2014
 *
 *  \description{ 
 *  
 *  Script to simulate the problem of the drop jet in
 *  crossflow.
 *
 *  Boundary conditions
 *  ===================
 *  
 *  + kinds:
 *
 *  outflow: top
 *  periodic: sides
 *  crossflow: bottom
 *  slip: crossflow-top/bottom
 *  
 *  + physical groups in setGenericBC():
 *  'outflow', 'wallLeft', 'wallRight',
 *  'wallInflowVTransverse','wallNormalW'.
 *
 *  + calls setGenericBC(velU,velV,velW) 
 *
 *  }
 */ 
#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "InOut.h"
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 int iter = 1;
 double Re = 50.0;
 double We = 6.0;
 double Fr = 100.0;
 double alpha = 1;
 double mu_out = 0.048;
 double mu_in = 0.15*mu_out;
 double rho_out = 1136.0;
 double rho_in = 1.18*rho_out;
 double betaGrad = 0.0;
 double velVCrossflow = 2.0; // i.e. V_crossflow = velVCrossflow x V_jet

 //const char* _frame = "fixed";
 const char* _frame = "moving";
 
 string _physGroup = "\"wallInflowVTransverse\"";

 // fixed
 double c1 = 0.0;      // lagrangian
 double c2 = 1.0;      // smooth vel 
 double c3 = 3.0;     // smooth coord (fujiwara)
 double d1 = 1.0;      // surface tangent vel = (u-ut)
 double d2 = 0.5;      // surface smooth coord (fujiwara)

 // moving
 if( strcmp( _frame,"moving") == 0 )
 {
  c1 = 0.0;      // lagrangian
  c2 = 1.0;      // smooth vel: OBS - different result with c1=0.0
  c3 = 3.0;      // smooth coord (fujiwara)
  d1 = 0.0;      // surface tangent velocity u_n=u-u_t 
  d2 = 0.5;      // surface smooth cord (fujiwara)
 }


 Solver *solverP = new PCGSolver(); 
 //Solver *solverP = new PetscSolver(KSPCG,PCICC);
 //Solver *solverP = new PetscSolver(KSPCG,PCILU); 
 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU); 
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();
 
 // moving
 //string meshFile = "crossflow.msh";
 string meshFile = "crossflow-holed-3d.msh";
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/crossflow/bin/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/crossflow/vtk/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/crossflow/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/crossflow/msh/";

 string meshDir = (string) getenv("MESH3D_DIR");

 if( strcmp( _frame,"moving") == 0 )
  meshDir += "/rising/movingFrame/" + meshFile;
 else
  meshDir += "/rising/" + meshFile;

 const char *mesh = meshDir.c_str();

 Model3D m1;

 const char *mesh1 = mesh;

  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  m1.readMSH(mesh1);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setCrossflowVVelocity(velVCrossflow); 
  m1.setGenericBCPBCNew(_physGroup);
  m1.setGenericBC();

  Periodic3D pbc(m1);
  pbc.MountPeriodicVectorsNew("noPrint");

  Simulator3D s1(pbc,m1);

  s1.setRe(Re);
  s1.setWe(We);
  s1.setFr(Fr);
  s1.setC1(c1);
  s1.setC2(c2);
  s1.setC3(c3);
  s1.setD1(d1);
  s1.setD2(d2);
  s1.setAlpha(alpha);
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.initJetVelocity(1.0);
  s1.setBetaPressureLiquid(betaGrad); 
  s1.setDtALETwoPhase();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
 

 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initRisingBubble();
 h1.assemble();
 h1.setk(0.2);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 //h1.saveVTK(vtkFolder,"edge");
 h1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTKPBC(vtkFolder,"initial",0,betaGrad);
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 double uinst=0;
 double vinst=0;
 double winst=0;
 double uref=0;
 double vref=0;
 double wref=0;
 double xref=0;
 double yref=0;
 double zref=0;
 double xinit=0;
 double yinit=0;
 double zinit=0;
 double dx=0;
 double dy=0;
 double dz=0;

 if( strcmp( _frame,"moving") == 0 )
 {
  s1.setCentroidVelPos();
  xinit = s1.getCentroidPosXAverage(); // initial x centroid
  yinit = s1.getCentroidPosYAverage(); // initial y centroid
  zinit = s1.getCentroidPosZAverage(); // initial z centroid
 }

 int nIter = 1000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl;
   cout << resetColor();

   if( strcmp( _frame,"moving") == 0 )
   {
	// moving frame
	dx = s1.getCentroidPosXAverage() - xinit; // mov. frame x displacement 
	dy = s1.getCentroidPosYAverage() - yinit; // mov. frame y displacement
	dz = s1.getCentroidPosZAverage() - zinit; // mov. frame z displacement
	uinst = s1.getCentroidVelXAverage() + dx/s1.getDt(); // uc + correction
	vinst = s1.getCentroidVelYAverage() + dy/s1.getDt(); // vc + correction
	winst = s1.getCentroidVelZAverage() + dz/s1.getDt(); // wc + correction
	uref += uinst; // sums to recover inertial vel
	vref += vinst;
	wref += winst;
	xref += uref*s1.getDt();
	yref += vref*s1.getDt();
	zref += wref*s1.getDt();
	cout << "uref: " << uref << " xref: " << xref << endl;
	cout << "vref: " << vref << " yref: " << yref << endl;
	cout << "wref: " << wref << " zref: " << zref << endl;
	cout << "uinst: " << uinst << " vinst: " << vinst << " winst: " << winst << endl;
	cout << "dx: " << dx << endl;
	cout << "dy: " << dy << endl;
	cout << "dz: " << dz << endl;
	s1.setUSol(uinst); // subtraction for inertial frame: u - u_MFR
	s1.setVSol(vinst); // subtraction for inertial frame: v - v_MFR
	s1.setWSol(winst); // subtraction for inertial frame: w - w_MFR
	m1.setCrossflowVVelocity(velVCrossflow); // set of crossflow velocity for BC
    m1.setGenericBCPBCNew(_physGroup);
	m1.setGenericBC(uref,vref,wref); // crossflow condition
    pbc.MountPeriodicVectorsNew("print");
	s1.setURef(uref); // sets to recover the inertial physics; ease prints in InOut
	s1.setVRef(vref);
	s1.setWRef(wref);
	s1.setXRef(xref);
	s1.setYRef(yref);
	s1.setZRef(zref);
   }

   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepALEPBC();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC(); // <<
   s1.setGravity("-Z");
   //s1.setBetaFlowLiq("+X");
   s1.setRHS();
   s1.setCopyDirectionPBC("RL");
   s1.setInterfaceGeo();
   //s1.setInterfaceLevelSet();
   s1.unCoupledPBCNew();

   save.saveVTKPBC(vtkFolder,"sim",iter,betaGrad);
   save.saveVTKSurfacePBC(vtkFolder,"sim",iter,betaGrad);
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.bubbleWallDistance(datFolder,"dist",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl;
   cout << resetColor();

   s1.timeStep();

   iter++;
  }

  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initRisingBubble();
  h2.assemble();
  h2.setk(0.2);
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.saveVTK(vtkFolder,"edge",iter-1);
  //h2.saveChordalEdge(datFolder,"edge",iter-1);
  h2.setModel3DEdgeSize();
  

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  
  m1.setNormalAndKappa();
  m1.initMeshParameters();
  // 3D mesh operations
  m1.insert3dMeshPointsByDiffusion(6.5); //<<
  m1.remove3dMeshPointsByDiffusion(1.5); //<<
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints(); //<<

  // surface mesh operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength("curvature"); //<<
  //m1.insertPointsByCurvature("flat");
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance("flat");
  m1.contractEdgesByLength("curvature"); //<<
  //m1.removePointsByLength();
  m1.flipTriangleEdges();
  
  m1.removePointsByNeighbourCheck();
  //m1.checkAngleBetweenPlanes();

  /* **************************************** */
  
  m1.mesh3DPoints();
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setInterfaceBC();

  if( strcmp( _frame,"moving") == 0 )
  {
	m1.setCrossflowVVelocity(velVCrossflow); 
    m1.setGenericBCPBCNew(_physGroup);
	m1.setGenericBC(uinst,vinst,winst);
    pbc.MountPeriodicVectorsNew("noPrint");
  }
  else
  { 
	m1.setCrossflowVVelocity(velVCrossflow); 
    m1.setGenericBCPBCNew(_physGroup);
    m1.setGenericBC();
    pbc.MountPeriodicVectorsNew("noPrint");
  }

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
  
 }

 PetscFinalize();
 return 0;
}


