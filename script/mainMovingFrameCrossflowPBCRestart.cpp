// =================================================================== //
// this is file mainRisingBubble.cpp, created at 10-Jun-2009           //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "TElement.h"
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

 int iter = 0;

 double velVCrossflow = 1.0; // i.e. V_crossflow = velVCrossflow x V_jet
 double velWCrossflow = velVCrossflow; 
 
 //const char* _frame = "fixed";
 const char* _frame = "moving";

 string _physGroup = "\"wallInflowVTransverse\"";
 double betaGrad = 0.0;
 
 Solver *solverP = new PetscSolver(KSPCG,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 string meshFile = "crossflow-holed-3d.msh";
 //string meshFile = "crossflow-holed-3d-Lp3.msh";

 const char* name = "wl";
 
 const char *binFolder = "null";
 const char *datFolder = "null";
 const char *mshFolder = "null";
 const char *vtkFolder = "null";

 if( strcmp( name,"wl") == 0 )
 {
 binFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/bin/";
 datFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/dat/";
 mshFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/msh/";
 vtkFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/vtk/";
 }
 else
 { 
 binFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/bin/";
 datFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/dat/";
 mshFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/msh/";
 vtkFolder  = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/vtk/";
 } 

 string meshDir = (string) getenv("MESH3D_DIR");

 if( strcmp( _frame,"moving") == 0 )
  meshDir += "/rising/movingFrame/" + meshFile;
 else
  meshDir += "/rising/" + meshFile;

 const char *mesh = meshDir.c_str();

 Model3D m1;


  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  string mshBase = "null"; 
 if( strcmp( name,"wl") == 0 )
  mshBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/msh/newMesh-";
 else
  mshBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/msh/newMesh-";

  // load surface mesh
  string aux = *(argv+1);
  string file = mshBase + *(argv+1) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();

  string vtkBase = "null";
 if( strcmp( name,"wl") == 0 )
  vtkBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/vtk/sim-";
 else
  vtkBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/vtk/sim-";
  
  // load 3D mesh
  file = vtkBase + *(argv+1) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.readVTKHeaviside(vtkFile);
  m1.setSurfaceConfig();
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setCrossflowVVelocity(velVCrossflow); 
  m1.setCrossflowWVelocity(velWCrossflow); 
  m1.setGenericBCPBCNew(_physGroup);
  m1.setGenericBC();

  Periodic3D pbc(m1);
  //pbc.MountPeriodicVectorsNew("noPrint");
  pbc.MountPeriodicVectors("noPrint");

  Simulator3D s1(pbc,m1);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

 const char *dirBase = "null";
 if( strcmp( name,"wl") == 0 )
  dirBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-longmire-Lp5/";
 else
  dirBase = "/work/gcpoliveira/post-processing/3d/crossflow-lambda-1.0-meister-Lp5/";
  
  iter = s1.loadSolution(dirBase,"sim",atoi(*(argv+1)));
  
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
 h1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveVTKPBC(vtkFolder,"initial",0,betaGrad);
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
  // moving
  uref = s1.getURef();
  xref = s1.getXRef();
  vref = s1.getVRef();
  yref = s1.getYRef();
  wref = s1.getWRef();
  zref = s1.getZRef();

  s1.setCentroidVelPos();
  xinit = s1.getCentroidPosXAverage(); // initial x centroid
  yinit = s1.getCentroidPosYAverage(); // initial y centroid
  zinit = s1.getCentroidPosZAverage(); // initial z centroid
 }

 int nIter = 50000;
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
	m1.setCrossflowWVelocity(velWCrossflow); // set of crossflow velocity for BC
    m1.setGenericBCPBCNew(_physGroup);
	m1.setGenericBC(uref,vref,wref); // crossflow condition
    //pbc.MountPeriodicVectorsNew("print");
    pbc.MountPeriodicVectors("print");
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
   s1.setUnCoupledBC();
   //s1.setGravity("-Z");
   //s1.setBetaFlowLiq("+X");
   s1.setRHS();
   s1.setCopyDirectionPBC("RL");
   s1.setInterfaceGeo();
   //s1.unCoupledPBCNew();
   s1.unCoupledPBC();

   if ( i%15 == 0 )
   {
   save.saveVTKPBC(vtkFolder,"sim",iter,betaGrad);
   save.saveVTKSurfacePBC(vtkFolder,"sim",iter,betaGrad);
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.bubbleWallDistance(datFolder,"dist",iter);
   save.saveBubbleShapeFactors(datFolder,"shapeFactors",iter);
   }
   
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
  m1.insert3dMeshPointsByDiffusion(6.5);
  m1.remove3dMeshPointsByDiffusion(1.5);
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface mesh operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength("curvature");
  //m1.insertPointsByCurvature("flat");
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance("flat");
  m1.contractEdgesByLength("curvature");
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
  

  if( strcmp( _frame,"moving") == 0 )
  {
	m1.setCrossflowVVelocity(velVCrossflow); 
	m1.setCrossflowWVelocity(velWCrossflow); 
    m1.setGenericBCPBCNew(_physGroup);
	m1.setGenericBC(uref,vref,wref);
    //pbc.MountPeriodicVectorsNew("noPrint");
    pbc.MountPeriodicVectors("noPrint");
  }
  else
  { 
	m1.setCrossflowVVelocity(velVCrossflow); 
	m1.setCrossflowWVelocity(velWCrossflow); 
    m1.setGenericBCPBCNew(_physGroup);
    m1.setGenericBC();
    //pbc.MountPeriodicVectorsNew("noPrint");
    pbc.MountPeriodicVectors("noPrint");
  }

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  if ( i%15 == 0 )
  {
  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
  }
 }
 
 PetscFinalize();
 return 0;
}


