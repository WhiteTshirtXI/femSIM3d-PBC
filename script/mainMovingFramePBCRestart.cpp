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

 //const char* _frame = "fixed";
 const char* _frame = "moving";

 string physGroup = "\"wallInflowZeroU\"";
 //string physGroup = "\"wallNoSlip\"";
 double betaGrad = 1.0;
 
 Solver *solverP = new PetscSolver(KSPCG,PCICC);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 // rising-pbc
 //string meshFile = "airWaterSugarPBC-wallLeftRight.msh";
 /*
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/rising-pbc/bin/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/rising-pbc/msh/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/rising-pbc/dat/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/rising-pbc/vtk/";
 */

 // rising-beta
 /*
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/bin/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/msh/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/dat/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/vtk/";
 */

 // cell-2D
 string meshFile = "unit-cell-s-2D-3d.msh";
 // 1
 /* 
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/cell-1/bin/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/cell-1/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/cell-1/msh/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/cell-1/vtk/";
 */ 
 // 2
 /*
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/cell-2/bin/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/cell-2/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/cell-2/msh/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/cell-2/vtk/";
 */
 // 3
 
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/cell-3/bin/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/cell-3/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/cell-3/msh/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/cell-3/vtk/";
 

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

  //string mshBase = "/work/gcpoliveira/post-processing/3d/cell-1/msh/newMesh-";
  //string mshBase = "/work/gcpoliveira/post-processing/3d/cell-2/msh/newMesh-";
  string mshBase = "/work/gcpoliveira/post-processing/3d/cell-3/msh/newMesh-";

  //string mshBase = "/work/gcpoliveira/post-processing/3d/rising-pbc/msh/newMesh-";
  //string mshBase = "/work/gcpoliveira/post-processing/3d/rising-beta/msh/newMesh-";

  // load surface mesh
  string aux = *(argv+1);
  string file = mshBase + *(argv+1) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();

  //string vtkBase = "/work/gcpoliveira/post-processing/3d/cell-1/vtk/sim-";
  //string vtkBase = "/work/gcpoliveira/post-processing/3d/cell-2/vtk/sim-";
  string vtkBase = "/work/gcpoliveira/post-processing/3d/cell-3/vtk/sim-";
 
  //string vtkBase = "/work/gcpoliveira/post-processing/3d/rising-pbc/vtk/sim-";
  //string vtkBase = "/work/gcpoliveira/post-processing/3d/rising-beta/vtk/sim-";
  
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
  m1.setGenericBCPBCNew(physGroup);
  m1.setGenericBC();

  Periodic3D pbc(m1);
  pbc.MountPeriodicVectorsNew("print");

  Simulator3D s1(pbc,m1);
  s1.setBetaPressureLiquid(betaGrad);

  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  //const char *dirBase = "/work/gcpoliveira/post-processing/3d/cell-1/";
  //const char *dirBase = "/work/gcpoliveira/post-processing/3d/cell-2/";
  const char *dirBase = "/work/gcpoliveira/post-processing/3d/cell-3/";
  
  //const char *dirBase = "/work/gcpoliveira/post-processing/3d/rising-pbc/";
  //const char *dirBase = "/work/gcpoliveira/post-processing/3d/rising-beta/";

  iter = s1.loadSolution(dirBase,"sim",atoi(*(argv+1)));
  
 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initRisingBubble();
 //h1.initThreeBubbles();
 h1.assemble();
 h1.setk(0.2);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 double vinst=0;
 double vref=0;
 double xref=0;
 double xinit=0;
 double dx=0;
 if( strcmp( _frame,"moving") == 0 )
 {
  // moving
  vref = s1.getURef();
  xref = s1.getXRef();
  s1.setCentroidVelPos();
  xinit = s1.getCentroidPosXAverage();
 }

 int nIter = 30000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   // moving
   if( strcmp( _frame,"moving") == 0 )
   {
	// moving frame
	dx = s1.getCentroidPosXAverage() - xinit;
	vinst = s1.getCentroidVelXAverage() + dx/s1.getDt();
    vref += vinst;
	xref += vref*s1.getDt();
	cout << "vref: " << vref << " xref: " << xref << endl;
	cout << "dx: " << dx << endl;
	s1.setUSol(vinst);
    m1.setGenericBCPBCNew(physGroup);
	m1.setGenericBC(vref);
    pbc.MountPeriodicVectorsNew("print");
	s1.setURef(vref);
	s1.setXRef(xref);
   }

   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepALEPBC();
   //s1.stepALE();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setGravity("-X");
   s1.setBetaFlowLiq("+X");
   s1.setRHS();
   s1.setCopyDirectionPBC("RL");
   s1.setInterfaceGeo();
   //s1.setInterfaceLevelSet();
   s1.unCoupledPBCNew();
   //s1.unCoupledBetaPBC(); // rising-beta

   if ( i%5 == 0 )
   {
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);
   save.saveBubbleShapeFactors(datFolder,"shapeFactors",iter);
   }
   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   s1.timeStep();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  //h2.initRisingBubble();
  h2.initThreeBubbles();
  h2.assemble();
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();
  m1.initMeshParameters();

  // 3D operations
  m1.insert3dMeshPointsByDiffusion(6.0);
  m1.remove3dMeshPointsByDiffusion(0.5);
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface operations
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

  //m1.mesh2Dto3DOriginal();
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
    m1.setGenericBCPBCNew(physGroup);
    m1.setGenericBC(vref);
    pbc.MountPeriodicVectorsNew("print");
  }
  else
  {
    m1.setGenericBCPBCNew(physGroup);
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


