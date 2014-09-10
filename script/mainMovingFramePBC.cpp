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

 /* 
  * bogdan's thesis 2010 - Film thickness
  * Air-Ethanol
  * Circular channels
  * Tube radius: 0.3 - 1.3mm
  *
  * Reference: Han and Shikazono
  *
  * */
 /*
 double Re = 100;
 double We = 10;
 double We = 0.1162;
 double Fr = 1;
 double mu_in = 0.01;
 double mu_out = 1.00;
 double rho_in = 0.001;
 double rho_out = 1.0;
 */

 // Rabello's thesis: sugar-syrup 1
 double Re = 33.0413; 
 double Sc = 1.0;
 double Fr = 1.0;
 double We = 115.662;
 double mu_in = 1.78E-5;
 double mu_out = 0.5396;
 double rho_in = 1.225;
 double rho_out = 1350.0;

 int iter = 0;
 double alpha = 1.0;
 double cfl = 0.8;

 double c1 = 0.0;  // lagrangian
 double c2 = 1.0;  // smooth vel
 double c3 = 10.0;  // smooth coord (fujiwara)
 double d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 double d2 = 0.1;  // surface smooth cord (fujiwara)

 //const char* _frame = "fixed";
 const char* _frame = "moving";

 //string meshFile = "circular.msh";
 string meshFile = "rising-x-moving.msh";
 
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "/home/gcpoliveira/post-processing/vtk/3d/rising-pbc-moving/bin/";
 const char *mshFolder  = "/home/gcpoliveira/post-processing/vtk/3d/rising-pbc-moving/msh/";
 const char *datFolder  = "/home/gcpoliveira/post-processing/vtk/3d/rising-pbc-moving/dat/";
 const char *vtkFolder  = "/home/gcpoliveira/post-processing/vtk/3d/rising-pbc-moving/";
 
 string meshDir = (string) getenv("MESH3D_DIR");

 if( strcmp( _frame,"moving") == 0 )
  meshDir += "/rising/movingFrame/" + meshFile;
 else
  meshDir += "/rising/" + meshFile;
 
 const char *mesh = meshDir.c_str();

 Model3D m1;

  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;

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
  m1.setGenericBCPBC();

  Periodic3D pbc(m1);
  pbc.MountPeriodicVectorsNew("print");

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
  s1.setCfl(cfl);
  s1.init();
  s1.setBetaPressureLiquid();
  s1.setDtALETwoPhase();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initThreeBubbles();
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
 if( strcmp( _frame,"moving") == 0 )
 {
  // moving
  vref = s1.getURef();
  s1.setCentroidVelPos();
 }

 int nIter = 3000;
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
    vinst = s1.getCentroidVelXAverage();
    vref += vinst;
    cout << vref << " " << vinst << endl;
    s1.setUSol(vinst);
    m1.setGenericBCPBC(vref);
	pbc.MountPeriodicVectorsNew("noprint");
    s1.setURef(vref);
   }

   //s1.stepLagrangian();
   s1.stepALEPBC();
   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.movePoints();
   s1.assemble();
   //s1.assembleBetaFlow(); //<<<
   s1.matMount();
   s1.setUnCoupledPBC();
   s1.setGravity("-X");
   s1.setBetaFlowLiq("+X");
   s1.setRHS_PBC();
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.setCopyDirectionPBC("RL");
   s1.unCoupledPBCNew();

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);
   save.saveBubbleShapeFactors(datFolder,"shapeFactors",iter);

   s1.saveOldData();

   s1.timeStep();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
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
  m1.insert3dMeshPointsByDiffusion(3.5);
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
   m1.setGenericBCPBC(vref);
   pbc.MountPeriodicVectorsNew("noPrint");
  }
  else
  {
   m1.setGenericBCPBC();
   pbc.MountPeriodicVectorsNew("print");
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


