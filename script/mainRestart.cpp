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
#include "Laplace3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 int iter = 1;

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";

 Model3D m1;
 Simulator3D s1;

 // load surface mesh
 string aux = *(argv+1);
 string file = (string) "./msh/newMesh-" + *(argv+1) + (string) ".msh";
 const char *mesh2 = file.c_str();
 m1.readMSH(mesh2);
 m1.setInterfaceBC();

 // load 3D mesh
 file = (string) "./vtk/sim-" + *(argv+1) + (string) ".vtk";
 const char *vtkFile = file.c_str();

 m1.readVTK(vtkFile);
 m1.setMiniElement();
 m1.readVTKHeaviside(vtkFile);
 m1.setOFace();
 m1.setSurfaceConfig();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.setInitSurfaceRadius();
 m1.setWallBC();

 s1(m1);

 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 iter = s1.loadSolution("sim",atoi(*(argv+1)));

 // Point's distribution
 Laplace3D d1(m1);
 d1.setk(0.1);
 d1.init();
 d1.assemble();
 d1.setBC();
 d1.matMountC();
 d1.setUnCoupledCBC(); 
 d1.setCRHS();
 d1.unCoupledC();
 //d1.saveVTK(vtkFolder,"edge");
 d1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 //save.saveInfo(datFolder,"info",mesh);
 //save.printInfo(meshFile.c_str());

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

   s1.stepALEVel();
   s1.setDtALETwoPhase();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   s1.setGravity("Z");
   s1.setInterfaceGeo();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKQuarter(vtkFolder,"simCutPlane",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.printSimulationReport();
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
  Laplace3D d2(m1,d1);
  d2.assemble();
  d2.setBC();
  d2.matMountC();
  d2.setUnCoupledCBC(); 
  d2.setCRHS();
  d2.unCoupledC();
  d2.saveVTK(vtkFolder,"edge",iter-1);
  d2.saveChordalEdge(datFolder,"edge",iter-1);
  d2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  //m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.delete3DPoints();

  // surface operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength();
  //m1.insertPointsByCurvature();
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance();
  m1.contractEdgeByLength();
  //m1.removePointsByLength();
  m1.flipTriangleEdge();

  m1.removePointByNeighbourCheck();
  m1.checkAngleBetweenPlanes();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveMSH(mshFolder,"newMesh",iter-1);
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKQuarter(vtkFolder,"simCutPlane",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  //saveEnd.saveVTU(vtkFolder,"sim",iter-1);
  //saveEnd.saveSolTXT(binFolder,"sim",iter-1);
  saveEnd.saveMeshInfo(datFolder);
  saveEnd.printMeshReport();
 }


 PetscFinalize();
 return 0;
}


