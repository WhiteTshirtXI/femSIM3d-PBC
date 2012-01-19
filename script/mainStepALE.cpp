// =================================================================== //
// this is file mainStep.cpp, created at 10-Jun-2007                   //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "Simulator3D.h"
#include "TElement.h"
#include "InOut.h"
#include "Laplace3D.h"
#include "PetscSolver.h"
#include "petscksp.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 vector< real > triEdge;
 triEdge.resize(2);
 triEdge[0] = 0.1; // none
 triEdge[1] = 0.3; // wall 

 int iter = 1;
 real Re = 10000;
 real Sc = 2000;
 real Fr = 10;
 //real alpha = 1;
 //real beta = 0;
 real c1 = 0.3; // lagrangian
 real c2 = 0.00; // smooth vel
 real c3 = 0.2; // smooth - fujiwara
 real c4 = 0.00; // smooth surface - fujiwara
 real cfl = 1;
 real mu_l = 1.0;
 real rho_l = 1.0;
 //Solver *solverP = new PCGSolver();
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU);
 //Solver *solverP = new PetscSolver(KSPLSQR,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 string meshFile = "stepSimple.msh";

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 //m1.setMeshStep(40,20,4);
 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge(triEdge);
 m1.checkTriangleOrientation();
 m1.mesh2Dto3D();
 //m1.setAdimenStep();
 //m1.setSingleElement();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setVertNeighbour();
 m1.setInOutVert();
 m1.setMapEdge();
 m1.setStepBC();
 m1.setCStepBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setFr(Fr);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setC4(c4);
 s1.setCfl(cfl);
 s1.setDtEulerian();
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 s1.init();

 if( (*(argv+1)) == NULL )
 {
  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;
 }
 else if( strcmp( *(argv+1),"restart") == 0 )
 {
  cout << endl;
  cout << "--------------> RE-STARTING..." << endl;
  cout << endl;

  string file = (string) "sim-" + *(argv+2);
  iter = s1.loadSolution("sim",atoi(*(argv+2)));
 }
 // Point's distribution
 Laplace3D d1(m1);
 d1.setk(0.7);
 d1.init();
 d1.assemble();
 d1.setBC();
 d1.matMountC();
 d1.setUnCoupledCBC(); 
 d1.setCRHS();
 d1.unCoupledC();
 //d1.saveVTK("./vtk/","edge");
 d1.setModel3DEdgeSize();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 1000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   //s1.stepSL();
   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.setDtEulerian();
   s1.movePoints();
   s1.assembleSlip();
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();
   s1.setUnCoupledCBC();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();

   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKHalf(vtkFolder,"simCutPlane",iter);
   save.saveVTU(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.printSimulationReport();

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
  d2.saveVTK("./vtk/","edge",iter-1);
  d2.setModel3DEdgeSize();

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume(0.005);
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.delete3DPoints();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setVertNeighbour();
  m1.setInOutVert();
  m1.setMapEdge();
  m1.setStepBC();
  m1.setCStepBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKHalf(vtkFolder,"simCutPlane",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  saveEnd.saveMeshInfo(datFolder);
  saveEnd.printMeshReport();
 }

 PetscFinalize();
 return 0;
}


