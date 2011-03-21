// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 0;
 real Re = 5;
 real Sc = 2;
 real We = 2;
 real alpha = 1;
 real beta = -20;
 real cfl = 10;
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *mesh = "../../db/gmsh/3D/bubble-tube1.msh";

 Model3D m1,mOld;

 if( *(argv+1) == NULL )     
 {
  const char *mesh1 = mesh;
  //const char *mesh = "../../db/gmsh/3D/cube-cube2D.msh";
  //const char *mesh = "../../db/gmsh/3D/3D-bubble-cube1.msh";
  m1.readMSH(mesh1);
 }
 else if( strcmp( *(argv+1),"restart") == 0 || 
          strcmp( *(argv+1),"remesh") == 0  ) 
 {
  string aux = *(argv+2);
  string file = (string) "./msh/newMesh-" + *(argv+2) + (string) ".msh";
  const char *mesh2 = file.c_str();
  m1.readMSH(mesh2);
 }

 m1.setInterfaceBC();
 m1.mesh2Dto3D();
 m1.setMiniElement();
 m1.setOFace();
 m1.setSurfaceConfig();
 m1.setWallBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setWe(We);
 s1.setAlpha(alpha);
 s1.setBeta(beta);
 //s1.setSigma(sigma);
 s1.setCflBubble(cfl);
 s1.init();
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

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

  string aux = *(argv+2);
  string file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
  m1.setMiniElement();
  m1.readVTKCC(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallBC();

  Simulator3D s2(m1,s1);
  s1 = s2; 
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  file = (string) "sim-" + *(argv+2);
  const char *sol = file.c_str();
  s1.loadSolution(binFolder,sol);
  iter = s1.loadIteration(vtkFolder,sol);
 }
 else if( strcmp( *(argv+1),"remesh") == 0 )  
 {
  cout << endl;
  cout << "--------------> RE-MESHING (NO ITERATION)..." << endl;
  cout << endl;

  string aux = *(argv+2);
  string file = (string) "./vtk/sim-" + *(argv+2) + (string) ".vtk";
  const char *vtkFile = file.c_str();

  m1.readVTK(vtkFile);
  m1.setMiniElement();
  m1.readVTKCC(vtkFile);
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallBC();
  mOld = m1; 
  m1.mesh2Dto3DOriginal();
  m1.setMiniElement();
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallBC();

  file = (string) "sim-" + *(argv+2);
  const char *sol = file.c_str();
  s1.loadSolution(binFolder,sol);
  iter = s1.loadIteration(vtkFolder,sol);

  Simulator3D s2(m1,s1);

  s2.applyLinearInterpolation(mOld);
  s1 = s2;

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveMSH(mshFolder,"newMesh",atoi(*(argv+2)));
  saveEnd.saveVTKSurface(vtkFolder,"sim",atoi(*(argv+2)));
  saveEnd.saveSol(binFolder,"sim",atoi(*(argv+2)));
  return 0;
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry",0);
 save.saveMeshInfo("./","meshingInfo" );
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);

 int nIter = 1;
 int nReMesh = 1;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << i*nReMesh+j+iter << endl;

   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravity();
   s1.setGravityBoussinesq();
   s1.setInterfaceGeo();
   //s1.setInterfaceGeoTest();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveMSH(mshFolder,"newMesh",i*nReMesh+j+iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"sim",i*nReMesh+j+iter);
   save.oscillating("./","oscillating",i*nReMesh+j+iter);
   save.oscillatingD("./","oscillatingD",i*nReMesh+j+iter);
   save.oscillatingKappa("./","oscillatingKappa",i*nReMesh+j+iter);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",i*nReMesh+j+iter);

   cout << "________________________________________ END of " 
	    << i*nReMesh+j+iter << endl;
  }
  mOld = m1; 
  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMiniElement();
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
  saveEnd.saveVTK(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveMSH(mshFolder,"newMesh",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveSol(binFolder,"sim",nReMesh+i*nReMesh+iter-1);
  saveEnd.saveSimTime(nReMesh+i*nReMesh+iter-1);
  saveEnd.saveMeshInfo("./","meshingInfo" );
 }

 PetscFinalize();
 return 0;
}


