// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "Simulator3D.h"
#include "Interface3D.h"
#include "InOut.h"
#include "Mumps_Petsc.h"
#include "PetscSolver.h"
#include "petscksp.h"

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 0;
 real Re = 1500;
 real Sc = 2;
 real We = 10;
 real alpha = 1;
 real beta = 0;
 real cfl = 50;
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 //const char *mesh = "../../db/gmsh/3D/cube-tube2D.msh";
 const char *mesh = "../../db/gmsh/3D/cube-cube2D.msh";
 const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
//--------------------------------------------------
//  const char *txtFolder  = "/scratch2/gustavo/cubeCube/txt/";
//  const char *binFolder  = "/scratch2/gustavo/cubeCube/bin/";
//  const char *vtkFolder  = "/scratch2/gustavo/cubeCube/vtk/";
//  const char *datFolder  = "/scratch2/gustavo/cubeCube/dat/";
//-------------------------------------------------- 

 Model3D m1,mOld;

 m1.readMSH(mesh);
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

  const char *mesh2 = "./vtk/sim-last.vtk";
  //const char *mesh2 = "/scratch2/gustavo/cubeCube/vtk/sim-last-0.vtk";
  //const char *mesh2 = "/scratch2/gustavo/cubeCube/vtk/sim-1740.vtk";
  //const char *mesh2 = "./vtk/sim-4970.vtk";

  m1.readVTK(mesh2);
  m1.setMiniElement();
  m1.readVTKCC(mesh2);
  m1.setWallBC();
  m1.setOFace();
  m1.setSurfaceConfig();

  Simulator3D s2(m1,s1);
  s1 = s2; 
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
  s1.loadSolution(binFolder,"sim-last");
  iter = s1.loadIteration(vtkFolder,"sim-last");
  //s1.loadSolution(binFolder,"UVWPC",4970); // set para velocidade no simulador
  //iter = s1.loadIteration(vtkFolder,"sim",4970);
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry",0);
 save.saveMeshInfo("./","meshingInfo" );
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);

 int nIter = 100;
 int nReMesh = 5;
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
   //s1.setGravityBoussinesq();
   s1.setInterfaceGeo();
   //s1.setInterfaceGeoTest();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   //save.saveVTU(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"UVWPC",i*nReMesh+j+iter);
   save.oscillating("./","oscillating",i*nReMesh+j+iter);
   save.oscillatingD("./","oscillatingD",i*nReMesh+j+iter);
   save.oscillatingKappa("./","oscillatingKappa",i*nReMesh+j+iter);

   cout << "________________________________________ END of " 
	    << i*nReMesh+j+iter << endl;
  }
  mOld = m1; 
  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMiniElement();
  m1.setWallBC();
  m1.setOFace();
  m1.setSurfaceConfig();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2; 
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim-remeshing",iter+nReMesh*nIter-1);
  saveEnd.saveSol(binFolder,"UVWPC-remeshing",iter+nReMesh*nIter-1);
  saveEnd.saveSimTime(iter+nReMesh*nIter-1);
  saveEnd.saveMeshInfo("./","meshingInfo" );
 }

 PetscFinalize();
 return 0;
}


