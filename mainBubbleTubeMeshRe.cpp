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

 const char *mesh = "../../db/mesh/3d/bubble-tube2D-2.vtk";
 const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 int iter = 0;

 Model3D m1,mNew,mOld,mOriginal;
 m1.readVTKSurface(mesh);
 m1.setSphere(1.5,1.5,3.5,0.5,1E-10); 
 m1.mesh2Dto3D();
 m1.setMiniElement();
 m1.setCubeBC2();
 m1.setOFace();
 m1.setSurfaceConfig();

 mOriginal = m1;

 Simulator3D s1(m1);
 s1.setRe(100);
 s1.setSc(2);
 s1.setWe(10);
 s1.setAlpha(1);
 s1.setBeta(-2.0);
 //s1.setSigma(0.0);
 s1.setCflBubble(30);
 s1.init();

 s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
 s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
 s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

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

  const char *mesh2 = "./vtk/sim-last-0.vtk";
  //const char *mesh2 = "./vtk/sim-6.vtk";

  m1.readVTK(mesh2);
  m1.setMiniElement();
  m1.readVTKCC(mesh2);
  m1.setCubeBC2();
  m1.setOFace();
  m1.setSurfaceConfig();

  //Simulator3D s2(m1,s1);
  Simulator3D s2(m1);

  s1 = s2; 
  s1.setRe(100);
  s1.setSc(2);
  s1.setWe(10);
  s1.setAlpha(1);
  s1.setBeta(-2.0);
  s1.setCflBubble(30);
  s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
  s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
  s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

  s1.loadSolution(binFolder,"sim-last",0);
  iter = s1.loadIteration();
  //s1.loadSolution(binFolder,"UVWPC",6);
  //iter = s1.loadIteration(vtkFolder,"sim",6);
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry",0);
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);

 int nIter = 1;
 int nReMesh = 5;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << "____________________________________ Iteration: " 
	<< i*nReMesh+j+iter << endl;
   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALE2();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravity();
   //s1.setGravityBoussinesq();
   s1.setInterfaceGeo();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"UVWPC",i*nReMesh+j+iter);
   save.oscillating("./","oscillating",i*nReMesh+j+iter);
   save.oscillatingD("./","oscillatingD",i*nReMesh+j+iter);
   save.oscillatingKappa("./","oscillatingKappa",i*nReMesh+j+iter);
  }
  mOld = m1; 
  m1.mesh2Dto3DOriginal();
  m1.setMiniElement();
  m1.setCubeBC2();
  m1.setOFace();
  m1.setSurfaceConfig();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2; 
  s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
  s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
  s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.saveVTK(vtkFolder,"sim-last",0);
  saveEnd.saveSol(binFolder,"sim-last",0);
  saveEnd.saveSimTime(iter+nReMesh*nIter);
 }

 PetscFinalize();
 return 0;
}


