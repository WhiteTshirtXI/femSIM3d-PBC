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

 const char *mesh = "../../db/mesh/3d/2bubble.vtk";
 //const char *mesh = "../../db/mesh/3d/2bubble2.vtk";
 const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";
 const char *simFolder  = "./sim/";

 int iter = 0;

 Model3D m1,mNew;
 m1.readVTK(mesh);
 m1.setMiniElement();
 m1.set2BubbleBC(); // malha do disco
 m1.setOFace();
 m1.setSurfaceConfig();

 Simulator3D s1(m1);

 s1.setRe(100);
 s1.setSc(2);
 s1.setWe(2);
 s1.setAlpha(1);
 s1.setBeta(-2.0);
 //s1.setSigma(0.0);
 s1.setCflBubble(5/2);
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

  m1.readVTK(mesh2);
  m1.setMiniElement();
  m1.readVTKCC(mesh2);

  Simulator3D s2(m1);
  s1 = s2; 

  s1.setRe(100);
  s1.setSc(2);
  s1.setWe(2);
  s1.setAlpha(1);
  s1.setBeta(-2.0);
  //s1.setSigma(0.0);
  s1.setCflBubble(5/2);
  //s1.init();
  s1.setSolverPressure(new PetscSolver(KSPGMRES,PCILU));
  s1.setSolverVelocity(new PetscSolver(KSPCG,PCICC));
  s1.setSolverConcentration(new PetscSolver(KSPCG,PCICC));

  s1.loadSolution(binFolder,"sim-last",0);
  iter = s1.loadIteration();
  //s1.loadSolution(binFolder,"UVWPC",22); // set para velocidade no simulador
  //iter = s1.loadIteration(vtkFolder,"sim",22);
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
   s1.stepALE2();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravityBoussinesq();
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveVTKTest(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"UVWPC",i*nReMesh+j+iter);
   save.saveDistance(simFolder,"dist",s1.getTime2());
  }
  mOld = m1;
  m1.mesh3DPoints();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation();
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


