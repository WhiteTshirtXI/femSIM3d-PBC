// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
//#include "GMRes.h"
#include "Simulator3D.h"
#include "InOut.h"
#include "PetscSolver.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

 int iter = 0;
 real Re = 100;
 real Sc = 1;
 real Fr = 2;
 int beta = 1;
 real cfl = 10;
 real mu_l = 1.0;
 real rho_l = 1.0;
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 const char *mesh = "../../db/mesh/3d/disk6-10-20.vtk";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *simFolder  = "./sim/";

 Model3D m1;
 m1.setMeshDisk(5,5,20);
 m1.setMiniElement();
 m1.setDiskFSBC();
 m1.setPerturbSurf();
 //m1.setPerturbSurf2();
 //m1.setPerturbSurfSquare();
 m1.setOFace();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setCflDisk(cfl);
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setFr(Fr);
 s1.setBeta(beta);
 s1.setSolverVelocity(solverV);
 s1.setSolverPressure(solverP);

 s1.init();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);

 for( int i=0;i<1000;i++ )
 {
  cout << "____________________________________ Iteration: " << i << endl;
  //s1.stepLagrangian();
  s1.stepALE();
  s1.matMount();
  s1.matMountC();
  s1.setUnCoupledBC();
  s1.setUnCoupledCBC();
  s1.setRHS();
  s1.setCRHS();
  s1.setGravity();
  s1.unCoupled();
  s1.unCoupledC();
  save.saveVTK(vtkFolder,"sim",i);
  save.saveSol(binFolder,"UVWPC",i);

  cout << "__________________________________________ End: " 
       << i << endl;
 }

 PetscFinalize();
 return 0;
}


