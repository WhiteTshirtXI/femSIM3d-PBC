/** 
 * \file    mainBetaFlowPBC.cpp
 * \author  Gustavo Peixoto de Oliveira 
 * \email   tavolesliv@gmail.com
 * \date    Created on August 19th, 2014
 *
 */

#include <cmath>
#include "Model3D.h"
#include "Simulator3D.h"
#include "CGSolver.h"
#include "PCGSolver.h"
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

 int iter = 1;
 double Re = 35.0;  
 /* NaCl (salt) diffusing in water at 20 oC and salinity 35 g/kg
  * k = 1.611E-9 m2/s (diffusivity coeff.) at 25 oC;
  * rho = 1035.0 kg/m3
  * mu = 1.08E-3 Pa.s;
  * nu = 1.05E-6 m2/s 
  * Sc = nu/k;
  */
 double Sc = 1.05E-6/1.611E-9; cout << "Sc = " << Sc << endl; 
 double Fr = 1.0;
 double alpha = 1.0;
 double cfl = 0.1; 
 double mu_l = 1.0;
 double rho_l = 2.0;

 string meshFile = "cuboid.msh";

 //string physGroup = "\"wallNoSlip\""; 
 string physGroup = "\"wallNormalW\""; 
 double betaGrad = 12.0/Re;

 string meshDir = (string) getenv("MESH3D_DIR");
 
 meshDir += "/cuboid/" + meshFile;
 const char *mesh = meshDir.c_str();

 Solver *solverP = new PetscSolver(KSPCG,PCILU);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PCGSolver(); 
 //Solver *solverV = new PetscSolver(KSPCG,PCILU);
 //Solver *solverV = new PCGSolver();
 //Solver *solverC = new PCGSolver();
 Solver *solverC = new PetscSolver(KSPCG,PCICC);
 
 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/bin/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/vtk/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/msh/";

 Model3D m1;

  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;

  m1.readMSH(mesh1);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D("QYYAzpaq1.414");
  m1.setMapping();
#if NUMGLEU == 5
  m1.setMiniElement();
#else
  m1.setQuadElement();
#endif
  m1.setNeighbour();
  //m1.setVertNeighbour();
  //m1.setInOutVert();
  m1.setGenericBCPBCNew(physGroup);
  m1.setGenericBC();

  Periodic3D pbc(m1);
  pbc.MountPeriodicVectorsNew("print");

  Simulator3D s1(pbc,m1);

  s1.setRe(Re);
  s1.setSc(Sc);
  s1.setFr(Fr);
  s1.setAlpha(alpha);
  s1.setMu(mu_l);
  s1.setRho(rho_l);
  s1.setCfl(cfl);
  s1.setDtEulerian(); //<<< for fixed mesh 
  s1.initTaylorVortex(); 
  s1.initCGaussian(0.8);
  s1.setBetaPressureLiquid(betaGrad);
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
 
 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"sim",0);
 save.saveVTK(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 50;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {
   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl;
   cout << resetColor();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepSLPBCFix(); 
   s1.assembleSlip(); 
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();  
   s1.setUnCoupledCBC();
   //s1.setGravity("+X");
   //s1.setBetaFlowLiq("+X"); 
   s1.setRHS();
   s1.setCRHS();
   s1.setCopyDirectionPBC("RL");
   s1.unCoupledPBCNew(); 
   s1.unCoupledCPBCNew();  

   save.saveVTKPBC(vtkFolder,"sim",iter,betaGrad);
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveSol(binFolder,"sim",iter);
   
   // error
   s1.calcTaylorVortexError();
   save.saveTaylorVortexError(datFolder); 

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl;
   cout << resetColor();

   s1.timeStep();

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}


