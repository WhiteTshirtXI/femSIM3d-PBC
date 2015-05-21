/* \file    mainTaylorVortexPBC.cpp
 * \author  Gustavo Peixoto de Oliveira 
 * \email   gustavo.oliveira@uerj.br
 * \date    Created on August 19th, 2014
 *
 * \Remark {Slow convergence without pressure point; faster with
 * one-point pressure setting.}
 *
 * Transport case:
 *
 * NaCl (salt) diffusing in water at 20 oC and salinity 35 g/kg
 * k = 1.611E-9 m2/s (diffusivity coeff.) at 25 oC;
 * rho = 1035.0 kg/m3
 * mu = 1.08E-3 Pa.s;
 * nu = 1.05E-6 m2/s 
 * Sc = nu/k;
 *
 */

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "Simulator3D.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "InOut.h"
#include "petscksp.h"
#include "PetscSolver.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 double Re = 600.0;  
 double Sc = 1.05E-6/1.611E-9; 
 double Fr = 1.0;
 double alpha = 1;
 double cfl = 0.5;
 double mu_l = 1.08E-3;
 double rho_l = 1035.0;

 /* Overlapping phys. group */
 string physGroup = "\"wallInflowU\""; // required for slip-kind condition

 /* \attention{ Certify that these environmental 
  * variables are correctly defined. Otherwise, 
  * SIGABORT runtime error will be launched. }
  */
 string meshDir = (string) getenv("MESH3D_DIR");
 string ppd = getenv("POST_PROCESSING3D_DIR");

 string meshFile = "taylorVortexPBC.msh";
 meshDir += "/singlePhase/pbc/" + meshFile;
 const char *mesh = meshDir.c_str();

 /* \remark PCGSolver family gives the best results for PBC.
  * Compare with PetscSolver to see error propagation over the
  * periodic boundaries. Scalar is unchangeed. */
 //Solver *solverP = new PCGSolver();
 Solver *solverP = new PetscSolver(KSPCG,PCICC);
 Solver *solverV = new PCGSolver(); 
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 string bin = "/taylorVortexPBC/bin/";
 string vtk = "/taylorVortexPBC/vtk/";
 string dat = "/taylorVortexPBC/dat/";
 string msh = "/taylorVortexPBC/msh/";
 bin = ppd + bin;
 vtk = ppd + vtk;
 dat = ppd + dat;
 msh = ppd + msh;
 const char *binFolder = bin.c_str();
 const char *vtkFolder = vtk.c_str();
 const char *datFolder = dat.c_str();
 const char *mshFolder = msh.c_str();

 Model3D m1;

  m1.readMSH(mesh);
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

 // mesh statistics info 
 m1.setInitSurfaceVolume();
 m1.setSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.tetMeshStats();

 // boundary conditions
 m1.setGenericBCPBCNew(physGroup); // set periodic vectors
 m1.setGenericBC(); 

 // Periodic call
 Periodic3D pbc(m1);
 pbc.MountPeriodicVectorsNew("print");

 Simulator3D s1(pbc,m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setAlpha(alpha);
 s1.setFr(Fr);
  s1.setMu(mu_l);
  s1.setRho(rho_l);
  s1.setCfl(cfl);
  s1.setDtEulerian();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  // initial condition vortex + scalar layer
  s1.initTaylorVortex(); 
  s1.initCTwoShearLayers(1.0,0.0); 

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
  iter = s1.loadSolution("./","sim",atoi(*(argv+2))); // change dir
  s1.setCfl(cfl);
  s1.setDtEulerian();
 }
 
 InOut save(m1,s1); 
 save.saveVTK(vtkFolder,"sim",0);
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 100;
 int nRe = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   save.printSimulationReport();

   s1.stepSLPBCFix(); 
   s1.assembleSlip(); // required for moving walls and scalar solving
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();  
   s1.setUnCoupledCBC();
   s1.setRHS();
   s1.setCRHS();
   s1.setCopyDirectionPBC("RL");
   s1.unCoupledPBCNew(); 
   s1.unCoupledCPBCNew();  

   save.saveVTK(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   
   // error
   s1.calcTaylorVortexError();
   save.saveTaylorVortexError(datFolder); 

   s1.saveOldData();
   s1.timeStep();

   cout << "________________________________________ END of " 
	    << iter << endl;

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}


