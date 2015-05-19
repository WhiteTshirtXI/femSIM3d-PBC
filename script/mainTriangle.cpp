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
#include "PetscSolver.h"
#include "petscksp.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;

 // R1234ze at saturation temperature (31.5 celcius)
 double D       = 100E-06; // m
 double G       = 695; // kg/m^2.s
 double mu_l    = 187.4E-06; // Pa.s
 double rho_l   = 1137.6; // kg/m^3
 double v       = G/rho_l; // m/s
 double sigma   = 7.7E-03; // N/m
 double g       = 9.81; // m/s^2
 double Re = rho_l*v*D/mu_l;
 double Fr = v/(sqrt(g*D));

 double cfl = 0.3;

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 string meshFile = "triangle.msh";

 //const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 //const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/singlePhase/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 //m1.setMeshStep(40,20,2);
 //m1.setAdimenStep();
 m1.readMSH(mesh);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.mesh2Dto3D("QYYApa0.0001q1.41q10");
 //m1.mesh2Dto3D();
 m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setNeighbour();
 m1.setVertNeighbour();
 m1.setInOutVert();
 m1.setGenericBC();
 //m1.setCStepBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setFr(Fr);
 s1.setCfl(cfl);
 s1.setDtEulerian();
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 s1.init();
 s1.assemble();
 s1.matMount();
 //s1.matMountC();
 s1.setUnCoupledBC(); 
 //s1.setUnCoupledCBC(); 

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
  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 int nIter = 1000;
 int nRe = 5;
 for( int i=0;i<nIter;i++ )
 {
  save.printSimulationReport();
  //save.printMeshReport();
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   s1.stepSL();
   s1.setRHS();
   //s1.setCRHS();
   s1.setGravity("-Z");
   s1.unCoupled();
   //s1.unCoupledC();

   save.saveVTK(vtkFolder,"sim",iter);
   //save.saveVTU(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.crossSectionSol("uSol","YZ",5.0/7.0,vtkFolder,"secU",iter );
   save.crossSectionSol("uSol","YZ",6.0/7.0,vtkFolder,"secU6-7",iter );

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
