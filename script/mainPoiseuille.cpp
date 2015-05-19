// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "Simulator3D.h"
#include "InOut.h"
#include "petscksp.h"
#include "PetscSolver.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 double Re = 10000;
 double Sc = 2000;
 double alpha = 1;
 double cfl = 1;

 /* \attention{ Certify that these environmental 
  * variables are correctly defined. Otherwise, 
  * SIGABORT runtime error will be launched. }
  */
 string meshDir = (string) getenv("MESH3D_DIR");
 string ppd = getenv("POST_PROCESSING3D_DIR");

 string meshFile = "poiseuille.msh";
 meshDir += "/singlePhase/internal/" + meshFile;
 const char *mesh = meshDir.c_str();

 Solver *solverP = new PetscSolver(KSPCG,PCICC);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 string bin = "/poiseuille/bin/";
 string vtk = "/poiseuille/vtk/";
 bin = ppd + bin;
 vtk = ppd + vtk;
 const char *binFolder = bin.c_str();
 const char *vtkFolder = vtk.c_str();

 Model3D m1;

  m1.readMSH(mesh);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D();
  m1.setMapping();

#if NUMGLEU == 4
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setSurfaceConfig();
 m1.setInitSurfaceArea();
 m1.setInitSurfaceVolume();

 // boundary conditions
 m1.setGenericBC();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc); // coeficiente do laplaciano de suavizacao
 s1.setAlpha(alpha);
 s1.setMu(1.0);
 s1.setRho(1.0);
 s1.setCp(1.0);
 s1.setCfl(cfl);
 s1.setDtEulerian();
 s1.init();
 s1.assemble();
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
  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
  s1.setCfl(cfl);
  s1.setDtEulerian();
 }

 InOut save(m1,s1); // objeto para gravacao de dados
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 int nIter = 1000;
 int nRe = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: "
	    << iter << endl;

    save.printSimulationReport();

	cout << "SL" << endl;
    s1.stepSL();
	cout << "ASSEMBLE" << endl;
	s1.assemble();
	cout << "MATMOUNT" << endl;
    s1.matMount();
	cout << "MATMOUNTC" << endl;
    s1.matMountC();
	cout << "SETUNCOUPLEDBC" << endl;
    s1.setUnCoupledBC();
	cout << "SETUNCOUPLEDCBC" << endl;
    s1.setUnCoupledCBC();
	cout << "SETRHS" << endl;
    s1.setRHS();
	cout << "SETCRHS" << endl;
    s1.setCRHS();
	cout << "UNCOUPLED" << endl;
    s1.unCoupled();
	cout << "UNCOUPLEDC" << endl;
    s1.unCoupledC();
	save.saveVTK(vtkFolder,"sim",iter);
	save.saveSol(binFolder,"sim",iter);

	s1.saveOldData();
	
	cout << "________________________________________ END of "
	     << iter << endl;

	iter++;
  }
 }
 PetscFinalize();
 return 0;
}


