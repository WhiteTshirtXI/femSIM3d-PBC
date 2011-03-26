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
#include "InOut.h"

int main(int argc, char **argv)
{
 int iter = 0;
 real Re = 10000;
 real Sc = 2000;
 real Fr = 10;
 real alpha = 1;
 real beta = 0;
 real cfl = 3;
 Solver *solverP = new PCGSolver();
 Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 const char *dir  = "./";
 const char *mesh = "../../db/mesh/3d/step40-20-2.vtk";
 const char *txtFolder  = "./txt/";
 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *datFolder  = "./dat/";

 Model3D m1;
 m1.setMeshStep(40,20,2);
 m1.setAdimenStep();
 m1.setMiniElement();
 //m1.setQuadElement();
 m1.setStepBC();
 m1.setCStepBC();
 m1.setOFace();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setFr(Fr);
 //s1.setDt(dt);
 s1.setCfl(cfl);
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 s1.init();
 s1.assembleSlip();
 s1.matMount();
 s1.matMountC();
 s1.setUnCoupledBC(); 
 s1.setUnCoupledCBC(); 

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
  const char *sol = file.c_str();
  s1.loadSolution(binFolder,sol);
  iter = s1.loadIteration(vtkFolder,sol);
 }

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);
 save.printInfo(mesh);

 int nIter = 100;
 int nRe = 5;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << i*nRe+j+iter << endl;

   s1.stepSL();
   s1.setRHS();
   s1.setCRHS();
   s1.unCoupled();
   s1.unCoupledC();
   save.saveVTK(vtkFolder,"sim",i*nRe+j+iter);
   save.saveVTU(vtkFolder,"sim",i*nRe+j+iter);
   save.saveSol(binFolder,"UVWPC",i*nRe+j+iter);

   cout << "________________________________________ END of "
	    << i*nRe+j+iter << endl;
  }
 }

 return 0;
}


