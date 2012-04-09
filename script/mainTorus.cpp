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
#include "TElement.h"
#include "GMRes.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 
 // static bubble test (Fabricio's thesis (2005))
 int iter = 1;
 real Re = 10;
 real Sc = 1;
 real We = 1;
 real Fr = 1.0;
 real c1 = 1.0; // lagrangian
 real c2 = 1.0; // smooth vel
 real c3 = 0.0; // smooth coord (fujiwara)
 real d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.0;  // surface smooth cord (fujiwara)
 real alpha = 1;
 real beta = 1;

 real sigma = 1.0;

 real mu_in = 1.0;
 real mu_out = 0.1;

 real rho_in = 1.0;
 real rho_out = 0.01;

 real cfl = 0.8;

 string meshFile = "torus.msh";

 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/torus/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 cout << endl;
 cout << "--------------> STARTING FROM 0" << endl;
 cout << endl;

 const char *mesh1 = mesh;
 m1.readMSH(mesh1);
 m1.setInterfaceBC();
 m1.setTriEdge();
 //m1.checkTriangleOrientation();
 m1.mesh2Dto3D();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setOFace();
 m1.setSurfaceConfig();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.setGenericBC();

 s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setWe(We);
 s1.setFr(Fr);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);
 s1.setAlpha(alpha);
 s1.setBeta(beta);
 s1.setSigma(sigma);
 //s1.setDtALETwoPhase(dt);
 s1.setMu(mu_in,mu_out);
 s1.setRho(rho_in,rho_out);
 s1.setCfl(cfl);
 s1.init();
 s1.setDtALETwoPhase();
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 3000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   //s1.stepLagrangian();
   //s1.stepALE();
   s1.stepALEVel();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setRHS();
   //s1.setGravity("Z");
   //s1.setInterface();
   s1.setInterfaceGeo();
   s1.unCoupled();

   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   save.chordalPressure(datFolder,"chordalPressure",iter);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}


