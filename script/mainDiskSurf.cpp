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

 int iter = 1;
 real Re = 500;
 real Sc = 1;
 real Fr = 1;
 int beta = 1;
 real c1 = 0.6; // lagrangian
 real c2 = 0.1; // smooth
 real c3 = 0.06; // smooth
 real c4 = 0.0; // surface
 real mu_l = 1.0;
 real rho_l = 1.0;
 real cfl = 0.3;
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 string meshFile = "disk6-10-20.vtk";

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 //const char *simFolder  = "./sim/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/mesh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 m1.setMeshDisk(8,8,7);
 m1.setMiniElement();
 m1.setDiskFSBC();
 m1.setDiskCFSBC();
 m1.setPerturbSurf();
 //m1.setPerturbSurf2();
 //m1.setPerturbSurfSquare();
 m1.setOFace();
 m1.setSurface();
 m1.setInOutVert();
 //m1.setInOutElem();
 m1.setVertNeighbour();
 m1.setInOutElem();
 m1.setNeighbourSurfaceElem();
 m1.setNeighbourSurfacePoint();
 //m1.setMapEdgeTri();
 //m1.computeNormalAndKappa();

 Simulator3D s1(m1);

 s1.setRe(Re);
 s1.setSc(Sc);
 s1.setCfl(cfl);
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setFr(Fr);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setC4(c4);
 s1.setBeta(beta);
 s1.setSolverVelocity(solverV);
 s1.setSolverPressure(solverP);
 s1.setSolverConcentration(solverC);

 s1.init();
 s1.setDtDisk();

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);
 save.printInfo(meshFile.c_str());

 int nIter = 3000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << i*nReMesh+j+iter << endl << endl;
   cout << resetColor();

   /* dt variable */
   //s1.setDtDisk();

   //s1.stepLagrangian();
   s1.stepALEVel();
   //s1.stepALE();
   s1.assembleSlip();
   s1.matMount();
   s1.matMountC();
   s1.setUnCoupledBC();
   s1.setUnCoupledCBC();
   s1.setRHS();
   s1.setCRHS();
   s1.setGravity("Z");
   s1.unCoupled();
   s1.unCoupledC();

   InOut save(m1,s1); // cria objeto de gravacao
   save.saveMSH(mshFolder,"newMesh",i*nReMesh+j+iter);
   save.saveVTK(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveVTKHalf(vtkFolder,"simCutPlane",i*nReMesh+j+iter);
   save.saveVTKSurface(vtkFolder,"sim",i*nReMesh+j+iter);
   save.saveSol(binFolder,"sim",i*nReMesh+j+iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",i*nReMesh+j+iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << i*nReMesh+j+iter << endl << endl;;
   cout << resetColor();
  }
//--------------------------------------------------
//   Model3D mOld = m1; 
//   //m1.mesh2Dto3DOriginal();
//   m1.mesh3DPoints();
//   m1.setMiniElement();
//   m1.setDiskFSBC();
//   m1.setOFace();
//   m1.setSurface();
//   m1.setInOutVert();
//   //m1.setInOutElem();
//   m1.setVertNeighbour();
//   m1.setInOutElem();
//   m1.setNeighbourSurface();
//   //m1.setMapEdgeTri();
//   //m1.computeNormalAndKappa();
// 
//   Simulator3D s2(m1,s1);
//   s2.applyLinearInterpolation(mOld);
//   s1 = s2;
//   s1.setSolverPressure(solverP);
//   s1.setSolverVelocity(solverV);
//   s1.setSolverConcentration(solverC);
// 
//   InOut saveEnd(m1,s1); // cria objeto de gravacao
//   saveEnd.saveMSH(mshFolder,"newMesh",nReMesh+i*nReMesh+iter-1);
//   saveEnd.saveVTK(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
//   saveEnd.saveVTKSurface(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
//   saveEnd.saveVTKTest(vtkFolder,"simCutPlane",nReMesh+i*nReMesh+iter-1);
//   saveEnd.saveSol(binFolder,"sim",nReMesh+i*nReMesh+iter-1);
//   //saveEnd.saveVTU(vtkFolder,"sim",nReMesh+i*nReMesh+iter-1);
//   //saveEnd.saveSolTXT(binFolder,"sim",nReMesh+i*nReMesh+iter-1);
//   saveEnd.saveMeshInfo(datFolder);
//-------------------------------------------------- 
 }

 PetscFinalize();
 return 0;
}


