// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
//                                                                     //
// Obs.: Annular flow setup:                                           //
//                                                                     //
// Simulator3D: setAnnularALEVelBC()                                   //
//              mu and rho not zero for Z.Max() and Z.Min()            //
//              turn off applyBubbleVolumeCorrection()                 //
//                                                                     //
// Model3D: insertPointsByLength()                                     //
//           - no insertion at the boundaries Z.Max() and Z.Min()      //
//          setPolyhedron()                                            //
//           - uncomment wSwap2 == node2 until break;                  //
//          uncomment setNormalAndKappa2D()                            //
//                                                                     //
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
 
 // set each bubble length
 vector< real > triEdge;
 triEdge.resize(2);
 triEdge[0] = 0.1; // wall
 triEdge[1] = 0.1; // bubble 1 

 int iter = 1;
 real Re = 100;
 real We = 5;
 real c1 = 0.00; // lagrangian
 real c2 = 1.00; // smooth vel
 real c3 = 1.00; // smooth coord (fujiwara)
 real d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 real d2 = 0.1;  // surface smooth cord (fujiwara)
 real alpha = 1;
 real beta = 1;

 real sigma = 1.0;

 real mu_in = 1.0;
 real mu_out = 0.01;

 real rho_in = 1.0;
 real rho_out = 0.001;

 real cfl = 0.5;

 string meshFile = "annular.msh";
 //string meshFile = "annularSquare.msh";

 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("DATA_DIR");
 meshDir += "/gmsh/3d/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 const char *mesh1 = mesh;
 m1.readMSH(mesh1);
 m1.setInterfaceBC();
 m1.setWallInterfaceBC();
 m1.setTriEdge(triEdge);
 m1.checkTriangleOrientation();
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
 m1.setWallAnnularBC();

 s1(m1);

 s1.setRe(Re);
 s1.setWe(We);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setC4(c4);
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
       << i << endl << endl;
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
  //s1.setInterface();
  s1.setInterfaceGeo();
  s1.unCoupled();

  save.saveMSH(mshFolder,"newMesh",i);
  save.saveVTK(vtkFolder,"sim",i);
  save.saveVTKSurface(vtkFolder,"sim",i);
  save.saveSol(binFolder,"sim",i);
  save.saveBubbleInfo(datFolder);
  save.chordalPressure(datFolder,"chordalPressure",i);

  s1.saveOldData();

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();

  iter++;
 }
  Model3D mOld = m1; 
  m1.setTriEdge(triEdge);

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  //m1.setNormalAndKappa();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  //m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume(0.005);
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  //m1.delete3DPoints();

  // surface operations
  m1.insertPointsByLength();
  //m1.insertPointsByCurvature();
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance();
  m1.contractEdgeByLength();
  //m1.removePointsByLength();
  //m1.flipTriangleEdge();
  //m1.removePointByNeighbourCheck();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setOFace();
  m1.setSurfaceConfig();
  m1.setWallAnnularBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMSH(mshFolder,"newMesh",iter-1);
  saveEnd.saveVTK(vtkFolder,"sim",iter-1);
  saveEnd.saveVTKSurface(vtkFolder,"sim",iter-1);
  saveEnd.saveSol(binFolder,"sim",iter-1);
  //saveEnd.saveVTU(vtkFolder,"sim",iter-1);
  //saveEnd.saveSolTXT(binFolder,"sim",iter-1);
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}


