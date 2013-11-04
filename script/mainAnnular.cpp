// =================================================================== //
// this is file main.cpp, created at 10-Jun-2007                       //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
//                                                                     //
// Obs.: Annular flow setup:                                           //
//                                                                     //
// Simulator3D: setAnnularALEBC()                                      //
//              mu and rho not zero for Z.Max() and Z.Min()            //
//              turn off applyBubbleVolumeCorrection()                 //
//                                                                     //
// Model3D: insertPointsByLength("flat")                                     //
//           - no insertion at the boundaries Z.Max() and Z.Min()      //
//           - uncomment double v2 = mapEdgeTri.Get(edge,2)              //
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
#include "Helmholtz3D.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "colors.h"
#include "Periodic3D.h"

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 
 int iter = 1;
 double Re = 91;
 double We = 0.23;
 double c1 = 0.0;  // lagrangian
 double c2 = 1.0;  // smooth vel
 double c3 = 10.0;  // smooth coord (fujiwara)
 double d1 = 1.0;  // surface tangent velocity u_n=u-u_t 
 double d2 = 0.1;  // surface smooth cord (fujiwara)
 double alpha = 1.0;

 double mu_in = 1.78E-5;
 double mu_out = 1.0E-3;

 double rho_in = 1.25;
 double rho_out = 1000;

 double cfl = 0.5;

 string meshFile = "annular-3d.msh";

 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "./bin/";
 const char *vtkFolder  = "./vtk/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 string meshDir = (string) getenv("MESH3D_DIR");
 //string meshDir = (string) getenv("MESH_DIR");
 meshDir += "/" + meshFile;
 //meshDir += "/gustavo-anjos-meshes/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;
 Simulator3D s1;

 const char *mesh1 = mesh;
 m1.readMSH(mesh1);
 m1.setInterfaceBC();
 m1.setTriEdge();
 m1.setWallInterfaceBC();
 m1.mesh2Dto3D();
 m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setSurfaceConfig();
 m1.setInitSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.setGenericBC();

 s1(m1);

 s1.setRe(Re);
 s1.setWe(We);
 s1.setC1(c1);
 s1.setC2(c2);
 s1.setC3(c3);
 s1.setD1(d1);
 s1.setD2(d2);
 s1.setAlpha(alpha);
 s1.setMu(mu_in,mu_out);
 s1.setRho(rho_in,rho_out);
 s1.setCfl(cfl);
 s1.initAnnular();
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
  s1.stepALE();
  s1.movePoints();
  s1.assemble();
  s1.matMount();
  s1.setUnCoupledBC();
  s1.setRHS();
  s1.setGravity("+Z");
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

  s1.timeStep();

  cout << color(none,magenta,black);
  cout << "________________________________________ END of " 
       << i << endl << endl;;
  cout << resetColor();

  iter++;
 }

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();
  m1.initMeshParameters();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion();
  m1.remove3dMeshPointsByDiffusion();
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  m1.remove3dMeshPointsByHeight();
  m1.delete3DPoints();

  // surface operations
  m1.smoothPointsByCurvature();

  m1.insertPointsByLength("flat");
  //m1.insertPointsByCurvature("flat");
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance("flat");
  m1.contractEdgesByLength("flat");
  //m1.removePointsByLength();
  //m1.flipTriangleEdges();

  //m1.removePointsByNeighbourCheck();
  //m1.checkAngleBetweenPlanes();
  /* **************************************** */

  //m1.mesh2Dto3DOriginal();
  m1.mesh3DPoints();
  m1.setMapping();
#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
  m1.setSurfaceConfig();
  m1.setGenericBC();

  Simulator3D s2(m1,s1);
  s2.applyLinearInterpolation(mOld);
  s1 = s2;
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);

  InOut saveEnd(m1,s1); // cria objeto de gravacao
  saveEnd.printMeshReport();
  saveEnd.saveMeshInfo(datFolder);
 }

 PetscFinalize();
 return 0;
}


