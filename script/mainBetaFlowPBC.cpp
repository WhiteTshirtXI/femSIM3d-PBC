/** 
 * \file    mainBetaFlowPBC.cpp
 * \author  Gustavo Peixoto de Oliveira 
 * \email   tavolesliv@gmail.com
 * \date    Created on September 10th, 2014
 *
 * \brief   Benchmark test for validation of single-phase flows impelled
 *          by the pressure gradient. (Starting flow.) 
 *
 * \details This benchmark test doesn't include the interface force and 
 *          sets the same properties for both fluids, thus reducing the 
 *          two-phase dynamics to single-phase.
 *
 *          \f$ \beta = \frac{ \Delta p }{ L } \f$  
 *
 * \remark  Test scope: single-phase::beta::PBC::fixed
 *               
 * DIAGRAM
 * =======
 * 
 *
 *                  top
 *           ---  ---  ---  ---      _
 *          |                  |     |
 *          |        g         |  
 *          |        |         |
 *          |        ¿         |
 *          |                  |
 *          |                  | 
 *          |                  |
 *          |                  |        
 *          |                  |     L 
 *          |                  |
 *          |                  |
 *          |                  |
 *          |                  |
 *          |                  |
 *          |             ^    |
 *          |             |    |
 *          |            beta  |     |
 *           ---  ---  ---  ---      _
 *                 bottom
 *   
 *          |-       D        -|
 *
 *
 * + Variables:
 *   - \f$ beta \f$: pressure gradient
 *   - g:            gravity
 *   - L:            length
 *   - D:            diameter
 *
 * + Boundary conditions:
 *   - according to the problem
 * 
 * + Physical forces:
 *   - gravity: 1.0, dimensionless
 *   - beta: average pressure drop/length
 *
 */
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

#define NUMPHASES 2

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 // bogdan's thesis 2010 (Bhaga and Weber, JFM 1980)
 int iter = 1;

 double Re = 30000.0; 
 double Sc = 1.0;
 double We = 115.662;
 double Fr = 1.0;
 double alpha = 1.0;

 double rho_in = 1.225;
 double rho_out =1.225; 

 double mu_out = 0.001;
 double mu_in = 0.001;

 double cfl = 0.8;

 // fixed
 double c1 = 0.0;      // lagrangian
 double c2 = 1.0;      // smooth vel 
 double c3 = 10.0;     // smooth coord (fujiwara)
 double d1 = 1.0;      // surface tangent vel = (u-ut)
 double d2 = 0.1;      // surface smooth coord (fujiwara)
 
 //string physGroup = "\"wallInflowU\"";
 string physGroup1 = "\"wallNormalV\"";
 string physGroup2 = "\"wallNormalW\"";

 double betaGrad = 32.0/Re;
 
 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverP = new PetscSolver(KSPCG,PCJACOBI);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 //Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/bin/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/vtk/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/taylor-vortex/msh/";
 
 string meshDir = (string) getenv("MESH3D_DIR");
 
 string meshFile = "cuboid.msh";
 
 meshDir += "/cuboid/" + meshFile;
 const char *mesh = meshDir.c_str();

 Model3D m1;

  cout << endl;
  cout << "--------------> STARTING FROM 0" << endl;
  cout << endl;

  const char *mesh1 = mesh;

  m1.readMSH(mesh1);
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
  m1.setInitSurfaceVolume();
  m1.setInitSurfaceArea();
  m1.setGenericBC();
  //m1.setGenericBCPBCNew(physGroup);
  m1.setGenericBCPBCNewDuo(physGroup1,physGroup2);
  
  Periodic3D pbc(m1);
  pbc.MountPeriodicVectorsNew("print");

  Simulator3D s1(pbc,m1);

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
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  //s1.init();
  s1.initTaylorVortex();
  s1.setBetaPressureLiquid(betaGrad);
  s1.setDtALETwoPhase();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
 
 /*
 // Point's distribution
 Helmholtz3D h1(m1);
 h1.setBC();
 h1.initRisingBubble();
 h1.assemble();
 h1.setk(0.2);
 h1.matMountC();
 h1.setUnCoupledCBC(); 
 h1.setCRHS();
 h1.unCoupledC();
 h1.setModel3DEdgeSize();
 */

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"sim",0);
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 200;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl;
   cout << resetColor();

   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepALEPBC();
   //s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   //s1.setGravity("+X");
   //s1.setBetaFlowLiq("+X");
   s1.setRHS();
   s1.setCopyDirectionPBC("RL");
   s1.unCoupledPBCNew();

   save.saveVTKPBC(vtkFolder,"sim",iter,betaGrad);
   save.saveMSH(mshFolder,"newMesh",iter);
   //save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurfacePBC(vtkFolder,"sim",iter,betaGrad);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);

   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   s1.timeStep();

   iter++;
  }
  /*
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initRisingBubble();
  h2.assemble();
  h2.setk(0.2);
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();
  h2.saveVTK(vtkFolder,"edge",iter-1);
  h2.saveChordalEdge(datFolder,"edge",iter-1);
  h2.setModel3DEdgeSize();
  */

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  //m1.setNormalAndKappa();
  //m1.initMeshParameters();

  // 3D operations
  //m1.insert3dMeshPointsByDiffusion(6.0);
  //m1.remove3dMeshPointsByDiffusion(0.5);
  //m1.removePointByVolume();
  //m1.removePointsByInterfaceDistance();
  //m1.remove3dMeshPointsByDistance();
  //m1.remove3dMeshPointsByHeight();
  //m1.delete3DPoints();

  // surface operations
  //m1.smoothPointsByCurvature();

  //m1.insertPointsByLength("flat");
  //m1.insertPointsByCurvature("flat");
  //m1.removePointsByCurvature();
  //m1.insertPointsByInterfaceDistance("flat");
  //m1.contractEdgesByLength("flat");
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
  //m1.setGenericBCPBCNew(physGroup);
  m1.setGenericBCPBCNewDuo(physGroup1,physGroup2);
  pbc.MountPeriodicVectorsNew("noPrint");

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


