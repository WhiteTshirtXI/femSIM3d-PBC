/** 
 * \file    mainRisingBubbleBetaFlowPBC.cpp
 * \author  Gustavo Peixoto de Oliveira 
 * \email   tavolesliv@gmail.com
 * \date    Created on August 19th, 2014
 *
 * \brief   Benchmark test for validation of the rising bubble flow. 
 *
 * \details This benchmark test considers the problem of a rising bubble
 *          which goes upward through buoyancy, but considering open
 *          walls on the streamwise direction. When opening the top and
 *          bottom walls, the flow is pushed down due to the gravity
 *          force. Then, a body force is placed to act against the
 *          gravitational field so as to stop the flow. Such force is
 *          given by the pressure gradient given by 
 *
 *          \f$ \beta = \frac{ \Delta p }{ L } \f$  
 *               
 * \remark  Test scope: two-phase::beta::gravity::NBC::fixed
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
 *          |                  |
 *          |                  |
 *          |                  | 
 *          |                  |
 *          |        E         |        
 *          |        ^         |     L 
 *          |        |         |
 *          |       ___        |
 *          |      /   \       |
 *          |     |  b  |      |
 *          |      \___/       |
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
 *   - E:            buoyancy force
 *   - \f$ beta \f$: pressure gradient
 *   - g:            gravity
 *   - L:            length
 *   - b:            bubble
 *   - D:            diameter
 *
 * + Boundary conditions:
 *   - top, bottom: \f$ \textbf{n} \cdot \nabla \textbf{v} = 0, p = 0 \f$
 *   - side:        \f$ \textbf{n} \cdot \nabla \textbf{p} = 0, \textbf{v} = 0 \f$
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

 double Re = 100; 
 double Sc = 1;
 double We = 115.662;
 double Fr = 1.0;
 double alpha = 1.0;


 double rho_in = 1.225;
 double rho_out =1350; 

 double mu_out = 1;
 double mu_in = 0.0000178;

 const char* _case = "3";

 // case 1
 if( strcmp( _case,"1") == 0 )
 {
  Re = sqrt(42.895); 
  mu_out = 2.73;
 }

 else if( strcmp( _case,"2") == 0 )
 {
  Re = 13.8487; // case 2
  mu_out = 1.28; 
 }

 else if( strcmp( _case,"3") == 0 )
 {
  Re = 33.0413; // case 3
  mu_out = 0.5396; // case 3
 }

 else if( strcmp( _case,"6") == 0 )
 {
  Re = sqrt(3892.856); // case 6
  mu_out = 0.2857; // case 6
 }

 else if( strcmp( _case,"7") == 0 )
 {
  Re = sqrt(18124.092); // case 7
  mu_out = 0.1324; // case 7
 }

 else if( strcmp( _case,"8") == 0 )
 {
  Re = sqrt(41505.729); // case 8 (extream)
  mu_out = 0.0875134907735; // extream
 }
 else
 {
  cerr << "test case " << _case << " not available!" << endl;
  exit(1);
 }

 double cfl = 0.1;
 
 //string meshFile = "rising-x.msh";
 string meshFile = "airWaterSugarPBC-wallLeftRight.msh";
 
 string physGroup = "\"wallNoSlip\"";
 double betaGrad = 1.0;

 const char* _frame = "fixed";
 //const char* _frame = "moving";

 // fixed
 double c1 = 0.0;      // lagrangian
 double c2 = 1.0;      // smooth vel 
 double c3 = 10.0;     // smooth coord (fujiwara)
 double d1 = 1.0;      // surface tangent vel = (u-ut)
 double d2 = 0.2;      // surface smooth coord (fujiwara)

 // moving
 if( strcmp( _frame,"moving") == 0 )
 {
  c1 = 0.0;      // lagrangian
  c2 = 1.0;      // smooth vel: OBS - different result with c1=0.0
  c3 = 10.0;      // smooth coord (fujiwara)
  d1 = 0.0;      // surface tangent velocity u_n=u-u_t 
  d2 = 0.1;      // surface smooth cord (fujiwara)
 }

 Solver *solverP = new PetscSolver(KSPCG,PCJACOBI);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 //Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 const char *binFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/bin/";
 const char *vtkFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/vtk/";
 const char *datFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/dat/";
 const char *mshFolder  = "/work/gcpoliveira/post-processing/3d/rising-beta/msh/";
 
 string meshDir = (string) getenv("MESH3D_DIR");
 if( strcmp( _frame,"moving") == 0 )
  meshDir += "/rising/movingFrame/" + meshFile;
 else
  meshDir += "/rising/" + meshFile;
  //meshDir += "/test/" + meshFile;
 
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
  m1.setGenericBCPBCNew(physGroup);
  m1.setGenericBC();

  Simulator3D s1(m1);
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
  s1.setMu(mu_in,mu_out);
  s1.setRho(rho_in,rho_out);
  s1.setCfl(cfl);
  s1.init();
  s1.setBetaPressureLiquid(betaGrad);
  s1.setDtALETwoPhase();
  s1.setSolverPressure(solverP);
  s1.setSolverVelocity(solverV);
  s1.setSolverConcentration(solverC);
 

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

 InOut save(m1,s1); // cria objeto de gravacao
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 double vinst=0;
 double vref=0;
 double xref=0;
 double xinit=0;
 double dx=0;
 if( strcmp( _frame,"moving") == 0 )
 {
  // moving
  xref = s1.getURef();
  xref = s1.getXRef();
  s1.setCentroidVelPos();
  xinit = s1.getCentroidPosXAverage();
 }

 int nIter = 30000;
 int nReMesh = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nReMesh;j++ )
  {

   cout << color(none,magenta,black);
   cout << "____________________________________ Iteration: " 
	    << iter << endl << endl;
   cout << resetColor();

   // moving
   if( strcmp( _frame,"moving") == 0 )
   {
	// moving frame
	dx = s1.getCentroidPosXAverage() - xinit;
	vinst = s1.getCentroidVelXAverage() + dx/s1.getDt();
    vref += vinst;
	xref += vref*s1.getDt();
	cout << "vref: " << vref << " xref: " << xref << endl;
	cout << "dx: " << dx << endl;
    s1.setUSol(vinst);
    m1.setGenericBC(vref);
    s1.setURef(vref);
	s1.setXRef(xref);
   }

   s1.setDtALETwoPhase();

   InOut save(m1,s1); // cria objeto de gravacao
   save.printSimulationReport();

   s1.stepALE();
   s1.movePoints();
   s1.assemble();
   s1.matMount();
   s1.setUnCoupledBC();
   s1.setGravity("-X");
   s1.setBetaFlowLiq("+X");
   s1.setRHS(); //<<
   s1.setInterfaceGeo();
   s1.unCoupledBeta();

   if ( i%15 == 0 )
   {
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTK(vtkFolder,"sim",iter);
   save.saveVTKSurface(vtkFolder,"sim",iter);
   save.saveSol(binFolder,"sim",iter);
   save.saveBubbleInfo(datFolder);
   }
   s1.saveOldData();

   cout << color(none,magenta,black);
   cout << "________________________________________ END of " 
	    << iter << endl << endl;;
   cout << resetColor();

   s1.timeStep();

   iter++;
  }
  Helmholtz3D h2(m1,h1);
  h2.setBC();
  h2.initRisingBubble();
  h2.assemble();
  h2.setk(0.2);
  h2.matMountC();
  h2.setUnCoupledCBC(); 
  h2.setCRHS();
  h2.unCoupledC();

  if ( i%15 == 0 )
  {
   h2.saveVTK(vtkFolder,"edge",iter-1);
   h2.saveChordalEdge(datFolder,"edge",iter-1);
   h2.setModel3DEdgeSize();
  }

  Model3D mOld = m1; 

  /* *********** MESH TREATMENT ************* */
  // set normal and kappa values
  m1.setNormalAndKappa();
  m1.initMeshParameters();

  // 3D operations
  m1.insert3dMeshPointsByDiffusion(6.0);
  m1.remove3dMeshPointsByDiffusion(0.5);
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
  m1.flipTriangleEdges();

  m1.removePointsByNeighbourCheck();
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
  if( strcmp( _frame,"moving") == 0 )
  {
    m1.setGenericBCPBCNew(physGroup);
    m1.setGenericBC(vref);
  }
  else
  {
    m1.setGenericBCPBCNew(physGroup);
    m1.setGenericBC();
  }

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


