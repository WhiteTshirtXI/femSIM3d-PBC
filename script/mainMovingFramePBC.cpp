/* File: mainMovingFramePBC.cpp                                            
 * Created on February 17th, 2014                                      
 * Author: Gustavo Charles P. de Oliveira                              
 * e-mail: tavolesliv@gmail.com
 * Maintainance: Gustavo Rabello dos Anjos                            
 * E-mail: gustavo.rabello@gmail.com                                   
 *
 * Description: version adapted from mainMovingFrame.cpp to include periodic 
 * boundary conditions. Suitable to run two-phase simulations for 
 * longitudinal channels. PBC is implemented to work with meshes whose
 * x-axis is the streamwise direction of spatial periodicity.
 *
 * \remark: PBC working for velocity and pressure. 
 *         To be implemented to further scalar fields, such as
 *         temperature.  
 * 
 *
 *  ============
 *   Guidelines
 *  ============
 * 
 *
 * # PRE-PROCESSING
 *
 * - Choose mesh file extension: vtk or msh;
 * - Introduce path to file;
 * - Call method to check spatial periodicity of the mesh;
 * 
 * # PBC
 *
 * - Impose pressure gradient based on input parameters;
 *
 *
 * # ITERATIVE STEP
 *
 * - Allow gravity effects or not;
 * - Allow pressure gradient or not;
 *
 */

// C++ Header Files
#include <cmath>

// Code Header Files
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

// Number of phases (Do not change!)
#define NUMPHASES 2

// Main function
int main(int argc, char **argv)
{
 /* PRE-PROCESSING SECTION */

 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
 //PetscInitializeNoArguments();

 //** Numerical Parameters
 // bogdan's thesis 2010 (Bhaga and Weber, JFM 1980)
 int iter = 1;
 double vel = 0;
 double c1 = 0.0;      // lagrangian
 double c2 = 0.0;      // smooth vel
 double c3 = 3.0;      // smooth coord (fujiwara)
 double d1 = 1.0;      // surface tangent velocity u_n=u-u_t 
 double d2 = 0.1;      // surface smooth cord (fujiwara)
 double alpha = 1.0;   // time discrete method
 double cfl = 1.0;

 //** Physical Parameters
 double Re = 14.5767; 
 double We = 0.2007;
 double Fr = 13.2160;
 double mu_in = 1.820E-05;
 double mu_out = 1.002E-3;
 double rho_in = 1.205;
 double rho_out = 0.0728;

 //*** File
 string meshFile = "bubble-elongated-3d-nb3.msh";
 
 
 //** Solver and Pre-Conditioner Choice - pressure, velocity, scalar
 Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 //Solver *solverV = new PetscSolver(KSPCG,PCJACOBI);

 //** Data Saving Directories
 const char *binFolder  = "./bin/";
 const char *mshFolder  = "./msh/";
 const char *datFolder  = "./dat/";
 //const char *vtkFolder  = "/home/gcpoliveira/post-processing/vtk/3d/slug-nb3/";
 const char *vtkFolder  = "./vtk/";
 
 
 //*** Peixoto's tree
 string meshDir = (string) getenv("MESH3D_DIR");
 meshDir += "/" + meshFile;
 const char *mesh = meshDir.c_str();

 //** Model Constructor
 Model3D m1;

 if( *(argv+1) == NULL )     
 {
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

 //*** B.C.
 //m1.setGenericBC();
 m1.setBubbleArrayPeriodicBC(); // <<<
  	
	
 //*** Periodic Constructor
 Periodic3D pbc(m1);
 pbc.MountPeriodicVectors(m1);
 pbc.SetGeometricalShape("default");
 
 //*** Simulator Constructor
 Simulator3D s2(pbc,m1);

 //*** Physics
 s2.setRe(Re);
 s2.setWe(We);
 s2.setFr(Fr);
 s2.setMu(mu_in,mu_out);
 s2.setRho(rho_in,rho_out);
  	
 //*** ALE parameters
 s2.setC1(c1);
 s2.setC2(c2);
 s2.setC3(c3);
 s2.setD1(d1);
 s2.setD2(d2);
 
 //*** Numerics
 s2.setCfl(cfl);
 s2.setDtALETwoPhase();
 s2.setAlpha(alpha);
  	
 //*** Solver
 s2.setSolverPressure(solverP);
 s2.setSolverVelocity(solverV);
 s2.setSolverConcentration(solverC);

 s2.init();

 /* Saving mesh info */
 InOut save(m1,s2);
 save.saveVTK(vtkFolder,"geometry");
 save.saveVTKSurface(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 //*** Output (Initial Condition)
 save.saveVTK(vtkFolder,"initial",0);

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

 double vinst=0;
 double vref=0;
 int nIter = 30;
 int nReMesh = 1;

 for( int i=1;i<=nIter;i++ )
 {
   for( int j=0;j<nReMesh;j++ )
   {
      cout << color(none,magenta,black);
      cout << "____________________________________ Iteration: " 
	       << iter << endl << endl;
      cout << resetColor();

      vinst = s2.getCentroidVelXAverage();
      vref += vinst;
      cout << vref << " " << vinst << endl;
      s2.setUSol(vinst);
      //m1.setGenericBC(vref);
	  m1.setBubbleArrayPeriodicBC(vref);
      s2.setURef(vref);

      s2.setDtALETwoPhase();

      save.printSimulationReport();

      //s2.stepLagrangian();
      //s2.stepALE();
      s2.stepSLPBCFix(); // <<< 

      s2.movePoints();

      //s2.assemble();
      s2.assemblePBC(); // <<<
      s2.matMount();

      //s2.setUnCoupledBC();
      s2.setUnCoupledPBC(); // <<<
      
	  //s2.setRHS();
      s2.setRHS_PBC(); // <<<
      
	  //s2.setGravity("-Z");

      //s2.setInterface();
      s2.setInterfaceGeo();
      
	  //s2.unCoupled();
      s2.unCoupledPBC(); // <<<

	  //*** Solution Saving 
      save.saveMSH(mshFolder,"newMesh",iter);
      save.saveVTK(vtkFolder,"sim",iter);
      save.saveVTKSurface(vtkFolder,"sim",iter);
      save.saveSol(binFolder,"sim",iter);
      save.saveBubbleInfo(datFolder);
      //save.crossSectionalVoidFraction(datFolder,"voidFraction",iter);
      s2.saveOldData();

      s2.timeStep();

      cout << color(none,magenta,black);
      cout << "________________________________________ END of " 
	       << iter << endl << endl;;
      cout << resetColor();

      iter++;
     }
     /* Points distribution - Helmholtz eq. */

     Helmholtz3D h2(m1,h1);
     h2.setBC();
     h2.initRisingBubble();
     h2.assemble();
     h2.setk(0.2);
     h2.matMountC();
     h2.setUnCoupledCBC(); 
     h2.setCRHS();
     h2.unCoupledC();
     h2.saveChordalEdge(datFolder,"edge",iter-1);
     h2.setModel3DEdgeSize();

     Model3D mOld = m1; 

     /* *********** MESH TREATMENT ************* */
     // set normal and kappa values
     m1.setNormalAndKappa();
     m1.initMeshParameters();

     /* 3D operations */
     
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
    m1.setInterfaceBC();
    //m1.setGenericBC(vref);
	m1.setBubbleArrayPeriodicBC(vref);

    Simulator3D s3(m1,s2);
    s3.applyLinearInterpolation(mOld);
    s2 = s3;
    s2.setSolverPressure(solverP);
    s2.setSolverVelocity(solverV);
    s2.setSolverConcentration(solverC);

    InOut saveEnd(m1,s2); // cria objeto de gravacao
    saveEnd.printMeshReport();
    saveEnd.saveMeshInfo(datFolder);
    
   } // end for

  } // end if 
 
  PetscFinalize();
  return 0;
}


