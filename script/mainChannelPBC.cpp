/* File: mainChannelPBC.cpp                                            
 * Created on February 17th, 2014                                      
 * Author: Gustavo Charles P. de Oliveira                              
 * e-mail: tavolesliv@gmail.com
 * Maintainance: Gustavo Rabello dos Anjos                            
 * E-mail: gustavo.rabello@gmail.com                                   
 *
 * Description: version adapted from mainStep.cpp to include periodic 
 * boundary conditions. Suitable to run single-phase simulations for 
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
#include "CGSolver.h"
#include "PCGSolver.h"
#include "GMRes.h"
#include "Simulator3D.h"
#include "TElement.h"
#include "InOut.h"
#include "PetscSolver.h"
#include "petscksp.h"
#include "Periodic3D.h"

// Number of Phases (Do not change!)
#define NUMPHASES 1

// Main function
int main(int argc, char **argv)
{
 /* PRE-PROCESSING SECTION */

 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);


 int iter = 1;
 //double alpha = 1;
 double cfl = 1.0;

 double Re = 1.0;
 double Sc = 200;
 double Fr = 100;
 double mu_l = 1.0;
 double rho_l = 1.0;
 
 //** Solver and Pre-Conditioner Choice - pressure, velocity, scalar
 //Solver *solverP = new PCGSolver();
 //Solver *solverP = new PetscSolver(KSPGMRES,PCJACOBI);
 Solver *solverP = new PetscSolver(KSPCG,PCJACOBI);
 //Solver *solverP = new PetscSolver(KSPPREONLY,PCLU);
 //Solver *solverP = new PetscSolver(KSPLSQR,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 //Solver *solverV = new PCGSolver();
 Solver *solverC = new PCGSolver();

 //** Data Saving Directories
 //const char *binFolder  = "./bin/";
 //const char *datFolder  = "./dat/";
 //const char *datFolder  = "./sol/";
 //const char *txtFolder  = "./txt/";
 const char *vtkFolder = "/home/gcpoliveira/post-processing/vtk/3d/poiseuille-pbc/";
 //const char *vtkFolder = "/home/gcpoliveira/post-processing/vtk/3d/midwall-pbc/";
 //const char *vtkFolder = "/home/gcpoliveira/post-processing/vtk/3d/taylor-vortex/";
 //const char *vtkFolder = "/home/gcpoliveira/post-processing/vtk/3d/taylor-green-vortex/";

 //** Model Constructor
 Model3D m1;

 //** Mesh Loading: Selection of File and Path
 
 string selectionExtension = "msh";

 const char *mesh = NULL;
 if ( selectionExtension == "msh")
 {
  	//*** File
 	//string meshFile = "cuboid-3d.msh";
 	//string meshFile = "thesis-jet.msh";
 	//string meshFile = "cuboid-3d-w0.1.msh";
 	//string meshFile = "cylinder-3d.msh";
 	string meshFile = "cylinder-3d-L0.5D-w0.01.msh";

 	string meshDir = (string) getenv("MESH3D_DIR");
 	meshDir += "/" + meshFile;
 	const char *aux = meshDir.c_str();
	mesh = aux;

	//** Model Objects Call
 	m1.readMSH(mesh);
 	m1.setInterfaceBC();
 	m1.setTriEdge();
 	m1.mesh2Dto3D();
	m1.setJetMesh(0.05); // mesh geometrical transform
 	m1.setMapping();
 
	#if NUMGLEU == 5
 		m1.setMiniElement();
	#else
 		m1.setQuadElement();
	#endif
 	
 	m1.setVertNeighbour();
 	m1.setInOutVert(); // set of boundaryVert
	//m1.setGenericBC(); // only useful for PBC with phys. groups defined
	m1.setWallNormalVWBC();
	m1.setWallMovingPBC(0.0,0.0);
 	//m1.setOnePointPressureBC();
 }
 else if ( selectionExtension == "vtk" )
 {
	//*** File
 	string meshFile = "cube-h0.05.vtk";

 	string meshDir = (string) getenv("MESH3D_DIR");
 	meshDir += "/" + meshFile;
 	const char *aux = meshDir.c_str();
	mesh = aux;
	
	//** Model Objects Call
	m1.readVTK(mesh);
	#if NUMGLEU == 5
		m1.setMiniElement();
	#else
		m1.setQuadElement();
	#endif
	m1.setMapping();
	m1.setInOutVert();
	m1.setCubeVortexBC();
	m1.setOnePointPressureBC();
 }
 else
 {
	cerr << "Error Check mesh file: vtk/msh." << endl; 
	exit(1);
 }

 //* Periodic Objets Call
 Periodic3D pbc(m1);
 //pbc.MountPeriodicVectorsNew(m1);
 pbc.MountPeriodicVectors(m1);

 //** Simulator Objects Call
 Simulator3D sp(pbc,m1);
 
 //*** Further Setting
 //**** Physics
 sp.setRe(Re);
 sp.setSc(Sc);
 sp.setFr(Fr);
 sp.setMu(mu_l);
 sp.setRho(rho_l);

 //**** Numerics
 sp.setCfl(cfl);
 sp.setDtEulerian();

 //**** Solver
 sp.setSolverPressure(solverP);
 sp.setSolverVelocity(solverV);
 sp.setSolverConcentration(solverC);

 //*** Initial Conditions
 sp.init(); // default
 //sp.initTaylorVortex(); // Taylor vortex
 //sp.initTaylorGreenVortex();

 //*** Matrix Mounting and b.c. Setting - Velocity, Pressure
 sp.assemble();
 //sp.assembleSlip();
 sp.matMount();
 //sp.matMountC();

 //*** Starting Flow: Pressure Gradient Setting
 sp.setBetaPressureLiquid();


 /* PROCESSING / POST-PROCESSING SECTION */

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
  iter = sp.loadSolution("./","sim",atoi(*(argv+2)));
 }

 //** Save Objects Call
 InOut save(m1,sp); // cria objeto de gravacao
 
 //*** Mesh Information
 save.saveVTK(vtkFolder,"geometry");
 save.saveInfo("./","info",mesh);

 //*** Output (Initial Condition)
 save.saveVTK(vtkFolder,"initial",0);

 //*** Iterative Process (Temporal Loop)
 int nIter = 100;
 int nRe = 1;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   //**** Advective Term
   sp.stepSLPBCFix();
   //sp.stepNoConvection();
   //sp.stepSL();
   
   //**** B.C. update
   sp.setUnCoupledPBC(); 
   //sp.setUnCoupledBC(); 

   //**** Physical Effects
   //sp.setGravity("+X");
   sp.setBetaFlowLiq("+X");

   //**** r.h.s Vector
   sp.setRHS_PBC();

   //**** Periodic Copy
   sp.setCopyDirectionPBC("RL");
   
   //**** Matricial System Solution
   //sp.unCoupledPBCVector();
   sp.unCoupledPBC();
   //sp.unCoupled();
   
   //**** Solution Saving
   save.saveVTK(vtkFolder,"sim",iter);
   //save.saveVTU(vtkFolder,"sim",iter);
   //save.saveSol(binFolder,"sim",iter);

   //**** Updating Quantities
   sp.saveOldData();

   //**** Updating Time
   sp.timeStep();

   cout << "________________________________________ END of "
	    << iter << endl;

   iter++;
  }
 }

 PetscFinalize();
 return 0;
}
