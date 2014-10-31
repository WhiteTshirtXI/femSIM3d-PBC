/* File: mainChannelPBC.cpp                                            
 * Created on February 17th, 2014                                      
 * Author: Gustavo Charles P. de Oliveira                              
 * e-mail: tavolesliv@gmail.com
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

 /* PETSC Call */
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 /* PRE-PROCESSING SECTION */

 //** Numerical Parameters
 int iter = 1;
 //double alpha = 1;
 double cfl = 1.0;

 //** Physical Parameters
 double Re = 1000;
 double Sc = 200;
 double Fr = 10;
 double mu_l = 1.0002E-3;
 double rho_l = 998.63;
 double betaGrad = 1.0;

 // double cp_l = 1.0; // check if necessary for scalar
 
 string physGroup1 = "\"wallNoSlip\"";
 string physGroup2 = "\"wallNoSlip\"";
 
 //** Solver and Pre-Conditioner Choice - pressure, velocity, scalar
 //Solver *solverP = new PCGSolver();
 Solver *solverP = new PetscSolver(KSPCG,PCILU);
 //Solver *solverP = new PetscSolver(KSPGMRES,PCILU);
 Solver *solverV = new PetscSolver(KSPCG,PCILU);
 Solver *solverC = new PCGSolver();

 //** Data Saving Directories
 const char *vtkFolder = "/work/gcpoliveira/post-processing/3d/taylor-vortex/vtk/";
 const char *datFolder = "/work/gcpoliveira/post-processing/3d/taylor-vortex/dat/";
 const char *binFolder = "/work/gcpoliveira/post-processing/3d/taylor-vortex/bin/";
 const char *mshFolder = "/work/gcpoliveira/post-processing/3d/taylor-vortex/msh/";

 //** Model Constructor
 Model3D m1;

 /* Mesh Loading: Selection of File and Path
 * 
 * \remark: Gmsh seems not generating the volume mesh for vtk files.
 * Consequently, .msh should be used. Check it out!
 */

 string selectionExtension = "msh";

 const char *mesh = NULL;
 if ( selectionExtension == "msh")
 {
  	//*** File
 	string meshFile = "cuboid.msh";

 	string meshDir = (string) getenv("MESH3D_DIR");
 	meshDir += "/cuboid/" + meshFile;
 	const char *aux = meshDir.c_str();
	mesh = aux;

	//** Model Objects Call
 	m1.readMSH(mesh);
 	m1.setInterfaceBC();
 	m1.setTriEdge();
 	m1.mesh2Dto3D();
 	m1.setMapping();
 
	#if NUMGLEU == 5
 		m1.setMiniElement();
	#else
 		m1.setQuadElement();
	#endif
 	
 	m1.setNeighbour();
	m1.setVertNeighbour();
 	m1.setInOutVert();
	m1.setGenericBCPBCNewDuo(physGroup1,physGroup2);
	m1.setGenericBC();
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
	//m1.setNeumannPressureBC(); // change to setOnePointPressureBC()
 }
 else
 {
	cerr << "Error. Check mesh file: vtk/msh." << endl; 
	exit(1);
 }

 //* Periodic Objets Call
 Periodic3D pbc(m1);
 pbc.MountPeriodicVectorsNew("print");

 //** Simulator Objects Call
 Simulator3D sp(pbc,m1);
 
 //*** Further Setting
 //**** Physics
 sp.setRe(Re);
 sp.setSc(Sc);
 sp.setFr(Fr);
 sp.setMu(mu_l);
 sp.setRho(rho_l);
 //sp.setCp(cp_l); // check for scalar

 //**** Numerics
 sp.setCfl(cfl);
 sp.setDtEulerian();

 //**** Solver
 sp.setSolverPressure(solverP);
 sp.setSolverVelocity(solverV);
 sp.setSolverConcentration(solverC);

 //*** Initial Conditions
 //sp.init(); // default
 sp.initTaylorVortex(); // Taylor vortex

 //*** Starting Flow: Pressure Gradient Setting
 sp.setBetaPressureLiquid(betaGrad);

 //**** Mounting
 sp.assemble(); 
 //sp.assembleSlip(); 
 sp.matMount();
 //sp.matMountC();

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
 save.saveInfo(datFolder,"info",mesh);
 save.saveMeshInfo(datFolder);

 //*** Output (Initial Condition)
 save.saveVTK(vtkFolder,"initial",0);

 //*** Iterative Process (Temporal Loop)
 int nIter = 10;
 int nRe = 1;
 for( int i=0;i<nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: " 
	    << iter << endl;

   //**** Advective Term
   sp.stepSLPBCFix(); // semi-lagrangian repair

   sp.setUnCoupledBC(); 
   //sp.setUnCoupledCBC();

   //**** Physical Effects
   //sp.setGravity("+X");
   sp.setBetaFlowLiq("+X");

   //**** r.h.s Vector
   sp.setRHS();

   //**** Periodic Copy
   sp.setCopyDirectionPBC("RL");
   
   //**** Matricial System Solution
   sp.unCoupledPBCNew();
   //sp.unCoupledCPBCNew();
   
   //**** Solution Saving
   save.saveMSH(mshFolder,"newMesh",iter);
   save.saveVTKPBC(vtkFolder,"sim",iter,betaGrad);
   save.saveSol(binFolder,"sim",iter);
   save.saveMeshInfo(datFolder);

   //**** Updating Quantities
   sp.saveOldData();


   cout << "________________________________________ END of "
	    << iter << endl;

   iter++;
   
   //**** Updating Time
   sp.timeStep();
  
  }
 }

 PetscFinalize();
 return 0;
}
