/* \file mainPoiseuilleBeta.cpp                                            
 * \author Gustavo P. de Oliveira                              
 * \email: tavolesliv@gmail.com
 * \date Created on November 22nd, 2013                                      
 *
 * \description Poiseuille Starting Flow
 * 
 *  walls: no-slip
 *  inlet/outlet: OBC 
 *  flow initialized by pressure gradient $ Eu_{\beta} $
 *
 * \remark{ The pressure gradient term is included in 
 *  Simulator3D::unCoupledBeta() }
 */

#include <cmath>
#include "Model3D.h"
#include "CGSolver.h"
#include "Simulator3D.h"
#include "InOut.h"
#include "petscksp.h"
#include "PetscSolver.h"

#define NUMPHASES 1

int main(int argc, char **argv)
{
 PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

 int iter = 1;
 double Re = 1.0;
 double Sc = 2000;
 double alpha = 1;
 double cfl = 1.0;
 double mu_l = 1.0;
 double rho_l = 1.0;
 double EuBeta = 32.0; ///< modified Euler (pressure gradient)

 /* Flag to select pressure term placement: 
  * 1: RHS; 2: unCoupled() (requires smaller dt on setDt(dt) )
  */
 string method = "1"; 

 /* Overlapping phys. group */
 string physGroup = "\"wallNoSlip\"";

 /* \attention{ Certify that these environment 
  * variables are correctly defined. Otherwise, 
  * SIGABORT runtime error will be launched. }
  */
 string meshDir = (string) getenv("MESH3D_DIR");
 string ppd = getenv("POST_PROCESSING3D_DIR");

 string meshFile = "poiseuilleBeta.msh"; 
 meshDir += "/singlePhase/internal/" + meshFile;
 const char *mesh = meshDir.c_str();

 Solver *solverP = new PetscSolver(KSPCG,PCICC);
 Solver *solverV = new PetscSolver(KSPCG,PCICC);
 Solver *solverC = new PetscSolver(KSPCG,PCICC);

 string bin = "/poiseuilleBeta/bin/";
 string vtk = "/poiseuilleBeta/vtk/";
 string dat = "/poiseuilleBeta/dat/";
 bin = ppd + bin;
 vtk = ppd + vtk;
 dat = ppd + dat;
 const char *binFolder = bin.c_str();
 const char *vtkFolder = vtk.c_str();
 const char *datFolder = dat.c_str();

 Model3D m1;

  m1.readMSH(mesh);
  m1.setInterfaceBC();
  m1.setTriEdge();
  m1.mesh2Dto3D("QDYYAzpqa0.1");
  m1.setMapping();

#if NUMGLEU == 5
 m1.setMiniElement();
#else
 m1.setQuadElement();
#endif
 m1.setSurfaceConfig();
 
 // mesh statistics info
 m1.setInitSurfaceVolume();
 m1.setSurfaceVolume();
 m1.setInitSurfaceArea();
 m1.setSurfaceArea();
 m1.tetMeshStats();

 // boundary conditions
 m1.setGenericBC();

 Simulator3D s1(m1);

 s1.setBetaPressureLiquid(EuBeta); // pressure gradient

 s1.setRe(Re);
 s1.setSc(Sc); 
 s1.setAlpha(alpha);
 s1.setMu(mu_l);
 s1.setRho(rho_l);
 s1.setCfl(cfl);
 s1.setDt(0.01); 
 s1.init();
 s1.assemble();
 s1.setSolverPressure(solverP);
 s1.setSolverVelocity(solverV);
 s1.setSolverConcentration(solverC);

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
  iter = s1.loadSolution("./","sim",atoi(*(argv+2)));
  s1.setCfl(cfl);
  s1.setDtEulerian();
 }

 InOut save(m1,s1); // objeto para gravacao de dados
 save.saveVTK(vtkFolder,"geometry");
 save.saveMeshInfo(datFolder);
 save.saveInfo(datFolder,"info",mesh);

 int nIter = 40;
 int nRe = 1;
 for( int i=1;i<=nIter;i++ )
 {
  for( int j=0;j<nRe;j++ )
  {
   cout << "____________________________________ Iteration: "
	    << iter << endl;

    save.printSimulationReport();

	s1.stepSL();
	s1.assemble();
    s1.matMount();
    s1.matMountC();
    s1.setUnCoupledBC();

	s1.setBetaFlowLiq("+X"); 
	if ( method == "1" ) // pressure gradient term goes on rhs
	{
	  s1.setRHSBeta(); 
	  s1.unCoupled(); 
	}
	else 
	{
      s1.setRHS(); 
	  s1.unCoupledBeta();  
	}

	save.saveVTK(vtkFolder,"sim",iter);
	save.saveSol(binFolder,"sim",iter);

	s1.saveOldData();
	s1.timeStep();
	
	cout << "________________________________________ END of "
	     << iter << endl;

	iter++;
  }
 }
 PetscFinalize();
 return 0;
}


