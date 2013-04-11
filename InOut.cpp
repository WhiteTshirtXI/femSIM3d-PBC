// =================================================================== // 
// this is file InOut, created at 23-Ago-2007                         //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "InOut.h"

//using namespace std;

InOut::InOut( Model3D &_m )
{
 m = &_m;
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 heaviside = m->getHeaviside();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 outflow = m->getOutflow();

 IEN = m->getIEN();
 surface = m->getSurface();
 surfMesh = m->getSurfMesh();
 interfaceDistance = m->getInterfaceDistance();
 triEdge = m->getTriEdge();
 inElem = m->getInElem();
 outElem = m->getOutElem();
 edgeSize = m->getEdgeSize();
 vertIdRegion = m->getVertIdRegion();
 elemIdRegion = m->getElemIdRegion();
 neighbourPoint = m->getNeighbourPoint();
 averageTriLength = m->getAverageTriLength();
 initSurfaceArea = m->getInitSurfaceArea();
 surfaceArea = m->getSurfaceArea();
 errorArea = m->getErrorArea();
 dArea = m->getDArea();
 initSurfaceVolume = m->getInitSurfaceVolume();
 surfaceVolume = m->getSurfaceVolume();
 errorVolume = m->getErrorVolume();
 dVolume = m->getDVolume();

 // surface mesh indexes:
 isp = m->getISP();
 ispc = m->getISPC();
 rsp = m->getRSP();
 rspn = m->getRSPN();
 rspc = m->getRSPC();
 flip = m->getFLIP();
 spc = m->getSPC();
 spp = m->getSPP();
 intet = m->getINTET();
 maxArea = m->getMaxArea();
 minArea = m->getMinArea();
 idMaxArea = m->getIdMaxArea();
 idMinArea = m->getIdMinArea();
 maxLength = m->getMaxLength();
 minLength = m->getMinLength();
 numSurfElems = m->getNumSurfElems();
 numSurfVerts = m->getNumSurfVerts();

 // volumetric mesh indexes:
 ip = m->getIP();
 ipd = m->getIPD();
 rp = m->getRP();
 rpi = m->getRPI();
 rpd = m->getRPD();
 rpdist = m->getRPDist();
 rph = m->getRPH();
 rpv = m->getRPV();
 csp = m->getCSP();
 maxVolume = m->getMaxVolume();
 minVolume = m->getMinVolume();
 idMaxVolume = m->getIdMaxVolume();
 idMinVolume = m->getIdMinVolume();
}

InOut::InOut( Model3D &_m, Simulator3D &_s )
{
 m = &_m;
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 heaviside = m->getHeaviside();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 outflow = m->getOutflow();
 IEN = m->getIEN();
 surface = m->getSurface();
 surfMesh = m->getSurfMesh();
 interfaceDistance = m->getInterfaceDistance();
 triEdge = m->getTriEdge();
 inElem = m->getInElem();
 outElem = m->getOutElem();
 edgeSize = m->getEdgeSize();
 vertIdRegion = m->getVertIdRegion();
 elemIdRegion = m->getElemIdRegion();
 neighbourPoint = m->getNeighbourPoint();
 averageTriLength = m->getAverageTriLength();
 initSurfaceArea = m->getInitSurfaceArea();
 surfaceArea = m->getSurfaceArea();
 errorArea = m->getErrorArea();
 dArea = m->getDArea();
 initSurfaceVolume = m->getInitSurfaceVolume();
 surfaceVolume = m->getSurfaceVolume();
 errorVolume = m->getErrorVolume();
 dVolume = m->getDVolume();

 // surface mesh indexes:
 isp = m->getISP();
 ispc = m->getISPC();
 rsp = m->getRSP();
 rspn = m->getRSPN();
 rspc = m->getRSPC();
 flip = m->getFLIP();
 spc = m->getSPC();
 spp = m->getSPP();
 intet = m->getINTET();
 maxArea = m->getMaxArea();
 minArea = m->getMinArea();
 idMaxArea = m->getIdMaxArea();
 idMinArea = m->getIdMinArea();
 maxLength = m->getMaxLength();
 minLength = m->getMinLength();

 // volumetric mesh indexes:
 ip = m->getIP();
 ipd = m->getIPD();
 rp = m->getRP();
 rpi = m->getRPI();
 rpd = m->getRPD();
 rpdist = m->getRPDist();
 rph = m->getRPH();
 rpv = m->getRPV();
 csp = m->getCSP();
 maxVolume = m->getMaxVolume();
 minVolume = m->getMinVolume();
 idMaxVolume = m->getIdMaxVolume();
 idMinVolume = m->getIdMinVolume();

 s = &_s;
 Re = s->getRe();
 Sc = s->getSc();
 We = s->getWe();
 Fr = s->getFr();
 dt = s->getDt();
 dtLagrangian = s->getDtLagrangian();
 dtSurfaceTension = s->getDtSurfaceTension();
 dtSemiLagrangian = s->getDtSemiLagrangian();
 dtGravity = s->getDtGravity();
 cfl = s->getCfl();
 alpha = s->getAlpha();
 beta = s->getBeta();
 simTime = s->getTime();
 mu_in = s->getMu_in();
 mu_out = s->getMu_out();
 rho_in = s->getRho_in();
 rho_out = s->getRho_out();
 sigma = s->getSigma();
 iter = s->getIter();
 c1 = s->getC1();
 c2 = s->getC2();
 c3 = s->getC3();
 d1 = s->getD1();
 d2 = s->getD2();
 uRef = s->getURef();
 vRef = s->getVRef();
 wRef = s->getWRef();
 xRef = s->getXRef();
 yRef = s->getYRef();
 zRef = s->getZRef();

 K = s->getK();
 M = s->getM();
 G = s->getG();
 D = s->getD();
 gx = s->getGx();
 gy = s->getGy();
 gz = s->getGz();
 uSol = s->getUSol();
 vSol = s->getVSol();
 wSol = s->getWSol();
 uALE = s->getUALE();
 vALE = s->getVALE();
 wALE = s->getWALE();
 pSol = s->getPSol();
 cSol = s->getCSol();
 kappa = s->getKappa();
 fint = s->getFint();
 heatFlux = s->getHeatFlux();
 mu = s->getMu();
 rho = s->getRho();
 //hSmooth = s->getHSmooth();
 gravity = s->getGravity();
 uSolOld = s->getUSolOld();
 vSolOld = s->getVSolOld();
 wSolOld = s->getWSolOld();
 pSolOld = s->getPSolOld();
 cSolOld = s->getCSolOld();
}

InOut::~InOut(){}

void InOut::saveSol( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 clVector vec = *uSol;
 vec.Append(*vSol);
 vec.Append(*wSol);
 vec.Append(*pSol);
 vec.Append(*cSol); // adicionando o vetor concentracao no vetor Vel+Pressao
 vec.Append(*uALE);
 vec.Append(*vALE);
 vec.Append(*wALE);

 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "-" + str + ".bin";
 const char* filenameUVWPC = fileUVWPC.c_str();

 ofstream UVWPC_file( filenameUVWPC,ios::binary | ios::trunc ); 

 UVWPC_file.write( (const char*) vec.GetVec(),vec.Dim()*sizeof(real) );

 UVWPC_file.close();

 copyLastFile(_dir,filenameUVWPC,_filename);

 cout << "solution No. " << _iter << " saved in binary" << endl;
 
} // fecha metodo saveSol 

void InOut::loadSol( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "UVWPC" + "-" + str + ".bin";
 const char* filenameUVWPC = fileUVWPC.c_str();

 clVector aux2(6*numNodes+2*numVerts); // vetor tambem carrega a concentracao

 ifstream UVWPC_file( filenameUVWPC,ios::in | ios::binary ); 

 UVWPC_file.read( (char*) aux2.GetVec(),aux2.Dim()*sizeof(real) );

 UVWPC_file.close();

 aux2.CopyTo(0,*uSol);
 aux2.CopyTo(numNodes,*vSol);
 aux2.CopyTo(2*numNodes,*wSol);
 aux2.CopyTo(3*numNodes,*pSol);
 aux2.CopyTo(3*numNodes+numVerts,*cSol);
 aux2.CopyTo(3*numNodes+2*numVerts,*uALE);
 aux2.CopyTo(4*numNodes+2*numVerts,*vALE);
 aux2.CopyTo(5*numNodes+2*numVerts,*wALE);

 cout << "solution No.  " << _iter << " read in binary" << endl;
 
} // fecha metodo loadSol 

void InOut::saveSolTXT( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 clVector vec = *uSol;
 vec.Append(*vSol);
 vec.Append(*wSol);
 vec.Append(*pSol);
 vec.Append(*cSol); // adicionando o vetor concentracao no vetor Vel+Pressao
 vec.Append(*uALE);
 vec.Append(*vALE);
 vec.Append(*wALE);

 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "-" + str + ".txt";
 const char* filenameUVWPC = fileUVWPC.c_str();

 ofstream UVWPC_file( filenameUVWPC, ios::trunc ); 

 UVWPC_file << setprecision(10) << scientific; 
 for( int i=0;i<vec.Dim();i++ )
  UVWPC_file << vec.Get(i) << endl;

 UVWPC_file.close();

 copyLastFile(_dir,filenameUVWPC,_filename);

 cout << "solution No. " << _iter << " saved in ascii" << endl;
 
} // fecha metodo saveSolTXT 

void InOut::saveTXT( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 clVector vec = *uSol;
 vec.Append(*vSol);
 vec.Append(*wSol);
 vec.Append(*pSol);
 vec.Append(*cSol); // adicionando o vetor concentracao no vetor Vel+Pressao

 int i,j,k;

 string fileUVWPC   = _dir;
 string fileKM      = _dir;
 string fileG       = _dir;
 string fileD       = _dir;
 string fileOutflow = _dir;
 string aux = _filename;

 fileUVWPC   += aux + "UVWPC" + "-" + str + ".dat";
 fileKM      += aux + "KM" + "-" + str + ".dat";
 fileG       += aux + "G" + "-" + str + ".dat";
 fileD       += aux + "D" + "-" + str + ".dat";
 fileOutflow += aux + "outflow" + "-" + str + ".dat";

 const char* filenameUVWPC = fileUVWPC.c_str();
 const char* filenameKM = fileKM.c_str();
 const char* filenameG = fileG.c_str();
 const char* filenameD = fileD.c_str();
 const char* filenameOutflow = fileOutflow.c_str();

 ofstream UVWPC_file( filenameUVWPC ); 
 ofstream KM_file( filenameKM ); 
 ofstream G_file( filenameG ); 
 ofstream D_file( filenameD ); 
 ofstream outflow_file( filenameOutflow ); 

 UVWPC_file << numVerts << numNodes << endl << endl;

 for( i=0;i<vec.Dim();i++ )
   UVWPC_file << vec.Get(i) << endl;

 UVWPC_file << endl;
 UVWPC_file.close();
 
 for( i=0;i<K->DimI();i++ )
 {
  for( j=0;j<K->DimJ();j++ )
  {
   KM_file << K->Get(i,j) << " " << M->Get(i,j) << endl;
  }
  for( k=0;k<G->DimJ();k++ )
  {
   G_file << G->Get(i,k) << endl;
  }
 }

 for( i=0;i<D->DimI();i++ )
 {
  outflow_file << outflow->Get(i) << endl;
  for( j=0;j<D->DimJ();j++ )
  {
   D_file << D->Get(i,j) << endl;
  }
 }

 KM_file.close();
 G_file.close();
 D_file.close();
 outflow_file.close();

 cout << "sistema no.  " << _iter << " gravado em TXT" << endl;

} // fecha metodo saveTXT 

void InOut::saveMatrix( clMatrix &_matrix,const char* _filename,string &_filename2 )
{
 // concatenando nomes para o nome do arquivo final
 string file = _filename;
 file += _filename2 + ".dat";
 const char* filename = file.c_str();

 ofstream matrixFile( filename ); 

 for( int i=0;i<_matrix.DimI();i++ )
 {
  for( int j=0;j<_matrix.DimJ();j++ )
   matrixFile << _matrix.Get(i,j) << " " ;
  matrixFile << endl;
 }

 matrixFile << endl;
 matrixFile.close();

 cout << "matriz oFace gravada em ASCII" << endl;

} // fecha metodo saveMatrix 

void InOut::saveVTK( const char* _dir,const char* _filename )
{
 IEN = m->getIEN();
 numElems = m->getNumElems();

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkHeader(vtkFile);
 vtkCoords(vtkFile);
 vtkCellArray(vtkFile);
 vtkCellType(vtkFile);
 vtkScalarHeader(vtkFile);
 vtkScalar(vtkFile,"pressure",*pc);
 vtkVector(vtkFile,"boundary_velocity",*uc,*vc,*wc);

 // este if existe pois nem todos os metodos tem cc
 if( cc->Dim() > 0 )
  vtkScalar(vtkFile,"concentration",*cc);

 if( heaviside->Dim() > 0 )
  vtkScalar(vtkFile,"heaviside",*heaviside);

 if( edgeSize->Dim() > 0 )
  vtkScalar(vtkFile,"edgeSize",*edgeSize);

 vtkFile.close();

 cout << "mesh saved in VTK" << endl;

} // fecha metodo saveVTK 

void InOut::saveVTK( const char* _dir,const char* _filename, int _iter )
{
 IEN = m->getIEN();
 numElems = m->getNumElems();

 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkHeader(vtkFile,_iter);
 vtkCoords(vtkFile);
 vtkCellArray(vtkFile);
 vtkCellType(vtkFile);

 vtkScalarCellHeader(vtkFile);
 if( elemIdRegion->Dim() > 0 )
  vtkScalarCell(vtkFile,"elemId",*elemIdRegion);

 vtkScalarHeader(vtkFile);
 vtkScalar(vtkFile,"pressure",*pSol);
 vtkVector(vtkFile,"velocity",*uSol,*vSol,*wSol);

 // version 3.98 of Paraview has Crinkle Slice feature
 setCutPlane(vtkFile); // set cut plane functions

 if( uALE->Dim() > 0 )
  vtkVector(vtkFile,"ALE_velocity",*uALE,*vALE,*wALE);

 // este if existe pois nem todos os metodos tem cc
 if( cSol->Dim() > 0 )
  vtkScalar(vtkFile,"concentration",*cSol);

 if( heaviside->Dim() > 0 )
  vtkScalar(vtkFile,"heaviside",*heaviside);

 if( kappa->Dim() > 0 )
  vtkScalar(vtkFile,"kappa",*kappa);

 if( interfaceDistance->Dim() > 0 )
  vtkScalar(vtkFile,"distance",*interfaceDistance);

 if( vertIdRegion->Dim() > 0 )
  vtkScalar(vtkFile,"vertId",*vertIdRegion);

 if( gravity->Dim() > 0 )
  vtkVector(vtkFile,"gravity",*gravity);

 if( fint->Dim() > 0 )
  vtkVector(vtkFile,"surface_force",*fint);

 if( heatFlux->Dim() > 0 )
  vtkScalar(vtkFile,"heatFlux",*heatFlux);

 if( edgeSize->Dim() > 0 )
  vtkScalar(vtkFile,"edgeSize",*edgeSize);

 vtkScalar(vtkFile,"viscosity",*mu);
 vtkScalar(vtkFile,"density",*rho);

 vtkFile.close();

 copyLastFile(_dir,filename,_filename);

 cout << "solution No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVtk

void InOut::saveVTK( const char* _coord1, const char* _coord2,
                     const char* _dir, const char* _filename, 
					 int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = (string) _dir + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkHeader(vtkFile,_iter);
 vtkCoords(vtkFile);

 clVector *coordPlane1,*coordPlane2;
 if( strcmp( _coord1,"x") == 0 || strcmp( _coord1,"X") == 0 )
  coordPlane1 = X;
 if( strcmp( _coord1,"y") == 0 || strcmp( _coord1,"Y") == 0 )
  coordPlane1 = Y;
 if( strcmp( _coord1,"z") == 0 || strcmp( _coord1,"Z") == 0 )
  coordPlane1 = Z;
 if( strcmp( _coord2,"x") == 0 || strcmp( _coord2,"X") == 0 )
  coordPlane2 = X;
 if( strcmp( _coord2,"y") == 0 || strcmp( _coord2,"Y") == 0 )
  coordPlane2 = Y;
 if( strcmp( _coord2,"z") == 0 || strcmp( _coord2,"Z") == 0 )
  coordPlane2 = Z;

 real plane1 = ( coordPlane1->Max()+coordPlane1->Min() )/2.0;
 real plane2 = ( coordPlane2->Max()+coordPlane2->Min() )/2.0;

 // conta numero de elementos
 int count = 0;
 for( int i=0;i<numElems;i++ )
 {
  int v1 = IEN->Get(i,0);
  int v2 = IEN->Get(i,1);
  int v3 = IEN->Get(i,2);
  int v4 = IEN->Get(i,3);
  bool planeTest = (coordPlane1->Get( v1 ) <  plane1) && 
                   (coordPlane1->Get( v2 ) <  plane1) && 
 	               (coordPlane1->Get( v3 ) <  plane1) && 
				   (coordPlane1->Get( v4 ) <  plane1) &&
                   (coordPlane2->Get( v1 ) <  plane2) && 
                   (coordPlane2->Get( v2 ) <  plane2) && 
 	               (coordPlane2->Get( v3 ) <  plane2) && 
				   (coordPlane2->Get( v4 ) <  plane2);

  real hTest = heaviside->Get(v1)+heaviside->Get(v2)+
	           heaviside->Get(v3)+heaviside->Get(v4) > 1.5;

  if( planeTest || hTest )
   count++;
 }
 
 vtkFile << "CELLS " << count << " " << 5*count << endl;
 vtkFile << setprecision(0) << fixed; 
 for( int i=0;i<numElems;i++ )
 {
  int v1 = IEN->Get(i,0);
  int v2 = IEN->Get(i,1);
  int v3 = IEN->Get(i,2);
  int v4 = IEN->Get(i,3);
  bool planeTest = (coordPlane1->Get( v1 ) <  plane1) && 
                   (coordPlane1->Get( v2 ) <  plane1) && 
 	               (coordPlane1->Get( v3 ) <  plane1) && 
				   (coordPlane1->Get( v4 ) <  plane1) &&
                   (coordPlane2->Get( v1 ) <  plane2) && 
                   (coordPlane2->Get( v2 ) <  plane2) && 
 	               (coordPlane2->Get( v3 ) <  plane2) && 
				   (coordPlane2->Get( v4 ) <  plane2);

  real hTest = heaviside->Get(v1)+heaviside->Get(v2)+
	           heaviside->Get(v3)+heaviside->Get(v4) > 1.5;

  if( planeTest || hTest )
  {
   vtkFile << "4 " << IEN->Get(i,0) << " "  
            	   << IEN->Get(i,1) << " " 
				   << IEN->Get(i,2) << " " 
				   << IEN->Get(i,3) << endl;
  }
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << count << endl;
 for( int i=0;i<numElems;i++ )
 {
  int v1 = IEN->Get(i,0);
  int v2 = IEN->Get(i,1);
  int v3 = IEN->Get(i,2);
  int v4 = IEN->Get(i,3);
  bool planeTest = (coordPlane1->Get( v1 ) <  plane1) && 
                   (coordPlane1->Get( v2 ) <  plane1) && 
 	               (coordPlane1->Get( v3 ) <  plane1) && 
				   (coordPlane1->Get( v4 ) <  plane1) &&
                   (coordPlane2->Get( v1 ) <  plane2) && 
                   (coordPlane2->Get( v2 ) <  plane2) && 
 	               (coordPlane2->Get( v3 ) <  plane2) && 
				   (coordPlane2->Get( v4 ) <  plane2);

  real hTest = heaviside->Get(v1)+heaviside->Get(v2)+
	           heaviside->Get(v3)+heaviside->Get(v4) > 1.5;

  if( planeTest || hTest )
   vtkFile << "10 ";
 }

 vtkFile << endl;
 vtkFile << endl;

 vtkScalarHeader(vtkFile);
 vtkScalar(vtkFile,"pressure",*pSol);
 vtkScalar(vtkFile,"concentration",*cSol);
 vtkVector(vtkFile,"velocity",*uSol,*vSol,*wSol);
 vtkVector(vtkFile,"ALE_velocity",*uALE,*vALE,*wALE);

 // este if existe pois nem todos os metodos tem cc
 if( cSol->Dim() > 0 )
  vtkScalar(vtkFile,"concentration",*cSol);

 if( heaviside->Dim() > 0 )
  vtkScalar(vtkFile,"heaviside",*heaviside);

 if( kappa->Dim() > 0 )
  vtkScalar(vtkFile,"kappa",*kappa);

 if( interfaceDistance->Dim() > 0 )
  vtkScalar(vtkFile,"distance",*interfaceDistance);

 if( vertIdRegion->Dim() > 0 )
  vtkScalar(vtkFile,"vertId",*vertIdRegion);

 if( gravity->Dim() > 0 )
  vtkVector(vtkFile,"gravity",*gravity);

 if( fint->Dim() > 0 )
  vtkVector(vtkFile,"surface_force",*fint);

 if( edgeSize->Dim() > 0 )
  vtkScalar(vtkFile,"edgeSize",*edgeSize);

 vtkScalar(vtkFile,"viscosity",*mu);
 vtkScalar(vtkFile,"density",*rho);

 vtkFile.close();

 copyLastFile(_dir,filename,_filename);

 cout << "solution Cut-Plane No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVtkHalf

void InOut::matrixPrint( clMatrix &_m, const char* _filename )
{
 int numRows = _m.DimI();
 int numColunms = _m.DimJ();
 char matrixOut[numRows][numColunms];
 
 string file = _filename;
 file += ".dat";
 const char* filename = file.c_str();

 ofstream matrixFile;
 matrixFile.open( filename );


 int zeros = 0;
 int nonZeros = 0;
 for( int i=0;i<numRows;i++ )
 {
  for( int j=0;j<numColunms;j++ )
  {
   if( _m.Get(i,j) != 0 || fabs(_m.Get(i,j)) > 10e-10 )
   {
	matrixOut[i][j] = 'X';
	nonZeros++;
   }
   else
   {
	matrixOut[i][j] = '-';
	zeros++;
   }
  }
 }

 cout << "total de elementos zeros = " << zeros << endl;
 cout << "total de elementos nao zeros = " << nonZeros << endl;

 for( int i=0;i<numRows;i++ )
 {
  for( int j=0;j<numColunms;j++ )
   matrixFile << matrixOut[i][j] << " " ;
  matrixFile << endl;
 }

 matrixFile.close();
} // fecha metodo matrixPrint


// fazer aparecer todos os campos com '-' e a diagonal com 'X' para
// matrizes diagonais
void InOut::matrixPrint( clDMatrix &_m,const char* _filename )
{
 int numColunms = _m.Dim();
 char matrixOut[numColunms][numColunms];
 
 string file = _filename;
 file += ".dat";
 const char* filename = file.c_str();

 ofstream matrixFile;
 matrixFile.open( filename );


 int zeros = 0;
 int nonZeros = 0;
 for( int i=0;i<numColunms;i++ )
 {
  if( _m.Get(i) != 0 || fabs(_m.Get(i)) > 10e-10 )
  {
   matrixOut[i][i] = 'X';
   nonZeros++;
  }
  else
  {
   matrixOut[i][i] = '-';
   zeros++;
  }
 }

 cout << "total of zero elements = " << zeros << endl;
 cout << "total of non-zero elements = " << nonZeros << endl;

 for( int i=0;i<numColunms;i++ )
 {
  for( int j=0;j<numColunms;j++ )
   matrixFile << matrixOut[i][j] << " " ;
  matrixFile << endl;
 }

 matrixFile.close();
} // fecha metodo matrixPrint para matrizes diagonais

void InOut::saveVonKarman(const char* _dir,const char* _filename)
{
 int count = 0;
 int iter = s->getIter();
 for( int i=0;i<numVerts;i++ )
 {
  if( Z->Get(i) == Z->Min() && Y->Get(i) == 0 )
  {
   stringstream ss1,ss2;  //convertendo int --> string
   string str1,str2;
   ss1 << iter;
   ss1 >> str1;
   ss2 << count;
   ss2 >> str2;

   string file = _dir;
   string aux = _filename;
   file += aux + str2 + "." + str1;
   const char* filename = file.c_str();

   ofstream vonKarmanFile;
   vonKarmanFile.open( filename );

   vonKarmanFile << setprecision(10) << scientific; 
   vonKarmanFile << "#z" << setw(17) 
	             << "F"  << setw(18)
	             << "G"  << setw(18)
	             << "H"  << setw(18)
	             << "c"  << setw(18)
	             << "p"  << setw(19)
	             << "mu"  << endl;

   for( int j=0;j<numVerts;j++ )
   {
	if( X->Get(j) == X->Get(i) && Y->Get(j) == Y->Get(i) )
	{
	 real radius = sqrt( X->Get(i)*X->Get(i) +
	                     Y->Get(i)*Y->Get(i) );
	 int vert = j;

	 //      X*u      Y*v
	 // F = ------ + ------
	 //      R^2      R^2
	 real F = X->Get(vert)*uSol->Get(vert)/(radius*radius)+  
	          Y->Get(vert)*vSol->Get(vert)/(radius*radius);
	 //      X*v      Y*u
	 // G = ------ - ------
	 //      R^2      R^2
	 real G = X->Get(vert)*vSol->Get(vert)/(radius*radius)-  
	          Y->Get(vert)*uSol->Get(vert)/(radius*radius);
	 //      
	 // H = -1*w
	 //     
	 real H = (-1)*wSol->Get(vert); 
	 real c = cSol->Get(vert);
	 real p = pSol->Get(vert);
	 real muValue = mu->Get(vert);

     vonKarmanFile << setw(16) << Z->Get(vert)  
	               << setw(18) << F
				   << setw(18) << G
				   << setw(18) << H
				   << setw(18) << c  
				   << setw(18) << p 
				   << setw(18) << muValue << endl;
	}
   }

   vonKarmanFile << endl;
   vonKarmanFile << fixed; 
   vonKarmanFile << "Radius = " << X->Get(i) << endl;

   vonKarmanFile.close();

   /* ----------- copying to file vk?.last ----------- */
   ifstream inFile( filename,ios::binary ); 

   string fileCopy = (string) _dir + (string) _filename + str2 + "." + "last";
   const char* filenameCopy = fileCopy.c_str();
   ofstream outFile( filenameCopy,ios::binary ); 

   outFile << inFile.rdbuf();
   inFile.close();
   outFile.close();
   /* ------------------------------------------------ */

   count++;
  }
 }
 cout << "von Karman num. " << iter << " saved in ASCII" << endl;
}

void InOut::saveDiskRadiusError(const char* _dir,
                                const char* _filename,
								const char* _fileAnalytic)
{
 int iter = s->getIter();

 stringstream ss1;  //convertendo int --> string
 string str1;
 ss1 << iter;
 ss1 >> str1;

 string file = _dir;
 string aux = _filename;
 file += aux + "." + str1;
 const char* filename = file.c_str();

 ofstream vonKarmanFile;
 vonKarmanFile.open( filename );
 vonKarmanFile << "#radius" << setw(10)
               << "U_error" << setw(18)
               << "V_error" << setw(18)
               << "W_error" << setw(18)
               << "c_error" << setw(20)
               << "UVW_error" << setw(19)
               << "UVWc_error" << endl;


 real auxA;
 real dist1,dist2;
 clMatrix solFile(2401,5); 
 clVector uExact(numVerts);
 clVector vExact(numVerts);
 clVector wExact(numVerts);
 clVector cExact(numVerts);

 ifstream fileA( _fileAnalytic,ios::in );

 if( !fileA )
 {
  cerr << "Esta faltando o arquivo de perfis!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !fileA.eof() )
 {
  for( int i=0;i<solFile.DimI();i++ )
  {
   fileA >> auxA;
   solFile.Set(i,0,auxA);
   fileA >> auxA;
   solFile.Set(i,1,auxA);
   fileA >> auxA;
   solFile.Set(i,2,auxA);
   fileA >> auxA;
   solFile.Set(i,3,auxA);
   fileA >> auxA;
   solFile.Set(i,4,auxA);
  }
 }

 int j;
 real omega = 1.0;
 for( int i=0;i<numVerts;i++ )
 {
  for( j=0;j<solFile.DimI()-1;j++ )
  {
   dist1 = fabs( Z->Get(i) - solFile(j,0) );
   dist2 = fabs( Z->Get(i) - solFile(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  auxA = ( solFile(j,1)*X->Get(i)-solFile(j,2)*Y->Get(i) )*omega; // F
  //aux = solFile(j,1); // F
  uExact.Set(i,auxA); 
  auxA = ( solFile(j,2)*X->Get(i)+solFile(j,1)*Y->Get(i) )*omega; // G
  //aux = solFile(j,2); // G
  vExact.Set(i,auxA);
  auxA = (-1)*solFile(j,3); // H (positive on file)
  wExact.Set(i,auxA);
  auxA = solFile(j,4); // C
  cExact.Set(i,auxA);
 }
 for( int i=0;i<numVerts;i++ )
 {
  if( Z->Get(i) == Z->Min() && 
	  Y->Get(i) == 0 )
  {
   real radius = sqrt( X->Get(i)*X->Get(i) +
                       Y->Get(i)*Y->Get(i) );
   if( radius > 0 )
   {
    real sumUDiff = 0.0;
    real sumVDiff = 0.0;
    real sumWDiff = 0.0;
    real sumcDiff = 0.0;
    real sumUVWDiff = 0.0;
    real sumUVWcDiff = 0.0;
    real sumU = 0.0;
    real sumV = 0.0;
    real sumW = 0.0;
    real sumc = 0.0;
    real sumUVW = 0.0;
    real sumUVWc = 0.0;
    for( int j=0;j<numVerts;j++ )
    {
     if( X->Get(j) == X->Get(i) && 
         Y->Get(j) == Y->Get(i) )
     {
      real UVW = uSol->Get(j)+
                 vSol->Get(j)+
				 wSol->Get(j); 
   
      real UVWc = UVW+cSol->Get(j); 
     
      real UVWExact = uExact.Get(j)+
                      vExact.Get(j)+
					  wExact.Get(j); 
   
      real UVWcExact = UVWExact+cExact.Get(j); 
     
      sumUDiff += (uSol->Get(j)-uExact.Get(j))*
                  (uSol->Get(j)-uExact.Get(j));
      sumVDiff += (vSol->Get(j)-vExact.Get(j))*
                  (vSol->Get(j)-vExact.Get(j));
      sumWDiff += (wSol->Get(j)-wExact.Get(j))*
                  (wSol->Get(j)-wExact.Get(j));
      sumcDiff += (cSol->Get(j)-cExact.Get(j))*
                  (cSol->Get(j)-cExact.Get(j));
     
      sumUVWDiff += (UVW-UVWExact)*
                    (UVW-UVWExact);
      sumUVWcDiff += (UVWc-UVWcExact)*
                     (UVWc-UVWcExact);
     
      sumU += uSol->Get(j)*
              uSol->Get(j); 
      sumV += vSol->Get(j)*
              vSol->Get(j); 
      sumW += wSol->Get(j)*
              wSol->Get(j); 
      sumc += cSol->Get(j)*
              cSol->Get(j); 
     
      sumUVW += UVW*
                UVW; 
      sumUVWc += UVWc*
                 UVWc; 
     }
    }
    /*  
     *           (  sum( sol[i] - sol_a )^2    )
     *  _e = sqrt( --------------------------- )
     *           (      sum( sol[i]^2 )        )
     * */
    real errorU = sqrt( sumUDiff/(sumU+EPS) );
    real errorV = sqrt( sumVDiff/(sumV+EPS) );
    real errorW = sqrt( sumWDiff/(sumW+EPS) );
    real errorc = sqrt( sumcDiff/(sumc+EPS) );
    real errorUVW = sqrt( sumUVWDiff/(sumUVW+EPS) );
    real errorUVWc = sqrt( sumUVWcDiff/(sumUVWc+EPS) );
   
    vonKarmanFile << setprecision(4) << fixed; 
    vonKarmanFile << setw(8) << X->Get(i);
    vonKarmanFile << setprecision(10) << scientific; 
    vonKarmanFile << setw(18) << errorU 
                  << setw(18) << errorV 
                  << setw(18) << errorW 
                  << setw(18) << errorc 
                  << setw(18) << errorUVW
                  << setw(18) << errorUVWc << endl;
   
   }
  }
 }
  vonKarmanFile.close();
  /* ----------- copying to file vk?.last ----------- */
  ifstream inFile( filename,ios::binary ); 

  string fileCopy = (string) _dir + (string) _filename + "." + "last";
  const char* filenameCopy = fileCopy.c_str();
  ofstream outFile( filenameCopy,ios::binary ); 

  outFile << inFile.rdbuf();
  inFile.close();
  outFile.close();
  /* ------------------------------------------------ */
 cout << "disk radius num. " << iter << " saved in ASCII" << endl;
}

void InOut::savePert(const char* _dir,const char* _filename,int _iter,
                     int vertice)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".dat";
 const char* filename = file.c_str();

 ofstream vonKarmanFile;
 vonKarmanFile.open( filename );

 clVector Yaux = *Y==0.0; // todos os pontos que estao no eixo X
 clVector Y0 = Yaux.Find(); // retorna o vertice dos pontos em X
 clVector Xaux = *X==X->Get( (int) Y0.Get(vertice) ); // escolhe um X qualquer
 clVector X0 = Xaux.Find(); // retorna os vertices com o mesmo X em Z variavel
 real raio = X->Get( (int) X0.Get(0)); // raio do ponto

 clVector aux2;
 for( int i=0;i<numVerts;i++ )
 {
  if( (Z->Get(i) == Z->Min()) &&
	(X->Get(i)*X->Get(i)+Y->Get(i)*Y->Get(i)>=raio*raio - 0.01) &&
	(X->Get(i)*X->Get(i)+Y->Get(i)*Y->Get(i)<=raio*raio + 0.01) )
  {
   aux2.AddItem(i);
  }
 }

 for( int ntheta=0;ntheta<aux2.Dim();ntheta++ )
 {
  Yaux = *Y==Y->Get( (int) aux2.Get(ntheta) );
  Y0 = Yaux.Find(); // retorna o vertice dos pontos em X
  Xaux = *X==X->Get( (int) Y0.Get(ntheta) );
  X0 = Xaux.Find(); // retorna os vertices com o mesmo X em Z variavel

  int v1 = (int) X0.Get(0);
  int v2;
  for( int i=0;i<X0.Dim();i++ )
  {
   v2 = (int) X0.Get(i);
   vonKarmanFile << Z->Get( v2 ) << " " 
	<< uSol->Get( v2 )/X->Get( v1 ) << " " 
	<< vSol->Get( v2 )/X->Get( v1 ) << " " 
	<< (-1)*wSol->Get( v2 ) << " " 
	<< cSol->Get( v2 ) << endl;
  }
  vonKarmanFile << endl;
  vonKarmanFile << "# theta = " 
   << atan(Y->Get( (int) aux2.Get(ntheta))/Y->Max())*180.0/3.141592
   << " graus" << endl;
  vonKarmanFile << "# raio = " << raio << endl;
  vonKarmanFile << "# X = " << X->Get( (int) Y0.Get(ntheta) ) << endl;
  vonKarmanFile << "# Y = " << Y->Get( (int) aux2.Get(ntheta) ) << endl;
  vonKarmanFile << endl;
 }

 vonKarmanFile.close();

 cout << "von Karman perturb no. " << _iter << " saved in ASCII" << endl;
}

void InOut::saveVortX(const char* _dir,const char* _filename,int _iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile;
 vtkFile.open( filename );

 clMatrix Gy = *gy;
 clMatrix Gz = *gz;
 clMatrix Mx(numNodes,numNodes);
 clVector vVert(numVerts);
 clVector wVert(numVerts);
 M->CopyTo(0,0,Mx);
 vSol->CopyTo(0,vVert);
 wSol->CopyTo(0,wVert);

 clVector IntVort(numNodes);
 IntVort = (Gy*wVert) - (Gz*vVert);
 clVector vort(numNodes);

 PCGSolver pcg;

 cout << " ----> calculando vorticidade em X -- " << endl;
 pcg.solve(0.000001,Mx,vort,IntVort);
 cout << " ------------------------------------ " << endl;

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Vorticidade em X 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << s->getTime() << endl;
 vtkFile << "ITERATION 1 1 int" << endl;
 vtkFile << _iter << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 vtkFile << endl;
 
 vtkFile << "CELLS " << numElems << " " << 5*numElems << endl;
 for( int i=0;i<numElems;i++ )
 {
  vtkFile << "4 " << IEN->Get(i,0) << " "  
                  << IEN->Get(i,1) << " " 
                  << IEN->Get(i,2) << " " 
                  << IEN->Get(i,3) << endl;
 };
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
  vtkFile << "10 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS vort double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em X no. " << _iter << " gravada em VTK" << endl;

}

void InOut::saveVortY(const char* _dir,const char* _filename,int _iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile;
 vtkFile.open( filename );

 clMatrix Gx = *gx;
 clMatrix Gz = *gz;
 clMatrix Mx(numNodes,numNodes);
 clVector uVert(numVerts);
 clVector wVert(numVerts);
 M->CopyTo(0,0,Mx);
 uSol->CopyTo(0,uVert);
 wSol->CopyTo(0,wVert);

 clVector IntVort(numNodes);
 IntVort = (Gx*wVert) - (Gz*uVert);
 clVector vort(numNodes);

 PCGSolver pcg;
 cout << " ----> calculando vorticidade em Y -- " << endl;
 pcg.solve(0.000001,Mx,vort,IntVort);
 cout << " ------------------------------------ " << endl;

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Vorticidade 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << s->getTime() << endl;
 vtkFile << "ITERATION 1 1 int" << endl;
 vtkFile << _iter << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 vtkFile << endl;
 
 vtkFile << "CELLS " << numElems << " " << 5*numElems << endl;
 for( int i=0;i<numElems;i++ )
 {
  vtkFile << "4 " << IEN->Get(i,0) << " "  
                  << IEN->Get(i,1) << " " 
                  << IEN->Get(i,2) << " " 
                  << IEN->Get(i,3) << endl;
 };
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
  vtkFile << "10 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS vort double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em Y no. " << _iter << " gravada em VTK" << endl;

}

void InOut::saveVortZ(const char* _dir,const char* _filename,int _iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile;
 vtkFile.open( filename );

 clMatrix Mx(numNodes,numNodes);
 clMatrix Gx = *gx;
 clMatrix Gy = *gy;
 clVector uVert(numVerts);
 clVector vVert(numVerts);
 M->CopyTo(0,0,Mx);
 uSol->CopyTo(0,uVert);
 vSol->CopyTo(0,vVert);

 clVector IntVort(numNodes);
 IntVort = (Gx*vVert) - (Gy*uVert);
 clVector vort(numNodes);

 PCGSolver pcg;
 cout << " ----> calculando vorticidade em Z -- " << endl;
 pcg.solve(0.000001,Mx,vort,IntVort);
 cout << " ------------------------------------ " << endl;

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Vorticidade 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << s->getTime() << endl;
 vtkFile << "ITERATION 1 1 int" << endl;
 vtkFile << _iter << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 vtkFile << endl;
 
 vtkFile << "CELLS " << numElems << " " << 5*numElems << endl;
 for( int i=0;i<numElems;i++ )
 {
  vtkFile << "4 " << IEN->Get(i,0) << " "  
                  << IEN->Get(i,1) << " " 
                  << IEN->Get(i,2) << " " 
                  << IEN->Get(i,3) << endl;
 };
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
  vtkFile << "10 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS vort double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em W no. " << _iter << " gravada em VTK" << endl;

}

void InOut::saveTime(const char* _comment)
{
 ofstream file( "relatorio.dat",ios::app );
 if( !file )
 {
  cerr << "error to open report file";
  exit(1);
 }

 time_t currentTime;

 time( &currentTime );
 file << asctime( localtime( &currentTime ) ) << "   - " << _comment << endl;
 file << endl;
}

void InOut::saveSimTime(int _iter)
{
 ofstream file( "./sim/simTime.dat",ios::trunc );
 if( !file )
 {
  cerr << "error to open simTime file";
  exit(1);
 }

 file << _iter << " " << s->getTime() << endl;
 file << endl;
}

void InOut::saveSimTime( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".dat";
 const char* filename = file.c_str();

 ofstream datFile( filename ); 
 datFile << _iter << " " << s->getTime() << endl;
 datFile << endl;
}

int InOut::loadIter()
{
 ifstream simTime( "./sim/simTime.dat",ios::in ); 

 real time;
 int iter;

 simTime >> iter;
 simTime >> time;
 s->setTime(time); 

 simTime.close();

 cout << "time = " << time << " " << "itereracao: " << iter << endl;

 return iter;
} // fecha metodo loadIter

int InOut::loadIter( const char* filename )
{
 char auxstr[255];
 real time;
 int iter;

 ifstream vtkFile( filename,ios::in );

 if( !vtkFile )
 {
  cerr << "VTK Mesh file is missing for TIME reading!" << endl;
  exit(1);
 }

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"TIME") != 0) )
  vtkFile >> auxstr;

 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> time;

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"ITERATION") != 0) )
  vtkFile >> auxstr;

 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> iter;

 vtkFile.close();

 s->setTime(time); 

 return iter;
} // fecha metodo loadIter

void InOut::saveInfo(const char* _dir,const char* _filename,const char* _mesh)
{
 string path = (string) _dir + (string) _filename + ".dat";
 const char* filePath = path.c_str();
 ofstream file( filePath,ios::app );

 time_t currentTime;

 time( &currentTime );
 file << endl;
 file << endl;
 file << "----------------------- BEGIN ----------------------" << endl; 
 file << "start process:    " << asctime( localtime( &currentTime ) );
 file << "mesh:             " << _mesh << endl;
 file << "numVerts:         " << numVerts << endl;
 file << "numNodes:         " << numNodes << endl;
 file << "numElems:         " << numElems << endl;
 file << "Lin Sys dim UVW:  " << 3*numNodes << " x " << 3*numNodes << endl;
 file << "Lin Sys dim P:    " << numVerts << " x " << numVerts << endl;
 file << "Lin Sys dim C:    " << numVerts << " x " << numVerts << endl;
 file << "Reynolds number:  " << Re << endl;
 file << "Schmidt number:   " << Sc << endl;
 file << "Froud number:     " << Fr << endl;
 file << "Webber number:    " << We << endl;
 file << "alpha number:     " << alpha << endl;
 file << "beta number:      " << beta << endl;
 file << "parameter c1:     " << c1 << endl;
 file << "parameter c2:     " << c2 << endl;
 file << "parameter c3:     " << c3 << endl;
 file << "parameter d1:     " << d1 << endl;
 file << "parameter d2:     " << d2 << endl;
 file << "liquid viscosity: " << mu_in << endl;
 file << "gas viscosity:    " << mu_out << endl;
 file << "liquid density:   " << rho_in << endl;
 file << "gas density:      " << rho_out << endl;
 file << "CFL number:       " << cfl << endl;
 file << "dt:               " << dt << endl;
 file << "----------------------------------------------------" << endl; 
 file << endl;
 file << endl;

 file.close();
 cout << "simulation INFO file saved!" << endl;
}

void InOut::oscillating(const char* _dir,const char* _filename, int _iter)
{
 // concatenando nomes para o nome do arquivo final
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(22) << "velocity U" 
                  << setw(17) << "velocity V" 
				  << setw(17) << "velocity W"
				  << setw(17) << "iteration" 
				  << endl;
 }

 real aux;

 int dim = surface->Dim();
 clVector xSurface(dim);
 clVector ySurface(dim);
 clVector zSurface(dim);
 clVector uSolSurface(dim);
 clVector vSolSurface(dim);
 clVector wSolSurface(dim);
 for( int i=0;i<dim;i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get(aux));
  ySurface.Set(i,Y->Get(aux));
  zSurface.Set(i,Z->Get(aux));
  uSolSurface.Set(i,uSol->Get(aux));
  vSolSurface.Set(i,vSol->Get(aux));
  wSolSurface.Set(i,wSol->Get(aux));
 }
 clVector xSurfaceAux = xSurface==xSurface.Max();
 clVector xSurfaceMax = xSurfaceAux.Find(); // retorna o vertice de maior X
 clVector ySurfaceAux = ySurface==ySurface.Max();
 clVector ySurfaceMax = ySurfaceAux.Find(); // retorna o vertice de maior Y
 clVector zSurfaceAux = zSurface==zSurface.Max();
 clVector zSurfaceMax = zSurfaceAux.Find(); // retorna o vertice de maior Z

 // retorna o valor de maior Y da interface 
 real pointX = uSolSurface.Get((int) xSurfaceMax.Get(0));
 real pointY = vSolSurface.Get((int) ySurfaceMax.Get(0));
 real pointZ = wSolSurface.Get((int) zSurfaceMax.Get(0));

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << pointX << " " 
                  << pointY << " " 
                  << pointZ << " " 
				  << setprecision(0) << fixed << _iter 
				  << endl;

 file.close();
}

void InOut::oscillating(int point1,int point2,int point3,const char* _filename)
{
 ifstream testFile( _filename );
 ofstream file( _filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(22) << "velocity U" 
                  << setw(17) << "velocity V" 
				  << setw(17) << "velocity W"
				  << endl;
 }

 // retorna o valor de maior Y da interface 
 real pointX = uSol->Get(point1);
 real pointY = vSol->Get(point2);
 real pointZ = wSol->Get(point3);

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << pointX << " " 
                  << pointY << " " 
                  << pointZ << " " 
				  << endl;
}

void InOut::oscillatingD(int point1,int point2,int point3,int point4,
						 int point5,int point6,const char* _filename)
{
 ifstream testFile( _filename );
 ofstream file( _filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(22) << "diameter X" 
                  << setw(17) << "diameter Y" 
				  << setw(17) << "diameter Z"
				  << endl;
 }

 // retorna o valor de maior Y da interface 
 real pointX1 = X->Get(point1);
 real pointX2 = X->Get(point2);
 real pointY1 = Y->Get(point3);
 real pointY2 = Y->Get(point4);
 real pointZ1 = Z->Get(point5);
 real pointZ2 = Z->Get(point6);

 // retorna o valor do maior diametro em Z na interface 
 real diameterX = fabs( pointX1-pointX2 ); 
 // retorna o valor do maior diametro em Y na interface 
 real diameterY = fabs( pointY1-pointY2 );
 // retorna o valor do maior diametro em Y na interface 
 real diameterZ = fabs( pointZ1-pointZ2 );
                 
 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << diameterX << " " 
                  << diameterY << " " 
                  << diameterZ << " " 
				  << endl;
 file.close();
}

void InOut::oscillatingD(const char* _dir,const char* _filename, int _iter)
{
 // concatenando nomes para o nome do arquivo final
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(22) << "diameter X" 
                 << setw(17) << "diameter Y" 
				 << setw(17) << "diameter Z"
				 << setw(16) << "iteration" 
				 << endl;
 }
 
 real aux;

 int dim = surface->Dim();
 clVector xSurface(dim);
 clVector ySurface(dim);
 clVector zSurface(dim);
 for( int i=0;i<dim;i++ )
 {
  aux = surface->Get(i);
  xSurface.Set(i,X->Get(aux));
  ySurface.Set(i,Y->Get(aux));
  zSurface.Set(i,Z->Get(aux));
 }

 // pega o 1o. valor da interface
 real xMax = xSurface.Get(0); 
 real yMax = ySurface.Get(0); 
 real zMax = zSurface.Get(0); 
 real xMin = xSurface.Get(0); 
 real yMin = ySurface.Get(0); 
 real zMin = zSurface.Get(0); 
 for( int i=1;i<dim;i++ )
 {
  if( xSurface.Get(i) > xMax )
   xMax = xSurface.Get(i);
  if( ySurface.Get(i) > yMax )
   yMax = ySurface.Get(i);
  if( zSurface.Get(i) > zMax )
   zMax = zSurface.Get(i);
  if( xSurface.Get(i) < xMin )
   xMin = xSurface.Get(i);
  if( ySurface.Get(i) < yMin )
   yMin = ySurface.Get(i);
  if( zSurface.Get(i) < zMin )
   zMin = zSurface.Get(i);
 }

 // retorna o valor do maior diametro em X na interface 
 real diameterX = xMax - xMin; 
 // retorna o valor do maior diametro em Y na interface 
 real diameterY = yMax - yMin; 
 // retorna o valor do maior diametro em Z na interface 
 real diameterZ = zMax - zMin; 

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << diameterX << " " 
                  << diameterY << " " 
                  << diameterZ << " " 
				  << setprecision(0) << fixed << _iter 
				  << endl;
 file.close();
}

void InOut::oscillatingKappa(const char* _dir,const char* _filename, int _iter)
{
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(17) << "kappa" 
				  << setw(21) << "iteration" 
				  << endl;
 }

 real aux = 0;
 for( int i=0;i<kappa->Dim();i++ )
  aux += kappa->Get(i);

 real kappaAverage = aux/kappa->Dim();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << kappaAverage << " " 
				  << setprecision(0) << fixed << _iter 
				  << endl;
 file.close();
}

void InOut::saveVTKSurface( const char* _dir,const char* _filename )
{
 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "TRI" + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkHeader(vtkFile);
 vtkSurfaceCoords(vtkFile);

 int numTri = 0;
 for( int i=0;i<surfMesh->numElems;i++ )
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   numTri++;

 vtkFile << setprecision(0) << fixed; 
 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   vtkFile << "3 " << surfMesh->IEN.Get(i,0) << " "  
                   << surfMesh->IEN.Get(i,1) << " " 
                   << surfMesh->IEN.Get(i,2) << endl;
 };
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << setprecision(0) << fixed; 
 vtkFile << "CELL_DATA " << numTri << endl;
 vtkFile << "NORMALS " << "cell_normal_test" << " double" << endl;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   vtkFile << m->getNormalElem(i).Get(0) << " " 
           << m->getNormalElem(i).Get(1) << " " 
	   	   << m->getNormalElem(i).Get(2) << endl;
 };
 vtkFile << endl;

 //vtkSurfaceCellHeader(vtkFile);
 //vtkSurfaceCellNormalVector(vtkFile,"cell_normal_test");

 vtkSurfaceScalarHeader(vtkFile);
 vtkSurfaceScalar(vtkFile,"pressure",*pc);
 vtkSurfaceScalar(vtkFile,"concentration",*cc);
 vtkSurfaceVector(vtkFile,"boundary_velocity",*uc,*vc,*wc);

 // este if existe pois nem todos os metodos tem cc
 if( cc->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"concentration",*cc);

 if( heaviside->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"heaviside",*heaviside);

 if( edgeSize->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"edgeSize",*edgeSize);

  vtkSurfaceNormalVector(vtkFile,"normal",surfMesh->xNormal,
	                                      surfMesh->yNormal,
						     			  surfMesh->zNormal);

 vtkSurfaceScalar(vtkFile,"distance",*interfaceDistance);

 vtkFile.close();

 string aux = (string) _filename + "TRI";
 const char* filenameAux = aux.c_str();
 copyLastFile(_dir,filename,filenameAux);

 cout << "surface mesh saved in VTK" << endl;

} // fecha metodo saveVTKSurface

void InOut::saveVTKSurface( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "TRI" + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkHeader(vtkFile,_iter);
 vtkSurfaceCoords(vtkFile);

 int numTri = 0;
 for( int i=0;i<surfMesh->numElems;i++ )
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   numTri++;

 vtkFile << setprecision(0) << fixed; 
 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   vtkFile << "3 " << surfMesh->IEN.Get(i,0) << " "  
                   << surfMesh->IEN.Get(i,1) << " " 
                   << surfMesh->IEN.Get(i,2) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << setprecision(0) << fixed; 
 vtkFile << "CELL_DATA " << numTri << endl;
 vtkFile << "NORMALS " << "cell_normal_test" << " double" << endl;
 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  if( surfMesh->elemIdRegion.Get(i) > 0 )
   vtkFile << m->getNormalElem(i).Get(0) << " " 
           << m->getNormalElem(i).Get(1) << " " 
	   	   << m->getNormalElem(i).Get(2) << endl;
 };
 vtkFile << endl;

 //vtkSurfaceCellHeader(vtkFile);
 //vtkSurfaceCellNormalVector(vtkFile,"cell_normal_test");

 vtkSurfaceScalarHeader(vtkFile);
 vtkSurfaceScalar(vtkFile,"pressure",*pSol);
 vtkSurfaceScalar(vtkFile,"concentration",*cSol);
 vtkSurfaceVector(vtkFile,"velocity",*uSol,*vSol,*wSol);
 vtkSurfaceVector(vtkFile,"ALE_velocity",*uALE,*vALE,*wALE);

 // este if existe pois nem todos os metodos tem cc
 if( cc->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"concentration",*cSol);

 if( heaviside->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"heaviside",*heaviside);

 if( edgeSize->Dim() > 0 )
  vtkSurfaceScalar(vtkFile,"edgeSize",*edgeSize);

 if( surfMesh->numInterfaces > 0 )
 {
  vtkSurfaceScalar(vtkFile,"kappa",*kappa);
  vtkSurfaceScalar(vtkFile,"distance",*interfaceDistance);
  vtkSurfaceVector(vtkFile,"gravity",*gravity);
  vtkSurfaceVector(vtkFile,"surface_force",*fint);
  vtkSurfaceNormalVector(vtkFile,"normal",surfMesh->xNormal,
	                                      surfMesh->yNormal,
						     			  surfMesh->zNormal);
 }

 vtkSurfaceScalar(vtkFile,"viscosity",*mu);
 vtkSurfaceScalar(vtkFile,"density",*rho);

 vtkFile.close();

 string aux = (string) _filename + "TRI";
 const char* filenameAux = aux.c_str();
 copyLastFile(_dir,filename,filenameAux);

 cout << "surface mesh No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVTKSurface

void InOut::bubblesDistance(const char* _dir,const char* _filename, int _iter)
{
 // concatenando nomes para o nome do arquivo final
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(17) << "Ymin1" 
				  << setw(18) << "Ymax1"
                  << setw(18) << "Ymin2" 
				  << setw(17) << "Ymax2"
				  << setw(17) << "dist1" 
				  << setw(17) << "dist2" 
				  << setw(15) << "iter" 
				  << endl;
 }

 real Ymax1=100;
 real Ymin1=-100;
 real Ymax2=-100;
 real Ymin2=100;
 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  // bubble 1 (Y<0)
  if( surfMesh->Y.Get(i) < 0 && 
	  heaviside->Get(i)==0.5 )
  {
   if(surfMesh->Y.Get(i)>Ymin1) Ymin1=surfMesh->Y.Get(i);
   if(surfMesh->Y.Get(i)<Ymax1) Ymax1=surfMesh->Y.Get(i);
  }
  // bubble 2 (Y>0)
  if( surfMesh->Y.Get(i) > 0 && heaviside->Get(i)==0.5 )
  {
   if(surfMesh->Y.Get(i)<Ymin2) Ymin2=surfMesh->Y.Get(i);
   if(surfMesh->Y.Get(i)>Ymax2) Ymax2=surfMesh->Y.Get(i);
  }
 }

 real dist1 = Ymin2-Ymin1;
 real dist2 = Ymax2-Ymax1;

 file << setprecision(10) << scientific; 
 file << setw(9) <<  s->getTime() << " " << Ymin1 << " " << Ymax1 << " " 
                                     << Ymin2 << " " << Ymax2 << " " 
							    	 << dist1 << " " << dist2 << " "
									 << _iter << endl;

 file.close();

 cout << "2 bubbles distances saved in ASCII " << dist1 << endl;
}

void InOut::saveMeshInfo(const char* _dir)
{
 real time = s->getTime();
 string file = (string) _dir + "meshingInfo.dat";
 const char* filename = file.c_str();

 ifstream testFile( filename );
 ofstream mesh( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file meshingInfo.dat" << endl;
 }
 else
 {
  cout << "Creating file meshingInfo.dat" << endl;
  mesh << "#time" << setw(21) << "numVerts" 
                  << setw(10) << "numNodes" 
				  << setw(10) << "numElems"
                  << setw(11) << "surfVerts" 
				  << setw(11) << "surfElems"
				  << setw(6) << "iter"
				  << endl;
 }


 mesh << setprecision(10) << scientific; 
 mesh << setw(9) <<  time << " " << setw(9) << numVerts << " " 
                                 << setw(9) <<numNodes << " " 
						   	     << setw(9) <<numElems << " "
							     << setw(10) <<surfMesh->numVerts << " "
							     << setw(10) <<surfMesh->numElems << " "
								 << setw(5) << iter << " "
								 << endl; 

 mesh.close();

 cout << "meshing info saved" << endl;

 for(int nb=0;nb<=elemIdRegion->Max();nb++ )
 {
  stringstream ss;  //convertendo int --> string
  string str;
  ss << nb;
  ss >> str;

  file = (string) _dir + "meshReport" + str + ".dat";
  const char* filename2 = file.c_str();

  ifstream testFile2( filename2 );
  ofstream mesh2( filename2,ios::app );
  if( testFile2 )
  {
   testFile2.close();
   cout << "appending on file meshReport" << nb << ".dat" << endl;
  }
  else
  {
   cout << "Creating file meshReport.dat" << endl;
   mesh2 << "#time" << setw(16) << "isp" 
					<< setw(7) << "ispc" 
					<< setw(6) << "rsp" 
					<< setw(7) << "rspn" 
					<< setw(7) << "rspc" 
					<< setw(5) << "ip" 
					<< setw(6) << "ipd" 
					<< setw(5) << "rp" 
					<< setw(6) << "rpi" 
					<< setw(6) << "rpd" 
					<< setw(6) << "rpdist" 
					<< setw(6) << "rph" 
					<< setw(6) << "rpv" 
					<< setw(6) << "csp" 
					<< setw(7) << "flip" 
					<< setw(6) << "spc" 
					<< setw(6) << "spp" 
					<< setw(8) << "intet" 
					<< setw(20) << "maxLength"
					<< setw(20) << "minLength"
					<< setw(20) << "maxArea"
					<< setw(20) << "minArea"
					<< setw(14) << "idMaxArea"
					<< setw(14) << "idMinArea"
					<< setw(20) << "maxVolume"
					<< setw(20) << "minVolume"
					<< setw(14) << "idMaxVolume"
					<< setw(14) << "idMinVolume"
					<< setw(7) << "iter"
					<< endl;
  }


  mesh2 << setprecision(10) << scientific; 
  mesh2 << setw(9) <<  time << " " << setw(4)  << isp[nb] << " " 
 								   << setw(6)  << ispc[nb] << " " 
								   << setw(5)  << rsp[nb] << " "
								   << setw(6)  << rspn[nb] << " "
								   << setw(6)  << rspc[nb] << " "
								   << setw(4)  << ip[nb] << " "
								   << setw(5)  << ipd[nb] << " "
								   << setw(4)  << rp[nb] << " "
								   << setw(5)  << rpi[nb] << " "
								   << setw(5)  << rpd[nb] << " "
								   << setw(5)  << rpdist[nb] << " "
								   << setw(5)  << rph[nb] << " "
								   << setw(5)  << rpv[nb] << " "
								   << setw(5)  << csp[nb] << " "
								   << setw(6)  << flip[nb] << " "
								   << setw(5)  << spc[nb] << " "
								   << setw(5)  << spp[nb] << " "
								   << setw(7)  << intet[nb] << " "
								   << setw(19) << maxLength[nb] << " "
								   << setw(19) << minLength[nb] << " "
								   << setw(19) << maxArea[nb] << " "
								   << setw(19) << minArea[nb] << " "
								   << setw(13) << idMaxArea[nb] << " "
								   << setw(13) << idMinArea[nb] << " "
								   << setw(19) << maxVolume[nb] << " "
								   << setw(19) << minVolume[nb] << " "
								   << setw(13) << idMaxVolume[nb] << " "
								   << setw(13) << idMinVolume[nb] << " "
								   << setw(6)  << iter << " "
								   << endl; 

  mesh2.close();
  cout << "mesh report surface (" << nb << ") saved" << endl;
 }

}

void InOut::saveVTU( const char* _dir,const char* _filename, int _iter )
{
 IEN = m->getIEN();
 numElems = m->getNumElems();

 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".vtu";
 const char* filename = file.c_str();

 ofstream vtuFile( filename ); 

 vtuFile << "<?xml version=\"1.0\"?> " << endl;
 vtuFile << endl;

 /* comments */
 vtuFile << "<!-- " << endl;
 vtuFile << "3D Simulation C++ " << endl;
 vtuFile << "ASCII " << endl;
 vtuFile << "DATASET UNSTRUCTURED_GRID " << endl;
 vtuFile << "FIELD FieldData 2 " << endl;
 vtuFile << "TIME 1 1 double " << endl;
 vtuFile << setprecision(10) << scientific;
 vtuFile << s->getTime() << endl;
 vtuFile << "ITERATION 1 1 int " << endl;
 vtuFile << setprecision(0) << fixed;
 vtuFile << _iter << endl;
 vtuFile << "--> " << endl;
 vtuFile << endl;

 /* header */
 vtuFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"> " << endl;
 vtuFile << " <UnstructuredGrid> " << endl;
 vtuFile << "  <Piece NumberOfPoints=\""<< numVerts 
         << "\" NumberOfCells=\"" << numElems << "\"> " << endl;
 vtuFile << "   <Points> " << endl;

 /* points */
 vtuFile << "    <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"> " << endl;

 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 vtuFile << "    </DataArray> " << endl;
 vtuFile << "   </Points> " << endl;

 /* cells */
 vtuFile << "   <Cells> " << endl;
 vtuFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"> " << endl;
 
 vtuFile << setprecision(0) << fixed;
 for( int i=0;i<numElems;i++ )
  vtuFile << "     " << IEN->Get(i,0) << " "  
                     << IEN->Get(i,1) << " " 
				     << IEN->Get(i,2) << " " 
				     << IEN->Get(i,3) << endl;
 vtuFile << "    </DataArray>" << endl;

 vtuFile << "    <DataArray type=\"Int32\"  Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;

 vtuFile << "     ";
 for( int i=1;i<=numElems;i++ )
  vtuFile << 4*i << " ";
 vtuFile << endl;
 vtuFile << "    </DataArray> " << endl;

 vtuFile << "    <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"> " << endl;
 vtuFile << "     ";
 for( int i=0;i<numElems;i++ )
  vtuFile << "10 ";
 vtuFile << endl;
 vtuFile << "    </DataArray> " << endl;
 vtuFile << "   </Cells> " << endl;

 /* data point header */
 vtuFile << "   <PointData  Scalars=\"Scalars\"> " << endl;

 /* data points */
 vtuFile << "    <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << pSol->Get(i) << endl;
 vtuFile << "    </DataArray>" << endl;

 vtuFile << "    <DataArray type=\"Float32\" Name=\"Concentration\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << cSol->Get(i) << endl;
 vtuFile << "    </DataArray>" << endl;

 /* vector data points */
 vtuFile << "    <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << uSol->Get(i) << " " 
                     << vSol->Get(i) << " " 
				     << wSol->Get(i) << endl;
 vtuFile << "    </DataArray>" << endl;

 vtuFile << "    <DataArray type=\"Float32\" Name=\"ALE Velocity\" NumberOfComponents=\"3\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << uALE->Get(i) << " " 
                     << vALE->Get(i) << " " 
				     << wALE->Get(i) << endl;
 vtuFile << "    </DataArray>" << endl;

 // este if existe pois nem todos os metodos tem cc
 if( heaviside->Dim() > 0 )
 {
  vtuFile << "    <DataArray type=\"Float32\" Name=\"Heaviside\" format=\"ascii\"> " << endl;
  vtuFile << setprecision(10) << scientific;
  for( int i=0;i<numVerts;i++ )
   vtuFile << "     " << heaviside->Get(i) << endl;
  vtuFile << "    </DataArray>" << endl;
 }

 if( kappa->Dim() > 0 )
 {
 vtuFile << "    <DataArray type=\"Float32\" Name=\"kappa\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << kappa->Get(i) << endl;
 vtuFile << "    </DataArray>" << endl;

 vtuFile << "    <DataArray type=\"Float32\" Name=\"surface force\" NumberOfComponents=\"3\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << fint->Get(i) << " " 
                     << fint->Get(i+numVerts) << " " 
				     << fint->Get(i+numVerts*2) << endl;
 vtuFile << "    </DataArray>" << endl;
 }

 /* end of file */
 vtuFile << "   </PointData> " << endl;
 vtuFile << "  </Piece> " << endl;
 vtuFile << " </UnstructuredGrid> " << endl;
 vtuFile << "</VTKFile> " << endl;

 vtuFile.close();

 cout << "solution No. " << _iter << " saved in VTU" << endl;

} // fecha metodo saveVTU


void InOut::saveConvergence(const char* _dir,const char* _filename)
{
 // concatenando nomes para o nome do arquivo final
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(18) << "uError" 
				  << setw(17) << "vError"
				  << setw(17) << "wError" 
				  << setw(17) << "pError" 
				  << setw(17) << "cError" 
				  << setw(19) << "uvwError" 
				  << setw(18) << "uvwpError" 
				  << setw(18) << "uvwpcError" 
				  << setw(12) << "iter" 
				  << endl;
 }

 // EPS to avoid division by 0
 real uDiff = ( (*uSol - *uSolOld).Abs() ).Sum();
 real uSum = ( uSol->Abs() ).Sum();
 real uError = (uDiff/(uSum+EPS))/dt;

 real vDiff = ( (*vSol - *vSolOld).Abs() ).Sum();
 real vSum = ( vSol->Abs() ).Sum();
 real vError = (vDiff/(vSum+EPS))/dt;

 real wDiff = ( (*wSol - *wSolOld).Abs() ).Sum();
 real wSum = ( wSol->Abs() ).Sum();
 real wError = (wDiff/(wSum+EPS))/dt;

 real pDiff = ( (*pSol - *pSolOld).Abs() ).Sum();
 real pSum = ( pSol->Abs() ).Sum();
 real pError = (pDiff/(pSum+EPS))/dt;

 real cDiff = ( (*cSol - *cSolOld).Abs() ).Sum();
 real cSum = ( cSol->Abs() ).Sum();
 real cError = (cDiff/(cSum+EPS))/dt;
 
 real uvwError = ( (uDiff+vDiff+wDiff) / (uSum+vSum+wSum+EPS) ) / dt;
 real uvwpError = ( (uDiff+vDiff+wDiff+pDiff) / (uSum+vSum+wSum+pSum+EPS) )/dt;
 real uvwpcError = ( (uDiff+vDiff+wDiff+pDiff+cDiff) /
                     (uSum+vSum+wSum+pSum+cSum+EPS) ) / dt;

 iter = s->getIter();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << uError << " " 
                  << vError << " " 
                  << wError << " " 
                  << pError << " " 
                  << cError << " " 
                  << uvwError << " " 
                  << uvwpError << " " 
                  << uvwpcError  
                  << setw(6) << iter 
				  << endl;

 file.close();

 cout << endl;
 cout << "relative error: " << uvwpcError << endl;
 cout << endl;

} // fecha metodo saveConvergence

//--------------------------------------------------
// void InOut::saveDiskError2(const char* _dir,const char* _filename )
// {
//  real aux;
//  real dist1,dist2;
//  clMatrix solFile(2401,5); 
//  clVector solF(numVerts);
//  clVector solG(numVerts);
//  clVector solH(numVerts);
//  clVector solC(numVerts);
// 
//  ifstream file( _filename,ios::in );
// 
//  if( !file )
//  {
//   cerr << "Esta faltando o arquivo de perfis!" << endl;
//   exit(1);
//  }
// 
//  // leitura do arquivo e transferencia para matriz
//  if( !file.eof() )
//  {
//   for( int i=0;i<solFile.DimI();i++ )
//   {
//    file >> aux;
//    solFile.Set(i,0,aux);
//    file >> aux;
//    solFile.Set(i,1,aux);
//    file >> aux;
//    solFile.Set(i,2,aux);
//    file >> aux;
//    solFile.Set(i,3,aux);
//    file >> aux;
//    solFile.Set(i,4,aux);
//   }
//  }
// 
//  int j;
//  for( int i=0;i<numVerts;i++ )
//  {
//   for( j=0;j<solFile.DimI()-1;j++ )
//   {
//    dist1 = fabs( Z->Get(i) - solFile(j,0) );
//    dist2 = fabs( Z->Get(i) - solFile(j+1,0) );
//    if( dist2 > dist1 ) break;
//   }
//   aux = solFile(j,1); // F
//   solF.Set(i,aux); 
//   aux = solFile(j,2); // G
//   solG.Set(i,aux);
//   aux = solFile(j,3); // H
//   solH.Set(i,aux);
//   aux = solFile(j,4); // C
//   solC.Set(i,aux);
//  }
// 
//  // this loop retrives all the points with Y=0 and Z varying from 0 to
//  // Z.Max for all radius X 
//  // count starts at 1 because 0 and the last radius are not used
//  // (boundary nodes)
//    string fileAux = (string) _dir + "diskError2" + ".dat";
//    const char* filename = fileAux.c_str();
// 
//    ofstream errorFile;
//    errorFile.open( filename );
// 
//    errorFile << setprecision(10) << scientific; 
//    errorFile << "#F_Error" 
// 	         << setw(17) << "G_Error"
// 			 << setw(18) << "H_Error" 
// 			 << setw(18) << "C_Error" 
// 			 << setw(20) << "FGH_Error" 
// 			 << setw(19) << "FGHC_Error" 
// 			 << setw(16) << "numVerts" 
// 			 << setw(10) << "numElems" 
// 			 << endl;
// 
//    real sumFDiff = 0.0;
//    real sumGDiff = 0.0;
//    real sumHDiff = 0.0;
//    real sumcDiff = 0.0;
//    real sumFGHDiff = 0.0;
//    real sumFGHcDiff = 0.0;
//    real sumF = 0.0;
//    real sumG = 0.0;
//    real sumH = 0.0;
//    real sumc = 0.0;
//    real sumFGH = 0.0;
//    real sumFGHc = 0.0;
//    for( int i=0;i<numVerts;i++ )
//    {
// 	 real radius = sqrt( X->Get(i)*X->Get(i) +
// 	                     Y->Get(i)*Y->Get(i) )+EPS;
// 
// 	 //      X*u      Y*v
// 	 // F = ------ + ------
// 	 //      R^2      R^2
// 	 real F = X->Get(i)*uSol->Get(i)/(radius*radius)+  
// 	          Y->Get(i)*vSol->Get(i)/(radius*radius);
// 	 //      X*v      Y*u
// 	 // G = ------ - ------
// 	 //      R^2      R^2
// 	 real G = X->Get(i)*vSol->Get(i)/(radius*radius)-  
// 	          Y->Get(i)*uSol->Get(i)/(radius*radius);
// 	 //      
// 	 // H = -1*w
// 	 //     
// 	 real H = (-1)*wSol->Get(i); 
// 
// 	 real c = cSol->Get(i);
// 
// 	 real FGH = F+G+H; 
// 	 real FGHc = F+G+H+c;
// 
// 	 real FExact = solF.Get(i);
// 	 real GExact = solG.Get(i);
// 	 real HExact = solH.Get(i);
// 	 real cExact = solC.Get(i);
// 
// 	 real FGHExact = FExact+GExact+HExact;
// 	 real FGHcExact = FExact+GExact+HExact+cExact;
// 
// 	 sumFDiff += (F-FExact)*(F-FExact);
// 	 sumGDiff += (G-GExact)*(G-GExact);
// 	 sumHDiff += (H-HExact)*(H-HExact);
// 	 sumcDiff += (c-cExact)*(c-cExact);
// 
// 	 sumFGHDiff += (FGH-FGHExact)*(FGH-FGHExact);
// 	 sumFGHcDiff += (FGHc-FGHcExact)*(FGHc-FGHcExact);
// 
// 	 sumF += F*F; 
// 	 sumG += G*G; 
// 	 sumH += H*H; 
// 	 sumc += c*c; 
// 
// 	 sumFGH += FGH*FGH; 
// 	 sumFGHc += FGHc*FGHc; 
//    }
//    /*  
// 	*           (  sum( sol[i] - sol_a )^2    )
// 	*  _e = sqrt( --------------------------- )
// 	*           (      sum( sol[i]^2 )        )
// 	* */
//    real errorF = sqrt( sumFDiff/(sumF+EPS) );
//    real errorG = sqrt( sumGDiff/(sumG+EPS) );
//    real errorH = sqrt( sumHDiff/(sumH+EPS) );
//    real errorc = sqrt( sumcDiff/(sumc+EPS) );
//    real errorFGH = sqrt( sumFGHDiff/(sumFGH+EPS) );
//    real errorFGHc = sqrt( sumFGHcDiff/(sumFGHc+EPS) );
// 
//    errorFile << errorF 
//              << setw(18) << errorG 
// 			 << setw(18) << errorH 
// 			 << setw(18) << errorc
// 			 << setw(18) << errorFGH
// 			 << setw(18) << errorFGHc
// 			 << fixed
// 			 << setw(10) << numVerts 
// 			 << setw(10) << numElems
// 			 << endl;
// 
//    errorFile << endl;
//    errorFile.close();
// 
//  cout << "relative error for disk saved in dat" << endl;
// }
//-------------------------------------------------- 

void InOut::saveDiskError(const char* _dir,
                          const char* _filename,const
						  char* _fileAnalytic)
{
 int iter = s->getIter();
 real aux;
 real dist1,dist2;
 clMatrix solFile(2401,5); 
 clVector uExact(numVerts);
 clVector vExact(numVerts);
 clVector wExact(numVerts);
 //clVector pExact(numVerts);
 clVector cExact(numVerts);

 ifstream fileA( _fileAnalytic,ios::in );

 if( !fileA )
 {
  cerr << "Esta faltando o arquivo de perfis!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !fileA.eof() )
 {
  for( int i=0;i<solFile.DimI();i++ )
  {
   fileA >> aux;
   solFile.Set(i,0,aux);
   fileA >> aux;
   solFile.Set(i,1,aux);
   fileA >> aux;
   solFile.Set(i,2,aux);
   fileA >> aux;
   solFile.Set(i,3,aux);
   fileA >> aux;
   solFile.Set(i,4,aux);
  }
 }

 int j;
 real omega = 1.0;
 for( int i=0;i<numVerts;i++ )
 {
  for( j=0;j<solFile.DimI()-1;j++ )
  {
   dist1 = fabs( Z->Get(i) - solFile(j,0) );
   dist2 = fabs( Z->Get(i) - solFile(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = ( solFile(j,1)*X->Get(i)-solFile(j,2)*Y->Get(i) )*omega; // F
  //aux = solFile(j,1); // F
  uExact.Set(i,aux); 
  aux = ( solFile(j,2)*X->Get(i)+solFile(j,1)*Y->Get(i) )*omega; // G
  //aux = solFile(j,2); // G
  vExact.Set(i,aux);
  aux = (-1)*solFile(j,3); // H (positive on file)
  wExact.Set(i,aux);
  aux = solFile(j,4); // C
  cExact.Set(i,aux);
 }

 // this loop retrives all the points with Y=0 and Z varying from 0 to
 // Z.Max for all radius X 
 // count starts at 1 because 0 and the last radius are not used
 // (boundary nodes)

 // concatenando nomes para o nome do arquivo final
 string fileAux = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  file << "#time" << setw(18) << "uError" 
				  << setw(17) << "vError"
				  << setw(17) << "wError" 
				  //<< setw(17) << "pError" 
				  << setw(17) << "cError" 
				  << setw(19) << "uvwError" 
				  //<< setw(19) << "uvwpError" 
				  << setw(18) << "uvwcError" 
				  //<< setw(18) << "uvwpcError" 
				  << setw(17) << "numVerts" 
				  << setw(10) << "numNodes" 
				  << setw(10) << "numElems" 
				  << setw(8) << "iter" 
				  << endl;
 }

 real sumUDiff = 0.0;
 real sumVDiff = 0.0;
 real sumWDiff = 0.0;
 //real sumpDiff = 0.0;
 real sumcDiff = 0.0;
 real sumUVWDiff = 0.0;
 //real sumUVWpDiff = 0.0;
 real sumUVWcDiff = 0.0;
 //real sumUVWpcDiff = 0.0;
 real sumU = 0.0;
 real sumV = 0.0;
 real sumW = 0.0;
 //real sump = 0.0;
 real sumc = 0.0;
 real sumUVW = 0.0;
 //real sumUVWp = 0.0;
 real sumUVWc = 0.0;
 //real sumUVWpc = 0.0;
 for( int i=0;i<numVerts;i++ )
 {
  real UVW = uSol->Get(i)+
             vSol->Get(i)+
			 wSol->Get(i); 
  //real UVWp = UVW + pSol->Get(i); 
  real UVWc = UVW + cSol->Get(i); 
  //real UVWpc = UVWp + cSol->Get(i); 

  real UVWExact = uExact.Get(i)+
                  vExact.Get(i)+
				  wExact.Get(i); 
  //real UVWpExact = UVWExact+pExact.Get(i); 
  real UVWcExact = UVWExact+cExact.Get(i); 
  //real UVWpcExact = UVWpExact+cExact.Get(i); 

  sumUDiff += (uSol->Get(i)-uExact.Get(i))*
              (uSol->Get(i)-uExact.Get(i));
  sumVDiff += (vSol->Get(i)-vExact.Get(i))*
              (vSol->Get(i)-vExact.Get(i));
  sumWDiff += (wSol->Get(i)-wExact.Get(i))*
              (wSol->Get(i)-wExact.Get(i));
  //sumpDiff += (pSol->Get(i)-pExact.Get(i))*
  //            (pSol->Get(i)-pExact.Get(i));
  sumcDiff += (cSol->Get(i)-cExact.Get(i))*
              (cSol->Get(i)-cExact.Get(i));

  sumUVWDiff += (UVW-UVWExact)*
                (UVW-UVWExact);
  //sumUVWpDiff += (UVWp-UVWpExact)*
  //               (UVWp-UVWpExact);
  sumUVWcDiff += (UVWc-UVWcExact)*
                 (UVWc-UVWcExact);
  //sumUVWpcDiff += (UVWpc-UVWpcExact)*
  //                (UVWpc-UVWpcExact);

  sumU += uSol->Get(i)*
          uSol->Get(i); 
  sumV += vSol->Get(i)*
          vSol->Get(i); 
  sumW += wSol->Get(i)*
          wSol->Get(i); 
  //sump += pSol->Get(i)*
  //        pSol->Get(i); 
  sumc += cSol->Get(i)*
          cSol->Get(i); 

  sumUVW += UVW*
            UVW; 
  //sumUVWp += UVWp*
  //           UVWp; 
  sumUVWc += UVWc*
             UVWc; 
  //sumUVWpc += UVWpc*
  //            UVWpc; 
 }

 /*  
  *           (  sum( sol[i] - sol_a )^2    )
  *  _e = sqrt( --------------------------- )
  *           (      sum( sol[i]^2 )        )
  * */
 real uError = sqrt( sumUDiff/(sumU+EPS) );
 real vError = sqrt( sumVDiff/(sumV+EPS) );
 real wError = sqrt( sumWDiff/(sumW+EPS) );
 //real pError = sqrt( sumpDiff/(sump+EPS) );
 real cError = sqrt( sumcDiff/(sumc+EPS) );
 real uvwError = sqrt( sumUVWDiff/(sumUVW+EPS) );
 //real uvwpError = sqrt( sumUVWpDiff/(sumUVWp+EPS) );
 real uvwcError = sqrt( sumUVWcDiff/(sumUVWc+EPS) );
 //real uvwpcError = sqrt( sumUVWpcDiff/(sumUVWpc+EPS) );

 iter = s->getIter();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
                  << uError << " " 
                  << vError << " " 
                  << wError << " " 
 //                 << pError << " " 
                  << cError << " " 
                  << uvwError << " " 
 //                 << uvwpError << " " 
                  << uvwcError << " " 
 //                 << uvwpcError << " " 
				  << fixed
				  << setw(9) << numVerts 
				  << setw(10) << numNodes
				  << setw(10) << numElems
				  << setw(8) << iter  
				  << endl;

 file.close();

 cout << "relative error for disk saved in dat" << endl;
}

void InOut::chordalPressure( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "-" + str + ".dat";
 const char* filename = file.c_str();

 ofstream vtkFile( filename,ios::app );

 vtkFile << "x-position" 
         << setw(17) << "y-position" 
		 << setw(17) << "z-position"
		 << setw(17) << "pressure"
		 << setw(17) << "kappa" 
		 << endl;

 // get iteration number
 iter = s->getIter();

 // xVert da malha nova
 real nPoints = 1000;
 clVector xVert(nPoints);
 clVector yVert(nPoints);
 clVector zVert(nPoints);

 for( int i=0;i<nPoints;i++ )
 {
  real dx = i * ( (X->Max()-X->Min()) )/(nPoints-1);
  real pos = X->Min()+dx;
  xVert.Set(i,pos);
 }
 yVert.SetAll( (Y->Max()+Y->Min())/2.0 );
 zVert.SetAll( (Z->Max()+Z->Min())/2.0 );

 // interpolacao linear em numVerts
 clMatrix interpLin = meshInterp(*m,xVert,yVert,zVert);
 clVector pLin = interpLin*(*pSol);

 clVector kappaVec(numVerts);
 for( int i=0;i<numVerts;i++ )
  kappaVec.Set(i,kappa->Get(i));

 clDMatrix kappaLin = interpLin*(kappaVec);

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<nPoints;i++ )
  vtkFile << setw(10) << xVert.Get(i) << " " 
          << setw(17) << yVert.Get(i) <<  " "
          << setw(17) << zVert.Get(i) <<  " "
          << setw(17) << pLin.Get(i) <<  " "
          << setw(17) << kappaLin.Get(i) << endl;

 vtkFile.close();

 copyLastFile(_dir,filename,_filename);

 cout << "chordal pressure No. " << _iter << " saved in dat" << endl;

} // fecha metodo chordalPressure

void InOut::crossSectionalPressure( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "-" + str + ".dat";
 const char* filename = file.c_str();

 ofstream pFile( filename ); 

 // xVert da malha nova
 real nPoints = 100;
 clVector xVert(nPoints*nPoints);
 clVector yVert(nPoints*nPoints);
 clVector zVert(nPoints*nPoints);

 int count = 0;
 for( int i=0;i<nPoints;i++ )
 {
  for( int j=0;j<nPoints;j++ )
  {
   real dx = i * (3.0-0.0)/(nPoints-1);
   xVert.Set(count,dx);
   real dy = j * (3.0-0.0)/(nPoints-1);
   yVert.Set(count,dy);
   count++;
  }
 }
 zVert.SetAll(1.5);

 // interpolacao linear em numVerts
 clMatrix interpLin = meshInterp(*m,xVert,yVert,zVert);
 clVector pLin(nPoints);
 pLin = interpLin*(*pSol);

 pFile << setprecision(10) << scientific;
 for( int i=0;i<nPoints*nPoints;i++ )
  pFile << setw(10) << xVert.Get(i) << " " 
                    << yVert.Get(i) << " " 
					<< pLin.Get(i) << " "
                    << setprecision(0) << fixed
                    << _iter << endl;

 pFile.close();

 cout << "cross sectional pressure No. " << _iter << " saved in dat" << endl;

} // fecha metodo chordalPressure

void InOut::vtkHeader(ofstream& _file)
{
 _file << "# vtk DataFile Version 1.0" << endl;
 _file << "3D Simulation C++" << endl;
 _file << "ASCII" << endl;
 _file << "DATASET UNSTRUCTURED_GRID" << endl;
 _file << endl;
}

void InOut::vtkHeader(ofstream& _file,int _iter)
{
 _file << "# vtk DataFile Version 1.0" << endl;
 _file << "3D Simulation C++" << endl;
 _file << "ASCII" << endl;
 _file << "DATASET UNSTRUCTURED_GRID" << endl;
 _file << "FIELD FieldData 8" << endl;
 _file << "TIME 1 3 double" << endl;
 _file << dt << " " << cfl << " " << s->getTime() << endl;
 _file << "ITERATION 1 1 int" << endl;
 _file << _iter << endl;
 _file << "NODES 1 3 int" << endl;
 _file << numVerts << " " << numNodes << " " << numElems << endl;
 _file << "PARAMETERS 1 4 float" << endl;
 _file << Re << " " << Sc << " " << Fr << " " << We << endl;
 _file << "PROPERTIES 1 5 float" << endl;
 _file << mu_in << " " << mu_out << " " 
       << rho_in << " " << rho_out << " " 
	   << sigma << endl;
 _file << "REFERENCEVELOCITY 1 6 float" << endl;
 _file << uRef << " " << vRef << " " << wRef << " " 
       << xRef << " " << yRef << " " << zRef << endl;
 _file << "COEFFICIENTS 1 7 float" << endl;
 _file << c1 << " " << c2 << " " << c3 << " " 
       << d1  << " " << d2 << " "
       << alpha << " " << beta << endl;

 if( surfMesh->elemIdRegion.Dim() > 0 )
 {
  _file << "CHARACTERISTICLENGTH 1 " << surfMesh->elemIdRegion.Max()+1
        << " float" << endl;
  for( int nb=0;nb<=surfMesh->elemIdRegion.Max();nb++ )
   _file << triEdge[nb] << " ";
  _file << endl;
 }
 else 
 {
  _file << "CHARACTERISTICLENGTH 1 1 float" << endl;
  _file << 0 << endl;
 }
 _file << endl;
}

void InOut::vtkCoords(ofstream& _file)
{
 _file << "POINTS " << numVerts << " double" << endl;
 //_file << "POINTS " << numNodes << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
 //for( int i=0;i<numNodes;i++ )
  _file << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 _file << endl;
}

void InOut::vtkSurfaceCoords(ofstream& _file)
{
 _file << "POINTS " << surfMesh->numVerts << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << surfMesh->X.Get(i) << " " 
        << surfMesh->Y.Get(i) << " " 
		<< surfMesh->Z.Get(i) << endl;

 _file << endl;
}

void InOut::vtkCellArray(ofstream& _file)
{
 _file << "CELLS " << numElems << " " << 5*numElems << endl;
 //_file << "CELLS " << numElems << " " << 11*numElems << endl;
 _file << setprecision(0) << fixed;
 for( int i=0;i<numElems;i++ )
 {
  _file << "4 " << IEN->Get(i,0) << " "  
                << IEN->Get(i,1) << " " 
                << IEN->Get(i,2) << " " 
                << IEN->Get(i,3) << endl;
 }
 _file << endl;
}

void InOut::vtkCellType(ofstream& _file)
{
 _file <<  "CELL_TYPES " << numElems << endl;
 _file << setprecision(0) << fixed;
 for( int i=0;i<numElems;i++ )
  _file << "10 ";

 _file << endl;
 _file << endl;
}

void InOut::vtkScalarHeader(ofstream& _file)
{
 _file << "POINT_DATA " << numVerts << endl;
}

void InOut::vtkSurfaceScalarHeader(ofstream& _file)
{
 _file << "POINT_DATA " << surfMesh->numVerts << endl;
}

void InOut::vtkScalarCellHeader(ofstream& _file)
{
 _file << "CELL_DATA " << numElems << endl;
}

void InOut::vtkSurfaceCellHeader(ofstream& _file)
{
 _file << "CELL_DATA " << surfMesh->numElems << endl;
}

void InOut::vtkScalar(ofstream& _file,string _name,clVector &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkScalar(ofstream& _file,string _name,clDMatrix &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkScalarCell(ofstream& _file,string _name,clVector &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<numElems;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkScalarCell(ofstream& _file,string _name,clDMatrix &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<numElems;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkSurfaceScalar(ofstream& _file,string _name,clVector &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkSurfaceScalar(ofstream& _file,string _name,clDMatrix &_scalar)
{
 _file << "SCALARS " << _name << " double" << endl;
 _file << "LOOKUP_TABLE default"  << endl;

 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _scalar.Get(i) << endl;

 _file << endl;
}

void InOut::vtkVector(ofstream& _file,string _name,clVector &_v)
{
 _file << "VECTORS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  _file << _v.Get(i) << " " 
        << _v.Get(i+numNodes) << " " 
		<< _v.Get(i+numNodes*2) << endl;
 _file << endl;
}

void InOut::vtkVector(ofstream& _file,string _name,
                      clVector &_vx,clVector &_vy,clVector &_vz)
{
 _file << "VECTORS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  _file << _vx.Get(i) << " " 
        << _vy.Get(i) << " " 
		<< _vz.Get(i) << endl;
 _file << endl;
}

void InOut::vtkSurfaceVector(ofstream& _file,string _name,clVector &_v)
{
 _file << "VECTORS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _v.Get(i) << " " 
        << _v.Get(i+numNodes) << " " 
		<< _v.Get(i+numNodes*2) << endl;
 _file << endl;
}

void InOut::vtkSurfaceVector(ofstream& _file,string _name,
                             clVector &_vx,clVector &_vy,clVector &_vz)
{
 _file << "VECTORS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _vx.Get(i) << " " 
        << _vy.Get(i) << " " 
		<< _vz.Get(i) << endl;
 _file << endl;
}

void InOut::vtkSurfaceNormalVector(ofstream& _file,string _name,clVector &_v)
{
 _file << "NORMALS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _v.Get(i) << " " 
        << _v.Get(i+numNodes) << " " 
		<< _v.Get(i+numNodes*2) << endl;
 _file << endl;
}

void InOut::vtkSurfaceNormalVector(ofstream& _file,string _name,
                             clVector &_vx,clVector &_vy,clVector &_vz)
{
 _file << "NORMALS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _vx.Get(i) << " " 
        << _vy.Get(i) << " " 
		<< _vz.Get(i) << endl;
 _file << endl;
}

void InOut::vtkSurfaceCellNormalVector(ofstream& _file,string _name)
                                  
{
 _file << "NORMALS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numElems;i++ )
  _file << m->getNormalElem(i).Get(0) << " " 
        << m->getNormalElem(i).Get(1) << " " 
		<< m->getNormalElem(i).Get(2) << endl;
 _file << endl;
}

void InOut::vtkSurfaceCellNormalVector(ofstream& _file,string _name,
                             clVector &_vx,clVector &_vy,clVector &_vz)
{
 _file << "NORMALS " << _name << " double" << endl;
 _file << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  _file << _vx.Get(i) << " " 
        << _vy.Get(i) << " " 
		<< _vz.Get(i) << endl;
 _file << endl;
}


void InOut::saveMSH( const char* _dir,const char* _filename )
{
 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + ".msh";
 const char* filename = file.c_str();

 ofstream mshFile( filename ); 

 mshFile << "$MeshFormat" << endl;
 mshFile << "2.1 0 8" << endl;
 mshFile << "$EndMeshFormat" << endl;
 mshFile << "$PhysicalNames" << endl;
 mshFile << surfMesh->phyNames.size() << endl;
 for( int nb=0;nb<(int)surfMesh->phyNames.size();nb++ )
   mshFile << "2 " << nb+1 << " " << surfMesh->phyNames.at(nb) << endl;
 mshFile << "$EndPhysicalNames" << endl;
 mshFile << "$Nodes" << endl;
 mshFile << surfMesh->numVerts << endl;

 mshFile << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  mshFile << i+1 << " " << X->Get(i) << " " 
                        << Y->Get(i) << " " 
						<< Z->Get(i) << endl;

 mshFile << "$EndNodes" << endl;
 mshFile << "$Elements" << endl;
 mshFile << surfMesh->numElems << endl;

 mshFile << setprecision(0) << fixed;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  int v1 = surfMesh->IEN.Get(i,0);
  int v2 = surfMesh->IEN.Get(i,1);
  int v3 = surfMesh->IEN.Get(i,2);

  mshFile << i+1 
          << " 2" 
		  << " 2" 
		  << " " << surfMesh->idRegion.Get(i) + 1 
		  << " " << surfMesh->idRegion.Get(i) + 1  // surface number Gmsh
		  << " " << v1+1 
		  << " " << v2+1 
		  << " " << v3+1 
		  << endl; 

 }
 mshFile << "$EndElements" << endl;

 mshFile.close();

 cout << "mesh saved in MSH" << endl;

} // fecha metodo saveMSH

void InOut::saveMSH( const char* _dir,const char* _filename, int _iter )
{
 SurfaceMesh *surfMesh = m->getSurfMesh();

 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "-" + str + ".msh";
 const char* filename = file.c_str();

 ofstream mshFile( filename ); 

 mshFile << "$MeshFormat" << endl;
 mshFile << "2.1 0 8" << endl;
 mshFile << "$EndMeshFormat" << endl;
 mshFile << "$PhysicalNames" << endl;
 mshFile << surfMesh->phyNames.size() << endl;
 for( int nb=0;nb<(int)surfMesh->phyNames.size();nb++ )
   mshFile << "2 " << nb+1 << " " << surfMesh->phyNames.at(nb) << endl;
 mshFile << "$EndPhysicalNames" << endl;
 mshFile << "$Nodes" << endl;
 mshFile << surfMesh->numVerts << endl;

 mshFile << setprecision(10) << scientific;
 for( int i=0;i<surfMesh->numVerts;i++ )
  mshFile << i+1 << " " << X->Get(i) << " " 
                        << Y->Get(i) << " " 
						<< Z->Get(i) << endl;

 mshFile << "$EndNodes" << endl;
 mshFile << "$Elements" << endl;
 mshFile << surfMesh->numElems << endl;

 mshFile << setprecision(0) << fixed;
 for( int i=0;i<surfMesh->numElems;i++ )
 {
  int v1 = surfMesh->IEN.Get(i,0);
  int v2 = surfMesh->IEN.Get(i,1);
  int v3 = surfMesh->IEN.Get(i,2);

  mshFile << i+1 
          << " 2" 
		  << " 2" 
		  << " " << surfMesh->idRegion.Get(i) + 1 
		  << " " << surfMesh->idRegion.Get(i) + 1  // surface number Gmsh
		  << " " << v1+1 
		  << " " << v2+1 
		  << " " << v3+1 
		  << endl; 

 }
 mshFile << "$EndElements" << endl;

 mshFile.close();

 copyLastFile(_dir,filename,_filename);

 cout << "mesh No. " << _iter << " saved in MSH" << endl;

} // fecha metodo saveMSH

void InOut::saveBubbleInfo(const char* _dir)
{
 saveOscillatingError(_dir);    // oscillating velocity and diameter
 saveKappaErrorSphere(_dir);       // kappa sphere
 //saveKappaErrorCylinder(_dir);  // kappa cylinder
 //saveKappaErrorTorus(_dir);     // kappa torus
 savePressureError(_dir);       // pressure 
 saveVolumeError(_dir);         // bubble volume
 saveVolumeCorrection(_dir);    // volume correction
 saveTimeError(_dir);           // time step
 //saveParasiticCurrent(_dir);   // velocity
 saveFilmThickness(_dir);    // oscillating velocity and diameter
}

/*
 * Zone (1): surfMesh wall
 *           mesh3D outside bubbles
 * Zone (2): surfMesh bubble1
 *           mesh3D inside bubble1
 * Zone (3): surfMesh bubble2
 *           mesh3D inside bubble2
 * Zone (4): surfMesh bubble3, etc.
 *           mesh3D inside bubble3, etc.
 * */
void InOut::printMeshReport()
{
 intet = m->getINTET();
 maxVolume = m->getMaxVolume();
 minVolume = m->getMinVolume();
 idMaxVolume = m->getIdMaxVolume();
 idMinVolume = m->getIdMinVolume();
 maxArea = m->getMaxArea();
 minArea = m->getMinArea();
 idMaxArea = m->getIdMaxArea();
 idMinArea = m->getIdMinArea();
 maxLength = m->getMaxLength();
 minLength = m->getMinLength();
 numSurfElems = m->getNumSurfElems();
 numSurfVerts = m->getNumSurfVerts();

 time_t currentTime;

 time( &currentTime );

 // ATTRIBUTES  :none,underscore,blink,reverse,concealed
 // COLORS      :black,red,green,yellow,blue,magenta,cyan,white
 cout << endl;
 cout << endl;
 cout << "   |------------------------------ Mesh Report ----------------------------|" 
      << endl;
 cout << "   |                                                                       |" 
      << endl;
 cout << "      number of 3D points          (numVerts):        " 
      << numVerts << endl;
 cout << "      number of 3D nodes           (numNodes):        " 
      << numNodes << endl;
 cout << "      number of tetrahedrons        (numEles):        " 
      << numElems << endl;
 cout << "      number of surface points     (surfMesh):        " 
      << surfMesh->numVerts << endl;
 cout << "      number of surface triangles    (numTri):        " 
      << surfMesh->numElems << endl;
 cout << "      min tetrahedron edge size:                      " 
      << m->getMinEdge() << endl;

 cout << endl;
 for(int nb=0;nb<=elemIdRegion->Max();nb++ )
 {
  cout << "      surface (" << nb << ")" << endl;
  cout << "       |area (initArea):                              "
       << surfaceArea[nb] << " (" << initSurfaceArea[nb] << ")" << endl;
  cout << "       |volume (initVolume):                          "
       << surfaceVolume[nb] << " (" << initSurfaceVolume[nb] << ")" << endl;
  cout << "       |average element edge length:                  " 
       << averageTriLength[nb] << endl;
  cout << "       |desired tetrahedron volume:                   "   
       << averageTriLength[nb]*
	      averageTriLength[nb]*
		  averageTriLength[nb]*sqrt(2)/12 << endl;
  cout << "       |triangle edge size:                           "  
       << triEdge[nb] << endl;
  cout << "       |number of surface triangle:                   "  
       << numSurfElems[nb] << endl;
  cout << "       |min triangle length:                          " 
       << minLength[nb] << endl;
  cout << "       |max triangle length:                          " 
       << maxLength[nb] << endl;
  cout << "       |min triangle area:                            " 
       << minArea[nb] << " (" << idMinArea[nb] << ")" << endl;
  cout << "       |max triangle area:                            " 
       << maxArea[nb] << " (" << idMaxArea[nb] << ")" << endl;
  cout << "       |min tetrahedron volume:                       " 
       << minVolume[nb] << " (" << idMinVolume[nb] << ")" << endl;
  cout << "       |max tetrahedron volume:                       " 
       << maxVolume[nb] <<  " (" << idMaxVolume[nb] << ")" << endl;
  cout << "       |number of tets with 4 verts on surface:       " 
       << intet[nb] << endl;

  cout << "       |" << color(none,yellow,black) 
	   << "inserted" << resetColor() 
	   << " surface points by length:            " 
	   << isp[nb] << endl;
  cout << "       |" << color(none,red,black) 
	   << "removed" << resetColor() 
   	   << " surface points by lenght:             " 
 	   << rsp[nb] << endl;
  cout << "       |" << color(none,yellow,black) 
       << "inserted" << resetColor() 
 	   << " surface points by curvature:         " 
 	   << ispc[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
 	   << " surface points by neighbour check:    " 
 	   << rspn[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
 	   << " surface points by curvature:          " 
 	   << rspc[nb] << endl;
  cout << "       |" << color(none,green,black) 
       << "flipped" << resetColor() 
 	   << " operations at surface:                " 
 	   << flip[nb] << endl; 
  cout << "       |" << color(none,magenta,black) 
       << "smoothed" << resetColor() 
 	   << " points by curvature:                 " 
 	   << spc[nb] << endl; 
  cout << "       |" << color(none,magenta,black) 
       << "smoothed" << resetColor() 
 	   << " points by angle between planes:      " 
 	   << spp[nb] << endl; 
  cout << "       |" << color(none,cyan,black) 
       << "contracted" << resetColor() 
 	   << " surface points by lenght:          " 
 	   << csp[nb] << endl; 
  cout << "       |" << color(none,yellow,black) 
       << "total number of " << color(none,yellow,black)
       << "inserted" << resetColor()
 	   << " surface mesh points: " << color(none,yellow,black)
 	   << isp[nb]+ispc[nb] << resetColor() << endl;
  cout << "       |" << color(none,red,black) 
       << "total number of " << color(none,red,black)
       << "removed" << resetColor() 
       << " surface mesh points:  " << color(none,red,black)
	   << rsp[nb]+rspn[nb]+rspc[nb]+csp[nb] << resetColor() << endl;
  cout << "       |" << color(none,yellow,black) 
       << "inserted" << resetColor()
 	   << " 3D mesh points by diffusion:         " << ipd[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
       << " 3D mesh points by diffusion:          " << rpd[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
       << " 3D mesh points by distance:           " << rpdist[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
       << " 3D mesh points by height:             " << rph[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
       << " 3D mesh points by volume:             " << rpv[nb] << endl;
  cout << "       |" << color(none,red,black) 
       << "removed" << resetColor() 
       << " 3D mesh points by interface distance: " << rpi[nb] << endl;
  cout << "       |" << color(none,yellow,black) 
       << "total number of " << color(none,yellow,black)
       << "inserted" << resetColor()
 	   << " 3D mesh points:      " << color(none,yellow,black)
 	   << ip[nb]+ipd[nb] << resetColor() << endl;
  cout << "       |" << color(none,red,black) 
       << "total number of " << color(none,red,black)
       << "removed" << resetColor() 
       << " 3D mesh points:       " << color(none,red,black)
	   << rp[nb]+rpd[nb]+rph[nb]+rpi[nb]+rpv[nb]+rpdist[nb] 
	   << resetColor() << endl;
 }
 cout << endl;

 cout << "   |                                                                       |" 
      << endl;
 cout << "   |-----------------------------------------------------------------------|" 
      << endl;
 cout << endl;
}

void InOut::printSimulationReport()
{
//--------------------------------------------------
//  time_t currentTime;
//  time( &currentTime );
//  cout << color(none,magenta,black)
//       << "          start process" << resetColor()
// 	  << ":                       " 
// 	  << asctime( localtime( &currentTime ) ) << endl;
//-------------------------------------------------- 

 real rho_ratio,mu_ratio;
 if( rho_in >= rho_out )
 {
  rho_ratio = rho_in/rho_out;
  mu_ratio = mu_in/mu_out;
 }
 else
 {
  rho_ratio = rho_out/rho_in;
  mu_ratio = mu_out/mu_in;
 }

 cout << endl;
 cout << "   |-------------------------- Simulation Report --------------------------|" 
      << endl;
 cout << "   |                                                                       |" 
      << endl;
 cout << color(none,magenta,black)
      << "          linear system dimension " << resetColor()
	  << " UVW:                  " 
	  << 3*numNodes << " x " << 3*numNodes << endl;
 cout << color(none,magenta,black)
      << "          linear system dimension " << resetColor()
	  << " P:                    " 
	  << numVerts << " x " << numVerts << endl;
 cout << color(none,magenta,black)
      << "          linear system dimension " << resetColor()
	  << " C/T:                  " 
	  << numVerts << " x " << numVerts << endl;
 cout << endl;
 cout << color(none,magenta,black)
      << "          Reynolds/Archimedes" << resetColor()
	  << " number:                    " << Re << endl;
 cout << color(none,magenta,black)
      << "          Froud" << resetColor()
	  << " number:                                  " 
	  << Fr << endl;
 cout << color(none,magenta,black)
      << "          Schmidt" << resetColor()
	  << " number:                                " << Sc << endl;
 cout << color(none,magenta,black)
      << "          Webber/Eotvos" << resetColor()
	  << " number:                          " << We << endl;
 cout << color(none,magenta,black)
      << "          alpha (time method)" << resetColor()
	  << " number:                    " << alpha << endl;
 cout << "          parameter " << color(none,magenta,black)
      << "c1 (Lagrangian velocity)" 
	  << resetColor() << ":            " << c1 << endl;
 cout << "          parameter " << color(none,magenta,black)
      << "c2 (smooth by velocity)" 
	  << resetColor() << ":             " << c2 << endl;
 cout << "          parameter " << color(none,magenta,black)
      << "c3 (smooth by coordinates)" 
	  << resetColor() << ":          " << c3 << endl;
 cout << "          parameter " << color(none,magenta,black)
      << "d1 (remove tangent velocity)" 
	  << resetColor() << ":        " << d1 << endl;
 cout << "          parameter " << color(none,magenta,black)
      << "d2 (smooth surface by coordinates)" 
	  << resetColor() << ":  " << d2 << endl;
 cout << color(none,magenta,black)
      << "          liquid viscosity" << resetColor() 
	  << ":                              " 
	  << mu_out << endl;
 cout << color(none,magenta,black)
      << "          gas viscosity" << resetColor() 
	  << ":                                 " 
	  << mu_in << endl;
 cout << color(none,magenta,black)
      << "          liquid density" << resetColor() 
	  << ":                                " 
	  << rho_out << endl;
 cout << color(none,magenta,black)
      << "          gas density" << resetColor() 
	  << ":                                   " 
	  << rho_in << endl;
 cout << endl;
 cout << color(none,magenta,black)
      << "          density ratio" << resetColor() 
	  << ":                                 " 
	  << rho_ratio << endl;
 cout << color(none,magenta,black)
      << "          viscosity ratio" << resetColor() 
	  << ":                               " 
	  << mu_ratio << endl;
 cout << endl;
 cout << "          time step:"  << endl;
 cout << "              |" << color(none,magenta,black) 
      << "Lagrangian" << resetColor()
	  << ":                               " 
	  << dtLagrangian << endl;
 cout << "              |" << color(none,magenta,black) 
      << "Semi-Lagrangian" << resetColor()
	  << ":                          " 
	  << dtSemiLagrangian << endl;
 cout << "              |" << color(none,magenta,black) 
      << "Surface tension" << resetColor()
	  << ":                          " 
	  << dtSurfaceTension << endl;
 cout << "              |" << color(none,magenta,black) 
      << "Gravity" << resetColor()
	  << ":                                  " 
	  << dtGravity << endl;
 cout << endl;
 cout << color(none,magenta,black)
      << "          iteration" << resetColor() 
	  << ":                                     " 
	  << iter << endl;
 cout << color(none,magenta,black)
      << "          CFL" << resetColor()
	  << " number:                                    " 
	  << cfl << endl;
 cout << color(none,magenta,black)
      << "          dt" << resetColor() 
	  << ":                                            " 
	  << dt << endl;
 cout << color(none,magenta,black)
      << "          time" << resetColor() 
	  << ":                                          " 
	  << s->getTime() << endl;
 cout << "   |                                                                       |" 
      << endl;
 cout << "   |-----------------------------------------------------------------------|" 
      << endl;
 cout << endl;
}

/*
 * This method saves into the file the information about the volume
 * occupied by 2 phases. It uses a 2D structured cartesian grid, defined by
 * np1 and np2. Thus, the 3D mesh data is interpolated to the 2D one
 * using the meshInterp function. Depending on the _filename, one of the
 * 3 planes are used (XY, XZ and YZ). 
 *
 * */
void InOut::crossSectionalPlane( const char* _dir,const char* _filename, int _iter )
{
 // file output
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 real _pos = (Y->Max()+Y->Min())/2.0;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + ".dat";
 const char* filename = file.c_str();

 // xVert da malha nova
 int np1 = 80;
 int np2 = 80;
 int nTotal = np1*np2;
 clVector xVert(nTotal);
 clVector yVert(nTotal);
 clVector zVert(nTotal);

 real plane1_i,plane1_f,plane2_i,plane2_f;
 real dp1,dp2; // mesh space in 2 directions: plane1 and plane2
 if( (strcmp(_filename,"XY") == 0) ||
     (strcmp(_filename,"YX") == 0)  )
 {
  // structured mesh points generator
  real xi = X->Min();
  real xf = X->Max();
  real yi = Y->Min();
  real yf = Y->Max();
  real zi = _pos;
  real dx = (xf-xi)/(np1-1);
  real dy = (yf-yi)/(np2-1);
  int count = 0;
  for( int i=0;i<np1;i++ )
  {
   for( int j=0;j<np2;j++ )
   {
	real x = xi + i * dx ;
	xVert.Set(count,x);
	real y = yi + j * dy;
	yVert.Set(count,y);
	count++;
   }
  }
  zVert.SetAll(zi);
  plane1_i = xi;
  plane1_f = xf;
  plane2_i = yi;
  plane2_f = yf;
  dp1=dx;
  dp2=dy;
 }
 else if( (strcmp(_filename,"XZ") == 0) ||  
          (strcmp(_filename,"ZX") == 0) )
 {
  // structured mesh points generator
  real xi = X->Min();
  real xf = X->Max();
  real yi = _pos;
  real zi = Z->Min();
  real zf = Z->Max();
  real dx = (xf-xi)/(np1-1);
  real dz = (zf-zi)/(np2-1);
  int count = 0;
  for( int i=0;i<np1;i++ )
  {
   for( int j=0;j<np2;j++ )
   {
	real x = xi + i * dx ;
	xVert.Set(count,x);
	real z = zi + j * dz;
	zVert.Set(count,z);
	count++;
   }
  }
  yVert.SetAll(yi);
  plane1_i = xi;
  plane1_f = xf;
  plane2_i = zi;
  plane2_f = zf;
  dp1=dx;
  dp2=dz;
 }
 else if( (strcmp(_filename,"YZ") == 0) || 
          (strcmp(_filename,"ZY") == 0) )
 {
  // structured mesh points generator
  real xi = _pos;
  real yi = Y->Min();
  real yf = Y->Max();
  real zi = Z->Min();
  real zf = Z->Max();
  real dy = (yf-yi)/(np1-1);
  real dz = (zf-zi)/(np2-1);
  int count = 0;
  for( int i=0;i<np1;i++ )
  {
   for( int j=0;j<np2;j++ )
   {
	real y = yi + i * dy;
	yVert.Set(count,y);
	real z = zi + j * dz;
	zVert.Set(count,z);
	count++;
   }
  }
  xVert.SetAll(xi);
  plane1_i = yi;
  plane1_f = yf;
  plane2_i = zi;
  plane2_f = zf;
  dp1=dy;
  dp2=dz;
 }
 else
 {
  cout << endl;
  cout << "       * ************************************** *" << endl;
  cout << "       *  You should define a plane:            *" << endl;
  cout << "       *  Ex. XY,XZ,YZ,YX,YZ,ZX                 *" << endl;
  cout << "       *             Using plane XZ             *" << endl;
  cout << "       * ************************************** *" << endl;
  cout << endl;
  // structured mesh points generator
  real xi = X->Min();
  real xf = X->Max();
  real yi = _pos;
  real zi = Z->Min();
  real zf = Z->Max();
  real dx = (xf-xi)/(np1-1);
  real dz = (zf-zi)/(np2-1);
  int count = 0;
  for( int i=0;i<np1;i++ )
  {
   for( int j=0;j<np2;j++ )
   {
	real x = xi + i * dx ;
	xVert.Set(count,x);
	real z = zi + j * dz;
	zVert.Set(count,z);
	count++;
   }
  }
  yVert.SetAll(yi);
  plane1_i = xi;
  plane1_f = xf;
  plane2_i = zi;
  plane2_f = zf;
  dp1=dx;
  dp2=dz;
 }

 // linear interpolation on nTotal using the 2D generated grid.
 clMatrix interpLin = meshInterp(*m,xVert,yVert,zVert);
 clVector cLin(nTotal);
 cLin = interpLin*(*heaviside);

 // gas is the number of vertices inside the bubble/drop.
 int gas=0;
 for( int i=0;i<cLin.Dim();i++ )
 {
  if( cLin.Get(i) >= 0.5 )
   gas++;
 }

 // saving in DAT format
 ifstream testFile( filename );
 ofstream voidFile( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  voidFile << "#time" << setw(20) << "total c-s area" 
				      << setw(15) << "bubble c-s area"
					  << setw(15) << "fluid c-s area" 
					  << setw(15) << "void fraction" 
					  << setw(15) << "iteration" 
					  << endl;
 }

 real totalArea = (plane1_f-plane1_i)*(plane2_f-plane2_i);
 real bubbleArea = gas*dp1*dp2;
 real fluidArea = totalArea - bubbleArea;
 real voidFraction = bubbleArea/totalArea;
 voidFile << setprecision(10) << scientific;
 voidFile << setw(10) << s->getTime() << " " 
                      << totalArea << " " 
					  << bubbleArea << " " 
					  << fluidArea << " " 
					  << voidFraction << " " 
                      << setprecision(0) << fixed
                      << _iter << endl;

 voidFile.close();

 cout << "cross sectional void fraction No. " << _iter 
      << " saved in dat" << endl;
 
 // concatenando nomes para o nome do arquivo final
 file = (string) _dir + (string) _filename + "pressure-" + str + ".dat";
 const char* filenameP = file.c_str();
 
 // saving in DAT format
 ifstream testFileP( filenameP );
 ofstream pFile( filenameP,ios::app );
 if( testFileP )
 {
  testFileP.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  pFile << "#x-position" << setw(17) << "z-position" 
				         << setw(16) << "pressure"
				         << setw(16) << "iteration"
				         << endl;
 }

 // interpolacao linear em numVerts
 clVector pLin(nTotal);
 pLin = interpLin*(*pSol);

 pFile << setprecision(10) << scientific;
 for( int i=0;i<nTotal;i++ )
 {
  pFile << setw(10) << setprecision(10) << scientific 
        << xVert.Get(i) << " " 
        << zVert.Get(i) << " " 
		<< pLin.Get(i) << " "
		<< setprecision(0) << fixed
		<< _iter << endl;
 }

 pFile.close();

 cout << "cross sectional pressure No. " << _iter << " saved in dat" << endl;

} // fecha metodo crossSectionalPlane

void InOut::bubbleWallDistance( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 //string file = (string) _dir + (string) _filename + "-" + str + ".dat";
 string file  = (string) _dir + (string) _filename + "." + str;
 const char* filename = file.c_str();

 ofstream dist ( filename ); 

 dist << "#P1x" << setw(18)
      << "P1y"  << setw(20) 
      << "P1z"  << setw(20) 
	  << "zMin" << setw(21)
	  << "length" << setw(18)
	  << "kappa" << setw(22)
	  << "pressure" << endl; 

 real zMin = surfMesh->Z.Min();
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->Marker.Get(i) == 0.5 )
  {
   real P1x = X->Get(i);
   real P1y = Y->Get(i);
   real P1z = Z->Get(i);
   real length = fabs(P1z-zMin);

   dist << setprecision(10) << scientific; 
   dist << setw(11) << P1x << " " 
	    << setw(18) << P1y << " " 
	    << setw(18) << P1z << " " 
		<< setw(18) << zMin << " " 
	    << setw(18) << length << " " 
		<< setw(18) << kappa->Get(i) << " " 
		<< setw(18) << pSol->Get(i) << endl;
  }
 }

 dist.close();

 copyLastFile(_dir,filename,_filename);

 cout << "Drop-Wall distance No. " << _iter << " saved in dat" << endl;
} // fecha metodo bubbleWallDistance

void InOut::saveKappaErrorSphere(const char* _dir)
{
 // kappa
 string fileAux = (string) _dir + "kappa" + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename);
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file kappa.dat" << endl;
 }
 else
 {
  cout << "Creating file kappa.dat" << endl;
  file << "#time" << setw(29) << "kappa" 
		    	  << setw(18) << "analytic" 
		    	  << setw(18) << "error" 
		    	  << setw(18) << "errorRel" 
		    	  << setw(18) << "stand deviat" 
		    	  << setw(14) << "averag edge" 
				  << setw(14) << "edge/radius" 
				  << setw(14) << "area" 
				  << setw(14) << "volume" 
				  << setw(15) << "averag neigh" 
		    	  << setw(13) << "num points" 
		    	  << setw(6) << "iter" 
				  << endl;
 }
 
 real xMax = -1E-10; 
 real yMax = -1E-10; 
 real zMax = -1E-10; 
 real xMin = 1E10; 
 real yMin = 1E10; 
 real zMin = 1E10; 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == 1 )
  {
   if( surfMesh->X.Get(i) > xMax )
	xMax = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) > yMax )
	yMax = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) > zMax )
	zMax = surfMesh->Z.Get(i);
   if( surfMesh->X.Get(i) < xMin )
	xMin = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) < yMin )
	yMin = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) < zMin )
	zMin = surfMesh->Z.Get(i);
  }
 }

 // retorna o valor do maior diametro em X na interface 
 real diameterX = xMax - xMin; 
 // retorna o valor do maior diametro em Y na interface 
 real diameterY = yMax - yMin; 
 // retorna o valor do maior diametro em Z na interface 
 real diameterZ = zMax - zMin; 

 real radius = (diameterX+diameterY+diameterZ)/6.0;

 // sphere
 real kappaAnalytic = 2.0/radius; 

 real surfacePoints = surface->Dim();
 real sumKappa = 0;
 real sumKappaSquare = 0;
 for( int i=0;i<surfacePoints;i++ )
 {
  int node = surface->Get(i);
  sumKappa += surfMesh->curvature.Get(node);
 }
 real kappaAverage = sumKappa/surfacePoints;

 real sumKappaSD = 0;
 real sumKappaError = 0;
 real sumNeighbours = 0;
 int countK = 0;
 for( int i=0;i<surfacePoints;i++ )
 {
  int node = surface->Get(i);
  sumKappaError += (surfMesh->curvature.Get(node)-kappaAnalytic)*
                   (surfMesh->curvature.Get(node)-kappaAnalytic);
  sumKappaSD += (surfMesh->curvature.Get(node)-kappaAverage)*
                (surfMesh->curvature.Get(node)-kappaAverage);
  sumNeighbours += neighbourPoint->at(i).size();
  sumKappaSquare += surfMesh->curvature.Get(node)*
                    surfMesh->curvature.Get(node);
  countK++;
 }

 /*  
  *            (  sum( kappa[i] - kappa_a )^2    )
  *  k_e = sqrt( -----------------------------   )
  *            (        sum( kappa[i]^2 )        )
  * */
 real kappaError = sqrt( sumKappaError/sumKappaSquare );

 /*
  *         kappaAverage-kappaAnalytic
  *  k_r = ----------------------------
  *               kappaAverage
  * */
 real kappaErrorRel = (kappaAverage-kappaAnalytic)/kappaAverage;

 /* 
  *             (  sum( kappa[i] - kappaAverage )  )
  *  SD =  sqrt ( -------------------------------- )
  *             (           number of i            )
  * */
 real kappaSD = sqrt( sumKappaSD/surfacePoints );

 real averageNeigh = sumNeighbours/surfacePoints;

 averageTriLength = m->getAverageTriLength();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
      << setw(17) << kappaAverage << " " 
      << setw(17) << kappaAnalytic << " " 
      << setw(17) << kappaError << " " 
      << setw(17) << kappaErrorRel << " " 
      << setw(17) << kappaSD << " " 
	  << setprecision(3) << fixed
      << setw(13) << averageTriLength[1] << " " 
	  << setprecision(4) << fixed
      << setw(13) << averageTriLength[1]/radius << " " 
	  << setprecision(3) << fixed
      << setw(14) << surfaceArea[1] << " " 
      << setw(14) << surfaceVolume[1] << " " 
      << setw(14) << averageNeigh << " " 
	  << setprecision(0) << fixed
	  << setw(12) << surfacePoints << " " 
	  << setw(5) << iter << endl;
 file.close();
}

void InOut::saveKappaErrorCylinder(const char* _dir)
{
 // kappa
 string fileAux = (string) _dir + "kappa" + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename);
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file kappa.dat" << endl;
 }
 else
 {
  cout << "Creating file kappa.dat" << endl;
  file << "#time" << setw(29) << "kappa" 
		    	  << setw(18) << "analytic" 
		    	  << setw(18) << "error" 
		    	  << setw(18) << "errorRel" 
		    	  << setw(18) << "stand deviat" 
		    	  << setw(14) << "averag edge" 
				  << setw(14) << "edge/radius" 
				  << setw(14) << "area" 
				  << setw(14) << "volume" 
				  << setw(15) << "averag neigh" 
		    	  << setw(13) << "num points" 
		    	  << setw(6) << "iter" 
				  << endl;
 }

 /* ------------ Radius calculation ------------ */
 real xMax = -1E-10; 
 real zMax = -1E-10; 
 real xMin = 1E10; 
 real zMin = 1E10; 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == 1 )
  {
   if( surfMesh->X.Get(i) > xMax )
	xMax = surfMesh->X.Get(i);
   if( surfMesh->Z.Get(i) > zMax )
	zMax = surfMesh->Z.Get(i);
   if( surfMesh->X.Get(i) < xMin )
	xMin = surfMesh->X.Get(i);
   if( surfMesh->Z.Get(i) < zMin )
	zMin = surfMesh->Z.Get(i);
  }
 }
 real diameterX = xMax - xMin; 
 real diameterZ = zMax - zMin; 
 real radius = (diameterX+diameterZ)/4.0;
 /* -------------------------------------------- */

 /* -------------- Kappa Analytic -------------- */
 real kappaAnalytic = 1.0/radius;
 /* -------------------------------------------- */

 int countK = 0;
 real sumKappa = 0;
 real sumKappaSquare = 0;
 for( int i=0;i<surface->Dim();i++ )
 {
  int node = surface->Get(i);
  
  // k=1/R
  if( (Y->Get(node)<0.75) && 
	  (Y->Get(node)>-0.75) )
  {
   sumKappa += surfMesh->curvature.Get(node);
   sumKappaSquare += surfMesh->curvature.Get(node)*
	                 surfMesh->curvature.Get(node);
   countK++;
  }
 }
 real kappaAverage = sumKappa/countK;
 

 real sumKappaSD = 0;
 real sumKappaError = 0;
 real sumNeighbours = 0;
 countK = 0;
 for( int i=0;i<surface->Dim();i++ )
 {
  int node = surface->Get(i);
  
  // k=1/R
  if( (Y->Get(node)<0.75) && 
	  (Y->Get(node)>-0.75) )
  {
   sumKappaError += (surfMesh->curvature.Get(node)-kappaAnalytic)*
                    (surfMesh->curvature.Get(node)-kappaAnalytic);
   sumKappaSD += (surfMesh->curvature.Get(node)-kappaAverage)*
                 (surfMesh->curvature.Get(node)-kappaAverage);
   sumNeighbours += neighbourPoint->at(i).size();
   countK++;
  }
 }

 real kappaError = sqrt( sumKappaError/sumKappaSquare );

 /*
  *         kappaAverage-kappaAnalytic
  *  k_r = ----------------------------
  *               kappaAverage
  * */
 real kappaErrorRel = (kappaAverage-kappaAnalytic)/kappaAverage;

 real kappaSD = sqrt( sumKappaSD/countK );

 real averageNeigh = sumNeighbours/countK;

 averageTriLength = m->getAverageTriLength();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
      << setw(17) << kappaAverage << " " 
      << setw(17) << kappaAnalytic << " " 
      << setw(17) << kappaError << " " 
      << setw(17) << kappaErrorRel << " " 
      << setw(17) << kappaSD << " " 
	  << setprecision(3) << fixed
      << setw(13) << averageTriLength[1] << " " 
	  << setprecision(4) << fixed
      << setw(13) << averageTriLength[1]/radius << " " 
	  << setprecision(3) << fixed
      << setw(14) << surfaceArea[1] << " " 
      << setw(14) << surfaceVolume[1] << " " 
      << setw(14) << averageNeigh << " " 
	  << setprecision(0) << fixed
	  << setw(12) << countK << " " 
	  << setw(5) << iter << endl;
 file.close();
}

void InOut::saveKappaErrorTorus(const char* _dir)
{
 // kappa
 string fileAux = (string) _dir + "kappa" + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename);
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file kappa.dat" << endl;
 }
 else
 {
  cout << "Creating file kappa.dat" << endl;
  file << "#time" << setw(29) << "kappa" 
		    	  << setw(18) << "analytic" 
		    	  << setw(18) << "error" 
		    	  << setw(18) << "stand deviat" 
		    	  << setw(14) << "averag edge" 
				  << setw(14) << "edge/radius" 
				  << setw(14) << "area" 
				  << setw(14) << "volume" 
				  << setw(15) << "averag neigh" 
		    	  << setw(13) << "num points" 
		    	  << setw(6) << "iter" 
				  << endl;
 }

 real surfacePoints = surface->Dim();
 real sumKappa = 0;
 real sumKappaSquare = 0;
 for( int i=0;i<surfacePoints;i++ )
 {
  int node = surface->Get(i);
  sumKappa += surfMesh->curvature.Get(node);
  sumKappaSquare += surfMesh->curvature.Get(node)*
                    surfMesh->curvature.Get(node);
 }
 real kappaAverage = sumKappa/surfacePoints;


/*  TORUS
 *
 *
 *
 *
 *
 *                           @@@@@@@                            
 *                  @@@@$##*!!!===!!!*##$@@@@                   
 *             @@@$$#*!=;;:~~-~~-~~~:~:;==!*#$$@@@              
 *          $@$$$#*!=;:~~~-~-----~~~~~~:::;=!!*#$$$@$           
 *        $$$$$##*=;~~~-----~~~~:::::::::~::;=!**#$$$$$         
 *      *$$$$$#*!=:~~------~~~~~::::::::::::::;!**#$$$$$#       
 *     *#$$$$$#*=:~--,,,--             ~:::::::=!*##$$$$#*     -----
 *    !*#$$$$$#*=~-,,,                     ~:~::!*#$$$$$#*!       |
 *   ;!*#$$$$$#*=,.                           ~~!*#$$$$$#*!;      |
 *  -;!*#$$$$$$#!,              x              ~*#$$$$$$#*!;,     |
 *  ~:;==**#$$$@@@$          (center)         #@@@$$$#**!=;~,    D1
 *  ~~~:;=!*###$$@@@@                       @@@@$$###*!=;:~-,     |
 *   ~~~::;=!***#$$$@@@@@@             @@@@@@$$$#***!=;:~~--      |
 *   :~~~~~::;=!!**###$$$$@@@@@@@@@@@@@$$$$###**!!=;;:~~---,      |
 *    ::~~~~~~~:;;==!!****#############****!!==;;;:~~-----,    -----     
 *      ::~~~~-~~~~::;;=====!!!!!!!!!=====;;:::~---------      
 *       :::~~~~~~~-----~~:::::::::::::~~~-------------,        
 *         :::::~~~~~-~------------------------------,          
 *            :::::~~~~~~~~-----------------------,             
 *               :::::::~~~~~~~~~~~~~-----------                
 *                    :::::::~~~~~~~~~~~~~-   
 *
 *
 *   |------------------------- D2 -------------------------|
 *
 *
 * */

 
 real xMax = -1E-10; 
 real yMax = -1E-10; 
 real zMax = -1E-10; 
 real xMin = 1E10; 
 real yMin = 1E10; 
 real zMin = 1E10; 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == 1 )
  {
   if( surfMesh->X.Get(i) > xMax )
	xMax = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) > yMax )
	yMax = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) > zMax )
	zMax = surfMesh->Z.Get(i);
   if( surfMesh->X.Get(i) < xMin )
	xMin = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) < yMin )
	yMin = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) < zMin )
	zMin = surfMesh->Z.Get(i);
  }
 }

 // radius 1 = D1/2
 real D1 = zMax - zMin; 
 real radius1 = D1/2.0;
 //real xCenter1 = (xMax+xMin)/2.0;
 //real yCenter1 = (yMax+yMin)/2.0; 
 //real zCenter1 = (zMax+zMin)/2.0; 
 
 // radius 2 = D2/2
 real xCenter2 = (xMax+xMin)/2.0;
 real yCenter2 = (yMax+yMin)/2.0; 

 real sumKappaSD = 0;
 real sumKappaError = 0;
 real sumNeighbours = 0;
 int countK = 0;
 for( int i=0;i<surfacePoints;i++ )
 {
  int node = surface->Get(i);
  real kappaAnalytic = (1.0/radius1) + 
                        distance(xCenter2,yCenter2,
	                             surfMesh->X.Get(node),surfMesh->Y.Get(node));
  sumKappaError += (surfMesh->curvature.Get(node)-kappaAnalytic)*
                   (surfMesh->curvature.Get(node)-kappaAnalytic);
  sumKappaSD += (surfMesh->curvature.Get(node)-kappaAverage)*
                (surfMesh->curvature.Get(node)-kappaAverage);
  sumNeighbours += neighbourPoint->at(i).size();
  countK++;
 }

 real kappaError = sqrt( sumKappaError/sumKappaSquare );
 real kappaSD = sqrt( sumKappaSD/surfacePoints );

 real averageNeigh = sumNeighbours/surfacePoints;

 averageTriLength = m->getAverageTriLength();

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
      << setw(17) << kappaAverage << " " 
      << setw(17) << "spec for each point" << " " 
      << setw(17) << kappaError << " " 
      << setw(17) << kappaSD << " " 
	  << setprecision(3) << fixed
      << setw(13) << averageTriLength[1] << " " 
	  << setprecision(4) << fixed
      << setw(13) << averageTriLength[1]/radius1 << " " 
	  << setprecision(3) << fixed
      << setw(14) << surfaceArea[1] << " " 
      << setw(14) << surfaceVolume[1] << " " 
      << setw(14) << averageNeigh << " " 
	  << setprecision(0) << fixed
	  << setw(12) << surfacePoints << " " 
	  << setw(5) << iter << endl;
 file.close();
}

void InOut::savePressureError(const char* _dir)
{
 string fileAux = (string) _dir + "pressure" + ".dat";
 const char* filenamePressure = fileAux.c_str();

 ifstream testFilePressure( filenamePressure);
 ofstream filePressure( filenamePressure,ios::app );
 if( testFilePressure )
 {
  testFilePressure.close();
  cout << "appending on file pressure.dat" << endl;
 }
 else
 {
  cout << "Creating file pressure.dat" << endl;
  filePressure << "#time" << setw(29) << "pressure-in" 
			       	      << setw(18) << "pressure-out" 
			       	      << setw(18) << "dpressure" 
			       	      << setw(18) << "analytic" 
						  << setw(14) << "area" 
						  << setw(14) << "volume" 
						  << setw(13) << "num points" 
			    	      << setw(6) << "iter" 
					      << endl;
 }

 real xMax = -1E-10; 
 real yMax = -1E-10; 
 real zMax = -1E-10; 
 real xMin = 1E10; 
 real yMin = 1E10; 
 real zMin = 1E10; 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == 1 )
  {
   if( surfMesh->X.Get(i) > xMax )
	xMax = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) > yMax )
	yMax = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) > zMax )
	zMax = surfMesh->Z.Get(i);
   if( surfMesh->X.Get(i) < xMin )
	xMin = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) < yMin )
	yMin = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) < zMin )
	zMin = surfMesh->Z.Get(i);
  }
 }

 // retorna o valor do maior diametro em X na interface 
 real diameterX = xMax - xMin; 
 // retorna o valor do maior diametro em Y na interface 
 real diameterY = yMax - yMin; 
 // retorna o valor do maior diametro em Z na interface 
 real diameterZ = zMax - zMin; 

 int surfacePoints = surface->Dim();

 real radius = (diameterX+diameterY+diameterZ)/6.0;

 real pressureAnalytic = 2*sigma/radius; 

 real sumPressureIn = 0;
 //real sumPressureSquareIn = 0;
 int countIn = 0;
 real sumPressureOut = 0;
 //real sumPressureSquareOut = 0;
 int countOut = 0;

 for( int i=0;i<numVerts;i++ )
 {
  if( heaviside->Get(i) > 0.5 )
  {
   sumPressureIn += pSol->Get(i);
   //sumPressureSquareIn += pSol->Get(i)*pSol->Get(i);
   countIn++;
  }
  if( heaviside->Get(i) < 0.5 )
  {
   sumPressureOut += pSol->Get(i);
   //sumPressureSquareOut += pSol->Get(i)*pSol->Get(i);
   countOut++;
  }
 }

 real averagePIn = sumPressureIn/countIn;
 real averagePOut = sumPressureOut/countOut;

 real dp = averagePIn-averagePOut;

 filePressure << setprecision(10) << scientific; 
 filePressure << setw(10) << s->getTime() << " " 
              << setw(17) << averagePIn << " " 
              << setw(17) << averagePOut << " " 
              << setw(17) << dp << " " 
              << setw(17) << pressureAnalytic << " " 
              << setw(14) << surfaceArea[1] << " " 
              << setw(14) << surfaceVolume[1] << " " 
			  << setprecision(0) << fixed
			  << setw(12) << surfacePoints << " " 
			  << setw(5) << iter << endl;
 filePressure.close();
}

void InOut::saveVolumeError(const char* _dir)
{
 for(int nb=0;nb<=elemIdRegion->Max();nb++ )
 {
  stringstream ss;  //convertendo int --> string
  string str;
  ss << nb;
  ss >> str;

  string fileAux = (string) _dir + "volume" + str + ".dat";
  const char* filename = fileAux.c_str();
 
  ifstream testFile( filename );
  ofstream file( filename,ios::app );
  if( testFile )
  {
   testFile.close();
   cout << "appending on file volume" << nb << ".dat" << endl;
  }
  else
  {
   cout << "Creating file volume" << nb << ".dat" << endl;
   file << "#time" << setw(29) << "volume" 
                   << setw(18) << "vel centroid X" 
                   << setw(18) << "vel centroid Y" 
                   << setw(18) << "vel centroid Z" 
                   << setw(18) << "vel refenrece X" 
                   << setw(18) << "vel refenrece Y" 
                   << setw(18) << "vel refenrece Z" 
                   << setw(18) << "centroid X" 
                   << setw(18) << "centroid Y" 
                   << setw(18) << "centroid Z" 
                   << setw(18) << "centroid ref X" 
                   << setw(18) << "centroid ref Y" 
                   << setw(18) << "centroid ref Z" 
              	   << setw(6)  << "iter" 
     			   << endl;
  }
 
  file << setprecision(10) << scientific; 
  file << setw(10) << s->getTime() << " " 
       << setw(17) << m->getSurfaceVolume()[nb] << " " 
       << setw(17) << s->getCentroidVelX()[nb] << " " 
       << setw(17) << s->getCentroidVelY()[nb] << " " 
       << setw(17) << s->getCentroidVelZ()[nb] << " " 
       << setw(17) << s->getURef() << " " 
       << setw(17) << s->getVRef() << " " 
       << setw(17) << s->getWRef() << " " 
       << setw(17) << s->getCentroidPosX()[nb] << " " 
       << setw(17) << s->getCentroidPosY()[nb] << " " 
       << setw(17) << s->getCentroidPosZ()[nb] << " " 
       << setw(17) << s->getXRef() << " " 
       << setw(17) << s->getYRef() << " " 
       << setw(17) << s->getZRef() << " " 
       << setw(5) << setprecision(0) << fixed << iter 
       << endl;
  file.close();
 }
}

void InOut::saveVolumeCorrection(const char* _dir)
{
 for(int nb=0;nb<=elemIdRegion->Max();nb++ )
 {
  stringstream ss;  //convertendo int --> string
  string str;
  ss << nb;
  ss >> str;

  string fileAux = (string) _dir + "correctedVolume" + str + ".dat";
  const char* filename = fileAux.c_str();
 
  ifstream testFile( filename );
  ofstream file( filename,ios::app );
  if( testFile )
  {
   testFile.close();
   cout << "appending on file correctedVolume" << nb << ".dat" << endl;
  }
  else
  {
   cout << "Creating file volume" << nb << ".dat" << endl;
   file << "#time" << setw(29) << "volume" 
                   << setw(18) << "errorVolume" 
                   << setw(18) << "area" 
                   << setw(18) << "errorArea" 
              	   << setw(6)  << "iter" 
     			   << endl;
  }
 
  file << setprecision(10) << scientific; 
  file << setw(10) << s->getTime() << " " 
       << setw(17) << m->getSurfaceVolume()[nb] << " " 
       << setw(17) << m->getErrorVolume()[nb] << " " 
       << setw(17) << m->getSurfaceArea()[nb] << " " 
       << setw(17) << m->getErrorArea()[nb] << " " 
       << setw(5) << setprecision(0) << fixed << iter 
       << endl;
  file.close();
 }
}

void InOut::saveOscillatingError(const char* _dir)
{
 // oscillating velocity
 string fileAux = (string) _dir + "velocity" + ".dat";
 const char* filenameVel = fileAux.c_str();

 ifstream testFileVel( filenameVel );
 ofstream fileVel( filenameVel,ios::app );
 if( testFileVel )
 {
  testFileVel.close();
  cout << "appending on file velocity.dat" << endl;
 }
 else
 {
  cout << "Creating file velocity.dat" << endl;
  fileVel << "#time" << setw(29) << "vel Umax xMax" 
                     << setw(18) << "vel Umax xMin" 
                     << setw(18) << "vel Vmax vMax" 
                     << setw(18) << "vel Vmax vMin" 
                     << setw(18) << "vel Wmax zMax" 
                     << setw(18) << "vel Wmax zMin" 
					 << setw(6) << "iter" 
					 << endl;
 }
 
 real xMax,uMax = -1E-10;
 real yMax,vMax = -1E-10; 
 real zMax,wMax = -1E-10; 
 real xMin,uMin = 1E10; 
 real yMin,vMin = 1E10; 
 real zMin,wMin = 1E10; 

 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == 1 )
  {
   if( surfMesh->X.Get(i) > xMax )
   {
	xMax = surfMesh->X.Get(i);
	uMax = uSol->Get(xMax);
   }
   if( surfMesh->Y.Get(i) > yMax )
   {
	yMax = surfMesh->Y.Get(i);
	vMax = vSol->Get(yMax);
   }
   if( surfMesh->Z.Get(i) > zMax )
   {
	zMax = surfMesh->Z.Get(i);
	wMax = wSol->Get(zMax);
   }
   if( surfMesh->X.Get(i) < xMin )
   {
	xMin = surfMesh->X.Get(i);
	uMin = uSol->Get(xMin);
   }
   if( surfMesh->Y.Get(i) < yMin )
   {
	yMin = surfMesh->Y.Get(i);
	vMin = vSol->Get(yMin);
   }
   if( surfMesh->Z.Get(i) < zMin )
   {
	zMin = surfMesh->Z.Get(i);
	wMin = wSol->Get(zMin);
   }
  }
 }

 fileVel << setprecision(10) << scientific; 
 fileVel << setw(10) << s->getTime() << " " 
         << setw(17) << uMax << " " 
         << setw(17) << uMin << " " 
		 << setw(17) << vMax << " " 
		 << setw(17) << vMin << " " 
		 << setw(17) << wMax << " " 
		 << setw(17) << wMin << " " 
		 << setw(5) << setprecision(0) << fixed << iter 
		 << endl;

 fileVel.close();

 // diameter
 fileAux = (string) _dir + "diameter" + ".dat";
 const char* filenameD = fileAux.c_str();

 ifstream testFileD( filenameD );
 ofstream fileD( filenameD,ios::app );
 if( testFileD )
 {
  testFileD.close();
  cout << "appending on file diameter.dat" << endl;
 }
 else
 {
  cout << "Creating file diameter.dat" << endl;
  fileD << "#time" << setw(29) << "diameter X" 
                   << setw(18) << "diameter Y" 
				   << setw(18) << "diameter Z"
				   << setw(18) << "X max"
				   << setw(18) << "X min"
				   << setw(18) << "X mid"
				   << setw(18) << "Y max"
				   << setw(18) << "Y min"
				   << setw(18) << "Y mid"
				   << setw(18) << "Z max"
				   << setw(18) << "Z min"
				   << setw(18) << "Z mid"
				   << setw(6) << "iter" 
				   << endl;
 }

 // retorna o valor do maior diametro em X na interface 
 real diameterX = xMax - xMin; 
 // retorna o valor do maior diametro em Y na interface 
 real diameterY = yMax - yMin; 
 // retorna o valor do maior diametro em Z na interface 
 real diameterZ = zMax - zMin; 

 fileD << setprecision(10) << scientific; 
 fileD << setw(10) << s->getTime() << " " 
       << setw(17) << diameterX << " " 
	   << setw(17) << diameterY << " " 
	   << setw(17) << diameterZ << " " 
	   << setw(17) << xMin << " " 
	   << setw(17) << xMax << " " 
	   << setw(17) << (xMax-xMin)*0.5 << " " 
	   << setw(17) << yMin << " " 
	   << setw(17) << yMax << " " 
	   << setw(17) << (yMax-yMin)*0.5 << " " 
	   << setw(17) << zMin << " " 
	   << setw(17) << zMax << " " 
	   << setw(17) << (zMax-zMin)*0.5 << " " 
	   << setw(5) << setprecision(0) << fixed << iter 
	   << endl;
 fileD.close();
}

void InOut::saveTimeError(const char* _dir)
{
 string fileAux = (string) _dir + "time" + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file time.dat" << endl;
 }
 else
 {
  cout << "Creating file time.dat" << endl;
  file << "#time" << setw(29) << "lagrangian" 
                  << setw(18) << "semi-lagrangian" 
                  << setw(18) << "gravity"
                  << setw(18) << "surface tension"
                  << setw(18) << "dt"
                  << setw(6) << "cfl" 
		   		  << setw(7)  << "iter" 
		  		  << endl;
 }

 file << setprecision(10) << scientific; 
 file << setw(16) << s->getTime() << " " 
          << setw(17) << s->getDtLagrangian() << " " 
          << setw(17) << s->getDtSemiLagrangian() << " " 
          << setw(17) << s->getDtGravity() << " " 
          << setw(17) << s->getDtSurfaceTension() << " " 
          << setw(17) << s->getDt() << " " 
          << setw(5) << setprecision(2) << fixed << s->getCfl() << " " 
	      << setw(6) << setprecision(0) << fixed << iter 
	 	  << endl;
 file.close();
}

void InOut::saveParasiticCurrent(const char* _dir)
{
 string fileAux = (string) _dir + "parasitic" + ".dat";
 const char* filename = fileAux.c_str();

 ifstream testFile( filename );
 ofstream file( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file parasitic.dat" << endl;
 }
 else
 {
  cout << "Creating file parasitic.dat" << endl;
  file << "#time" << setw(29) << "vel X" 
                  << setw(18) << "vel Y" 
                  << setw(18) << "vel Z" 
                  << setw(18) << "vel XY" 
                  << setw(18) << "vel XZ" 
                  << setw(18) << "vel YZ" 
                  << setw(18) << "vel XYZ" 
	         	  << setw(6)  << "iter" 
				  << endl;
 }

 real uMax = uSol->Abs().Max();
 real vMax = vSol->Abs().Max();
 real wMax = wSol->Abs().Max();
 real uvMax = max(uMax,vMax);
 real uwMax = max(uMax,wMax);
 real vwMax = max(vMax,wMax);
 real uvwMax = max(uMax,vMax);
 uvwMax = max(uvwMax,wMax);

 file << setprecision(10) << scientific; 
 file << setw(10) << s->getTime() << " " 
      << setw(17) << uMax << " " 
	  << setw(17) << vMax << " " 
	  << setw(17) << wMax << " " 
	  << setw(17) << uvMax << " " 
	  << setw(17) << uwMax << " " 
	  << setw(17) << vwMax << " " 
	  << setw(17) << uvwMax << " " 
	  << setw(5) << setprecision(0) << fixed << iter 
	  << endl;
 file.close();
}

/* 
 * copying to file filename-last.ext 
 * */
void InOut::copyLastFile(const char* _dir,
                         const char* _filename,
						 const char* _name)
{
 string aux = (string) _filename;
 string ext(aux.end()-3,aux.end()); // get file extension
 
 ifstream inFile( _filename,ios::binary ); 

 string last = (string) _dir + 
               (string) _name + "-last." + 
			   (string) ext;

 const char* filename = last.c_str();
 ofstream outFile( filename,ios::binary ); 

 outFile << inFile.rdbuf();
 inFile.close();
 outFile.close();
}

void InOut::setCutPlane(ofstream& _file)
{
 real plane1a = X->Min() +  ( X->Max()-X->Min() )/4.0;
 real plane1b = X->Min() +2*( X->Max()-X->Min() )/4.0;
 real plane1c = X->Min() +3*( X->Max()-X->Min() )/4.0;
 real plane2a = Y->Min() +  ( Y->Max()-Y->Min() )/4.0;
 real plane2b = Y->Min() +2*( Y->Max()-Y->Min() )/4.0;
 real plane2c = Y->Min() +3*( Y->Max()-Y->Min() )/4.0;
 real plane3a = Z->Min() +  ( Z->Max()-Z->Min() )/4.0;
 real plane3b = Z->Min() +2*( Z->Max()-Z->Min() )/4.0;
 real plane3c = Z->Min() +3*( Z->Max()-Z->Min() )/4.0;
 clVector cutPlaneX(numVerts);
 clVector cutPlaneY(numVerts);
 clVector cutPlaneZ(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  bool planeTestXa = (X->Get( i ) >  plane1a);
  bool planeTestXb = (X->Get( i ) >  plane1b);
  bool planeTestXc = (X->Get( i ) >  plane1c);

  bool planeTestYa = (Y->Get( i ) >  plane2a); 
  bool planeTestYb = (Y->Get( i ) >  plane2b); 
  bool planeTestYc = (Y->Get( i ) >  plane2c); 

  bool planeTestZa = (Z->Get( i ) >  plane3a); 
  bool planeTestZb = (Z->Get( i ) >  plane3b); 
  bool planeTestZc = (Z->Get( i ) >  plane3c); 

  // allows only the outter mesh for two-phase flows, 
  // however it works for single phase flows too

  real hTest;
  if( surfMesh->numInterfaces > 0 )
   hTest = heaviside->Get(i) > 0.0;
  else 
   hTest = false;

  if( planeTestXa || hTest )
   cutPlaneX.Set(i,1.0);
  if( planeTestXb || hTest )
   cutPlaneX.Set(i,2.0);
  if( planeTestXc || hTest )
   cutPlaneX.Set(i,3.0);

  if( planeTestYa || hTest )
   cutPlaneY.Set(i,1.0);
  if( planeTestYb || hTest )
   cutPlaneY.Set(i,2.0);
  if( planeTestYc || hTest )
   cutPlaneY.Set(i,3.0);

  if( planeTestZa || hTest )
   cutPlaneZ.Set(i,1.0);
  if( planeTestZb || hTest )
   cutPlaneZ.Set(i,2.0);
  if( planeTestZc || hTest )
   cutPlaneZ.Set(i,3.0);
 }

  vtkScalar(_file,"cutPlaneX",cutPlaneX);
  vtkScalar(_file,"cutPlaneY",cutPlaneY);
  vtkScalar(_file,"cutPlaneZ",cutPlaneZ);
}

void InOut::saveFilmThickness(const char* _dir)
{
 // max and min of domain
 real xdMin = surfMesh->X.Min();  
 real ydMin = surfMesh->Y.Min();
 real zdMin = surfMesh->Z.Min(); 
 real xdMax = surfMesh->X.Max(); 
 real ydMax = surfMesh->Y.Max(); 
 real zdMax = surfMesh->Z.Max(); 

 for(int nb=0;nb<=elemIdRegion->Max();nb++ )
 {
  stringstream ss;  //convertendo int --> string
  string str;
  ss << nb;
  ss >> str;

  string fileAux = (string) _dir + "film" + str + ".dat";
  const char* filename = fileAux.c_str();
 
  ifstream testFile( filename );
  ofstream file( filename,ios::app );
  if( testFile )
  {
   testFile.close();
   cout << "appending on file film" << nb << ".dat" << endl;
  }
  else
  {
   cout << "Creating file film" << nb << ".dat" << endl;
   file << "#time" << setw(29) << "minFilmX" 
                   << setw(18) << "minFilmY" 
                   << setw(18) << "minFilmZ" 
                   << setw(18) << "maxFilmX" 
                   << setw(18) << "maxFilmY" 
                   << setw(18) << "maxFilmZ" 
                   << setw(18) << "min domain X" 
                   << setw(18) << "min domain Y" 
                   << setw(18) << "min domain Z" 
                   << setw(18) << "max domain X" 
                   << setw(18) << "max domain Y" 
                   << setw(18) << "max domain Z" 
                   << setw(18) << "min drop X" 
                   << setw(18) << "min drop Y" 
                   << setw(18) << "min drop Z" 
              	   << setw(6)  << "iter" 
     			   << endl;
  }

 real xMax = -1E-10;
 real yMax = -1E-10; 
 real zMax = -1E-10; 
 real xMin = 1E10; 
 real yMin = 1E10; 
 real zMin = 1E10; 
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->vertIdRegion.Get(i) == nb )
  {
   if( surfMesh->X.Get(i) > xMax )
	xMax = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) > yMax )
	yMax = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) > zMax )
	zMax = surfMesh->Z.Get(i);
   if( surfMesh->X.Get(i) < xMin )
	xMin = surfMesh->X.Get(i);
   if( surfMesh->Y.Get(i) < yMin )
	yMin = surfMesh->Y.Get(i);
   if( surfMesh->Z.Get(i) < zMin )
	zMin = surfMesh->Z.Get(i);
  }
 }

 real minFilmX = fabs(xMin-xdMin);
 real minFilmY = fabs(yMin-ydMin);
 real minFilmZ = fabs(zMin-zdMin);
 real maxFilmX = fabs(xMax-xdMax);
 real maxFilmY = fabs(yMax-ydMax);
 real maxFilmZ = fabs(zMax-zdMax);

  file << setprecision(10) << scientific; 
  file << setw(10) << s->getTime() << " " 
       << setw(17) << minFilmX << " " 
       << setw(17) << minFilmY << " " 
       << setw(17) << minFilmZ << " " 
       << setw(17) << maxFilmX << " " 
       << setw(17) << maxFilmY << " " 
       << setw(17) << maxFilmZ << " " 
       << setw(17) << xdMin << " " 
       << setw(17) << ydMin << " " 
       << setw(17) << zdMin << " " 
       << setw(17) << xdMax << " " 
       << setw(17) << ydMax << " " 
       << setw(17) << zdMax << " " 
       << setw(17) << xMin << " " 
       << setw(17) << yMin << " " 
       << setw(17) << zMin << " " 
       << setw(17) << xMax << " " 
       << setw(17) << yMax << " " 
       << setw(17) << zMax << " " 
       << setw(5) << setprecision(0) << fixed << iter 
       << endl;
  file.close();
 }
}

void InOut::savePoint( const char* _dir,int _point )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _point;
 ss >> str;

 // oscillating velocity
 string fileAux = (string) _dir + "point-" + str + ".dat";
 const char* filenamePoint = fileAux.c_str();

 ifstream testFilePoint( filenamePoint );
 ofstream filePoint( filenamePoint,ios::app );
 if( testFilePoint )
 {
  testFilePoint.close();
  cout << "appending on file point-" << _point << ".dat" << endl;
 }
 else
 {
  cout << "Creating file point-" << _point << ".dat" << endl;
  filePoint << "#time" << setw(29) << "u" 
                       << setw(18) << "v" 
					   << setw(18) << "w" 
					   << setw(18) << "p"
					   << setw(18) << "rho"
					   << setw(18) << "mu"
					   << setw(18) << "x" 
					   << setw(18) << "y" 
					   << setw(18) << "z" 
					   << setw(6) << "iter" 
					   << endl;
 }
 
 filePoint << setprecision(10) << scientific; 
 filePoint << setw(10) << s->getTime() << " " 
           << setw(17) << uSol->Get(_point) << " " 
		   << setw(17) << vSol->Get(_point) << " " 
		   << setw(17) << wSol->Get(_point) << " " 
		   << setw(17) << pSol->Get(_point) << " " 
		   << setw(17) << rho->Get(_point) << " " 
		   << setw(17) << mu->Get(_point) << " " 
		   << setw(17) << X->Get(_point) << " " 
		   << setw(17) << Y->Get(_point) << " " 
		   << setw(17) << Z->Get(_point) << " " 
		   << setw(5) << setprecision(0) << fixed << iter 
		   << endl;

 filePoint.close();
}
