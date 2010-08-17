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
 numGLEP = m->getNumGLEP();
 numGLEU = m->getNumGLEU();
 numGLEC = m->getNumGLEC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 outflow = m->getOutflow();
 IEN = m->getIEN();
 surface = m->getSurface();
}

InOut::InOut( Model3D &_m, Simulator3D &_s )
{
 m = &_m;
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 numGLEP = m->getNumGLEP();
 numGLEU = m->getNumGLEU();
 numGLEC = m->getNumGLEC();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 uc = m->getUC();
 vc = m->getVC();
 wc = m->getWC();
 pc = m->getPC();
 cc = m->getCC();
 idbcu = m->getIdbcu();
 idbcv = m->getIdbcv();
 idbcw = m->getIdbcw();
 idbcp = m->getIdbcp();
 idbcc = m->getIdbcc();
 outflow = m->getOutflow();
 IEN = m->getIEN();
 surface = m->getSurface();

 s = &_s;
 Re = s->getRe();
 Sc = s->getSc();
 We = s->getWe();
 Fr = s->getFr();
 dt = s->getDt();
 cfl = s->getCfl();
 alpha = s->getAlpha();
 beta = s->getBeta();
 simTime = s->getTime();
 uAnt = s->getUAnt();
 cAnt = s->getCAnt();
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
 distance = s->getDistance();
 nu = s->getNu();
 rho = s->getRho();
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

 // copiando para arquivo sim-last.bin
 ifstream inFile( filenameUVWPC,ios::binary ); 

 string last = (string) _dir + "sim-last" + ".bin";
 const char* filenameCopy = last.c_str();
 ofstream outFile( filenameCopy,ios::binary ); 

 outFile << inFile.rdbuf();
 inFile.close();
 outFile.close();

 cout << "solution No.  " << _iter << " saved in binary" << endl;
 
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

 clVector aux2(3*numNodes+2*numVerts); // vetor tambem carrega a concentracao

 ifstream UVWPC_file( filenameUVWPC,ios::in | ios::binary ); 

 UVWPC_file.read( (char*) aux2.GetVec(),aux2.Dim()*sizeof(real) );

 UVWPC_file.close();

 //aux2.CopyTo(0,*uAnt);
 aux2.CopyTo(0,*uSol);
 aux2.CopyTo(numNodes,*vSol);
 aux2.CopyTo(2*numNodes,*wSol);
 aux2.CopyTo(3*numNodes,*pSol);
 aux2.CopyTo(3*numNodes+numVerts,*cSol);

//--------------------------------------------------
//  clVector solVel(3*numNodes+numVerts);
//  aux2.CopyTo(0,solVel);
//  clVector solC(numVerts);
//  aux2.CopyTo(3*numNodes+numVerts,solC);
// 
//  s->setUAnt(solVel); // impondo a velocidade no simulador
//  s->setCSol(solC); // impondo a concentracao no simulador
//-------------------------------------------------- 

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

 int i;
 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "-" + str + ".dat";
 const char* filenameUVWPC = fileUVWPC.c_str();

 ofstream UVWPC_file( filenameUVWPC ); 
 UVWPC_file << numVerts << " "  << numNodes << endl << endl;

 for( i=0;i<vec.Dim();i++ )
   UVWPC_file << vec.Get(i) << endl;

 UVWPC_file << endl;
 UVWPC_file.close();

 cout << "solucao no.  " << _iter << " gravada em TXT" << endl;
 
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
 simTime = s->getTime();

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 
 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "3D Simulation C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 vtkFile << setprecision(10) << scientific; 
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << X->Get(i) << " " 
          << setw(10) << Y->Get(i) << " " 
		  << setw(10) << Z->Get(i) << endl;

 vtkFile << endl;
 
 vtkFile << "CELLS " << numElems << " " << 5*numElems << endl;
 vtkFile << setprecision(0) << fixed; 
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

 vtkFile << "SCALARS pressure double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific; 
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << pc->Get(i) << endl;

 vtkFile << endl;

 // este if existe pois nem todos os metodos tem cc
 if( cc->Dim() > 0 )
 {
  vtkFile << "SCALARS concentration double" << endl;
  vtkFile << "LOOKUP_TABLE default"  << endl;

  for( int i=0;i<numVerts;i++ )
   vtkFile << setw(10) << cc->Get(i) << endl;

  vtkFile << endl;
 }

 vtkFile << "VECTORS vectors double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << uc->Get(i) << " " 
          << setw(10) << vc->Get(i) << " " 
		  << setw(10) << wc->Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "mesh saved in VTK" << endl;

} // fecha metodo saveVTK 

void InOut::saveVTK( const char* _dir,const char* _filename, int _iter )
{
 IEN = m->getIEN();
 numElems = m->getNumElems();
 simTime = s->getTime();

 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "3D Simulation C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << *simTime << endl;
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

 vtkFile << "SCALARS pressure double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << pSol->Get(i) << endl;

 vtkFile << endl;

 // este if existe pois nem todos os metodos tem cc
 if( cSol->Dim() > 0 )
 {
  vtkFile << "SCALARS concentration double" << endl;
  vtkFile << "LOOKUP_TABLE default"  << endl;

  for( int i=0;i<numVerts;i++ )
   vtkFile << cSol->Get(i) << endl;

  vtkFile << endl;
 }

 vtkFile << "SCALARS kappa double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << kappa->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS distance double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << distance->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uSol->Get(i) << " " 
          << vSol->Get(i) << " " 
		  << wSol->Get(i) << endl;
 vtkFile << endl;

 vtkFile << "VECTORS ALE-velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uALE->Get(i) << " " 
          << vALE->Get(i) << " " 
		  << wALE->Get(i) << endl;
 vtkFile << endl;

 vtkFile << "VECTORS surface-force double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << fint->Get(i) << " " 
          << fint->Get(i+numVerts) << " " 
		  << fint->Get(i+numVerts*2) << endl;

 vtkFile.close();

 // copiando para arquivo sim-last.bin
 ifstream inFile( filename,ios::binary ); 

 string last = (string) _dir + "sim-last" + ".vtk";
 const char* filenameCopy = last.c_str();
 ofstream outFile( filenameCopy,ios::binary ); 

 outFile << inFile.rdbuf();
 inFile.close();
 outFile.close();

 cout << "solution No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVtk

void InOut::saveVTKTest( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = (string) _dir + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Simulacao 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << *simTime << endl;
 vtkFile << "ITERATION 1 1 int" << endl;
 vtkFile << _iter << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << X->Get(i) << " " << Y->Get(i) << " " << Z->Get(i) << endl;

 vtkFile << endl;
 
 // conta numero de elementos
 real plane1 = ( X->Max()-X->Min() )/2.0;
 real plane2 = ( Y->Max()-Y->Min() )/2.0;
 int count = 0;
 for( int i=0;i<numElems;i++ )
 {
  if( (X->Get( IEN->Get(i,0) ) <  plane1) && 
	  (X->Get( IEN->Get(i,1) ) <  plane1) && 
	  (X->Get( IEN->Get(i,2) ) <  plane1) && 
	  (X->Get( IEN->Get(i,3) ) <  plane1) &&
      (Y->Get( IEN->Get(i,0) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,1) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,2) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,3) ) <  plane2) ) 
   count++;
 }
 
 vtkFile << "CELLS " << count << " " << 5*count << endl;
 vtkFile << setprecision(0) << fixed; 
 for( int i=0;i<numElems;i++ )
 {
  if( (X->Get( IEN->Get(i,0) ) <  plane1) && 
	  (X->Get( IEN->Get(i,1) ) <  plane1) && 
	  (X->Get( IEN->Get(i,2) ) <  plane1) && 
	  (X->Get( IEN->Get(i,3) ) <  plane1) &&
      (Y->Get( IEN->Get(i,0) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,1) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,2) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,3) ) <  plane2) ) 
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
  if( (X->Get( IEN->Get(i,0) ) <  plane1) && 
	  (X->Get( IEN->Get(i,1) ) <  plane1) && 
	  (X->Get( IEN->Get(i,2) ) <  plane1) && 
	  (X->Get( IEN->Get(i,3) ) <  plane1) &&
      (Y->Get( IEN->Get(i,0) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,1) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,2) ) <  plane2) && 
	  (Y->Get( IEN->Get(i,3) ) <  plane2) ) 
   vtkFile << "10 ";
 }

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;

 vtkFile << "SCALARS pressure double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << pSol->Get(i) << endl;

 vtkFile << endl;

 // este if existe pois nem todos os metodos tem cc
 if( cSol->Dim() > 0 )
 {
  vtkFile << "SCALARS concentration double" << endl;
  vtkFile << "LOOKUP_TABLE default"  << endl;

  for( int i=0;i<numVerts;i++ )
   vtkFile << cSol->Get(i) << endl;

  vtkFile << endl;
 }

 vtkFile << "SCALARS kappa double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << kappa->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS distance double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << distance->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uSol->Get(i) << " " 
          << vSol->Get(i) << " " 
		  << wSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS ALE-Velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uALE->Get(i) << " " 
          << vALE->Get(i) << " " 
		  << wALE->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS surface-force double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << fint->Get(i) << " " 
          << fint->Get(i+numNodes) << " " 
		  << fint->Get(i+numNodes*2) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "solution No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVtk

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

void InOut::saveVonKarman(const char* _dir,const char* _filename,int _iter,
                          int vertice )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = _dir;
 string aux = _filename;
 file += aux + "." + str;
 const char* filename = file.c_str();

 ofstream vonKarmanFile;
 vonKarmanFile.open( filename );

 clVector Xaux(0);
 clVector Yaux = *Y==0.0;
 clVector Y0 = Yaux.Find();
 real x1 = X->Get( (int) Y0.Get(vertice) );

 for( int i=0;i<X->Dim();i++ )
 {
  if( Y->Get(i) == 0.0 && X->Get(i) == x1 )
   Xaux.AddItem(i);
 }

 int v1 = (int) Xaux.Get(0);
 int v2;
 vonKarmanFile << setprecision(10) << scientific; 
 vonKarmanFile << "#z     F      G      H      c      p     nu" << endl; 
 for( int i=0;i<Xaux.Dim();i++ )
 {
  v2 = (int) Xaux.Get(i);
  vonKarmanFile << setw(9) <<  Z->Get( v2 ) << " " 
	            << setw(9) <<  uSol->Get( v2 )/X->Get( v1 ) << " " 
	            << setw(9) <<  vSol->Get( v2 )/X->Get( v1 ) << " " 
		        << setw(9) <<  (-1)*wSol->Get( v2 ) << " " 
				<< setw(9) <<  cSol->Get(v2) << " " 
				<< setw(9) <<  pSol->Get(v2) << " "
				<< setw(9) <<  nu->Get(v2) << endl;
 }
 vonKarmanFile << endl;
 vonKarmanFile << "X = " << X->Get( (int) Xaux.Get(0) ) << endl;

 vonKarmanFile.close();

 cout << "von Karman num. " << _iter << " saved in ASCII" << endl;
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
 vtkFile << *simTime << endl;
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
 vtkFile << *simTime << endl;
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
 vtkFile << *simTime << endl;
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

 file << _iter << " " << *simTime << endl;
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
 datFile << _iter << " " << *simTime << endl;
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
 file << "----------------------- BEGIN ----------------------" << endl; 
 file << "start process:   " << asctime( localtime( &currentTime ) ) << endl;
 file << "mesh:            " << _mesh << endl;
 file << "numVerts:        " << numVerts << endl;
 file << "numNodes:        " << numNodes << endl;
 file << "numElems:        " << numElems << endl;
 file << "Lin Sys dim UVW: " << 3*numNodes << " x " << 3*numNodes << endl;
 file << "Lin Sys dim P:   " << numVerts << " x " << numVerts << endl;
 file << "Lin Sys dim C:   " << numVerts << " x " << numVerts << endl;
 file << "Reynolds number: " << Re << endl;
 file << "Schmidt number:  " << Sc << endl;
 file << "Froud number:    " << Fr << endl;
 file << "Webber number:   " << We << endl;
 file << "alpha number:    " << alpha << endl;
 file << "beta number:     " << beta << endl;
 file << "CFL number:      " << cfl << endl;
 file << "dt:              " << dt << endl;
 file << "----------------------------------------------------" << endl; 
 file << endl;
 file << endl;

 file.close();
 cout << "simulation INFO file saved!" << endl;
}

void InOut::printInfo(const char* _mesh)
{
 time_t currentTime;

 time( &currentTime );
 cout << endl;
 cout << endl;
 cout << "              ";
 cout << "----------------------- INFO -----------------------" << endl; 
 cout << "               ";
 cout << "start process:   " << asctime( localtime( &currentTime ) ) << endl;
 cout << "               ";
 cout << "mesh:            " << _mesh << endl;
 cout << "               ";
 cout << "numVerts:        " << numVerts << endl;
 cout << "               ";
 cout << "numNodes:        " << numNodes << endl;
 cout << "               ";
 cout << "numElems:        " << numElems << endl;
 cout << "               ";
 cout << "Lin Sys dim UVW: " << 3*numNodes << " x " << 3*numNodes << endl;
 cout << "               ";
 cout << "Lin Sys dim P:   " << numVerts << " x " << numVerts << endl;
 cout << "               ";
 cout << "Lin Sys dim C:   " << numVerts << " x " << numVerts << endl;
 cout << "               ";
 cout << "Reynolds number: " << Re << endl;
 cout << "               ";
 cout << "Schmidt number:  " << Sc << endl;
 cout << "               ";
 cout << "Froud number:    " << Fr << endl;
 cout << "               ";
 cout << "Webber number:   " << We << endl;
 cout << "               ";
 cout << "alpha number:    " << alpha << endl;
 cout << "               ";
 cout << "beta number:     " << beta << endl;
 cout << "               ";
 cout << "CFL number:      " << cfl << endl;
 cout << "               ";
 cout << "dt:              " << dt << endl;
 cout << "              ";
 cout << "----------------------------------------------------" << endl; 
 cout << endl;
 cout << endl;

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
 file << setw(10) << *simTime << " " 
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
 file << setw(10) << *simTime << " " 
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
 file << setw(10) << *simTime << " " 
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
 file << setw(10) << *simTime << " " 
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
 file << setw(10) << *simTime << " " 
                  << kappaAverage << " " 
				  << setprecision(0) << fixed << _iter 
				  << endl;
 file.close();
}

void InOut::saveVTKSurface( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "TRI" + "-" + str + ".vtk";
 const char* filename = file.c_str();
 clMatrix *IENTri = m->getIENTri();
 int numTri = IENTri->DimI();

 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Modelo 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "FIELD FieldData 2" << endl;
 vtkFile << "TIME 1 1 double" << endl;
 vtkFile << *simTime << endl;
 vtkFile << "ITERATION 1 1 int" << endl;
 vtkFile << _iter << endl;
 vtkFile << endl;
 vtkFile << "POINTS " << numVerts << " double" << endl;

 vtkFile << setprecision(10) << scientific; 
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << X->Get(i) << " " 
          << setw(10) << Y->Get(i) << " " 
		  << setw(10) << Z->Get(i) << endl;

 vtkFile << endl;
 
 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 vtkFile << setprecision(0) << fixed; 
 for( int i=0;i<numTri;i++ )
 {
  vtkFile << "3 " << IENTri->Get(i,0) << " "  
                  << IENTri->Get(i,1) << " " 
                  << IENTri->Get(i,2) << endl;
 };
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;
 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;

 vtkFile << "SCALARS pressure double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << pSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS kappa double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << kappa->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS distance double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << distance->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uSol->Get(i) << " " 
          << vSol->Get(i) << " " 
		  << wSol->Get(i) << endl;
 vtkFile << endl;

 vtkFile << "VECTORS ALE-velocity double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uALE->Get(i) << " " 
          << vALE->Get(i) << " " 
		  << wALE->Get(i) << endl;
 vtkFile << endl;

 vtkFile << "VECTORS surface-force double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << fint->Get(i) << " " 
          << fint->Get(i+numVerts) << " " 
		  << fint->Get(i+numVerts*2) << endl;

 vtkFile.close();

 cout << "surface mesh No. " << _iter << " saved in VTK" << endl;

} // fecha metodo saveVTKSurface

void InOut::saveDistance(const char* _dir,const char* _filename,real time )
{
 string file = _dir;
 string aux = _filename;
 file += aux;
 const char* filename = file.c_str();

 ofstream dist( filename,ios::app );

 double Ymax1=-100;
 double Ymin1=100;
 double Ymax2=-100;
 double Ymin2=100;

 for( int i=0;i<surface->Dim();i++ )
 {
  if( Y->Get(i) > 0 )
  {
   if(Y->Get(i)>Ymax1) Ymax1=Y->Get(i);
   if(Y->Get(i)<Ymin1) Ymin1=Y->Get(i);
  }
  if( Y->Get(i) < 0 )
  {
   if(Y->Get(i)>Ymax2) Ymax2=Y->Get(i);
   if(Y->Get(i)<Ymin2) Ymin2=Y->Get(i);
  }
 }
 real dist1 = Ymin1-Ymax2;
 real dist2 = Ymax1-Ymin2;

 dist << setprecision(10) << scientific; 
 dist << setw(9) <<  time << " " << Ymin1 << " " << Ymax2 << " " 
                                 << Ymax1 << " " << Ymin2 << " " 
								 << dist1 << " " << dist2 << endl;

 dist.close();

 cout << "2 bubbles distance iter saved in ASCII " << dist1 << endl;
}

void InOut::saveMeshInfo(const char* _dir,const char* _filename )
{
 real *time = s->getTime();
 string file = (string) _dir + (string) _filename + ".dat";
 const char* filename = file.c_str();

 ifstream testFile( filename );
 ofstream mesh( filename,ios::app );
 if( testFile )
 {
  testFile.close();
  cout << "appending on file " << _filename << ".dat" << endl;
 }
 else
 {
  cout << "Creating file " << _filename << ".dat" << endl;
  mesh << "#time" << setw(20) << "numVerts" 
                  << setw(9) << "numNodes" 
				  << setw(9) << "numElems"
				  << endl;
 }


 mesh << setprecision(10) << scientific; 
 mesh << setw(9) <<  *time << " " << numVerts << " " 
                                  << numNodes << " " 
								  << numElems << endl; 

 mesh.close();

 cout << "meshing info saved" << endl;
}

void InOut::saveVTU( const char* _dir,const char* _filename, int _iter )
{
 IEN = m->getIEN();
 numElems = m->getNumElems();
 simTime = s->getTime();

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
 vtuFile << *simTime << endl;
 vtuFile << "ITERATION 1 1 int " << endl;
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

 // este if existe pois nem todos os metodos tem cc
 if( cSol->Dim() > 0 )
 {
  vtuFile << "    <DataArray type=\"Float32\" Name=\"Concentration\" format=\"ascii\"> " << endl;
  vtuFile << setprecision(10) << scientific;
  for( int i=0;i<numVerts;i++ )
   vtuFile << "     " << cSol->Get(i) << endl;
  vtuFile << "    </DataArray>" << endl;
 }

 vtuFile << "    <DataArray type=\"Float32\" Name=\"kappa\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << kappa->Get(i) << endl;
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

 vtuFile << "    <DataArray type=\"Float32\" Name=\"surface force\" NumberOfComponents=\"3\" format=\"ascii\"> " << endl;
 vtuFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtuFile << "     " << fint->Get(i) << " " 
                     << fint->Get(i+numVerts) << " " 
				     << fint->Get(i+numVerts*2) << endl;
 vtuFile << "    </DataArray>" << endl;

 /* end of file */
 vtuFile << "   </PointData> " << endl;
 vtuFile << "  </Piece> " << endl;
 vtuFile << " </UnstructuredGrid> " << endl;
 vtuFile << "</VTKFile> " << endl;

 vtuFile.close();

 cout << "solution No. " << _iter << " saved in VTU" << endl;

} // fecha metodo saveVTU
