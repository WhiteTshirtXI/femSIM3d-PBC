// =================================================================== // 
// this is file InOut, created at 23-Ago-2007                         //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail gustavo.rabello@gmail.com                                    //
// =================================================================== //

#include "InOut.h"

//using namespace std;

InOut::InOut( Model3D &_m, Simulator3D &_s )
{
 m = &_m;
 numVerts = m->getNumVerts();
 numNodes = m->getNumNodes();
 numElems = m->getNumElems();
 numGLEP = m->getNumGLEP();
 numGLEU = m->getNumGLEU();
 numGLEC = m->getNumGLEC();
 X = m->getPointerX();
 Y = m->getPointerY();
 Z = m->getPointerZ();
 uc = m->getPointerUC();
 vc = m->getPointerVC();
 wc = m->getPointerWC();
 pc = m->getPointerPC();
 cc = m->getPointerCC();
 idbcu = m->getPointerIdbcu();
 idbcv = m->getPointerIdbcv();
 idbcw = m->getPointerIdbcw();
 idbcp = m->getPointerIdbcp();
 idbcc = m->getPointerIdbcc();
 outflow = m->getPointerOutflow();
 IEN = m->getPointerIEN();
 surface = m->getPointerSurface();

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
 uAnt = s->getPointerUAnt();
 cAnt = s->getPointerCAnt();
 K = s->getPointerK();
 M = s->getPointerM();
 G = s->getPointerG();
 D = s->getPointerD();
 gx = s->getPointerGx();
 gy = s->getPointerGy();
 gz = s->getPointerGz();
 uSol = s->getPointerUSol();
 vSol = s->getPointerVSol();
 wSol = s->getPointerWSol();
 pSol = s->getPointerPSol();
 cSol = s->getPointerCSol();
 kappa = s->getPointerKappa();
 distance = s->getPointerDistance();


}

InOut::~InOut(){}

void InOut::saveSol( const char* _dir,const char* _filename, int iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
 ss >> str;

 clVector vec = *uAnt;
 vec.Append(*cSol); // adicionando o vetor concentracao no vetor Vel+Pressao

 
 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "UVWPC" + "-" + str + ".bin";
 const char* filenameUVWPC = fileUVWPC.c_str();

 ofstream UVWPC_file( filenameUVWPC,ios::binary ); 

 UVWPC_file.write( (const char*) vec.GetVec(),vec.Dim()*sizeof(real) );

 UVWPC_file.close();

 cout << "solucao no.  " << iter << " gravada em binario" << endl;
 
} // fecha metodo saveSol 

void InOut::loadSol( const char* _dir,const char* _filename, int iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
 ss >> str;

 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "UVWPC" + "-" + str + ".bin";
 const char* filenameUVWPC = fileUVWPC.c_str();

 clVector aux2(3*numNodes+2*numVerts); // vetor tambem carrega a concentracao

 ifstream UVWPC_file( filenameUVWPC,ios::in | ios::binary ); 

 UVWPC_file.read( (char*) aux2.GetVec(),aux2.Dim()*sizeof(real) );

 UVWPC_file.close();
 
 aux2.CopyTo(0,*uAnt);  // copyFrom nao funciona para o caso disk6-10-20
 aux2.CopyTo(3*numNodes+numVerts,*cSol);

 s->setUAnt(*uAnt); // impondo a velocidade no simulador
 s->setCSol(*cSol); // impondo a concentracao no simulador

 cout << "solucao no.  " << iter << " lida em binario" << endl;
 
} // fecha metodo loadSol 

void InOut::saveSolTXT( const char* _dir,const char* _filename, int iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
 ss >> str;

 clVector vec = *uAnt;
 vec.Append(*cSol); // adicionando a concentracao no vetor Vel+Pressao

 int i;
 string fileUVWPC = _dir;
 string aux = _filename;
 fileUVWPC += aux + "UVWPC" + "-" + str + ".dat";
 const char* filenameUVWPC = fileUVWPC.c_str();

 ofstream UVWPC_file( filenameUVWPC ); 
 UVWPC_file << numVerts << numNodes << endl << endl;

 for( i=0;i<vec.Dim();i++ )
   UVWPC_file << vec.Get(i) << endl;

 UVWPC_file << endl;
 UVWPC_file.close();

 cout << "solucao no.  " << iter << " gravada em TXT" << endl;
 
} // fecha metodo saveSolTXT 

void InOut::saveTXT( const char* _dir,const char* _filename, int iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
 ss >> str;

 clVector vec = *uAnt;
 vec.Append(*cSol);

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

 cout << "sistema no.  " << iter << " gravado em TXT" << endl;

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
 // concatenando nomes para o nome do arquivo final
 string file = _dir;
 string aux = _filename;
 file += aux + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Modelo 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "POINTS " << numVerts << " float" << endl;

 vtkFile << setprecision(6) << fixed; 
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

 vtkFile << "SCALARS pressure float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(6) << fixed; 
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << pc->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS concentration float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << cc->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS u float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << uc->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS v float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << vc->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS w float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << wc->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS vectors float" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << setw(10) << uc->Get(i) << " " 
          << setw(10) << vc->Get(i) << " " 
		  << setw(10) << wc->Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "mesh saved in VTK" << endl;

} // fecha metodo saveVTK 

void InOut::saveVTK( const char* _dir,const char* _filename, int iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = _dir;
 string aux = _filename;
 file += aux + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename ); 

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "Simulacao 3D C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << "POINTS " << numVerts << " float" << endl;

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

 vtkFile << "SCALARS pressure float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << fixed;
 for( int i=0;i<numVerts;i++ )
  vtkFile << pSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS concentration float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << cSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS u float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << uSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS v float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS w float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << wSol->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS kappa float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << kappa->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "SCALARS distance float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << distance->Get(i) << endl;

 vtkFile << endl;

 vtkFile << "VECTORS vectors float" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << uSol->Get(i) << " " 
          << vSol->Get(i) << " " 
		  << wSol->Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "solution num. " << iter << " saved in VTK" << endl;

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

void InOut::saveVonKarman(const char* _dir,const char* _filename,int iter,
                          int vertice )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
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
 //vonKarmanFile << setprecision(4) << setw(11) << setfill(' ') << scientific; 
 vonKarmanFile << setprecision(6) << fixed; 
 vonKarmanFile << "#z         F         G         H         c         p" << endl; 
 for( int i=0;i<Xaux.Dim();i++ )
 {
  v2 = (int) Xaux.Get(i);
  vonKarmanFile << setw(9) <<  Z->Get( v2 ) << " " 
	            << setw(9) <<  uSol->Get( v2 )/X->Get( v1 ) << " " 
	            << setw(9) <<  vSol->Get( v2 )/X->Get( v1 ) << " " 
		        << setw(9) <<  (-1)*wSol->Get( v2 ) << " " 
				<< setw(9) <<  cSol->Get(v2) << " " 
				<< setw(9) <<  pSol->Get(v2) 
				<< endl;
 }
 vonKarmanFile << endl;
 vonKarmanFile << "X = " << X->Get( (int) Xaux.Get(0) ) << endl;

 vonKarmanFile.close();

 cout << "von Karman num. " << iter << " saved in ASCII" << endl;
}

void InOut::savePert(const char* _dir,const char* _filename,int iter,
                     int vertice)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
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
   << atan(Y->Get( (int) aux2.Get(ntheta))/Y->Max())*180.0/PI 
   << " graus" << endl;
  vonKarmanFile << "# raio = " << raio << endl;
  vonKarmanFile << "# X = " << X->Get( (int) Y0.Get(ntheta) ) << endl;
  vonKarmanFile << "# Y = " << Y->Get( (int) aux2.Get(ntheta) ) << endl;
  vonKarmanFile << endl;
 }

 vonKarmanFile.close();

 cout << "von Karman perturb no. " << iter << " saved in ASCII" << endl;
}

void InOut::saveVortX(const char* _dir,const char* _filename,int iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
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
 vtkFile << "POINTS " << numVerts << " float" << endl;

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
 vtkFile << "SCALARS vort float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em X no. " << iter << " gravada em VTK" << endl;

}

void InOut::saveVortY(const char* _dir,const char* _filename,int iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
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
 vtkFile << "POINTS " << numVerts << " float" << endl;

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
 vtkFile << "SCALARS vort float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em Y no. " << iter << " gravada em VTK" << endl;

}

void InOut::saveVortZ(const char* _dir,const char* _filename,int iter)
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << iter;
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
 vtkFile << "POINTS " << numVerts << " float" << endl;

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
 vtkFile << "SCALARS vort float" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 for( int i=0;i<numVerts;i++ )
  vtkFile << vort.Get(i) << endl;

 vtkFile << endl;
 vtkFile.close();

 cout << "vorticidade em W no. " << iter << " gravada em VTK" << endl;

}

void InOut::saveTime(const char* _comment)
{
 ofstream file( "relatorio.dat",ios::app );
 if( !file )
 {
  cerr << "erro na abertura do arquivo de relatorio";
  exit(1);
 }

 time_t currentTime;

 time( &currentTime );
 file << asctime( localtime( &currentTime ) ) << "   - " << _comment << endl;
 file << endl;
}

void InOut::saveInfo(const char* _dir,const char* _mesh)
{
 string path = _dir;
 path += + "info.dat";
 const char* filePath = path.c_str();
 ofstream file( filePath,ios::app );

 time_t currentTime;

 time( &currentTime );
 file << "start process:   " << asctime( localtime( &currentTime ) ) << endl;
 file << "work directory:  " << _dir << endl;
 file << "mesh:            " << _mesh << endl;
 file << "Reynolds number: " << Re << endl;
 file << "Schmidt number:  " << Sc << endl;
 file << "Froud number:    " << Fr << endl;
 file << "Webber number:   " << We << endl;
 file << "alpha number:    " << alpha << endl;
 file << "beta number:     " << beta << endl;
 file << "CFL number:      " << cfl << endl;
 file << "dt:              " << dt << endl;
 file << endl;

 file.close();
 cout << "arquivo de informacoes da simulacao gravado!" << endl;
}

void InOut::printInfo(const char* _dir,const char* _mesh)
{
 string path = _dir;
 path += + "info.dat";
 const char* filePath = path.c_str();
 ofstream file( filePath,ios::app );

 time_t currentTime;

 time( &currentTime );
 cout << "-------------- INFO ----------------" << endl; 
 cout << "start process:   " << asctime( localtime( &currentTime ) ) << endl;
 cout << "work directory:  " << _dir << endl;
 cout << "mesh:            " << _mesh << endl;
 cout << "Reynolds number: " << Re << endl;
 cout << "Schmidt number:  " << Sc << endl;
 cout << "Froud number:    " << Fr << endl;
 cout << "Webber number:   " << We << endl;
 cout << "alpha number:    " << alpha << endl;
 cout << "beta number:     " << beta << endl;
 cout << "CFL number:      " << cfl << endl;
 cout << "dt:              " << dt << endl;
 cout << "------------------------------------" << endl; 

 file.close();
 cout << "arquivo de informacoes da simulacao gravado!" << endl;
}

void InOut::oscillating(const char* _file)
{
 ofstream file( _file,ios::app );

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

 //file << "#time xPoint yPoint" << endl;
 file << *simTime << " " << pointX << " " 
                         << pointY << " "
						 << pointZ << endl;
}

void InOut::oscillating(int point1,int point2,int point3,const char* _file)
{
 ofstream file( _file,ios::app );

 // retorna o valor de maior Y da interface 
 real pointX = uSol->Get(point1);
 real pointY = vSol->Get(point2);
 real pointZ = wSol->Get(point3);

 //file << "#time xPoint yPoint" << endl;
 file << *simTime << " " << pointX << " " << pointY << " " << pointZ << endl;
}

void InOut::oscillatingD(int point1,int point2,int point3,int point4,
						 int point5,int point6,const char* _file)
{
 ofstream file( _file,ios::app );

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
                 

 //file << "#time Xdiameter Ydiameter" << endl;
 file << *simTime << " " << diameterX << " " 
                         << diameterY << " " 
			    		 << diameterZ << endl;
}

void InOut::oscillatingD(const char* _file)
{
 ofstream file( _file,ios::app );

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

 //file << "#time Xdiameter Ydiameter Zdiameter" << endl;
 file << *simTime << " " << diameterX << " " 
                         << diameterY << " " 
						 << diameterZ << endl;
}

void InOut::oscillatingKappa(const char* _file)
{
 ofstream file( _file,ios::app );

 real norm = kappa->Norm();

 //file << "#time kappaNorm" << endl;
 file << *simTime << " " << norm << endl;
}

