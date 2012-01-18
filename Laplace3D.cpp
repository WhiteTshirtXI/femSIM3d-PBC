// =================================================================== //
// this is file Laplace3D.cpp, created at 21-Sep-2011                  //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Laplace3D.h"

Laplace3D::Laplace3D(){}

Laplace3D::Laplace3D( Model3D &_m )  
{
 getModel3DAttrib(_m);

 setSolver( new PCGSolver() );

 allocateMemoryToAttrib();
}

Laplace3D::Laplace3D( Model3D &_m,Laplace3D &_d )  
{
 getModel3DAttrib(_m);

 setSolver( new PCGSolver() );

 allocateMemoryToAttrib();

 //cSolOld = *_d.getCSol();
 cSolOld = *_m.getEdgeSize();
}

void Laplace3D::init()
{
 cSolOld.CopyFrom( 0,cc );
}

void Laplace3D::assemble()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 real aux;
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  real dif = 1.0;

  linElem.getM(*v); 

  for( i=0;i<NUMGLEC;i++ )
  {
   ii=v[i];
   for( j=0;j<NUMGLEC;j++ )
   {
	jj=v[j];
	aux = KcMat.Get(ii,jj) + dif*( linElem.kxxc[i][j] + 
	                               linElem.kyyc[i][j] + 
								   linElem.kzzc[i][j] );
	KcMat.Set(ii,jj,aux);

	aux = McMat.Get(ii,jj) + linElem.masselec[i][j];
	McMat.Set(ii,jj,aux);
   }
  }
 }
 
 Kc.CopyFrom(         0,          0,     KcMat );
 Mc.CopyFrom(         0,          0,     McMat );
 
}; // fecha metodo ASSEMBLENUC

void Laplace3D::setCRHS()
{
 vcc = McLumped * convC;
}

void Laplace3D::matMountC()
{
 //real k=5;     //   /\     more diffusion
 //real k=4;      //  /||\
 //real k=3;      //   ||     moderate diffusion 
 //real k=1;        //  \||/
 real k=0.7;      //   \/     less diffusion
 for( int i=0;i<numVerts;i++ )
 {
  real sumMc = Mc.SumLine(i);
  McLumped.Set( i,sumMc );
  invMcLumped.Set( i,1.0/sumMc );
 }
 matc = McLumped + (k * Kc);

}

void Laplace3D::setBC()
{
 cc.Dim(numVerts);
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
   idbcc.AddItem(i);

   real aux = triEdge[surfMesh->vertIdRegion.Get(i)];
   cc.Set(i,aux);

   convC.Set(i,aux);
 }

 clVector* vertIdRegion = m->getVertIdRegion();
 real minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  real radius = 0.5; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   real factor = triEdge[vertIdRegion->Get(i)]/minEdge;
   if( interfaceDistance->Get(i) < 3*radius )
   {
	real aux = triEdge[vertIdRegion->Get(i)]/factor;
	convC.Set(i,aux);
   }
   else
   {
	real aux = triEdge[vertIdRegion->Get(i)]/(factor*0.4);
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   real aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}


void Laplace3D::setUnCoupledCBC()
{
 int nbc,i,j;
 b1c.Dim(numVerts,0);  // zerando o vetore b1c
 ipc.Dim(numVerts,1); // inicializando vetor ipc com 1
 
 AcTilde = matc;

 nbc = idbcc.Dim();
 for( i=0;i<nbc;i++ )
 {
  j=(int) idbcc.Get(i);
  b1c.CopyMult(j,AcTilde,cc);
  AcTilde.Set(j,j,1);
  b1c.Set(j,cc.Get(j));
  ipc.Set(j,0);
 }
 cout << " boundary condition C --> SET " << endl;

} // fecha metodo setUnCoupledCBC 

void Laplace3D::unCoupledC()
{
 clVector vcIp(numVerts);
 clVector b1cTilde;

 vcIp = vcc.MultVec(ipc); // operacao vetor * vetor (elemento a elemento)

 b1cTilde = b1c + vcIp;

 // resolve sitema ATilde uTilde = b
 cout << " --------> solving scalar ----------- " << endl;
 solver->solve(1E-15,AcTilde,cTilde,b1cTilde);
 cout << " ------------------------------------ " << endl;

 cSol = cTilde;
 cSolOld = cSol;
}

void Laplace3D::getModel3DAttrib(Model3D &_m)
{
 m = &_m;

 // mesh information vectors
 numVerts = m->getNumVerts();
 numElems = m->getNumElems();
 numNodes = m->getNumNodes();
 X = m->getX();
 Y = m->getY();
 Z = m->getZ();
 IEN = m->getIEN();

 heaviside = m->getHeaviside();
 surfMesh = m->getSurfMesh();
 interfaceDistance = m->getInterfaceDistance();

 triEdge = m->getTriEdge();
 boundaryVert = m->getBoundaryVert();
 edgeSize = m->getEdgeSize();
}

void Laplace3D::allocateMemoryToAttrib()
{
 // assembly matrix
 Kc.Dim( numVerts,numVerts );
 Mc.Dim( numVerts,numVerts );
 McLumped.Dim( numVerts );

 // right hand side vectors
 vcc.Dim( numVerts );
 b1c.Dim( numVerts );

 // boundary condiction configured matrix
 AcTilde.Dim( numVerts,numVerts );

 // K + M matrix set
 matc.Dim( numVerts,numVerts );
 //invC.Dim( numVerts );
 invMcLumped.Dim( numVerts );

 // solution vectors 
 // vetores solucao
 cTilde.Dim( numVerts );
 cSol.Dim( numVerts );

 // auxiliar vectors
 ipc.Dim( numVerts,1 );

 // oldSol vectors
 cSolOld.Dim( numVerts );
}

clVector* Laplace3D::getCSol(){return &cSol;}
void Laplace3D::setSolver(Solver *s){solver = s;}
void Laplace3D::setModel3DEdgeSize(){ m->setEdgeSize(cSol); }

void Laplace3D::saveVTK( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = (string) _dir + (string) _filename + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename );

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "3D Simulation C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << endl;

 vtkFile << "POINTS " << numVerts << " double" << endl;
 for( int i=0;i<numVerts;i++ )
  vtkFile << X->Get(i) << " " 
          << Y->Get(i) << " " 
		  << Z->Get(i) << endl;

 vtkFile << endl;

 // conta numero de elementos
 real plane1 = ( X->Max()+X->Min() )/2.0;
 int count = 0;
 for( int i=0;i<numElems;i++ )
 {
  int v1 = IEN->Get(i,0);
  int v2 = IEN->Get(i,1);
  int v3 = IEN->Get(i,2);
  int v4 = IEN->Get(i,3);
  if( (heaviside->Get(v1)+heaviside->Get(v2)+
	   heaviside->Get(v3)+heaviside->Get(v4) > 1.5) || 
    ( (X->Get( v1 ) <  plane1) && (X->Get( v2 ) <  plane1) && 
	  (X->Get( v3 ) <  plane1) && (X->Get( v4 ) <  plane1) ) )
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
  if( (heaviside->Get(v1)+heaviside->Get(v2)+
	   heaviside->Get(v3)+heaviside->Get(v4) > 1.5) || 
    ( (X->Get( v1 ) <  plane1) && (X->Get( v2 ) <  plane1) && 
	  (X->Get( v3 ) <  plane1) && (X->Get( v4 ) <  plane1) ) )
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
  if( (heaviside->Get(v1)+heaviside->Get(v2)+
	   heaviside->Get(v3)+heaviside->Get(v4) > 1.5) || 
    ( (X->Get( v1 ) <  plane1) && (X->Get( v2 ) <  plane1) && 
	  (X->Get( v3 ) <  plane1) && (X->Get( v4 ) <  plane1) ) )
   vtkFile << "10 ";
 }

 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS edge-boundary double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << cc.Get(i) << endl;

 vtkFile << endl;
 vtkFile << "SCALARS edge-d double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << convC.Get(i) << endl;

 vtkFile << endl;
 vtkFile << "SCALARS edgeSol double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << cSol.Get(i) << endl;

 vtkFile << endl;

 vtkFile.close();
}

void Laplace3D::saveChordalEdge( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 // concatenando nomes para o nome do arquivo final
 string file = (string) _dir + (string) _filename + "." + str;
 const char* filename = file.c_str();

 ofstream vtkFile( filename,ios::app );

 vtkFile << "x-position" 
         << setw(17) << "y-position" 
		 << setw(17) << "z-position"
		 << setw(17) << "edge-boundary"
		 << setw(17) << "edge-init" 
		 << setw(17) << "edgeSol" 
		 << endl;

 // xVert da malha nova
 real nPoints = 1000;
 clVector xVert(nPoints);
 clVector yVert(nPoints);
 clVector zVert(nPoints);

 xVert.SetAll( (X->Max()+X->Min())/2.0 );
 yVert.SetAll( (Y->Max()+Y->Min())/2.0 );
 for( int i=0;i<nPoints;i++ )
 {
  real dz = i * ( (Z->Max()-Z->Min()) )/(nPoints-1);
  real pos = Z->Min()+dz;
  zVert.Set(i,pos);
 }

 // interpolacao linear em numVerts
 clMatrix interpLin = meshInterp(*m,xVert,yVert,zVert);
 clVector ccLin = interpLin*(cc);
 clVector convCLin = interpLin*(convC);
 clVector cSolLin = interpLin*(cSol);

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<nPoints;i++ )
  vtkFile << setw(10) << xVert.Get(i) << " " 
          << setw(17) << yVert.Get(i) <<  " "
          << setw(17) << zVert.Get(i) <<  " "
          << setw(17) << ccLin.Get(i) <<  " "
          << setw(17) << convCLin.Get(i) <<  " "
          << setw(17) << cSolLin.Get(i) << endl;

 vtkFile.close();

 /* --------- copying to file chordal.dat --------- */
 ifstream inFile( filename,ios::binary ); 

 string last = (string) _dir + (string) _filename + ".last";
 const char* filenameCopy = last.c_str();
 ofstream outFile( filenameCopy,ios::binary ); 

 outFile << inFile.rdbuf();
 inFile.close();
 outFile.close();
 /* ------------------------------------------------ */ 

 cout << "chordal edge No. " << _iter << " saved in dat" << endl;

} // fecha metodo chordalEdge
