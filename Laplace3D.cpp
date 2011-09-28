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
 cSolOld.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  if( heaviside->Get(i) < 0.5 )
   cSolOld.Set(i,triEdge[1]);
  else if( heaviside->Get(i) > 0.5 )
   cSolOld.Set(i,triEdge[2]);
  else
   cSolOld.Set(i,triEdge[2]);
 }
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
 //real k = 0.1;
 //vcc = ( (1.0/dt) * Mc - (1-alpha) * (1.0/(Re*Sc)) * Kc ) * convC;
 //vcc = ( (1.0/dt) * McLumped - (1-alpha) * (1.0/(Re*Sc)) * Kc ) * convC;
 //vcc = ( (1.0/dt) * McLumped ) * convC;
 //vcc = (1.0/k) * (cSolOld);
}

void Laplace3D::matMountC()
{
 real k=1.0;
 matc = k * Kc;
}

void Laplace3D::setBC()
{
 cc.Dim(numVerts);

 for (list<int>::iterator it=boundaryVert->begin(); 
                          it!=boundaryVert->end(); ++it)
 {
  idbcc.AddItem(*it);

  real aux = triEdge[1]; // edge length
  //real aux = triEdge[1]*triEdge[1]*triEdge[1]*sqrt(2.0)/12.0; // volume
  cc.Set(*it,aux);
 }

 for( int i=0;i<numVerts;i++ )
 {
  if( heaviside->Get(i) == 0.5 )
  {
   idbcc.AddItem(i);

   real aux = triEdge[2]; // edge length
   //real aux = triEdge[2]*triEdge[2]*triEdge[2]*sqrt(2.0)/12.0; // volume
   cc.Set(i,aux);
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

 triEdge = m->getTriEdge();
 boundaryVert = m->getBoundaryVert();
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
 //invMcLumped.Dim( numVerts );

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

 vtkFile << "CELLS " << numElems << " " << 5*numElems << endl;
 for( int i=0;i<numElems;i++ )
 {
   vtkFile << "4 " << IEN->Get(i,0) << " "
	               << IEN->Get(i,1) << " "
	               << IEN->Get(i,2) << " "
				   << IEN->Get(i,3) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
  vtkFile << "10 ";

 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS edge double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << cSol.Get(i) << endl;

 vtkFile << endl;

 vtkFile.close();
}
