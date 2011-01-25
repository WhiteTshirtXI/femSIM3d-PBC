// =================================================================== //
// this is file Model3D.cpp, created at 23-Ago-2007                    //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Model3D.h"
#include "tetgen.h"

// necessario para leitura do objeto triangle.o escrito em ANSI C
extern "C"
{
#include "triangle.h"
}

using namespace std;

Model3D::Model3D(){}
Model3D::~Model3D(){}

void Model3D::readVTK( const char* filename )
{
 char auxstr[255];
 real coords[3];
 int i,j,k,idnv;
 int auxvtx[4];

 ifstream vtkFile( filename,ios::in );

 if( !vtkFile )
 {
  cerr << "Esta faltando o arquivo de Malha!" << endl;
  exit(1);
 }

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"POINTS") != 0) )
  vtkFile >> auxstr;

 if( !vtkFile.eof())
 {
  vtkFile >> numVerts;
  vtkFile >> auxstr;
  
  X.Dim(numVerts);
  Y.Dim(numVerts);
  Z.Dim(numVerts);
  
  for (i=0; i < numVerts; i++)
  {
   for(j = 0; j < 3; j++)
	vtkFile >> coords[j];
   
   X.Set(i,coords[0]);
   Y.Set(i,coords[1]);
   Z.Set(i,coords[2]);
  }

  while( (! vtkFile.eof())&&(strcmp(auxstr,"CELLS") != 0) )
   vtkFile >> auxstr;

  if( !vtkFile.eof() )
  {
   vtkFile >> numElems;
   vtkFile >> auxstr;

   IEN.Dim(numElems,4);

   for( i=0; i < numElems; i++ )
   {
	vtkFile >> idnv;
	for( j=0; j < 4 ; j++ )
	{
	 vtkFile >> k;
	 IEN.Set(i,j,k);
	 auxvtx[j] = k;
	}
   }
  }
 }
 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readVTKCC( const char* filename )
{
 char auxstr[255];
 real fl;
 cc.Dim(numVerts);

 ifstream vtkFile( filename,ios::in );

 if( !vtkFile )
 {
  cerr << "Esta faltando o arquivo de leitura de CC!" << endl;
  exit(1);
 }

 while( (! vtkFile.eof())&&(strcmp(auxstr,"concentration") != 0) )
  vtkFile >> auxstr;

 if( !vtkFile.eof() )
 {
  vtkFile >> auxstr;
  vtkFile >> auxstr;
  vtkFile >> auxstr;

  for( int i=0; i < numVerts; i++ )
  {
   vtkFile >> fl;
   cc.Set(i,fl);
  }
 }
 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readVTKSurface( const char* filename )
{
 char auxstr[255];
 real coords[3];
 int i,j,k,idnv;
 int auxvtx[4];

 ifstream vtkFile( filename,ios::in );

 if( !vtkFile )
 {
  cerr << "Esta faltando o arquivo de Malha!" << endl;
  exit(1);
 }

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"POINTS") != 0) )
  vtkFile >> auxstr;

 if( !vtkFile.eof())
 {
  vtkFile >> numVerts;
  vtkFile >> auxstr;
  
  X.Dim(numVerts);
  Y.Dim(numVerts);
  Z.Dim(numVerts);
  
  for (i=0; i < numVerts; i++)
  {
   for(j = 0; j < 3; j++)
	vtkFile >> coords[j];
   
   X.Set(i,coords[0]);
   Y.Set(i,coords[1]);
   Z.Set(i,coords[2]);
  }

  while( (! vtkFile.eof())&&(strcmp(auxstr,"CELLS") != 0) )
   vtkFile >> auxstr;

  if( !vtkFile.eof() )
  {
   vtkFile >> numElems;
   vtkFile >> auxstr;

   IEN.Dim(numElems,4);

   for( i=0; i < numElems; i++ )
   {
	vtkFile >> idnv;
	for( j=0; j < 3 ; j++ )
	{
	 vtkFile >> k;
	 IEN.Set(i,j,k);
	 auxvtx[j] = k;
	}
   }
  }
 }
 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readMSH( const char* filename )
{
 char auxstr[255];
 real coords[3];
 int i,j,k,id;
 int auxvtx[4];
 int elemNumber,type,numberOfTags;

 ifstream mshFile( filename,ios::in );

 if( !mshFile )
 {
  cerr << "Esta faltando o arquivo de Malha!" << endl;
  exit(1);
 }

//--------------------------------------------------
//  int numberOfPhyNames;
//  while( ( !mshFile.eof())&&(strcmp(auxstr,"$PhysicalNames") != 0) )
//   mshFile >> auxstr;
// 
//  mshFile >> numberOfPhyNames;
// 
//  idRegion.resize(numberOfPhyNames);
//  if( ( !mshFile.eof())&&(strcmp(auxstr,"$EndPhysicalNames") != 0) )
//  {
//   for (i=0; i < numberOfPhyNames; i++)
//   {
//    mshFile >> auxstr;
//    mshFile >> auxstr;
//    mshFile >> auxstr;
//    idRegion.at(i)=auxstr;
//   }
//  }
//-------------------------------------------------- 

 while( ( !mshFile.eof())&&(strcmp(auxstr,"$Nodes") != 0) )
  mshFile >> auxstr;

 if( !mshFile.eof()&&(strcmp(auxstr,"$EndNodes") != 0)   )
 {
  mshFile >> surfMesh.numVerts;

  surfMesh.X.Dim(surfMesh.numVerts);
  surfMesh.Y.Dim(surfMesh.numVerts);
  surfMesh.Z.Dim(surfMesh.numVerts);

  for (i=0; i < surfMesh.numVerts; i++)
  {
   mshFile >> auxstr;
   for(j = 0; j < 3; j++)
	mshFile >> coords[j];

   surfMesh.X.Set(i,coords[0]);
   surfMesh.Y.Set(i,coords[1]);
   surfMesh.Z.Set(i,coords[2]);
  }
 }

 while( (!mshFile.eof())&&(strcmp(auxstr,"$Elements") != 0) )
  mshFile >> auxstr;

 if( !mshFile.eof()&&(strcmp(auxstr,"$EndElements") != 0)   )
 {
  mshFile >> surfMesh.numElems;

  surfMesh.IEN.Dim(surfMesh.numElems,3);
  idRegion.Dim(surfMesh.numElems);

  for( i=0; i < surfMesh.numElems; i++ )
  {
   mshFile >> elemNumber;
   mshFile >> type; // 2-2D or 3-3D
   mshFile >> numberOfTags;  
   if( numberOfTags == 3 ) // msh file version 2.1
   {
	// idRegion 1 = surface
	// idRegion 2 = wall
	// idRegion 3 = bubble
	mshFile >> id;
	idRegion.Set(i,id);
	mshFile >> auxstr;
	mshFile >> auxstr;
   }
   else // msh file vertion 2.2
   {
	mshFile >> id;
	idRegion.Set(i,id);
	mshFile >> auxstr;
   }

   for( j=0; j < type+1 ; j++ )
   {
	mshFile >> k;
	k=k-1; // elem .msh comecando com 1
	surfMesh.IEN.Set(i,j,k);
	auxvtx[j] = k;
   }
//--------------------------------------------------
//    if( region == 1 ) // 1 = interface 
//    {
// 	int v1 = IEN.Get(i,0);
// 	int v2 = IEN.Get(i,1);
// 	int v3 = IEN.Get(i,2);
// 	cc.Set(v1,0.5);
// 	cc.Set(v2,0.5);
// 	cc.Set(v3,0.5);
//    }
//-------------------------------------------------- 
  }
 }
 mshFile.close();
} // fim do metodo readMsh

void Model3D::setInterfaceBC()
{
 surfMesh.Marker.Dim(surfMesh.numVerts);
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  if( idRegion.Get(i) == 1 )
  {
   int v1 = surfMesh.IEN.Get(i,0);
   int v2 = surfMesh.IEN.Get(i,1);
   int v3 = surfMesh.IEN.Get(i,2);

   real aux = 0.5;
   surfMesh.Marker.Set(v1,aux);
   surfMesh.Marker.Set(v2,aux);
   surfMesh.Marker.Set(v3,aux);
  }
 }
}

/**
 * @brief metodo para leitura de arquivo do tipo BC para condicoes de
 * contorno. O arquivo a ser lido deve conter todos os nos de condicao
 * de contorno contendo colunas de vertices, colunas de u,v e p com
 * numeracao 1 para condicao do tipo Dirichlet e 2 para Newmann e seu
 * respectivo valor
 *
 * @return verdadeiro ou falso dependendo do sucesso da leitura
 **/
void Model3D::readBC( const char* filename )
{

 char auxstr[255];
 real coords[9];
 int i,j,numVertsBC;

 ifstream bcFile( filename,ios::in );

 if( bcFile == NULL )
 {
  cerr << "Esta faltando o arquivo de condicao de contorno!" << endl;
  exit(1);
 }

 while( ( !bcFile.eof())&&(strcmp(auxstr,"BC_DATA") != 0) )
  bcFile >> auxstr; // atribui BC_DATA
 
 if( !bcFile.eof())
 {
  bcFile >> auxstr;
  bcFile >> numVertsBC;

  for (i=0; i < numVertsBC; i++)
  {
   for(j = 0; j < 9; j++)
	bcFile >> coords[j];
   
   if( coords[1] == 1 )
   {
	idbcu.AddItem( (int) coords[0]);
	uc.Set( (int) coords[0],coords[2]);
   }

   if( coords[3] == 1 )
   {
	idbcv.AddItem(coords[0]);
	vc.Set( (int) coords[0],coords[4]);
   }

   if( coords[5] == 1 )
   {
	idbcw.AddItem(coords[0]);
	wc.Set( (int) coords[0],coords[6]);
   }

   if( coords[7] == 1 )
   {
	idbcp.AddItem(coords[0]);
	pc.Set( (int) coords[0],coords[8]);
	outflow.Set( (int) coords[0], 0);
   }

  }
 }
} // fim do metodo readBC

// este metodo cria os pontos de forma ordenada e igualmente espacada e
// depois utiliza a biblioteca tetgen para gerar a tetraedralizacao
// seguida pela atualizacao da matriz de mapeamento de elementos IEN
void Model3D::setMeshStep(int nX,int nY,int nZ)
{
 tetgenio in,out;
 in.mesh_dim = 3;
 in.numberofpoints = nX*nY*nZ;
 in.pointlist = new REAL[in.numberofpoints * 3];

 real dx = nX/(nX-1.0);
 real dy = nY/(nY-1.0);
 real dz = nZ/(nZ-1.0);
 real aux;
 numVerts = nX*nY*nZ;
 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 outflow.Dim(numVerts);

 int count = 0;
 for( int i=0;i<nX;i++ )
 {
  for( int j=0;j<nY;j++ )
  {
   for( int k=0;k<nZ;k++)
   {
	aux = i*dx;
	X.Set(count,aux);
	in.pointlist[3*count+0] = aux;
	aux = j*dy;
	Y.Set(count,aux);
	in.pointlist[3*count+1] = aux;
	aux = k*dz;
	Z.Set(count,aux);
	in.pointlist[3*count+2] = aux;
	count++;
   }
  }
 }

 cout << endl;
 cout << "----> meshing... ";
 tetrahedralize( (char*) "Q",&in,&out,NULL,NULL );
 cout << "finished <---- " << endl;;
 cout << endl;

 //out.save_elements("out");
 //out.save_nodes("out");
 numElems = out.numberoftetrahedra;
 IEN.Dim(numElems,4);

 // varre lista de elementos e passa para estrutura IEN
 for( int i=0;i<out.numberoftetrahedra;i++ )
 {
  for( int j=0;j<4;j++ )
  {
   int vertice = out.tetrahedronlist[i*4+j];
   IEN.Set(i,j,vertice);
  }
 }
}

void Model3D::setStepBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max()) )
   {

	uc.Set(i,1.0);
   }
  }
  if( (Z.Get(i)==Z.Min()) || (Z.Get(i) == Z.Max()) ||
	  (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   idbcw.AddItem(i);
   wc.Set(i,0.0);
  }
  if( X.Get(i)==X.Max())
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
   outflow.Set(i,0);
  }
 }
}

void Model3D::setCStepBC()
{
 cc.Dim(numVerts);
 idbcc.Dim(0);
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max()) )
   {
	idbcc.AddItem(i);

	cc.Set(i,1.0);
   }
  }
 }
}

void Model3D::setStepReservBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) 
	  || Z.Get(i)==Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max())  
	    && Z.Get(i)>Z.Min() )
	uc.Set(i,1.0);
  }
  if( Z.Get(i) == Z.Max() && X.Get(i)<X.Max() )
  {
   idbcw.AddItem(i);
   wc.Set(i,0.0);
  }
  if( X.Get(i)==X.Max() && Z.Get(i)>Z.Min() )
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
   outflow.Set(i,0);
  }
 }
}

void Model3D::setStepReservInvBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Max()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) 
	  || Z.Get(i)==Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
   if( (X.Get(i)==X.Max()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max())  
	    && Z.Get(i)>Z.Min() )
	uc.Set(i,-1.0);
  }
  if( Z.Get(i) == Z.Max() && X.Get(i)>X.Min() )
  {
   idbcw.AddItem(i);
   wc.Set(i,0.0);
  }
  if( X.Get(i)==X.Min() && Z.Get(i)>Z.Min() )
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
   outflow.Set(i,0);
  }
 }
}

void Model3D::setCouetteBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>Y.Min()) && (Y.Get(i)<Y.Max()) )
	uc.Set(i,1.0);
  }
  if( (Z.Get(i)==Z.Min()) || (Z.Get(i) == Z.Max()) )
  {
   idbcw.AddItem(i);
   wc.Set(i,0.0);
  }
  if( X.Get(i)==X.Max())
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
   outflow.Set(i,0);
  }
 }
}

void Model3D::setAdimenStep()
{
 real aux;
 real factor = 1.0/(Y.Max()-Y.Min());

 for( int i=0;i<numVerts;i++ )
 {
  aux = (X.Get(i)-X.Min())*factor;
  X.Set(i,aux);
  aux = (Y.Get(i)-Y.Min())*factor;
  Y.Set(i,aux);
  aux = (Z.Get(i)-Z.Min())*factor;
  Z.Set(i,aux);
 }
}

void Model3D::setMeshDisk(int nLados1Poli,int nCircMax,int nZ)
{
 real dr = 1;
 real  r = dr;
 int j = 0;
 real z = 0;
 real dl = ( (2*3.141592)/nLados1Poli)*dr;
 real theta,dTheta,dz;
 real aux;
 //clVector xCirc(1+nLados1Poli*nCircMax!-fatorial);
 //clVector yCirc(1+nLados1Poli*nCircMax!-fatorial);
 clVector xCirc;
 clVector yCirc;
 clVector xAux,yAux,zAux;
 X.Dim(1000);
 Y.Dim(1000);
 Z.Dim(1000);

 xCirc.AddItem(j,0);
 yCirc.AddItem(j,0);
 j++;

 for( int nCirc=1;nCirc<=nCircMax;nCirc++ )
 {
  theta = 0.0;
  dTheta = (dl/nCirc)*dr;
  for( int jj=1;jj<=(nLados1Poli*nCirc);jj++ )
  {
   aux = r*cos(theta);
   xCirc.AddItem(j,aux);
   //xCirc.Set(j,aux);
   aux = r*sin(theta);
   yCirc.AddItem(j,aux);
   //yCirc.Set(j,aux);
   theta = theta + dTheta;
   j++;
   if( theta >= 2*3.141592 ) break;
  }
  r=r+dr;
 }

 rMax = r - dr;
 j=0;
 z=0;
 ///////real factor = 1.1 + 1.0/nZ;
 //real factor = 1.02 + 1.0/nZ;
 for( int jz=1;jz<=nZ;jz++ )
 {
  //if( jz<=nZ/2 ) dz=0.1;
  //else 
  //if( exp(z/nZ) <= exp(2.0/3.0) ) dz=exp(z/nZ);
  //else dz = exp(2.0/3.0);
  //else dz=0.1;
  ////////dz = (factor*exp(z/nZ));
  dz = 0.1;
  //cout << z  << " " << nZ << " " << dz << endl;
  for( int jCirc=1;jCirc<=xCirc.Dim();jCirc++ )
  {
   if( j == X.Dim() ) 
   {
	xAux = X;
	yAux = Y;
	zAux = Z;
	X.Dim(j+X.Dim());
	Y.Dim(j+Y.Dim());
	Z.Dim(j+Z.Dim());
	X.CopyFrom(0,xAux);
	Y.CopyFrom(0,yAux);
	Z.CopyFrom(0,zAux);
   }
   aux = xCirc.Get(jCirc-1);
   X.Set(j,aux);
   aux = yCirc.Get(jCirc-1);
   Y.Set(j,aux);
   Z.Set(j,z);
   j++;
  }
  z=z+dz;
 }
 xAux = X;
 yAux = Y;
 zAux = Z;
 numVerts = j; 
 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 xAux.CopyTo(0,X);
 yAux.CopyTo(0,Y);
 zAux.CopyTo(0,Z);

 tetgenio in,out;
 in.mesh_dim = 3;
 in.numberofpoints = numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];

 for( int i=0;i<numVerts;i++ )
 {
  in.pointlist[3*i+0] = X.Get(i);
  in.pointlist[3*i+1] = Y.Get(i);
  in.pointlist[3*i+2] = Z.Get(i);
 }

 cout << endl;
 cout << "----> meshing... ";
 tetrahedralize( (char*) "Q",&in,&out );
 cout << "finished <---- " << endl;;
 cout << endl;

 //out.save_elements("out");
 //out.save_nodes("out");
 numElems = out.numberoftetrahedra;
 IEN.Dim(numElems,4);

 // varre lista de elementos e passa para estrutura IEN
 for( int i=0;i<out.numberoftetrahedra;i++ )
 {
  for( int j=0;j<4;j++ )
  {
   int vertice = out.tetrahedronlist[i*4+j];
   IEN.Set(i,j,vertice);
  }
 }
}

/*  This method makes re-meshing test and can be modified */
void Model3D::meshTest()
{
//--------------------------------------------------
//  clVector test = cc==0.5;
//  clVector t1 = test.Find();
// 
//  clVector x1,y1,z1;
//  for( int i=0;i<t1.Dim();i++ )
//  {
//   x1.AddItem(X.Get(t1.Get(i)));
//   y1.AddItem(Y.Get(t1.Get(i)));
//   z1.AddItem(Z.Get(t1.Get(i)));
//  }
//  numVerts = x1.Dim();
//-------------------------------------------------- 

 // cria objeto de malha do tetgen
 tetgenio in,out;
 in.mesh_dim = 3;
 in.numberofpoints = numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 // adiciona na estrutura tetgen as coordenadas dos pontos
 for( int i=0;i<numVerts;i++ )
 {
  in.pointlist[3*i+0] = X.Get(i);
  in.pointlist[3*i+1] = Y.Get(i);
  in.pointlist[3*i+2] = Z.Get(i);
  if( cc.Get(i) == 0.0 )
   in.pointmarkerlist[i] = 11; // fora
  if( cc.Get(i) == 0.5 )
   in.pointmarkerlist[i] = 22; // interface
  if( cc.Get(i) == 1.0 )
   in.pointmarkerlist[i] = 33; // dentro
 }

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");
 cout << "numElems IN = " << numElems << endl;
 cout << "numNodes IN = " << numNodes << endl;
 cout << "numVerts IN = " << in.numberofpoints << endl;

 cout << endl;
 cout << "----> meshing... ";
 //tetrahedralize( (char*) "QYYApq1.4241a0.05",&in,&out );
 tetrahedralize( (char*) "",&in,&out );
 cout << "finished <---- " << endl;;
 cout << endl;
 
 numElems = out.numberoftrifaces;
 numNodes = out.numberofpoints+numElems;
 numVerts = out.numberofpoints;
 cout << "numElems OUT = " << out.numberoftetrahedra << endl;
 cout << "numNodes OUT = " << out.numberofpoints+numElems << endl;
 cout << "numVerts OUT = " << out.numberofpoints << endl;
 cout << "numfacets OUT = " << out.numberoftrifaces << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 IEN.Dim(numElems,4);
 for( int i=0;i<out.numberoftrifaces;i++ )
 {
  for( int j=0;j<3;j++ )
  {
   int vertice = out.trifacelist[i*3+j];
   IEN.Set(i,j,vertice);
  }
  //cout << out.trifacemarkerlist[0] << endl;
 }

//--------------------------------------------------
//  // varre lista de elementos e passa para estrutura IEN
//  IEN.Dim(numElems,5);
//  cc.Dim(numVerts);
//  inElem.resize (0);
//  outElem.resize (0);
//  // varre lista de elementos e passa para estrutura IEN
//  for( int i=0;i<out.numberoftetrahedra;i++ )
//  {
//   for( int j=0;j<4;j++ )
//   {
//    int vertice = out.tetrahedronlist[i*4+j];
//    IEN.Set(i,j,vertice);
//   }
//  }
//-------------------------------------------------- 

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 X.Dim(numNodes);
 Y.Dim(numNodes);
 Z.Dim(numNodes);
 for( int i=0;i<numVerts;i++ )
 {
  X.Set(i,out.pointlist[3*i+0]);
  Y.Set(i,out.pointlist[3*i+1]);
  Z.Set(i,out.pointlist[3*i+2]);
 }
}

void Model3D::mesh2Dto3D()
{
 // cria objeto de malha do tetgen
 tetgenio in,out;
 in.mesh_dim = 3;
 in.numberofpoints = surfMesh.numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 // adiciona na estrutura tetgen as coordenadas dos pontos
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  in.pointlist[3*i+0] = surfMesh.X.Get(i);
  in.pointlist[3*i+1] = surfMesh.Y.Get(i);
  in.pointlist[3*i+2] = surfMesh.Z.Get(i);
  if( surfMesh.Marker.Get(i) == 0.0 )
   in.pointmarkerlist[i] = 11; // fora
  if( surfMesh.Marker.Get(i) == 0.5 )
   in.pointmarkerlist[i] = 22; // interface
 }

 /* ESTE PROCEDIMENTO DEFINE REGIOES NA MALHA E APOS A INSERCAO/RETIRADA
  * DE PONTOS PELO TETGEN, CONSEGUIMOS RECONHECER A LOCALIZACAO DOS
  * PONTOS E ASSIM PODEMOS DEFINIR NOVAMENTE A FUNCAO MARCADORA COMO
  * SENDO 1.0 DENTRO DA BOLHA, 0.5 NA SUPERFICIE E 0.0 FORA 
  * E NECESSARIO DEFINIR 1 PONTO EM CADA REGIAO */
 // fluido interior + fluido exterior + superficie
 in.numberofregions = 1; 
 in.regionlist = new REAL[in.numberofregions*4];

 // fora da bolha
 //in.regionlist[0] = X.Min();
 //in.regionlist[1] = Y.Min();
 //in.regionlist[2] = Z.Min();
//--------------------------------------------------
//  in.regionlist[0] = 0.1;
//  in.regionlist[1] = 0.1;
//  in.regionlist[2] = 0.1;
//-------------------------------------------------- 
 in.regionlist[0] = -5.8;
 in.regionlist[1] = 0.0;
 in.regionlist[2] = -2.8;
 in.regionlist[3] = 1;

 tetgenio::facet *f;   // Define a pointer of facet. 
 tetgenio::polygon *p; // Define a pointer of polygon.
 in.numberoffacets = surfMesh.numElems; 
 in.facetlist = new tetgenio::facet[in.numberoffacets]; 
 in.facetmarkerlist = new int[in.numberoffacets];
 //in.trifacemarkerlist = new int[in.numberoffacets];

 // definindo a superficie da bolha e convex-hull
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  f = &in.facetlist[i];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0; 
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 3;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = v1; 
  p->vertexlist[1] = v2; 
  p->vertexlist[2] = v3;
  // melhorar esta configuracao de facet para bolha e convex hull
  if( surfMesh.Marker.Get(v1) + 
	  surfMesh.Marker.Get(v2) + 
	  surfMesh.Marker.Get(v3) > 0 )
   in.facetmarkerlist[i] = 10;
  else
   in.facetmarkerlist[i] = 20;
 }

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");
 cout << "numElems IN = " << surfMesh.numElems << endl;
 cout << "numNodes IN = " << surfMesh.numNodes << endl;
 cout << "numVerts IN = " << surfMesh.numVerts << endl;

 cout << endl;
 cout << "----> meshing... ";
 tetrahedralize( (char*) "QYYCApq1.414q10a0.1",&in,&out );
 //tetrahedralize( (char*) "QYYApq1.4241a0.1",&in,&out );
 //tetrahedralize( (char*) "QYYApq1.4241a0.05",&in,&out );
 //tetrahedralize( (char*) "QpYY",&in,&out );
 cout << "finished <---- " << endl;;
 cout << endl;
 
 numElems = out.numberoftetrahedra;
 numNodes = out.numberofpoints+numElems;
 numVerts = out.numberofpoints;
 cout << "numElems OUT = " << out.numberoftetrahedra << endl;
 cout << "numNodes OUT = " << out.numberofpoints+numElems << endl;
 cout << "numVerts OUT = " << out.numberofpoints << endl;
 cout << "numfacets OUT = " << out.numberoftrifaces << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 // varre lista de elementos e passa para estrutura IEN
 IEN.Dim(numElems,5);
 cc.Dim(numVerts);
 for( int i=0;i<out.numberoftetrahedra;i++ )
 {
  // setting de cc = 0 para fora da bolha e cc = 0.5 para interface
  if( out.tetrahedronattributelist[i] == 1 )
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,0.0);
   }
  }
  // setting de cc = 1 para dentro da bolha e cc = 0.5 para interface
  else 
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,1.0);
   }
  }
 }

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 X.Dim(numNodes);
 Y.Dim(numNodes);
 Z.Dim(numNodes);
 for( int i=0;i<numVerts;i++ )
 {
  X.Set(i,out.pointlist[3*i+0]);
  Y.Set(i,out.pointlist[3*i+1]);
  Z.Set(i,out.pointlist[3*i+2]);
  if( out.pointmarkerlist[i] == 10 ||
   	  out.pointmarkerlist[i] == 22 )
   cc.Set(i,0.5);
 }
}

void Model3D::setTriEdge()
{
 int v1,v2,v3;
 int numFace = 3; // triangulo tem 3 arestas
 clVector faceaux(2);
 IFACE3DSurface *faces = NULL;
 int listSize = numFace*surfMesh.numElems;
 faces = new IFACE3DSurface[listSize];
 for( int mele=0;mele<surfMesh.numElems;mele++ )
 {
  v1 = (int) surfMesh.IEN.Get(mele,0);
  v2 = (int) surfMesh.IEN.Get(mele,1);
  v3 = (int) surfMesh.IEN.Get(mele,2);

  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+0].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+0].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+0].p3 = v3;
  faces[numFace*mele+0].p4 = mele;

  faceaux.Set(0,v1);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+1].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+1].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+1].p3 = v2;
  faces[numFace*mele+1].p4 = mele;

  faceaux.Set(0,v2);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+2].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+2].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+2].p3 = v1;
  faces[numFace*mele+2].p4 = mele;
 }

 // ordena uma estrutura (faces) em ordem crescente na linha e coluna
 // as faces continuam repetidas neste ponto, porem ordenadas e prontas
 // para serem excluidas.
 qsort(faces,listSize,sizeof(IFACE3DSurface),IFACE2DCompare);

//--------------------------------------------------
//  for( int i=0;i<listSize;i++ )
//   cout << faces[i].p1 << " " << faces[i].p2 << " " 
//        << faces[i].p3 << " " << faces[i].p4 << endl;
//-------------------------------------------------- 

 /*        - nome: mapEdgeTri
           - definicao: matrix com mapeamento de arestas da superficie e
		                convex hull, a dimensao da matrix eh o numero de 
						total de arestas.

  ---   +---+---+---+---+---+---+---+
   |    | a | b | c | d | e | f | g |   a = coordenada X do ponto 
   |    +---+---+---+---+---+---+---+       medio da aresta
   |    .   .   .   .   .   .   .   .   b = coordenada Y do ponto  
   |    .   .   .   .   .   .   .   .       medio da aresta
   h    .   .   .   .   .   .   .   .   c = coordenada Z do ponto 
   |    +---+---+---+---+---+---+---=       media da aresta                    
   |    | a | b | c | d | e | f | g |   d = 1o. vertice da aresta 
   |    +---+---+---+---+---+---+---+
   |    | a | b | c | d | e | f | g |   e = 2o. vertice da aresta
  ---   +---+---+---+---+---+---+---+
                                        f = 3o. vertice do 1o. elemento
        |____________ i ____________|
	    |                           |   g = 3o. vertice do 2o. elemento
		
		                                h = numero de arestas e numercao
      
		                                i = 7 colunas
 
 */

 int j=0;
 mapEdgeTri.Dim(listSize/2,7);

 // numeracao de arestas a partir de numVerts e associacao das arestas
 // aos elementos para inclusao na IEN
 for( int i=0;i<listSize/2;i++ )
 {
  real x1=surfMesh.X.Get(faces[j].p1);
  real y1=surfMesh.Y.Get(faces[j].p1);
  real z1=surfMesh.Z.Get(faces[j].p1);
  real x2=surfMesh.X.Get(faces[j].p2);
  real y2=surfMesh.Y.Get(faces[j].p2);
  real z2=surfMesh.Z.Get(faces[j].p2);
  real length = sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );

  mapEdgeTri.Set(i,0,length); // tamanho da aresta
  mapEdgeTri.Set(i,1,faces[j].p1 ); // numero do 1o. vertice da aresta
  mapEdgeTri.Set(i,2,faces[j].p2 ); // numero do 2o. vertice da areata
  mapEdgeTri.Set(i,3,faces[j].p3 );   // numero do 3o. vertice do 1o. elemento
  mapEdgeTri.Set(i,4,faces[j+1].p3 ); // numero do 3o. vertice do 2o. elemento
  mapEdgeTri.Set(i,5,faces[j].p4 ); // 1o. elemento
  mapEdgeTri.Set(i,6,faces[j+1].p4 ); // 2o. elemento 
  j=j+2; // pois cada aresta eh dividida com apenas 2 elementos
 }
}

// procura o tamanho da menor aresta na superficie
void Model3D::setTriangleMinEdge()
{
//--------------------------------------------------
//  minEdge = 1E10; // initil minimum edge length
//  for( int i=0;i<mapEdgeTri.DimI();i++ )
//   if( minEdge>mapEdgeTri.Get(i,0) && cc.Get(mapEdgeTri.Get(i,1) == 0.5) )
//    minEdge = mapEdgeTri.Get(i,0); // setting the min edge of 3d surface
//  //minEdge = 0.0600023; // tamanho minimo da malha bubble-tube-1
//  minEdge = 0.09; // tamanho minimo da malha bubble-tube-4
//-------------------------------------------------- 
 
 // norm
 int count = 0;
 real aux  = 0;
 for( int i=0;i<mapEdgeTri.DimI();i++ )
  if( surfMesh.Marker.Get(mapEdgeTri.Get(i,1)) == 0.5 )
  {
   aux += mapEdgeTri.Get(i,0); 
   count++;
  }
 minEdge = aux/count;
 
 cout << endl;
 cout << "       ****************** " << endl;
 cout << "       |    "    << minEdge << "    |" << endl;
 cout << "       ****************** " << endl;
 cout << endl;

 minEdge = 0.1;
}

int Model3D::findEdge(int _v1,int _v2)
{
 int aux=0;
 for( int i=0;mapEdgeTri.DimI();i++ )
  if( (mapEdgeTri.Get(i,1) == _v1 && mapEdgeTri.Get(i,2) == _v2) ||
	  (mapEdgeTri.Get(i,1) == _v2 && mapEdgeTri.Get(i,2) == _v1) )
  {
   aux = i;
   break;
  }

 return aux;
}

void Model3D::insertPointsByLength()
{
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge length
  real edgeLength = mapEdgeTri.Get(i,0);
  if( surfMesh.Marker.Get(mapEdgeTri.Get(i,1)) == 0.5 && 
	  edgeLength > 1.5*minEdge )//&&
	//--------------------------------------------------
	//   ( Y.Get(mapEdgeTri.Get(i,1) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,1) != Y.Min()) || 
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Min())  ) ) 
	//-------------------------------------------------- 
  //if( cc.Get(mapEdgeTri.Get(i,1)) == 0.5 && edgeLength > 0.15 ) 
  //if( cc.Get(mapEdgeTri.Get(i,1)) == 0.5 && edgeLength > 0.157 ) 
  //if( cc.Get(mapEdgeTri.Get(i,1)) == 0.5 && edgeLength > 0.16 ) 
   insertPoint(i);
 }
}

void Model3D::surfaceTriangulator(int _v)
{
 int listSize = neighbourPoint.at(_v).size();
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 int vert0 = *point;++point;

 for( int i=0;i<listSize-2;i++ )
 {
  int vert1 = *point;++point;
  int vert2 = *point;

  // adding new element
  surfMesh.IEN.AddRow();
  int elem = surfMesh.IEN.DimI()-1;
  surfMesh.IEN.Set(elem,0,vert0);
  surfMesh.IEN.Set(elem,1,vert1);
  surfMesh.IEN.Set(elem,2,vert2);
 }
}

void Model3D::surfaceTriangulatorEarClipping(int _v)
{
 // removing last element that is equal to the first one
 neighbourPoint.at(_v).pop_back();

 int listSize = neighbourPoint.at(_v).size();
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 int vert1,vert2,vert3;

//--------------------------------------------------
//  for( int j=0;j<listSize;j++ )
//  {
//   cout << *point << " ";
//   ++point;
//  }
//  cout << endl;
//  point=plist.begin();
//-------------------------------------------------- 

 while( listSize>2 )
 {
  vert1 = *point;
  point=neighbourPoint.at(_v).erase(point); // removing 2nd vertice
  plist.push_back(vert1); // moving to last position
  vert2 = *point;
  point=neighbourPoint.at(_v).erase(point); // removing 2nd vertice
  vert3 = *point;

  // adding new element
  surfMesh.IEN.AddRow();
  int elem = surfMesh.IEN.DimI()-1;
  surfMesh.IEN.Set(elem,0,vert1);
  surfMesh.IEN.Set(elem,1,vert2);
  surfMesh.IEN.Set(elem,2,vert3);
  //cout << vert1 << " " << vert2 << " " << vert3 << endl;

  listSize--;
 }
}

// input: vertex to be deleted
// output: vector with 3 nodes to create 1 element
// description: this method check the quality of all ear triangles that
// can be created from the resulting polyhedron after deletion of vertex
// _v.
clVector Model3D::triangleQuality(int _v)
{
 // adding 2nd element to end of the list, but comparing if it is
 // already there
 // 0 1 2 3 0 --> polyhedron
 // 0 1 2 3 0 1 --> after adding this point
 // Doing so we can get all the triangles combinations:
 // 0 1 2, 1 2 3, 2 3 0, 3 0 1
 list<int> plist2 = neighbourPoint.at(_v);
 list<int>::iterator mele2=plist2.begin();
 ++mele2;
 int v2=*mele2;
 mele2=plist2.end();--mele2;
 if( v2 != *mele2  )
  neighbourPoint.at(_v).push_back(v2);

 int listSize = neighbourPoint.at(_v).size();
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 int vert1,vert2,vert3;
 real qMax = 0; 
 clVector vertexMax(3);

//--------------------------------------------------
//  for( int j=0;j<listSize;j++ )
//  {
//   cout << *point << " ";
//   ++point;
//  }
//  cout << endl;
//  point=plist.begin();
//-------------------------------------------------- 

 for( int i=0;i<listSize-2;i++ )
 {
  vert1 = *point;++point;
  vert2 = *point;++point;
  vert3 = *point;--point;

  // Triangle quality measure;
  real length12 = getLength(vert1,vert2);
  real length13 = getLength(vert1,vert3);
  real length23 = getLength(vert2,vert3);

  real semiPerimeter = 0.5*(length12+length13+length23);
  real area = getArea(vert1,vert2,vert3);
  real inRadius = area/semiPerimeter;

  /* Triangle quality measure;
   *
   *        6          r(t)       r -> in-radius
   * q = -------- * --------- 
   *      sqrt(3)      h(t)       h -> longest edge length
   *
   * Frey,P.,Borouchaki,H.:Surfacemeshevaluation.In:Intl. Mesh-ing
   * Roundtable, pp. 363:374 (1997)
   * */

  real h = length12;
  if( h<length13 )
   h = length13;
  if( h<length23 )
   h = length23;

  real q = 3.4641*inRadius/h;
  //cout << "triangle: " << i << " " << "quality: " << q << endl;

  if( qMax<q )
  {
   qMax = q; 
   vertexMax.Set(0,vert1);
   vertexMax.Set(1,vert2);
   vertexMax.Set(2,vert3);
  }
 }

 return vertexMax;
}

/*  Triangulator for the interface points after the deletion process.
 *  This method works with ear methodology, i.e., using a sort numbered
 *  polyhedron and making a triangulation considering the point n and
 *  the nexts neighbours n+1 and n+2, then removing from the polyhedron
 *  list the middle point n+1, therefore reducing the polyhedron to its
 *  next smaller shape.
 *  After the 1st triangulation by ear methodoly, the IEN mesh matrix
 *  needs to be updated
 *  */
void Model3D::surfaceTriangulatorQualityEarClipping(int _v)
{
 int listSize = neighbourPoint.at(_v).size();
 int vert1,vert2,vert3;

 while( listSize>4 )
 {
  clVector vertex = triangleQuality(_v);

  vert1 = vertex.Get(0);
  vert2 = vertex.Get(1);
  vert3 = vertex.Get(2);
  //cout << "chosen: " << vert1 << " " << vert2 << " " << vert3 << endl;
  neighbourPoint.at(_v).remove(vert2); // removing 2nd vertice
  listSize--;

  // adding new element
  surfMesh.IEN.AddRow();
  int elem = surfMesh.IEN.DimI()-1;
  surfMesh.IEN.Set(elem,0,vert1);
  surfMesh.IEN.Set(elem,1,vert2);
  surfMesh.IEN.Set(elem,2,vert3);
  surfMesh.numElems++;

 }
 // adding last remaning element
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 vert1 = *point;++point;
 vert2 = *point;++point;
 vert3 = *point;
 //cout << "chosen: " << vert1 << " " << vert2 << " " << vert3 << endl;
 // adding new element
 surfMesh.IEN.AddRow();
 int elem = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem,0,vert1);
 surfMesh.IEN.Set(elem,1,vert2);
 surfMesh.IEN.Set(elem,2,vert3);
 surfMesh.numElems++;
}

void Model3D::deleteSurfacePoint(int _v)
{
 X.Delete(_v);
 Y.Delete(_v);
 Z.Delete(_v);
 cc.Delete(_v);
 numVerts--;

 surfMesh.X.Delete(_v);
 surfMesh.Y.Delete(_v);
 surfMesh.Z.Delete(_v);
 surfMesh.Marker.Delete(_v);
 surfMesh.numVerts--;

 // updating surfMesh.IEN
 for( int i=0;i<surfMesh.IEN.DimI();i++ )
  for( int j=0;j<surfMesh.IEN.DimJ();j++ )
   if( surfMesh.IEN.Get(i,j)>_v )
	surfMesh.IEN.Set(i,j,surfMesh.IEN.Get(i,j)-1);
}

void Model3D::deleteSurfaceElementByPoint(int _v)
{
 // marking the desired elements for deletion
 list<int> plist = neighbourSurfaceElem.at(_v);
 for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
 {
  surfMesh.IEN.Set(*mele,0,-1);
  surfMesh.IEN.Set(*mele,1,-1);
  surfMesh.IEN.Set(*mele,2,-1);
 }

 // deleting elements
 for( int i=0;i<surfMesh.IEN.DimI();i++ )
  if( surfMesh.IEN.Get(i,0) == -1 && 
	  surfMesh.IEN.Get(i,1) == -1 && 
	  surfMesh.IEN.Get(i,2) == -1 )
  {
   surfMesh.IEN.DelLine(i);
   surfMesh.numElems--;
   i--; // should go back to verify the next element as well
  }
}

void Model3D::setPolyhedron(int _v)
{
 /* esta rotina nao esta otimizada; criamos uma matriz clMatrix para
 // organizar e ordenar os vertices do poliedro resultante da retirada
 // de um ponto da malha de superficie. Para isso configuramos a tal
 // matriz com as linhas representando cada aresta do poliedro sendo que
 // as 2 colunas representam o 1o. vertice e o 2o. vertice da aresta.
 // Esta aresta por sinal nao esta ordenada em algum sentido (horario ou
 // anti-horario) e para isso eh preciso organiza-la e ordenar as linhas
 // da matriz. Isto eh feito buscando o 2 vertice da 1a linha e
 // verificando nas linhas remanescentes da matriz qual delas contem o
 // mesmo vertices. Se for a 1a coluna basta substituir a 2a. linha da
 // matriz pela linha em que o elemento se encontra. Caso esteja na 2a.
 // coluna fazemos o mesmo passo acima porem invertemos tambem o vertice
 // da 1a. coluna com o da 2o coluna. Fazemos este loop ate a ultima
 // linha e suprimimos os elementos repetidos. No final temos uma matriz
 // organizada por arestas e nos.
 // Exemplo:
 //
 // test = [ 0  1 ]   -->   [ 0  1 ]
 //        [ 3  1 ]   -->   [ 1  3 ]
 //        [ 2  0 ]   -->   [ 3  2 ]
 //        [ 3  2 ]   -->   [ 2  0 ]                               
 */
  int listSize = neighbourPoint.at(_v).size();
  clMatrix test(listSize/2,2);
  list<int> plist = neighbourPoint.at(_v);
  list<int>::iterator mele=plist.begin();
  for( int i=0;i<listSize/2;i++ )
  {
   test.Set(i,0,*mele);++mele;
   test.Set(i,1,*mele);++mele;
  }
//--------------------------------------------------
//   if( ii == 309 )
//   {
//    cout << "-------- " << ii << endl;
//    test.Display();
//    cout << endl;
//   }
//-------------------------------------------------- 

  for( int k=0;k<test.DimI()-1;k++ )
  {
   int node2 = test.Get(k,1);
   //--------------------------------------------------
   // if( ii == 325 )
   //  cout << "== " << k << endl;
   //-------------------------------------------------- 
   for( int z=k+1;z<test.DimI();z++ )
   {
	int v2 = test.Get(k+1,0);
	int v3 = test.Get(k+1,1);
	int vSwap1 = test.Get(z,0);
	int vSwap2 = test.Get(z,1);
	//--------------------------------------------------
	// if( ii == 325 )
	// {
	//  cout << "---- " << node2 << " ----" << endl;
	//  cout << k+1 << " " << v2 << " " << v3 << endl;
	//  cout << z << " " << vSwap1 << " " << vSwap2 << endl;
	//  cout << "---------" << endl;
	// }
	//-------------------------------------------------- 
	if( vSwap1 == node2 )
	{
	 test.Set(k+1,0,vSwap1);
	 test.Set(k+1,1,vSwap2);
	 test.Set(z,0,v2);
	 test.Set(z,1,v3);
	 break;
	}
	if( vSwap2 == node2 )
	{
	 test.Set(k+1,0,vSwap2);
	 test.Set(k+1,1,vSwap1);
	 test.Set(z,0,v3);
	 test.Set(z,1,v2);
	 break;
	}
   }
//--------------------------------------------------
//    if( ii == 309 )
//    {
//    cout << "-------- " << ii << endl;
//    test.Display();
//    cout << endl;
//    }
//-------------------------------------------------- 
  }
  neighbourPoint.at(_v).clear();
  for( int i=0;i<test.DimI();i++ )
   for( int j=0;j<test.DimJ();j++ )
	neighbourPoint.at(_v).push_back(test.Get(i,j));

  // removing duplicated elements
  neighbourPoint.at(_v).unique();
  //--------------------------------------------------
  //   cout << "---------" << ii << "------------" << endl;
  //   std::ostream_iterator< int > output( cout, " " );
  //   std::copy( neighbourPoint.at(ii).begin(), 
  //              neighbourPoint.at(ii).end(), output );
  //   cout << endl;
  //-------------------------------------------------- 
}

void Model3D::setNeighbourSurface()
{
 // list of element neighbours
 neighbourSurfaceElem.resize (0);
 neighbourSurfaceElem.resize (surfMesh.numVerts);
 for( int i=0;i<surfMesh.numElems;i++ )
  for( int j=0;j<3;j++ )
   neighbourSurfaceElem.at( (int) surfMesh.IEN.Get(i,j) ).push_back(i);

 // list of point neighbours
 neighbourPoint.resize (0);
 neighbourPoint.resize (surfMesh.numVerts);
 for( int ii=0;ii<surfMesh.numVerts;ii++ )
 {
  list<int> plist = neighbourSurfaceElem.at(ii);
  for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
  {
   int v1 = (int) surfMesh.IEN.Get(*mele,0);
   int v2 = (int) surfMesh.IEN.Get(*mele,1);
   int v3 = (int) surfMesh.IEN.Get(*mele,2);

   neighbourPoint.at( ii ).push_back(v1);
   neighbourPoint.at( ii ).push_back(v2);
   neighbourPoint.at( ii ).push_back(v3);
   neighbourPoint.at( ii ).remove(ii);
  }
 }
}

void Model3D::flipTriangleEdge( int _edge )
{
 /* Triangle quality measure;
  *
  *        6          r(t)       r -> in-radius
  * q = -------- * --------- 
  *      sqrt(3)      h(t)       h -> longest edge length
  *
  * Frey,P.,Borouchaki,H.:Surfacemeshevaluation.In:Intl. Mesh-ing
  * Roundtable, pp. 363:374 (1997)
  * */

 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  _edge = i;
  mapEdgeTri.Get(_edge,0); // length
  int v1 = mapEdgeTri.Get(_edge,1); // v1
  int v2 = mapEdgeTri.Get(_edge,2); // v2
  int v3elem1 = mapEdgeTri.Get(_edge,3); // v3elem1
  int v3elem2 = mapEdgeTri.Get(_edge,4); // v3elem2
  int elem1 = mapEdgeTri.Get(_edge,5); // elem1
  int elem2 = mapEdgeTri.Get(_edge,6); // elem2

  real length12 = getLength(v1,v2);
  real length13_1 = getLength(v1,v3elem1);
  real length13_2 = getLength(v1,v3elem2);
  real length23_1 = getLength(v2,v3elem1);
  real length23_2 = getLength(v2,v3elem2);
  real length3_1_3_2 = getLength(v3elem1,v3elem2);

  // elem1
  real semiPerimeter1 = 0.5*(length12+length13_1+length23_1);
  real area1 = getArea(v1,v2,v3elem1);
  real inRadius1 = area1/semiPerimeter1;
  real h1 = length12;
  if( h1<length13_1 )
   h1 = length13_1;
  if( h1<length23_1 )
   h1 = length23_1;
  real q1 = 3.4641*inRadius1/h1;

  // elem2
  real semiPerimeter2 = 0.5*(length12+length13_2+length23_2);
  real area2 = getArea(v1,v2,v3elem2);
  real inRadius2 = area2/semiPerimeter2;
  real h2 = length12;
  if( h2<length13_2 )
   h2 = length13_2;
  if( h2<length23_2 )
   h2 = length23_2;
  real q2 = 3.4641*inRadius2/h2;

  // elem3
  real semiPerimeter3 = 0.5*(length3_1_3_2+length13_1+length13_2);
  real area3 = getArea(v1,v3elem1,v3elem2);
  real inRadius3 = area3/semiPerimeter3;
  real h3 = length3_1_3_2;
  if( h3<length13_1 )
   h3 = length13_1;
  if( h3<length13_2 )
   h3 = length13_2;
  real q3 = 3.4641*inRadius3/h3;

  // elem4
  real semiPerimeter4 = 0.5*(length3_1_3_2+length23_1+length23_2);
  real area4 = getArea(v2,v3elem1,v3elem2);
  real inRadius4 = area4/semiPerimeter4;
  real h4 = length3_1_3_2;
  if( h4<length23_1 )
   h4 = length23_1;
  if( h4<length23_2 )
   h4 = length23_2;
  real q4 = 3.4641*inRadius4/h4;

  // this works, but is not consistent!!! CHANGE IT SOON!
  if( surfMesh.Marker.Get(v1)==0.5 && 
	  q1+q2 < q3+q4 && 
	  area1+area2>=area3+area4 ) //&&
	//--------------------------------------------------
	//   ( Y.Get(mapEdgeTri.Get(i,1) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,1) != Y.Min()) || 
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Min())  ) ) 
	//-------------------------------------------------- 
  {
   //cout << area1+area2 << " " << area3+area4 << endl;
   //cout << q1 << " " << q2 << " " << q3 << " " << q4 << endl;
   cout << "----------------- " << color(none,green,black) 
	    << "flipping edge: " << resetColor() 
		<< v1 << " " << v2 
		<< color(none,green,black) 
		<< " --> " << resetColor()
		<< v3elem1 << " " << v3elem2 << endl;

   surfMesh.IEN.Set(elem1,0,v1);
   surfMesh.IEN.Set(elem1,1,v3elem1);
   surfMesh.IEN.Set(elem1,2,v3elem2);

   surfMesh.IEN.Set(elem2,0,v2);
   surfMesh.IEN.Set(elem2,1,v3elem1);
   surfMesh.IEN.Set(elem2,2,v3elem2);
   setTriEdge();
   setNeighbourSurface();
  }
 }
}

void Model3D::insertPoint(int _edge)
{
 int vAdd = surfMesh.numVerts; // aditional vertice
 cout << "-------------- " << color(none,yellow,black) << "inserting vertex: "
      << resetColor() << vAdd << endl;
 //saveVTKSurface("./vtk/","insertBefore",vAdd);

 // edge vertices
 int v1 = mapEdgeTri.Get(_edge,1);
 int v2 = mapEdgeTri.Get(_edge,2);
 int v3elem1 = mapEdgeTri.Get(_edge,3);
 int v3elem2 = mapEdgeTri.Get(_edge,4);

 // elements
 int elem1 = mapEdgeTri.Get(_edge,5);
 int elem2 = mapEdgeTri.Get(_edge,6);

 // add point in the middle of a edge 
 real XvAdd = surfMesh.X.Get(v1)+( surfMesh.X.Get(v2)-surfMesh.X.Get(v1) )*0.5;
 real YvAdd = surfMesh.Y.Get(v1)+( surfMesh.Y.Get(v2)-surfMesh.Y.Get(v1) )*0.5;
 real ZvAdd = surfMesh.Z.Get(v1)+( surfMesh.Z.Get(v2)-surfMesh.Z.Get(v1) )*0.5;

 // insert aditional vertice coordinate
 X.AddItem(vAdd,XvAdd);
 Y.AddItem(vAdd,YvAdd);
 Z.AddItem(vAdd,ZvAdd);

 surfMesh.X.AddItem(XvAdd);
 surfMesh.Y.AddItem(YvAdd);
 surfMesh.Z.AddItem(ZvAdd);
 surfMesh.Marker.AddItem(0.5); // interface set up

 // incremeting the number of points
 surfMesh.numVerts++;
 numVerts++;

 /* by adding 1 point on the edge it is necessary to divide the
  * original element and also the oposite element by 2, becoming 4
  * elements in total. */

 // 1st. new element (v1 - vAdd - v3elem1) 
 // on the same position of the OLD 1st. element (v1 - v2 - v3elem1)
 surfMesh.IEN.Set(elem1,0,v1);
 surfMesh.IEN.Set(elem1,1,vAdd);
 surfMesh.IEN.Set(elem1,2,v3elem1);

 // 2nd. new element (v1 - vAdd - v3elem2) 
 // on the same position of the OLD 2nd. element (v1 - v2 - v3elem2)
 surfMesh.IEN.Set(elem2,0,v1);
 surfMesh.IEN.Set(elem2,1,vAdd);
 surfMesh.IEN.Set(elem2,2,v3elem2);

 // 3rd. new element (v2 - vAdd - v3elem1) on the last row
 surfMesh.IEN.AddRow();
 int elem3 = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem3,0,v2);
 surfMesh.IEN.Set(elem3,1,vAdd);
 surfMesh.IEN.Set(elem3,2,v3elem1);
 surfMesh.numElems++;

 // 4th. new element (v2 - vAdd - v3elem2) on the last row
 surfMesh.IEN.AddRow();
 int elem4 = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem4,0,v2);
 surfMesh.IEN.Set(elem4,1,vAdd);
 surfMesh.IEN.Set(elem4,2,v3elem2);
 surfMesh.numElems++;
 
 setTriEdge();
 setNeighbourSurface();

//--------------------------------------------------
//  /* ********************************************************************* */
//  /* updating mapEdgeTri */
//  /* eh necessario atualizar as 4 arestas dos 2 triangulos trabalhados */
//  int lastRow;
//  real length;
// 
//  // 1st. new edge on the same place as the old edge
//  length = getLength(v1,vAdd);
//  mapEdgeTri.Set(_edge,0,length);
//  mapEdgeTri.Set(_edge,1,v1);
//  mapEdgeTri.Set(_edge,2,vAdd);
//  mapEdgeTri.Set(_edge,3,v3elem1);
//  mapEdgeTri.Set(_edge,4,v3elem2);
//  mapEdgeTri.Set(_edge,5,elem1);
//  mapEdgeTri.Set(_edge,6,elem2);
// 
//  // 2nd. new edge on the end of the mapEdgeTri matrix
//  mapEdgeTri.AddRow();
//  lastRow = mapEdgeTri.DimI()-1;
//  length = getLength(v2,vAdd);
//  mapEdgeTri.Set(lastRow,0,length);
//  mapEdgeTri.Set(lastRow,1,v2);
//  mapEdgeTri.Set(lastRow,2,vAdd);
//  mapEdgeTri.Set(lastRow,3,v3elem1);
//  mapEdgeTri.Set(lastRow,4,v3elem2);
//  mapEdgeTri.Set(lastRow,5,elem3); 
//  mapEdgeTri.Set(lastRow,6,elem4); // last row element
// 
//  // 3rd. new edge on the end of the mapEdgeTri matrix
//  length = getLength(v3elem1,vAdd);
//  mapEdgeTri.AddRow();
//  lastRow = mapEdgeTri.DimI()-1;
//  mapEdgeTri.Set(lastRow,0,length);
//  mapEdgeTri.Set(lastRow,1,v3elem1);
//  mapEdgeTri.Set(lastRow,2,vAdd);
//  mapEdgeTri.Set(lastRow,3,v1);
//  mapEdgeTri.Set(lastRow,4,v2);
//  mapEdgeTri.Set(lastRow,5,elem1); 
//  mapEdgeTri.Set(lastRow,6,elem3); // last row element
// 
//  // 4th. new edge on the end of the mapEdgeTri matrix
//  length = getLength(v3elem2,vAdd);
//  mapEdgeTri.AddRow();
//  lastRow = mapEdgeTri.DimI()-1;
//  mapEdgeTri.Set(lastRow,0,length);
//  mapEdgeTri.Set(lastRow,1,v3elem2);
//  mapEdgeTri.Set(lastRow,2,vAdd);
//  mapEdgeTri.Set(lastRow,3,v1);
//  mapEdgeTri.Set(lastRow,4,v2);
//  mapEdgeTri.Set(lastRow,5,elem2); 
//  mapEdgeTri.Set(lastRow,6,elem4); // last row element
// 
//  // 1st. updated edge 
//  int edgeUpdate = findEdge(v1,v3elem1);
//  if( mapEdgeTri.Get(edgeUpdate,5) == elem1 )
//  {
//   mapEdgeTri.Set(edgeUpdate,3,vAdd);
//   mapEdgeTri.Set(edgeUpdate,5,elem1);
//  }
//  else
//  {
//   mapEdgeTri.Set(edgeUpdate,4,vAdd);
//   mapEdgeTri.Set(edgeUpdate,6,elem1);
//  }
// 
//  // 2nd. updated edge 
//  edgeUpdate = findEdge(v1,v3elem2);
//  if( mapEdgeTri.Get(edgeUpdate,5) == elem2 )
//  {
//   mapEdgeTri.Set(edgeUpdate,3,vAdd);
//   mapEdgeTri.Set(edgeUpdate,5,elem2);
//  }
//  else
//  {
//   mapEdgeTri.Set(edgeUpdate,4,vAdd);
//   mapEdgeTri.Set(edgeUpdate,6,elem2);
//  }
// 
//  // 3rd. updated edge 
//  edgeUpdate = findEdge(v2,v3elem1);
//  if( mapEdgeTri.Get(edgeUpdate,5) == elem1 )
//  {
//   mapEdgeTri.Set(edgeUpdate,3,vAdd);
//   mapEdgeTri.Set(edgeUpdate,5,elem3);
//  }
//  else
//  {
//   mapEdgeTri.Set(edgeUpdate,4,vAdd);
//   mapEdgeTri.Set(edgeUpdate,6,elem3);
//  }
// 
//  // 4th. updated edge 
//  edgeUpdate = findEdge(v2,v3elem2);
//  if( mapEdgeTri.Get(edgeUpdate,5) == elem2 )
//  {
//   mapEdgeTri.Set(edgeUpdate,3,vAdd);
//   mapEdgeTri.Set(edgeUpdate,5,elem4);
//  }
//  else
//  {
//   mapEdgeTri.Set(edgeUpdate,4,vAdd);
//   mapEdgeTri.Set(edgeUpdate,6,elem4);
//  }
//  /* _____________________________________________________________________ */
// 
//  /* ********************************************************************* */
//  /* updating neighbourSurfaceElem */
//  neighbourSurfaceElem.at(v1).push_back(elem2);
//  neighbourSurfaceElem.at(v2).remove(elem2);
//  neighbourSurfaceElem.at(v2).push_back(elem3);
//  neighbourSurfaceElem.at(v2).push_back(elem4);
//  neighbourSurfaceElem.at(v3elem1).remove(elem2);
//  neighbourSurfaceElem.at(v3elem1).push_back(elem3);
//  neighbourSurfaceElem.at(v3elem2).remove(elem1);
//  neighbourSurfaceElem.at(v3elem2).push_back(elem2);
//  neighbourSurfaceElem.at(v3elem2).push_back(elem4);
//  list<int> myNewList;myNewList.resize(0);
//  myNewList.push_back(elem1);
//  myNewList.push_back(elem2);
//  myNewList.push_back(elem3);
//  myNewList.push_back(elem4);
//  neighbourSurfaceElem.push_back(myNewList);
// 
//  /* updating neighbourPoint */
//  // o problema neste update eh que neighbourPoint precisa ter uma
//  // arrumacao especial para o uso seguido de setPolyhedron. Isso quer
//  // dizer que os vertices tem que estar duplicados e nao pode entao ser
//  // atualizado desta forma abaixo.
//  
// //--------------------------------------------------
// //  // point v1
// //  neighbourPoint.at(v1).push_back(vAdd);
// //  neighbourPoint.at(v1).remove(v2);
// //  // point v2
// //  neighbourPoint.at(v2).push_back(vAdd);
// //  neighbourPoint.at(v2).remove(v1);
// //  // point v3elem1
// //  neighbourPoint.at(v3elem1).push_back(vAdd);
// //  // point v3elem2
// //  neighbourPoint.at(v3elem2).push_back(vAdd);
// //  // add vAdd
// //  list<int> myNewList;myNewList.resize(0);
// //  myNewList.push_back(v1);
// //  myNewList.push_back(v3elem1);
// //  myNewList.push_back(v2);
// //  myNewList.push_back(v3elem2);
// //  neighbourPoint.push_back(myNewList);
// //-------------------------------------------------- 
// 
//  /* *********************************** */
//  // provisory update
//  int lastLine = surfMesh.numVerts-1;
//  neighbourPoint.resize(surfMesh.numVerts);
//  list<int> plist = neighbourSurfaceElem.at(lastLine);
//  for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
//  {
//   int v1 = (int) surfMesh.IEN.Get(*mele,0);
//   int v2 = (int) surfMesh.IEN.Get(*mele,1);
//   int v3 = (int) surfMesh.IEN.Get(*mele,2);
// 
//   neighbourPoint.at( lastLine ).push_back(v1);
//   neighbourPoint.at( lastLine ).push_back(v2);
//   neighbourPoint.at( lastLine ).push_back(v3);
//   neighbourPoint.at( lastLine ).remove(lastLine);
//  }
//  /* *********************************** */
//  /* _____________________________________________________________________ */
//-------------------------------------------------- 
}

void Model3D::deletePoint(int _v)
{
 cout << "--------------- " << color(none,red,black) << "removing vertex: "
      << resetColor() << _v << endl;
 //saveVTKSurface("./vtk/","deleteBefore",_v);

 // deleting elements
 deleteSurfaceElementByPoint(_v);

 // after the deletion process it's necessary to create new elements
 // to fullfill the space left by the deletion of elements
 setPolyhedron(_v);
 //surfaceTriangulator(_v);
 //surfaceTriangulatorEarClipping(_v);
 surfaceTriangulatorQualityEarClipping(_v);

 // deleting X,Y and Z coordinate; deleting the point maker funcition
 deleteSurfacePoint(_v);

 // updating edge matrix
 setTriEdge();
 // updating surface neighbours
 setNeighbourSurface();

//--------------------------------------------------
//  // 2nd. updated edge 
//  edgeUpdate = findEdge(v1,v3elem2);
//  if( mapEdgeTri.Get(edgeUpdate,5) == elem2 )
//  {
//   mapEdgeTri.Set(edgeUpdate,3,vAdd);
//   mapEdgeTri.Set(edgeUpdate,5,elem2);
//  }
//  else
//  {
//   mapEdgeTri.Set(edgeUpdate,4,vAdd);
//   mapEdgeTri.Set(edgeUpdate,6,elem2);
//  }
//-------------------------------------------------- 

 //saveVTKSurface("./vtk/","deleteAfter",_v);
}

void Model3D::removePointsByLength()
{
 real test = 0.5*minEdge; // 50% of minEdge
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge vertices
  int v1 = mapEdgeTri.Get(i,1);
  int v2 = mapEdgeTri.Get(i,2);

  // verifying the length of each surface edge
  if( surfMesh.Marker.Get(v1) == 0.5 && mapEdgeTri.Get(i,0) < test ) //&&
	//--------------------------------------------------
	//   ( Y.Get(mapEdgeTri.Get(i,1) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,1) != Y.Min()) || 
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Max()) ||
	//     Y.Get(mapEdgeTri.Get(i,2) != Y.Min())  ) ) 
	//-------------------------------------------------- 
  {
   // sum of all neighbour edge length of the 1st. point
   real sumLength1=0;
   list<int> plist = neighbourPoint.at(v1);
   for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
   {
	// node to be removed
	int vert1 = v1;
	// oposite node
	int vert2 = *mele;
	sumLength1 += sqrt( ( X.Get(vert1)-X.Get(vert2) )*
	                    ( X.Get(vert1)-X.Get(vert2) )+
						( Y.Get(vert1)-Y.Get(vert2) )*
						( Y.Get(vert1)-Y.Get(vert2) )+
						( Z.Get(vert1)-Z.Get(vert2) )*
						( Z.Get(vert1)-Z.Get(vert2) ) );
   }

   // sum of all neighbour edge length of the 1st. point
   real sumLength2=0;
   plist = neighbourPoint.at(v2);
   for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
   {
	// node to be removed
	int vert1 = v1;
	// oposite node
	int vert2 = *mele;
	sumLength2 += sqrt( ( X.Get(vert1)-X.Get(vert2) )*
	                    ( X.Get(vert1)-X.Get(vert2) )+
						( Y.Get(vert1)-Y.Get(vert2) )*
						( Y.Get(vert1)-Y.Get(vert2) )+
						( Z.Get(vert1)-Z.Get(vert2) )*
						( Z.Get(vert1)-Z.Get(vert2) ) );
   }

   // check which node has the smallest length sum and proceed
   if( sumLength1 < sumLength2 )
	deletePoint(v1);
   else // if the 2nd. node has the smallest edge length sum
	deletePoint(v2);
  }
 }
}

void Model3D::removePointsByInterfaceDistance()
{
 clVector surfaceAux = cc==0.5;
 clVector surface = surfaceAux.Find();

 clVector xSurface( surface.Dim() );
 clVector ySurface( surface.Dim() );
 clVector zSurface( surface.Dim() );
 for( int i=0;i<surface.Dim();i++ )
 {
  real aux = surface.Get(i);
  xSurface.Set(i,X.Get( aux ));
  ySurface.Set(i,Y.Get( aux ));
  zSurface.Set(i,Z.Get( aux ));
 }

 clVector closer = dsearchn(xSurface,ySurface,zSurface,X,Y,Z);
 clVector xCloser( closer.Dim() );
 clVector yCloser( closer.Dim() );
 clVector zCloser( closer.Dim() );
 for( int i=0;i<closer.Dim();i++ )
 {
  real aux = closer.Get(i);
  closer.Set(i,surface.Get(aux)); // alterando os valores de closer(i)
  xCloser.Set(i,X.Get( closer.Get(i) ));
  yCloser.Set(i,Y.Get( closer.Get(i) ));
  zCloser.Set(i,Z.Get( closer.Get(i) ));
 }

 clVector distance(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  real aux = sqrt( (X.Get(i)-xCloser.Get(i))*(X.Get(i)-xCloser.Get(i))+
                   (Y.Get(i)-yCloser.Get(i))*(Y.Get(i)-yCloser.Get(i))+
                   (Z.Get(i)-zCloser.Get(i))*(Z.Get(i)-zCloser.Get(i)) );
  distance.Set(i,fabs(aux));
 }
 for( int i=0;i<numVerts;i++ )
 {
  real d = distance.Get(i);
  //if( d>0 && d<0.4*minEdge ) // mainBubble.cpp
  if( d>0 && d<0.2*minEdge )
  {
   cout << "--- " << color(none,red,black) << "removing vertex by distance: "
	    << resetColor() << i << endl;

   X.Delete(i);
   Y.Delete(i);
   Z.Delete(i);
   cc.Delete(i);
   distance.Delete(i);
   numVerts--;
   i--;
  }
 }
}

void Model3D::breakup()
{
 for( int i=0;i<IEN.DimI();i++ )
 {
  int v1 = IEN.Get(i,0);
  int v2 = IEN.Get(i,1);
  int v3 = IEN.Get(i,2);
  int v4 = IEN.Get(i,3);
  if( cc.Get(v1) == 0.5 && cc.Get(v2) == 0.5 &&
	  cc.Get(v3) == 0.5 && cc.Get(v4) == 0.5 )
  {
   cc.Set(v1,0.0);
   cc.Set(v2,0.0);
   cc.Set(v3,0.0);
   cc.Set(v4,0.0);
  }
 }
}

void Model3D::insertPointsByArea()
{
 int lastRow;
 real test = 0.008;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  if( surfMesh.Marker.Get(v1) == 0.5 && getArea(i) > test )
  {
   int v4 = surfMesh.numVerts;

   // coords of triangle centroid
   real centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3) )*0.3334;
   real centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3) )*0.3334;
   real centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3) )*0.3334;

   X.AddItem(v4,centroidX);
   Y.AddItem(v4,centroidY);
   Z.AddItem(v4,centroidZ);
   cc.AddItem(v4,0.5); // interface set up

   surfMesh.X.AddItem(v4,centroidX);
   surfMesh.Y.AddItem(v4,centroidY);
   surfMesh.Z.AddItem(v4,centroidZ);
   surfMesh.Marker.AddItem(v4,0.5);
   surfMesh.numVerts++;

   /* substitutes the original element to 3 smaller elements adding 1 
	* point on the centroid */

   // 1st. element on the same position of previous one
   surfMesh.IEN.Set(i,0,v1);
   surfMesh.IEN.Set(i,1,v2);
   surfMesh.IEN.Set(i,2,v4);

   // 2nd. element
   surfMesh.IEN.AddRow();
   lastRow = surfMesh.IEN.DimI()-1;
   surfMesh.IEN.Set(lastRow,0,v1);
   surfMesh.IEN.Set(lastRow,1,v3);
   surfMesh.IEN.Set(lastRow,2,v4);
   surfMesh.numElems++;

   // 3rd. element
   surfMesh.IEN.AddRow();
   lastRow = surfMesh.IEN.DimI()-1;
   surfMesh.IEN.Set(lastRow,0,v2);
   surfMesh.IEN.Set(lastRow,1,v3);
   surfMesh.IEN.Set(lastRow,2,v4);
   surfMesh.numElems++;
  }
 }
}
 
/* This method re-mesh completly the domain preserving only the points
 * located at surface and convex-hull. To do so, surfMesh.numVerts and
 * surfMesh.IEN need to be set on the beginning of the running program,
 * usually when the mesh is created from the .MSH file */
void Model3D::mesh2Dto3DOriginal()
{
//--------------------------------------------------
//  saveVTKSurface("./vtk/","before",0);
//  insertPointsByLength();
//  saveVTKSurface("./vtk/","between",0);
//  removePointsByLength();
//  saveVTKSurface("./vtk/","flipBetween",0);
//  flipTriangleEdge(0);
//  saveVTKSurface("./vtk/","after",0);
//-------------------------------------------------- 

 int ny = 4;
 int nPoints = 20;

 // tetgen mesh object 
 tetgenio in,out;
 in.mesh_dim = 3;
 //in.numberofpoints = surfMesh.numVerts;
 in.numberofpoints = surfMesh.numVerts+(nPoints*nPoints*ny);
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 // add to tetgen struct the point coords
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  in.pointlist[3*i+0] = surfMesh.X.Get(i);
  in.pointlist[3*i+1] = surfMesh.Y.Get(i);
  in.pointlist[3*i+2] = surfMesh.Z.Get(i);
  if( surfMesh.Marker.Get(i) == 0.0 )
   in.pointmarkerlist[i] = 11;
  if( surfMesh.Marker.Get(i) == 0.5 )
   in.pointmarkerlist[i] = 22; // same id of facetmarker
 }

 // ******************************************** //
 // strategy to ADD points - to be implemented - //
 // ******************************************** //
 
 real Ymax1=100;
 real Ymin1=-100;
 real Ymax2=-100;
 real Ymin2=100;
 
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  // bubble 1 (Y<0)
  if( Y.Get(i) < 0 && cc.Get(i)==0.5 )
  {
   if(Y.Get(i)>Ymin1) Ymin1=Y.Get(i);
   if(Y.Get(i)<Ymax1) Ymax1=Y.Get(i);
  }
  // bubble 2 (Y>0)
  if( Y.Get(i) > 0 && cc.Get(i)==0.5 )
  {
   if(Y.Get(i)<Ymin2) Ymin2=Y.Get(i);
   if(Y.Get(i)>Ymax2) Ymax2=Y.Get(i);
  }
 }

 real xi = -0.80;
 real yi = Ymin1;
 real zi = -0.40;
 int count = surfMesh.numVerts;
 for( int i=0;i<nPoints;i++ )
 {
  for( int j=0;j<nPoints;j++ )
  {
   for( int k=1;k<(ny+1);k++ )
   {
	real dx = (-2.0*xi)/(nPoints-1);
	in.pointlist[3*count+0] = xi + dx*i;

	real dy = (Ymin2-Ymin1)/(ny+1);
	in.pointlist[3*count+1] = yi + dy*k;

	real dz = (-2.0*zi)/(nPoints-1);
	in.pointlist[3*count+2] = zi + dz*j;

	in.pointmarkerlist[count] = 11;

	count++;
   }
  }
 }

//--------------------------------------------------
//  for( int i=0;i<in.numberofpoints;i++ )
//  {
//   if( (in.pointlist[3*i+0]*in.pointlist[3*i+0] ) < 0.0 + EPS && 
//       (in.pointlist[3*i+1]*in.pointlist[3*i+1] ) < 0.0 + EPS && 
//       (in.pointlist[3*i+2]*in.pointlist[3*i+2] ) < 0.0 + EPS )
//    cout << i << " " 
// 	    << in.pointlist[3*i+0] << " " 
// 	    << in.pointlist[3*i+1] << " " 
// 	    << in.pointlist[3*i+2] << " "
// 	    << in.pointmarkerlist[i] << endl;
//  }
//-------------------------------------------------- 



//--------------------------------------------------
 //int j=1;
//  for( int i=surfMesh.numVerts;i<surfMesh.numVerts+n;i++ )
//  {
//   in.pointlist[3*i+0] = 0.0;
//   in.pointlist[3*i+1] = y0+j*deltaY;
//   in.pointlist[3*i+2] = 0.0 ;
//   in.pointmarkerlist[i] = 11;
//   j++;
//  }
//-------------------------------------------------- 


 /* This procedure defines regions on the mesh and after
  * insertion/deletion of points by the tetgen program we can easily
  * recognize the point locations by the regionlist array and then we
  * can define the points relied inside the bubble (1.0), on the surface
  * (0.5) and outside the bubble (0.0). To do so it is necessery to
  * define AT LEAST one region, preferably the one between the interface 
  * and convex hull. 
  * To recovery the information of different zones, after the meshing
  * procedure is necessary to look at tetrahedronattributelist 
  * */
 in.numberofregions = 1; 
 in.regionlist = new REAL[in.numberofregions*4];

 // definindo regiao fora da bolha 
 //in.regionlist[0] = X.Min();
 //in.regionlist[1] = Y.Min();
 //in.regionlist[2] = Z.Min();
//--------------------------------------------------
//  in.regionlist[0] = 0.1;
//  in.regionlist[1] = 0.1;
//  in.regionlist[2] = 0.1;
//-------------------------------------------------- 
 in.regionlist[0] = -5.8;
 in.regionlist[1] = -2.8;
 in.regionlist[2] = -2.8;
 in.regionlist[3] = 1;

 tetgenio::facet *f;   // Define a pointer of facet. 
 tetgenio::polygon *p; // Define a pointer of polygon.
 in.numberoffacets = surfMesh.numElems; 
 in.facetlist = new tetgenio::facet[in.numberoffacets]; 
 in.facetmarkerlist = new int[in.numberoffacets];
 //in.trifacemarkerlist = new int[in.numberoffacets];

 // defining the interface and convex-hull
 // definindo a superficie da bolha e convex-hull
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  f = &in.facetlist[i];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0; 
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 3;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = v1; 
  p->vertexlist[1] = v2; 
  p->vertexlist[2] = v3;
  // melhorar esta configuracao de facet para bolha e convex hull
  if( surfMesh.Marker.Get(v1) + 
	  surfMesh.Marker.Get(v2) + 
	  surfMesh.Marker.Get(v3) > 0 )
   in.facetmarkerlist[i] = 10;
  else
   in.facetmarkerlist[i] = 20;
 }

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");
 cout << "numElems IN = " << surfMesh.numElems << endl;
 cout << "numNodes IN = " << surfMesh.numNodes << endl;
 cout << "numVerts IN = " << in.numberofpoints << endl;

 saveVTKSurface("./vtk/","test1",0);

 cout << endl;
 cout << "----> complete re-meshing the domain... ";
 //tetrahedralize( (char*) "VYYApq1.4241",&in,&out );
 //tetrahedralize( (char*) "QYYCApq1.4241a0.05",&in,&out );
 tetrahedralize( (char*) "QYYCApq1.4241a0.1",&in,&out );
 //tetrahedralize( (char*) "QYYCApq1.414q10a0.001",&in,&out );
 cout << "finished <---- " << endl;;
 cout << endl;

 numElems = out.numberoftetrahedra;
 numNodes = out.numberofpoints+numElems;
 numVerts = out.numberofpoints;
 cout << "numElems OUT = " << out.numberoftetrahedra << endl;
 cout << "numNodes OUT = " << out.numberofpoints+numElems << endl;
 cout << "numVerts OUT = " << out.numberofpoints << endl;
 cout << "numfacets OUT = " << out.numberoftrifaces << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 // varre lista de elementos e passa para estrutura IEN
 IEN.Dim(numElems,5);
 cc.Dim(numVerts);
 for( int i=0;i<out.numberoftetrahedra;i++ )
 {
  // setting de cc = 0 para fora da bolha e cc = 0.5 para interface
  if( out.tetrahedronattributelist[i] == 1 )
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,0.0);
   }
  }
  // setting de cc = 1 para dentro da bolha e cc = 0.5 para interface
  else 
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,1.0);
   }
  }
 }

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 X.Dim(numNodes);
 Y.Dim(numNodes);
 Z.Dim(numNodes);
 for( int i=0;i<numVerts;i++ )
 {
  X.Set(i,out.pointlist[3*i+0]);
  Y.Set(i,out.pointlist[3*i+1]);
  Z.Set(i,out.pointlist[3*i+2]);
  if( out.pointmarkerlist[i] == 10 ||
	  out.pointmarkerlist[i] == 22 )
   cc.Set(i,0.5);
 }
//--------------------------------------------------
//  for( int i=0;i<out.numberofpoints;i++ )
//  {
//   if( (out.pointlist[3*i+0]*out.pointlist[3*i+0] ) < 0.0 + EPS && 
//       (out.pointlist[3*i+1]*out.pointlist[3*i+1] ) < 0.0 + EPS && 
//       (out.pointlist[3*i+2]*out.pointlist[3*i+2] ) < 0.0 + EPS )
//    cout << i << " " 
// 	    << out.pointlist[3*i+0] << " " 
// 	    << out.pointlist[3*i+1] << " " 
// 	    << out.pointlist[3*i+2] << endl;
//  }
//-------------------------------------------------- 
 saveVTKSurface("./vtk/","test2",0);
}

void Model3D::mesh3DPoints()
{
 saveVTKSurface("./vtk/","before",0);
 insertPointsByLength();
 removePointsByLength();
 saveVTKSurface("./vtk/","between",0);
 flipTriangleEdge(0);
 saveVTKSurface("./vtk/","after",0);
 removePointsByInterfaceDistance();

 // cria objeto de malha do tetgen
 tetgenio in,mid,out;
 in.mesh_dim = 3;
 in.numberofpoints = numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 /* ------------ pontos da malha separados em 2 loops ------------ */
 // adiciona na estrutura tetgen as coordenadas dos pontos da 
 // superficie e do convex-hull
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  in.pointlist[3*i+0] = surfMesh.X.Get(i); // surfMesh.X nao esta no ALE
  in.pointlist[3*i+1] = surfMesh.Y.Get(i); // surfMesh.Z nao esta no ALE
  in.pointlist[3*i+2] = surfMesh.Z.Get(i); // surfMesh.Y nao esta no ALE
  if( surfMesh.Marker.Get(i) == 0.0 )
   in.pointmarkerlist[i] = 11;
  if( surfMesh.Marker.Get(i) == 0.5 )
   in.pointmarkerlist[i] = 22;
 }
 
 // adicionando pontos que nao sao da interface e do convex-hull
 for( int i=surfMesh.numVerts;i<numVerts;i++ )
 {
  in.pointlist[3*i+0] = X.Get(i);
  in.pointlist[3*i+1] = Y.Get(i);
  in.pointlist[3*i+2] = Z.Get(i);
  if( cc.Get(i) == 0.0 ) // fora da bolha
   in.pointmarkerlist[i] = 11;
  if( cc.Get(i) == 0.5 ) // na interface
   in.pointmarkerlist[i] = 22;
  if( cc.Get(i) == 1.0 ) // dentro da bolha
   in.pointmarkerlist[i] = 33;
 }
 /* -------------------------------------------------------------- */

 /* as regioes do dominio (dentro e fora da bolha) sao definidas usando o 
  * flag AA em  tetrahedralize. O tetgen marca 1 para fora da bolha e 2
  * para dentro da bolha. Eh importante notar que esta marcacao eh em
  * nivel de elementos e NAO de pontos. Para definir os pontos da
  * superficie e com isso a funcao marcadora corretamente eh necessario
  * utilizar uma funcao marcadora de pontos do tetgen: pointmarkerlist.
  * --> este flag foi testado e nao pode ser usado sem a definicao
  *  explicita da regiao externa a bolha (regionlist). O resultado foi
  *  que de um passo para o outro ele inverte o tag da regiao, i.e.,
  *  define a regiao dentro da bolha com 1 e depois de 60 passos define
  *  a mesma regiao com 2.
  */
 /* ESTE PROCEDIMENTO DEFINE REGIOES NA MALHA E APOS A INSERCAO/RETIRADA
  * DE PONTOS PELO TETGEN, CONSEGUIMOS RECONHECER A LOCALIZACAO DOS
  * PONTOS E ASSIM PODEMOS DEFINIR NOVAMENTE A FUNCAO MARCADORA COMO
  * SENDO 1.0 DENTRO DA BOLHA, 0.5 NA SUPERFICIE E 0.0 FORA 
  * E NECESSARIO DEFINIR 1 PONTO EM CADA REGIAO */
 // fluido interior + fluido exterior + superficie
 in.numberofregions = 1; 
 in.regionlist = new REAL[in.numberofregions*4];

 // fora da bolha
//--------------------------------------------------
//  in.regionlist[0] = X.Min();
//  in.regionlist[1] = Y.Min();
//  in.regionlist[2] = Z.Min();
//-------------------------------------------------- 
//--------------------------------------------------
//  in.regionlist[0] = 0.1;
//  in.regionlist[1] = 0.1;
//  in.regionlist[2] = 0.1;
//-------------------------------------------------- 
 in.regionlist[0] = -5.5;
 in.regionlist[1] = 0.0;
 in.regionlist[2] = 2.5;
 in.regionlist[3] = 1;

 tetgenio::facet *f;   // Define a pointer of facet. 
 tetgenio::polygon *p; // Define a pointer of polygon.
 in.numberoffacets = surfMesh.numElems; 
 in.facetlist = new tetgenio::facet[in.numberoffacets]; 
 in.facetmarkerlist = new int[in.numberoffacets];
 //in.trifacemarkerlist = new int[in.numberoffacets];

 // definindo a superficie da bolha e convex-hull
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  f = &in.facetlist[i];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0; 
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 3;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = v1; 
  p->vertexlist[1] = v2; 
  p->vertexlist[2] = v3;
  // melhorar esta configuracao de facet para bolha e convex hull
  if( surfMesh.Marker.Get(v1) + 
	  surfMesh.Marker.Get(v2) + 
	  surfMesh.Marker.Get(v3) > 0 )
   in.facetmarkerlist[i] = 10;
  else
   in.facetmarkerlist[i] = 20;
 }

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");
 cout << endl;
 cout << "numElems IN = " << surfMesh.numElems << endl;
 cout << "numNodes IN = " << in.numberoffacets+in.numberofpoints << endl;
 cout << "numVerts IN = " << in.numberofpoints << endl;

 cout << endl;
 cout << "----> re-meshing 3D points... ";
 tetrahedralize( (char*) "QYYRCAipq1.414q5a0.5",&in,&out );
 //tetrahedralize( (char*) "VYYRCAipq1.414q5a1.5",&in,&out );
 cout << "finished <---- " << endl;;
 cout << endl;

 numElems = out.numberoftetrahedra;
 numNodes = out.numberofpoints+out.numberoftetrahedra;
 numVerts = out.numberofpoints;
 cout << "numElems OUT = " << out.numberoftetrahedra << endl;
 cout << "numNodes OUT = " << out.numberofpoints+numElems << endl;
 cout << "numVerts OUT = " << out.numberofpoints << endl;
 cout << "numfacets OUT = " << out.numberoftrifaces << endl;
 cout << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 // varre lista de elementos e passa para estrutura IEN
 IEN.Dim(numElems,5);
 cc.Dim(numVerts);
 cc.SetAll(0.0);
 for( int i=0;i<out.numberoftetrahedra;i++ )
 {
  // setting de cc = 0 para fora da bolha e cc = 0.5 para interface
  if( out.tetrahedronattributelist[i] == 1 )
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,0.0);
   }
  }
  // setting de cc = 1 para dentro da bolha e cc = 0.5 para interface
  else 
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = out.tetrahedronlist[i*4+j];
	IEN.Set(i,j,vertice);
	cc.Set(vertice,1.0);
   }
  }
 }

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 X.Dim(numNodes);
 Y.Dim(numNodes);
 Z.Dim(numNodes);
 for( int i=0;i<numVerts;i++ )
 {
  X.Set(i,out.pointlist[3*i+0]);
  Y.Set(i,out.pointlist[3*i+1]);
  Z.Set(i,out.pointlist[3*i+2]);
  if( out.pointmarkerlist[i] == 10 ||
	  out.pointmarkerlist[i] == 22 )
   cc.Set(i,0.5);
 }

 //breakup();
}


void Model3D::setDiskCouetteBC()
{
 real aux;
 rMax = Y.Max(); // CONFERIR! NAO SEI SE FUNCIONARA DIREITO!

 //for( int i=0;i<numNodes;i++ )
 //{
 // factor = 0.03;
 // aux = X.Get(i)*Z.Get(i)*factor;
 // wc.Set(i,aux);
 //}

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Max() &&
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)<rMax*rMax - 0.001) )
  {
   idbcp.AddItem(i);

   aux = 0.0;
   pc.Set(i,aux);
  }

  if( ( (Z.Get(i) == Z.Max()) || (Z.Get(i) == Z.Min()) ) &&
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = 0.0;
   wc.Set(i,aux);
   uc.Set(i,aux);
   vc.Set(i,aux);
  }

  if( Z.Get(i) == Z.Max() &&
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)<rMax*rMax - 0.001) )
  {
   idbcp.AddItem(i);

   aux = 0.0;
   pc.Set(i,aux);
  }


  if( Z.Get(i) == Z.Min() &&
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)<rMax*rMax - 0.001) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = 0.0;
   wc.Set(i,aux);
   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
  }

  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
 }
}

void Model3D::setNuCteDiskBC()
{
 real omega,aux;
 real radius;
 rMax = Y.Max();

 for( int i=0;i<numVerts;i++ )
 {
  radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

  if( Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcp.AddItem(i); // caso com c.c. livre em w
   //idbcw.AddItem(i);

   //uc.Set(i,radius*2.5666593e-02); // Z=4
   //vc.Set(i,radius*3.4943977e-02); // Z=4
   //wc.Set(i,-8.2505646e-01); // Z=4

   //uc.Set(i,radius*4.5487756e-03); // Z=6
   //vc.Set(i,radius*5.9598499e-03); // Z=6
   //wc.Set(i,-8.7414071e-01); // Z=6

   //uc.Set(i,radius*1.3326987e-04); // Z=10
   //vc.Set(i,radius*1.7327920e-04); // Z=10
   //wc.Set(i,-8.8416563E-01); // Z=10
   
   uc.Set(i,0.0); // Z=12
   vc.Set(i,0.0); // Z=12
   pc.Set(i,-0.391141); // caso com c.c. livre em w
  }

  if( Z.Get(i) == Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   omega=1.0;

   aux = (-1.0)*Y.Get(i)*omega;
   uc.Set(i,aux);
   aux = X.Get(i)*omega;
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);
  }

  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   //idbcp.AddItem(i);
   //aux = 0.0;
   //pc.Set(i,aux);
   //outflow.Set(i,aux);
  }
 }
}

void Model3D::setNuCDiskBC()
{
 real omega,aux,radius;
 rMax = Y.Max();

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Max() )
  {
   radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   ////idbcp.AddItem(i); // caso com c.c. livre em w

   //uc.Set(i,radius*2.6505291e-02); // Sc = 2000 Z = 4
   //vc.Set(i,radius*3.6178208e-02); // Sc = 2000 Z = 4
   //wc.Set(i,-8.2427145e-01); // Sc = 2000 Z = 4

   uc.Set(i,radius*1.3690760e-04); // Sc = 2000 Z = 10
   vc.Set(i,radius*1.7819422e-04); // Sc = 2000 Z = 10
   wc.Set(i,-8.8528405e-01); // Sc = 2000 Z = 10
   ////pc.Set(i,0.0); // caso com c.c. livre em w
   
   //wc.Set(i,-8.8559326E-01); // Sc = 2000
   //wc.Set(i,-9.1044679e-01); // Sc = 10
   //wc.Set(i,-9.2281563e-01); // Sc = 5
  }

  if( Z.Get(i) == Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   omega=1.0;

   aux = (-1.0)*Y.Get(i)*omega;
   uc.Set(i,aux);
   aux = X.Get(i)*omega;
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);

  }

  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   idbcp.AddItem(i);
   aux = 0.0;
   pc.Set(i,aux);
   outflow.Set(i,aux);
  }
 }
}

void Model3D::readAndSetVelocityDiskBC(const char* _dir,const char* _filename)
{
 int size = 2401;
 real aux;
 real dist1,dist2;
 clMatrix fghMatrix(size,4);

 string fileConcat = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileConcat.c_str();

 ifstream file( filename,ios::in );

 cout << endl;
 cout << "reading: " << filename  << " ...finished!" << endl;
 cout << endl;

 if( !file )
 {
  cerr << "Esta faltando o arquivo de perfil de velocidade!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !file.eof() )
 {
  for( int i=0;i<size;i++ )
  {
   file >> aux; // Z
   fghMatrix.Set(i,0,aux);
   file >> aux; // F
   fghMatrix.Set(i,1,aux);
   file >> aux; // G
   fghMatrix.Set(i,2,aux);
   file >> aux; // H
   fghMatrix.Set(i,3,aux);
  }
 }

 int j=0;
 real omega=1.0;
 for( int i=0;i<numNodes;i++ )
 {
  for( j=0;j<size-1;j++ )
  {
   dist1 = fabs( Z.Get(i) - fghMatrix.Get(j,0) );
   dist2 = fabs( Z.Get(i) - fghMatrix.Get(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = ( fghMatrix.Get(j,1)*X.Get(i)-fghMatrix.Get(j,2)*Y.Get(i) )*omega;
  uc.Set(i,aux);
  aux = ( fghMatrix.Get(j,2)*X.Get(i)-fghMatrix.Get(j,1)*Y.Get(i) )*omega;
  vc.Set(i,aux);
  aux = fghMatrix.Get(j,3);
  wc.Set(i,-1.0*aux);
 }
}

void Model3D::readAndSetPressureDiskBC(const char* _dir,const char* _filename)
{
 int size = 2401;
 real aux;
 real dist1,dist2;
 clMatrix pFile(size,2);

 string fileConcat = (string) _dir + (string) _filename + ".dat";
 const char* filename = fileConcat.c_str();

 ifstream file( filename,ios::in );

 cout << endl;
 cout << "reading: " << filename  << " ...finished!" << endl;
 cout << endl;

 if( !file )
 {
  cerr << "Esta faltando o arquivo de perfil da pressao!" << endl;
  exit(1);
 }

 // leitura do arquivo e transferencia para matriz
 if( !file.eof() )
 {
  for( int i=0;i<size;i++ )
  {
   file >> aux;
   pFile.Set(i,0,aux);
   file >> aux;
   pFile.Set(i,1,aux);
  }
 }

 int j;
 for( int i=0;i<numVerts;i++ )
 {
  for( j=0;j<size;j++ )
  {
   dist1 = fabs( Z.Get(i) - pFile(j,0) );
   dist2 = fabs( Z.Get(i) - pFile(j+1,0) );
   if( dist2 > dist1 ) break;
  }
//--------------------------------------------------
//   // applying b.c. only on the sidewall
//   if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
// 	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
//   {
//    aux = pFile(j,1);
//    pc.Set(i,aux);
//    outflow.Set(i,aux);
//   }
//-------------------------------------------------- 
   aux = pFile(j,1);
   pc.Set(i,aux);
   outflow.Set(i,aux);
 }
}

void Model3D::setCDiskBC()
{
 real aux;
 cc.Dim(numVerts);
 idbcc.Dim(0);

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Min() )
  {
   idbcc.AddItem(i);

   aux = 1.0;
   cc.Set(i,aux);
  }
 }
}

void Model3D::setBubbleBubbleBC()
{
 //real eps2 = 1E-2;
 real eps2 = 10E-3;
 //real eps2 = 60E-02;
 xCenter = 1.5;
 yCenter = 1.5;
 zCenter = 1.5;
 bubbleRadius = 0.5;
 cout << "xMax = " << X.Max() << endl;
 cout << "xCenter = " << xCenter << endl;
 cout << "yCenter = " << yCenter << endl;
 cout << "zCenter = " << zCenter << endl;
 cout << "bubbleRadius = " << bubbleRadius << endl;

 cc.Dim(numVerts);
 cc.SetAll(0.5);
 for( int i=0;i<numVerts;i++ )
 {
//--------------------------------------------------
//   cout << i << " " << ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+ 
// 	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter)+
// 	  (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter)) << " " << bubbleRadius*bubbleRadius << endl;
//-------------------------------------------------- 

  // dentro da bolha
  // [X-xCenter]^2 + [Y-yCenter]^2 < r^2
  if( ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+ 
	   (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter)+
	   (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter)+eps2<bubbleRadius*bubbleRadius))
  {
   cc.Set(i,1.0);
  }
  if( ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+ 
	   (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter)+
	   (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) -eps2>bubbleRadius*bubbleRadius))
  {
   cc.Set(i,0.0);
  }
  if( (X.Get(i)-xCenter)*(X.Get(i)-xCenter) +
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) +
	  (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) >
	  (X.Max()-xCenter)*(X.Max()-xCenter)-eps2 )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
//--------------------------------------------------
//   if( X.Get(i) == X.Max() )
//   {
//    idbcp.AddItem(i);
//    pc.Set(i,0.0);
//   }
//-------------------------------------------------- 
 }
} // fecham metodo setBubbleBubbleBC

void Model3D::set2BubbleBC()
{
 real aux;

 for( int i=0;i<numVerts;i++ )
 {
  // condicao de velocidade
  if( (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )  
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = X.Get(i);
   uc.Set(i,aux);
   aux = (-1.0)*Y.Get(i);
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);
  }
  // condicao de outflow
  if( (X.Get(i)==X.Min()) || X.Get(i)==X.Max() )  
  {
   idbcp.AddItem(i);

   aux = 0.0;
   pc.Set(i,aux);
  }
  if( (Z.Get(i)==Z.Min()) || Z.Get(i)==Z.Max() )  
  {
   idbcw.AddItem(i);

   aux = 0.0;
   wc.Set(i,aux);
  }
 }
}

void Model3D::setSphere(real _xC,real _yC,real _zC,real _r,real _eps)
{
 real eps2 = _eps;
 real r = _r;
 real xCenter = _xC;
 real yCenter = _yC;
 real zCenter = _zC;

 cc.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  /* bubble 1 */
  // dentro da bolha
  // [X-xCenter]^2 + [Y-yCenter]^2 < r^2
  if( (X.Get(i)-xCenter)*(X.Get(i)-xCenter) + 
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) +
	  (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) < r*r+eps2 )
  {
   cc.Set(i,1.0);
  }
  // na superficie da bolha
  // ( [X-xCenter]^2 + [Y-yCenter]^2 - r^2 ) < 10E-4
  if( fabs((X.Get(i)-xCenter)*(X.Get(i)-xCenter) +
	       (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) + 
	       (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) - r*r ) < eps2 ) 

  {
   cc.Set(i,0.5);
  }
 }
}

void Model3D::setCubeBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  // condicao de parede v=0
  if( (X.Get(i)==X.Max()) || (X.Get(i)==X.Min()) || 
      (Y.Get(i)==Y.Max()) || (Y.Get(i)==Y.Min()) || 
      (Z.Get(i)==Z.Max()) || (Z.Get(i)==Z.Min()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   real aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
 }
}

void Model3D::setWallBC()
{    
 for (list<int>::iterator it=outVert.begin(); it!=outVert.end(); ++it)
 {
  idbcu.AddItem(*it);
  idbcv.AddItem(*it);
  idbcw.AddItem(*it);

  real aux = 0.0;
  uc.Set(*it,aux);
  vc.Set(*it,aux);
  wc.Set(*it,aux);
 }
}

void Model3D::setWallAnnularBC()
{    
 for (list<int>::iterator it=outVert.begin(); it!=outVert.end(); ++it)
 {
  if(Y.Get(*it) > Y.Min() || Y.Get(*it) < Y.Max() )
  {
   idbcv.AddItem(*it);

   real aux = 0.0;
   vc.Set(*it,aux);
  }
  else
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   real aux = 0.0;
   uc.Set(*it,aux);
   vc.Set(*it,aux);
   wc.Set(*it,aux);
  }
 }
}

void Model3D::setCubeBC2()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (X.Get(i)==X.Max()) || 
	  (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
  if( (Z.Get(i)==Z.Max()) &&
      (X.Get(i)>X.Min()) && (X.Get(i)<X.Max()) && 
	  (Y.Get(i)>Y.Min()) && (Y.Get(i)<Y.Max()) )
  {
   idbcw.AddItem(i);
   idbcu.AddItem(i);
   idbcv.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,1.0);
  }
  if( (Z.Get(i)==Z.Min()) &&
      (X.Get(i)>X.Min()) && (X.Get(i)<X.Max()) && 
	  (Y.Get(i)>Y.Min()) && (Y.Get(i)<Y.Max()) )
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
  }
 }
}

void Model3D::setCube(real _lim1,real _lim2,real _eps)
{
 real eps2 = _eps;
 real lim1 = _lim1;
 real lim2 = _lim2;

 cc.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  // na interface
  if( (X.Get(i)<lim2+eps2) && (X.Get(i)>lim1-eps2) && 
      (Y.Get(i)<lim2+eps2) && (Y.Get(i)>lim1-eps2) && 
      (Z.Get(i)<lim2+eps2) && (Z.Get(i)>lim1-eps2) )
  {
   cc.Set(i,0.5);
  }
  // dentro da bolha
  if( (X.Get(i)<lim2-eps2) && (X.Get(i)>lim1+eps2) && 
      (Y.Get(i)<lim2-eps2) && (Y.Get(i)>lim1+eps2) && 
      (Z.Get(i)<lim2-eps2) && (Z.Get(i)>lim1+eps2) )
  {
   cc.Set(i,1.0);
  }
 }
}

void Model3D::setCube(real _xlimInf,real _xlimSup,
                      real _ylimInf,real _ylimSup,
					  real _zlimInf,real _zlimSup,real _eps)

{
 real eps2 = _eps;
 real xlimInf = _xlimInf;
 real ylimInf = _ylimInf;
 real zlimInf = _zlimInf;
 real xlimSup = _xlimSup;
 real ylimSup = _ylimSup;
 real zlimSup = _zlimSup;

 cc.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  // na interface
  if( (X.Get(i)<xlimSup+eps2) && (X.Get(i)>xlimInf-eps2) && 
      (Y.Get(i)<ylimSup+eps2) && (Y.Get(i)>ylimInf-eps2) && 
      (Z.Get(i)<zlimSup+eps2) && (Z.Get(i)>zlimInf-eps2) )
  {
   cc.Set(i,0.5);
  }
  // dentro da bolha
  if( (X.Get(i)<xlimSup-eps2) && (X.Get(i)>xlimInf+eps2) && 
      (Y.Get(i)<ylimSup-eps2) && (Y.Get(i)>ylimInf+eps2) && 
      (Z.Get(i)<zlimSup-eps2) && (Z.Get(i)>zlimInf+eps2) )
  {
   cc.Set(i,1.0);
  }
 }
}

void Model3D::setBubble3DBC()
{
 real eps2 = 10E-3;
 //real eps2 = 60E-02;
 xCenter = 1.5;
 yCenter = 1.5;
 zCenter = 1.5;
 bubbleRadius = 0.5;
 cout << "xCenter = " << xCenter << endl;
 cout << "yCenter = " << yCenter << endl;
 cout << "zCenter = " << yCenter << endl;
 cout << "bubbleRadius = " << bubbleRadius << endl;

 for( int i=0;i<numVerts;i++ )
 {
  //if( (X.Get(i)==X.Max()) && (Y.Get(i)==Y.Min()) )
  if( (X.Get(i)==X.Max()) && (Y.Get(i)==0) )
  {
//--------------------------------------------------
//    idbcp.AddItem(i);
//    pc.Set(i,0.0);
//    outflow.Set(i,0);
//-------------------------------------------------- 
  }
  // dentro da bolha
  // [X-xCenter]^2 + [Y-yCenter]^2 < r^2
  if( (X.Get(i)-xCenter)*(X.Get(i)-xCenter) + 
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) +
	  (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) <bubbleRadius*bubbleRadius+eps2 )
  {
   cc.Set(i,1.0);
  }
  // fora da bolha
  // [X-xCenter]^2 + [Y-yCenter]^2 < r^2
  if( (X.Get(i)-xCenter)*(X.Get(i)-xCenter) + 
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) +
	  (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) >bubbleRadius*bubbleRadius-eps2 )
  {
   cc.Set(i,0.0);
  }
  // na superficie da bolha
  // ( [X-xCenter]^2 + [Y-yCenter]^2 - r^2 ) < 10E-4
  if( fabs((X.Get(i)-xCenter)*(X.Get(i)-xCenter) +
	       (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) + 
	       (Z.Get(i)-zCenter)*(Z.Get(i)-zCenter) - 
	        bubbleRadius*bubbleRadius ) < eps2 ) 

  {
   cc.Set(i,0.5);
  }
  if( X.Get(i) == 0 || X.Get(i) == X.Max() ||
      Y.Get(i) == 0 || Y.Get(i) == Y.Max() ||
      Z.Get(i) == 0 || Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
 }
}

void Model3D::setBubbleBC2()
{
 //real eps2 = 10E-02;
 real eps2 = 1E-00;
 xCenter = (X.Max()+X.Min())/2.0;
 yCenter = (Y.Max()+Y.Min())/2.0;
 zCenter = (Z.Max()+Z.Min())/2.0;
 bubbleRadius = (X.Max()-xCenter)/6.0;
 int n = 20;
 clVector zz(n);
 real aux;
 for( int i=0;i<n;i++ )
 {
  aux = Z.Get(180*i+30);
  zz.Set(i,aux);
 }

 cc.SetAll(0.0);
 for( int i=0;i<numVerts;i++ )
 {
//--------------------------------------------------
//   //if( (X.Get(i)==X.Max()) && (Y.Get(i)==Y.Min()) )
//   if( (X.Get(i)==X.Max()) && (Y.Get(i)==0) )
//   {
//    idbcp.AddItem(i);
//    pc.Set(i,0.0);
//    outflow.Set(i,0);
//   }
//-------------------------------------------------- 
  if( ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+ 
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) - 3E-00<bubbleRadius*bubbleRadius)
	&& (Z.Get(i) > zz.Get(5)) && (Z.Get(i) < zz.Get(13))
	)
  {
   cc.Set(i,0.5);
  }
  // dentro da bolha
  // [X-xCenter]^2 + [Y-yCenter]^2 < r^2
  if( ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+ 
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) - eps2<bubbleRadius*bubbleRadius)
	&& (Z.Get(i) > zz.Get(6)) && (Z.Get(i) < zz.Get(12))
	)
  {
   cc.Set(i,1.0);
  }
//--------------------------------------------------
//   // fora da bolha
//   // [X-xCenter]^2 + [Y-yCenter]^2 > r^2
//   if( ((X.Get(i)-xCenter)*(X.Get(i)-xCenter)+
// 	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) + eps2>bubbleRadius*bubbleRadius)&&
// 	  (Z.Get(i) < zz.Get(7)) && (Z.Get(i) > zz.Get(12))
//     )
//   {
//    cc.Set(i,0.0);
//   }
//-------------------------------------------------- 
  // na superficie da bolha
  // ( [X-xCenter]^2 + [Y-yCenter]^2 - r^2 ) < 10E-4
  if( (Z.Get(i)==Z.Min()) || (Z.Get(i) == Z.Max()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() &&
	 (X.Get(i)-xCenter)*(X.Get(i)-xCenter) +
	 (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) + eps2>X.Max()*X.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
 }
}

void Model3D::readBaseStateNu( const char* _filename )
{
 real omega;
 real aux,dist1,dist2;
 int lines = 2401;
 int j;
 clMatrix fMatrix(lines,2);
 clMatrix gMatrix(lines,2);
 clMatrix hMatrix(lines,2);
 clMatrix dHdZMatrix(lines,2);
 clVector pVector(lines);


 string file = "../db/baseState/nuC/Sc2000/F_" + (string) _filename + ".dat";
 const char* filenameF = file.c_str();
 ifstream fFile( filenameF,ios::in );

 file = "../db/baseState/nuC/Sc2000/G_" + (string) _filename + ".dat"; 
 const char* filenameG = file.c_str();
 ifstream gFile( filenameG,ios::in );

 file = "../db/baseState/nuC/Sc2000/H_" + (string) _filename + ".dat"; 
 const char* filenameH =  file.c_str();
 ifstream hFile( filenameH,ios::in );

 file = "../db/baseState/nuC/Sc2000/dHdZ_" + (string) _filename + ".dat"; 
 const char* filenamedHdZ = file.c_str();
 ifstream dHdZFile( filenamedHdZ,ios::in );

 if( !fFile || !gFile || !hFile || !dHdZFile )
 {
  cerr << "Esta faltando algum arquivo do estado base para NuC!" << endl;
  exit(1);
 }

 for( int i=0;i<2401;i++ )
 {
  fFile >> aux;
  fMatrix.Set(i,0,aux);
  fFile >> aux;
  fMatrix.Set(i,1,aux);
  gFile >> aux;
  gMatrix.Set(i,0,aux);
  gFile >> aux;
  gMatrix.Set(i,1,aux);
  hFile >> aux;
  hMatrix.Set(i,0,aux);
  hFile >> aux;
  hMatrix.Set(i,1,aux);
  dHdZFile >> aux;
  dHdZMatrix.Set(i,0,aux);
  dHdZFile >> aux;
  dHdZMatrix.Set(i,1,aux);
  aux = dHdZMatrix.Get(i,1) - 0.5*hMatrix.Get(i,1)*hMatrix.Get(i,1);
  pVector.Set(i,aux);
 }

 j=0;
 omega=1.0;
 for( int i=0;i<numNodes;i++ )
 {
  for( j=0;j<lines-1;j++ )
  {
   dist1 = fabs( Z.Get(i) - fMatrix.Get(j,0) );
   dist2 = fabs( Z.Get(i) - fMatrix.Get(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = ( fMatrix.Get(j,1)*X.Get(i)-gMatrix.Get(j,1)*Y.Get(i) )*omega;
  uc.Set(i,aux);
  aux = ( gMatrix.Get(j,1)*X.Get(i)-fMatrix.Get(j,1)*Y.Get(i) )*omega;
  vc.Set(i,aux);
  aux = hMatrix.Get(j,1);
  wc.Set(i,aux);
 }

 j=0;
 for( int i=0;i<numVerts;i++ )
 {
  for( j=0;j<lines-1;j++ )
  {
   dist1 = fabs( Z.Get(i) - fMatrix.Get(j,0) );
   dist2 = fabs( Z.Get(i) - fMatrix.Get(j+1,0) );
   if( dist2 > dist1 ) break;
  }
  aux = pVector.Get(j);
  pc.Set(i,aux);
 }
}

void Model3D::setNuZDiskBC()
{
 real omega,aux,radius;
 rMax = Y.Max();

 for( int i=0;i<numVerts;i++ )
 {
  radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

  if( Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   ////idbcp.AddItem(i); // caso com c.c. livre em w

   //uc.Set(i,radius*4.8979239E-02); // Z=4
   //vc.Set(i,radius*7.9371092E-02); // Z=4
   //wc.Set(i,-9.1936908E-01); // Z=4

   //uc.Set(i,radius*6.8574756E-03); // Z=6
   //vc.Set(i,radius*1.0335366E-02); // Z=6
   //wc.Set(i,-0.10196142E+01); // Z=6

   uc.Set(i,radius*0.11735664E-03); // Z=10
   vc.Set(i,radius*0.17501519E-03); // Z=10
   wc.Set(i,-0.10193840E+01); // Z=10
   /////pc.Set(i,0.0); // Z=10 // caso com c.c. livre em w
  }

  if( Z.Get(i) == Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   omega=1.0;

   aux = (-1.0)*Y.Get(i)*omega;
   uc.Set(i,aux);
   aux = X.Get(i)*omega;
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);
  }

  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)> (rMax*rMax - 0.001) ) )
  {
   idbcp.AddItem(i);
   aux = 0.0;
   pc.Set(i,aux);
   outflow.Set(i,aux);
  }
 }
}

// FSBC - Free Surface Boundary Condition
void Model3D::setDiskFSBC()
{
 real aux;
 rMax = Y.Max();// CONFERIR! NAO SEI SE FUNCIONARA DIREITO!

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
  if( Z.Get(i)<=Z.Max() && Z.Get(i)>Z.Min() && 
	(X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) > (rMax*rMax - 0.001) ) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
  if( Z.Get(i)==Z.Max() && 
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMax*rMax - 0.001) ) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcp.AddItem(i);

   // cc = 0.5 -> noh da superficie
   //cc.Set(i,0.5); // para funcionamento do ALE

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   pc.Set(i,aux);
   outflow.Set(i,aux);
  }
 }
}

void Model3D::setAdimenDisk()
{
 real aux;
 real Red = 50;
 //real Red = 100;
 real factorz = 1.0/(Z.Max()-Z.Min());
 rMax = Y.Max();

 for( int i=0;i<numVerts;i++ )
 {
  aux = (X.Get(i)/rMax)*Red;
  X.Set(i,aux);
  aux = (Y.Get(i)/rMax)*Red;
  Y.Set(i,aux);
  //aux = Z.Get(i)*factorz*4;
  //aux = Z.Get(i)*factorz*6;
  aux = Z.Get(i)*factorz*10;
  Z.Set(i,aux);
 }
}

void Model3D::setAdimenDiskCouette()
{
 real aux;
 real Red = 100;
 real factorz = 1.0/(Z.Max()-Z.Min());

 for( int i=0;i<numVerts;i++ )
 {
  aux = (X.Get(i)/rMax)*Red;
  X.Set(i,aux);
  aux = (Y.Get(i)/rMax)*Red;
  Y.Set(i,aux);
  aux = Z.Get(i)*factorz*100;
  Z.Set(i,aux);
 }
}

void Model3D::setPerturbSurf()
{
 real aux;
 real Red = 1;
 real factorz = 1.0/(Z.Max()-Z.Min());

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) != Z.Min() &&
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMax*rMax - 0.1) ) )
  {
   aux = Z.Get(i)*( 1.0 + (X.Get(i)/Red)*factorz*0.4 );
   Z.Set(i,aux);
  }
 }
}

void Model3D::setPerturbSurf2()
{
 real aux;
 real Red = 1;
 real factorz = 1.0/(Z.Max()-Z.Min());
 real xMid = X.Min()+(X.Max()-X.Min())/2.0;

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) != Z.Min() && X.Get(i) > xMid  &&
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMax*rMax - 0.001) ) )
  {
   aux = Z.Get(i)*( 0.65 + (X.Get(i)/Red)*factorz*0.6 );
   Z.Set(i,aux);
  }
  if( Z.Get(i) != Z.Min() && X.Get(i) <= xMid  &&
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMax*rMax - 0.001) ) )
  {
   aux = Z.Get(i)*( 0.65 - (X.Get(i)/Red)*factorz*0.6 );
   Z.Set(i,aux);
  }
 }
}

void Model3D::setPerturbSurfSquare()
{
 real aux;
 real Red = 1;
 real factorz = 1.0/(Z.Max()-Z.Min());
 real rMid = rMax/10.0;

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) != Z.Min() && 
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMid*rMid - 0.001) ) )
  {
   aux = Z.Get(i)*( 1.0 + (Z.Get(i)/Red)*factorz*0.6 );
   Z.Set(i,aux);
  }
 }
}

void Model3D::setMiniElement()
{
 V.Dim(numElems);
 real centroidX,centroidY,centroidZ;
 int v1,v2,v3,v4,v5;

 numGLEP = 4; // triangulo linear
 numGLEU = 5; // elemento MINI
 numGLEC = 4; // elemento linear
 numNodes = numVerts + numElems;

 clearBC();
 reAllocStruct();

 for( int i=0;i<numElems;i++ )
 {
  v1 = (int) IEN.Get(i,0);
  v2 = (int) IEN.Get(i,1);
  v3 = (int) IEN.Get(i,2);
  v4 = (int) IEN.Get(i,3);

  real volume = getVolume(i);
  V.Set(i,volume);

  if( fabs(volume)<1.0E-10)
  {
   cout << i << endl;
   cout << volume << endl;
   cerr << "tetraedro singular, verificar a qualidade da malha!" << endl;
  }

  if( volume<0.0 )
  {
   IEN.Set(i,1,v3);
   IEN.Set(i,2,v2);
  };

  v5=numVerts+i;

  IEN.Set(i,4,v5);
  centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )*0.25;
  centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )*0.25;
  centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )*0.25;
  X.Set(v5,centroidX);
  Y.Set(v5,centroidY);
  Z.Set(v5,centroidZ);
 
 }
}

// criando os nos na metade de cada aresta. Para isso eh necessario
// fazer uma lista de arestas e numera-las de forma que uma aresta comum
// tenha apenas 1 numero e seja compartilhada em todos os elementos que
// tem a aresta. 
void Model3D::setQuadElement()
{
 int v1,v2,v3,v4;
 int numFace = 6; // teraedro tem 6 arestas
 V.Dim(numElems);
 clVector faceaux(2);
 IFACE2D *faces = NULL;
 int listSize = numFace*numElems;
 faces = new IFACE2D[listSize];
 for( int mele=0;mele<numElems;mele++ )
 {
  v1 = (int) IEN.Get(mele,0);
  v2 = (int) IEN.Get(mele,1);
  v3 = (int) IEN.Get(mele,2);
  v4 = (int) IEN.Get(mele,3);

  real volume = getVolume(mele);
  V.Set(mele,volume);

  if( fabs(volume)<1.0E-10)
   cerr << "tetraedro singular, verificar a qualidade da malha!" << endl;

  if( volume<0.0 )
  {
   IEN.Set(mele,1,v3);
   IEN.Set(mele,2,v2);
  }
  // -------------------------------------------- //

  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+0].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+0].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+0].p3 = mele;

  faceaux.Set(0,v1);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+1].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+1].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+1].p3 = mele;

  faceaux.Set(0,v1);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+2].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+2].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+2].p3 = mele;

  faceaux.Set(0,v2);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+3].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+3].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+3].p3 = mele;

  faceaux.Set(0,v2);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+4].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+4].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+4].p3 = mele;

  faceaux.Set(0,v3);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+5].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+5].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+5].p3 = mele;
 }

 // ordena uma estrutura (matriz) em ordem crescente na linha e coluna
 // as faces continuam repetidas neste ponto.
 qsort(faces,listSize,sizeof(IFACE2D),IFACE2DCompare);

 /*        - nome: mapEdge
           - definicao: matrix com mapeamento de arestas
		                Identificacao da aresta, coordenadas X e Y
						Ultima coluna representa o elemento que contem a aresta

     ---   +---+---+---+---+---+---+---+
      |    | a | b | c | d | e | f | g |  a = identificacao da aresta
      |    +---+---+---+---+---+---+---+
	  |    .   .   .   .   .   .   .   .  b = coordenada X do noh
	  |	   .   .   .   .   .   .   .   .
	  g    .   .   .   .   .   .   .   .  c = coordenada Y do noh
      |    +---+---+---+---+---+---+---+                               
      |    | a | b | c | c | e | f | g |  d = coordenada Z do noh
      |    +---+---+---+---+---+---+---+ 
      |    | a | b | c | c | e | f | g |  e = vertice vizinho v1
     ---   +---+---+---+---+---+---+---+                                           
                                          f = vertice vizinho v2
	       |____________ h ____________|                                           
		   |                           |  g = elemento que contem o noh da aresta
                                                                                   
		                                  h = numero de arestas (inclui repetidas) 

										  i = 7 colunas
 
 */

 int it=0;
 clMatrix mapEdge(listSize,7);

 // numeracao de arestas a partir de numVerts e associacao das arestas
 // aos elementos para inclusao na IEN
 for( int i=1;i<=listSize;i++ )
 {
  real x=X.Get(faces[i-1].p1)+(X.Get(faces[i-1].p2)-X.Get(faces[i-1].p1))*0.5;
  real y=Y.Get(faces[i-1].p1)+(Y.Get(faces[i-1].p2)-Y.Get(faces[i-1].p1))*0.5;
  real z=Z.Get(faces[i-1].p1)+(Z.Get(faces[i-1].p2)-Z.Get(faces[i-1].p1))*0.5;
  mapEdge.Set(i-1,0,numVerts+it); // numero da aresta
  mapEdge.Set(i-1,1,x ); // coordenada X da aresta
  mapEdge.Set(i-1,2,y ); // coordenada Y da aresta
  mapEdge.Set(i-1,3,z ); // coordenada Y da aresta
  mapEdge.Set(i-1,4,faces[i-1].p1 ); // 1o noh
  mapEdge.Set(i-1,5,faces[i-1].p2 ); // 2o noh
  mapEdge.Set(i-1,6,faces[i-1].p3 ); // elemento que contem a aresta
  while( (faces[i].p1 == faces[i-1].p1) &&
         (faces[i].p2 == faces[i-1].p2) )
  {
  mapEdge.Set(i,0,mapEdge.Get(i-1,0));
  mapEdge.Set(i,1,mapEdge.Get(i-1,1));
  mapEdge.Set(i,2,mapEdge.Get(i-1,2));
  mapEdge.Set(i,3,mapEdge.Get(i-1,3));
  mapEdge.Set(i,4,faces[i-1].p1 );
  mapEdge.Set(i,5,faces[i-1].p2 );
  mapEdge.Set(i,6,faces[i].p3 );
   i++;
  }
  it++; // numero total de arestas
 }

 // atualizado vetores com numero total de nos
 numGLEP = 4; // tetraedro linear
 numGLEC = 4; // tetraedro linear
 numGLEU = 10; // tetraedro quadratico
 numNodes = numVerts+it; // atualizando numNodes

 clearBC();
 reAllocStruct();

 // adicionando os vertices das arestas nas estruturas X,Y,Z e IEN
 vector< list<int> > edge;
 edge.resize(0);
 edge.resize(numElems);
 list<int> plist;
 list<int>::iterator ele;
 for( int i=0;i<listSize;i++ )
 {
  int node = mapEdge.Get(i,0);
  int elem = mapEdge.Get(i,6);
  X.Set(node,mapEdge.Get(i,1));
  Y.Set(node,mapEdge.Get(i,2));
  Z.Set(node,mapEdge.Get(i,3));
  // vertice v5 que fica entre v1 e v2
  if( (mapEdge.Get(i,4) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,1)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,1)) )
  {
   IEN.Set(elem,4,node);
  }
  // vertice v6 que fica entre v1 e v3
  if( (mapEdge.Get(i,4) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,2)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,2)) )
  {
   IEN.Set(elem,5,node);
  }
  // vertice v7 que fica entre v1 e v4
  if( (mapEdge.Get(i,4) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,3)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,0) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,3)) )
  {
   IEN.Set(elem,6,node);
  }
  // vertice v8 que fica entre v2 e v3
  if( (mapEdge.Get(i,4) == IEN.Get(elem,1) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,2)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,1) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,2)) )
  {
   IEN.Set(elem,7,node);
  }
  // vertice v9 que fica entre v3 e v4
  if( (mapEdge.Get(i,4) == IEN.Get(elem,2) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,3)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,2) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,3)) )
  {
   IEN.Set(elem,8,node);
  }
  // vertice v10 que fica entre v2 e v4
  if( (mapEdge.Get(i,4) == IEN.Get(elem,1) && 
	   mapEdge.Get(i,5) == IEN.Get(elem,3)) || 
	  (mapEdge.Get(i,5) == IEN.Get(elem,1) && 
	   mapEdge.Get(i,4) == IEN.Get(elem,3)) )
  {
   IEN.Set(elem,9,node);
  }
 }

 delete[] faces;
}

void Model3D::setNeighbour()
{
 neighbourElem.resize (0);
 neighbourElem.resize (numVerts);
 for( int i=0;i<numElems;i++ )
  for( int j= 0;j<numGLEP;j++ )
   neighbourElem.at( (int)IEN.Get(i,j) ).push_back(i);
}

void Model3D::setVertNeighbour()
{
 // cria lista de vizinhos para toda a malha
 int v1,v2,v3,v4;
 list<int> plist;
 list<int>::iterator mele;
 neighbourVert.resize (0);
 neighbourVert.resize (numVerts);

 for( int ii=0;ii<numVerts;ii++ )
 {
  clVector vert(0);
  plist = neighbourElem.at(ii);
  for( mele=plist.begin(); mele != plist.end();++mele )
  {
   v1 = (int) IEN.Get(*mele,0);
   v2 = (int) IEN.Get(*mele,1);
   v3 = (int) IEN.Get(*mele,2);
   v4 = (int) IEN.Get(*mele,3);
  // cout << v1 << " " << v2 << " " << v3 << " " << v4 << endl;
   if( v1==ii )
   {
	neighbourVert.at( ii ).push_back(v2);
	neighbourVert.at( ii ).push_back(v3);
	neighbourVert.at( ii ).push_back(v4);
   }
   if( v2==ii )
   {
	neighbourVert.at( ii ).push_back(v1);
	neighbourVert.at( ii ).push_back(v3);
	neighbourVert.at( ii ).push_back(v4);
   }
   if( v3==ii )
   {
	neighbourVert.at( ii ).push_back(v1);
	neighbourVert.at( ii ).push_back(v2);
	neighbourVert.at( ii ).push_back(v4);
   }
   if( v4==ii )
   {
	neighbourVert.at( ii ).push_back(v1);
	neighbourVert.at( ii ).push_back(v2);
	neighbourVert.at( ii ).push_back(v3);
   }
  }
  //cout << "---------" << ii << "------------" << endl;
  neighbourVert.at( ii ).sort();
  neighbourVert.at( ii ).unique();
//--------------------------------------------------
//   std::ostream_iterator< int > output( cout, " " );
//   std::copy( neighbourVert.at(ii).begin(), 
//              neighbourVert.at(ii).end(), output );
//   cout << endl;
//-------------------------------------------------- 
 }
}

// configura os vetores: 
// surfaceViz  -> vizinhos de cada no na superficie
void Model3D::setSurface()
{
 int surfaceNode;
 list<int> plist;
 list<int>::iterator vert;

 /*        - nome: surfaceViz
           - definicao: vetor de listas de vizinhos de vertices na interface 
		   - obs: no caso 2D n=2. no caso 3D n pode variar de vertice 
		          para vertice.

     ---   +---+----+-----+----+
      |    | b | a1 | ... | an |   a1...an = identificacao dos vizinhos (lista)
      |    +---+----+-----+----+
    c |    | b | a1 | ... | an |   b = identificacao do vertice de trabalho   
      |    +---+----+-----+----+
      |    | b | a1 | ... | an |   c = numero de vertices na interface
     ---   +---+----+-----+----+

 */

 // procurando vertices da bolha
 clVector surfaceAux = cc==0.5;
 surface = surfaceAux.Find();
 clVector nonSurfaceAux = cc!=0.5;
 nonSurface = nonSurfaceAux.Find();

 // dimensionando vetores
 surfaceViz.resize ( 0);
 surfaceViz.resize ( surface.Dim() );

 for( int i=0;i<surface.Dim();i++ )
 {
  surfaceNode = surface.Get(i);
  plist = neighbourVert.at(surfaceNode); 

  // adicionando no primeiro elemento da lista o valor do vertice de
  // trabalho
  surfaceViz.at( i ).push_back(surfaceNode); 
  // procura dos vertices adjacentes na interface (um de cada lado)
  for( vert=plist.begin();vert!=plist.end();++vert )
  {
   if( cc.Get(*vert) == 0.5 )   
   {
	surfaceViz.at( i ).push_back(*vert);
   }
  }
 }
} // fecha metodo setSurface

// as faces da superficie sao numeradas de 0..n-1 onde n eh o numero de
// faces na superficie. Para o caso 3D as faces sao triangulos. Junto a
// esta numeracao ha um mapeamento dos vertices da superficie vizinhos de 
// cada vertice na superficie. Os pontos vizinhos estao ordenados no
// sentido do relogio, com isso eh possivel descobrir facilmente o vetor
// normal a cada ponto da superficie utilizando as arestas dos elementos
// e produto escalar dos vetores localizados na aresta. 
void Model3D::setSurfaceFace()
{
 int v1,v2,v3,v4;
 list<int> plist;
 list<int>::iterator mele;
 int surfaceNode;
 
 // mapeamento de faces da interface com numeracao
 elemSurface.resize ( 0 ); 
 elemSurface.resize ( numVerts ); 
 // mapeamento de vertices das faces numeradas por elemSurface
 // super dimensionado!!!
 neighbourFaceVert.resize ( 0 );
 neighbourFaceVert.resize ( 5*numNodes );
 // contador de faces na superficie
 int count = 0; 

 for( int i=0;i<surface.Dim();i++ )
 {
  surfaceNode = surface.Get(i);
  // lista de elementos que contem surfaceNode
  plist = neighbourElem.at(surfaceNode);
  for( mele=plist.begin(); mele != plist.end();++mele )
  {
   v1 = (int) IEN.Get(*mele,0);
   v2 = (int) IEN.Get(*mele,1);
   v3 = (int) IEN.Get(*mele,2);
   v4 = (int) IEN.Get(*mele,3);

   // pegando o elemento com 3 vertices na superfice e 1 fora da bolha
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	neighbourFaceVert.at( count ).push_back(v1);
	neighbourFaceVert.at( count ).push_back(v2);
	neighbourFaceVert.at( count ).push_back(v3);
	neighbourFaceVert.at( count ).push_back(v4);
	count++;
   }
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	neighbourFaceVert.at( count ).push_back(v4);
	neighbourFaceVert.at( count ).push_back(v1);
	neighbourFaceVert.at( count ).push_back(v2);
	neighbourFaceVert.at( count ).push_back(v3);
	count++;
   }
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	neighbourFaceVert.at( count ).push_back(v3);
	neighbourFaceVert.at( count ).push_back(v4);
	neighbourFaceVert.at( count ).push_back(v1);
	neighbourFaceVert.at( count ).push_back(v2);
	count++;
   }
   if( cc.Get(v1)==0 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	neighbourFaceVert.at( count ).push_back(v2);
	neighbourFaceVert.at( count ).push_back(v3);
	neighbourFaceVert.at( count ).push_back(v4);
	neighbourFaceVert.at( count ).push_back(v1);
	count++;
   }

   // reordenando faces da superfice na ordem dos ponteiros do relogio 
   // (clockwise)
   list<int> plist = elemSurface.at (surfaceNode);
   for( list<int>::iterator face=plist.begin();face!=plist.end();++face )
   {
	list<int> plist2 = neighbourFaceVert.at (*face);
	list<int>::iterator vert=plist2.begin();
	int v1Old = *vert;++vert;
	int v2Old = *vert;++vert;
	int v3Old = *vert;++vert;
	int vIn = *vert; // vertice dentro da bolha
	vert=plist2.end();

	if( v1Old == surfaceNode )
	{
	 if( checkNormal(surfaceNode,v2Old,v3Old,vIn) == true )
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v2Old);
	  neighbourFaceVert.at( *face ).push_back(v3Old);
	 }
	 else
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v3Old);
	  neighbourFaceVert.at( *face ).push_back(v2Old);
	 }
	}
	if( v2Old == surfaceNode )
	{
	 if( checkNormal(surfaceNode,v1Old,v3Old,vIn) == true )
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v1Old);
	  neighbourFaceVert.at( *face ).push_back(v3Old);
	 }
	 else
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v3Old);
	  neighbourFaceVert.at( *face ).push_back(v1Old);
	 }
	}
	if( v3Old == surfaceNode )
	{
	 if( checkNormal(surfaceNode,v1Old,v2Old,vIn) == true )
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v1Old);
	  neighbourFaceVert.at( *face ).push_back(v2Old);
	 }
	 else
	 {
	  neighbourFaceVert.at( *face ).clear();
	  neighbourFaceVert.at( *face ).push_back(v2Old);
	  neighbourFaceVert.at( *face ).push_back(v1Old);
	 }
	}
   }
  }
//--------------------------------------------------
//   cout << "---------" << surfaceNode << "------------" << endl;
//   std::ostream_iterator< int > output( cout, " " );
//   std::copy( elemSurface.at(surfaceNode).begin(), elemSurface.at(surfaceNode).end(), output );
//   ///std::copy( neighbourFaceVert.at(surfaceNode).begin(), neighbourFaceVert.at(surfaceNode).end(), output );
//   cout << endl;
//-------------------------------------------------- 
  //
//--------------------------------------------------
//   int c1=0;
//   if( surfaceNode == 29 )
//   //if( surfaceNode == 144 )
//   {
//    list<int> plist = elemSurface.at (surfaceNode);
//    for( list<int>::iterator face=plist.begin();face!=plist.end();++face )
//    {
// 	cout << "Triangulo: ------------------------- " << c1 << endl;
// 	list<int> plist2 = neighbourFaceVert.at (*face);
// 	list<int>::iterator vert=plist2.begin();
// 	int v0 = *vert;++vert;
// 	int v1 = *vert;
// 	vert=plist2.end();
// 	cout << "v0 = " << v0 << " " << "v1 = " << v1 << " " << endl;
// 
// 	real testx = ( ( Y.Get(v1)-Y.Get(v0) )*( Z.Get(v2)-Z.Get(v0) ) )-( (Z.Get(v1)-Z.Get(v0))*(Y.Get(v2)-Y.Get(v0) ) );
// 	real testy = ( ( Z.Get(v1)-Z.Get(v0) )*( X.Get(v2)-X.Get(v0) ) )-( (X.Get(v1)-X.Get(v0))*(Z.Get(v2)-Z.Get(v0) ) );
// 	real testz = ( ( X.Get(v1)-X.Get(v0) )*( Y.Get(v2)-Y.Get(v0) ) )-( (Y.Get(v1)-Y.Get(v0))*(X.Get(v2)-X.Get(v0) ) );
// 	cout << testx << endl;
// 	cout << testy << endl;
// 	cout << testz << endl;
// 
// 	if( testz < 0.0 )
// 	{
// 	 int aux = v0;
// 	 v0=v1;
// 	 v1=aux;
// 	}
// 	//cout << "v0 = " << v0 << " " << "v1 = " << v1 << " " << "v2 = " << v2 << endl;
// 	c1++;
//    }
//   }
//-------------------------------------------------- 
 }
 neighbourFaceVert.resize (count); // trim vector para numero real de itens
}

/* cria matriz IEN para os elementos da superficie que no caso 3D sao
 * triagulos. Este metodo utiliza neighbourFaceVert e elemSurface 
 * como input e a funcao qsort para ordenacao de arestas.
 * Note que a triangulacao esta referenciada pelas coordenadas globais
 * X,Y e Z. Isto quer dizer que a malha da superfice tem 'buracos' em
 * sua numeracao. Para reordenar a malha com os vertices comecando do
 * indice '0' eh necessario passar no metodo arrangeIEN.
 * */
void Model3D::setSurfaceTri()
{
 clVector edgeaux(3);
 IFACE3D *edge = NULL;
 int listSize = numElems;
 edge = new IFACE3D[listSize];
 list<int> plist,plist2;
 list<int>::iterator face,vert;

 int count = 0;
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  plist = elemSurface.at (surfaceNode);

  for(face=plist.begin();face!=plist.end();++face)
  {
   plist2 = neighbourFaceVert.at(*face);
   vert=plist2.begin();
   int v1 = *vert;++vert;
   int v2 = *vert;
   vert=plist2.end();

   edgeaux.Set(0,v1);
   edgeaux.Set(1,v2);
   edgeaux.Set(2,surfaceNode);
   edgeaux.Sort(); // para ordenar os vertices do triangulo

   edge[count].p1 = (int) edgeaux.Get(0);
   edge[count].p2 = (int) edgeaux.Get(1);
   edge[count].p3 = (int) edgeaux.Get(2);
   count++;
  }
 }
 // ordena uma estrutura (matriz) em ordem crescente na linha e coluna
 // as faces continuam repetidas neste ponto.
 qsort(edge,count,sizeof(IFACE3D),IFACE3DCompare);

 // count / 3 pois a numeracao do triangulo pode ser repetida por 3
 // vertices que sao os vertices do proprio triangulo
 interfaceMesh.numElems = count/3;
 interfaceMesh.IEN.Dim(interfaceMesh.numElems,3);
 interfaceMesh.numVerts = surface.Max()+1;
 interfaceMesh.X.Dim(interfaceMesh.numVerts);
 interfaceMesh.Y.Dim(interfaceMesh.numVerts);
 interfaceMesh.Z.Dim(interfaceMesh.numVerts);
 interfaceMesh.Marker.Dim(interfaceMesh.numVerts);

 int it=0;
 for( int i=0;i<count/3;i++ )
 {
  int v1 = edge[it].p1;
  int v2 = edge[it].p2;
  int v3 = edge[it].p3;

  interfaceMesh.IEN.Set(i,0,v1);
  interfaceMesh.IEN.Set(i,1,v2);
  interfaceMesh.IEN.Set(i,2,v3);

  interfaceMesh.X.Set(v1,surfMesh.X.Get(v1));
  interfaceMesh.X.Set(v2,surfMesh.X.Get(v2));
  interfaceMesh.X.Set(v3,surfMesh.X.Get(v3));

  interfaceMesh.Y.Set(v1,surfMesh.Y.Get(v1));
  interfaceMesh.Y.Set(v2,surfMesh.Y.Get(v2));
  interfaceMesh.Y.Set(v3,surfMesh.Y.Get(v3));

  interfaceMesh.Z.Set(v1,surfMesh.Z.Get(v1));
  interfaceMesh.Z.Set(v2,surfMesh.Z.Get(v2));
  interfaceMesh.Z.Set(v3,surfMesh.Z.Get(v3));

  interfaceMesh.Marker.Set(v1,0.5);
  interfaceMesh.Marker.Set(v2,0.5);
  interfaceMesh.Marker.Set(v3,0.5);

  it=it+3;
 }
 delete[] edge;
}

/*
 * cria matriz IEN para os elementos do convex-hull que no caso 3D sao
 * triagulos. Este metodo utiliza freeFace. 
 * Note que a triangulacao esta referenciada pelas coordenadas globais
 * X,Y e Z. Isto quer dizer que a malha do convex-hull tem 'buracos' em
 * sua numeracao. Para reordenar a malha com os vertices comecando do
 * indice '0' eh necessario passar no metodo arrangeIEN.
 */
void Model3D::setConvexTri()
{
 convexMesh.numElems = freeFace.DimI();
 convexMesh.IEN.Dim(convexMesh.numElems,3);

 for (int i=0;i<freeFace.DimI();i++ )
 {
  int v1 = freeFace.Get(i,2);
  int v2 = freeFace.Get(i,3);
  int v3 = freeFace.Get(i,4);

  convexMesh.IEN.Set(i,0,v1);
  convexMesh.IEN.Set(i,1,v2);
  convexMesh.IEN.Set(i,2,v3);
 }

 convexMesh.numVerts = convexMesh.IEN.Max()+1;
 convexMesh.X.Dim(convexMesh.numVerts);
 convexMesh.Y.Dim(convexMesh.numVerts);
 convexMesh.Z.Dim(convexMesh.numVerts);
 convexMesh.Marker.Dim(convexMesh.numVerts);

 for (int i=0;i<convexMesh.numElems;i++ )
 {
  int v1 = convexMesh.IEN.Get(i,0);
  int v2 = convexMesh.IEN.Get(i,1);
  int v3 = convexMesh.IEN.Get(i,2);

  convexMesh.X.Set(v1,surfMesh.X.Get(v1));
  convexMesh.X.Set(v2,surfMesh.X.Get(v2));
  convexMesh.X.Set(v3,surfMesh.X.Get(v3));

  convexMesh.Y.Set(v1,surfMesh.Y.Get(v1));
  convexMesh.Y.Set(v2,surfMesh.Y.Get(v2));
  convexMesh.Y.Set(v3,surfMesh.Y.Get(v3));

  convexMesh.Z.Set(v1,surfMesh.Z.Get(v1));
  convexMesh.Z.Set(v2,surfMesh.Z.Get(v2));
  convexMesh.Z.Set(v3,surfMesh.Z.Get(v3));

  convexMesh.Marker.Set(v1,0.0);
  convexMesh.Marker.Set(v2,0.0);
  convexMesh.Marker.Set(v3,0.0);
 }
}

/* This method build the surfMesh struct, combining interface and convex
 * structs. This struct can only be created after the assembling of
 * interfaceMesh and convexMesh.
 * The method is useful to define the 3D triangle surface mesh 
 * that is sent straigth to the mesh generator program (TETGEN). This is
 * essential to define mesh limits and interfaces between fluids and
 * solids. 
 * */
void Model3D::buildSurfMesh()
{
 surfMesh.numElems = interfaceMesh.numElems+convexMesh.numElems;

 int maxInterface = interfaceMesh.IEN.Max();
 int maxConvex = convexMesh.IEN.Max();
 int max = maxInterface;
 if( max < maxConvex )
  max = maxConvex;
 surfMesh.numVerts = max+1;

 surfMesh.IEN = interfaceMesh.IEN;
 surfMesh.IEN.Append(convexMesh.IEN);
 surfMesh.X.Dim(surfMesh.numVerts);
 surfMesh.Y.Dim(surfMesh.numVerts);
 surfMesh.Z.Dim(surfMesh.numVerts);
 surfMesh.Marker.Dim(surfMesh.numVerts);

 for( int elem=0;elem<interfaceMesh.numElems;elem++ )
 {
  int v1 = interfaceMesh.IEN.Get(elem,0);
  int v2 = interfaceMesh.IEN.Get(elem,1);
  int v3 = interfaceMesh.IEN.Get(elem,2);

  real x1 = interfaceMesh.X.Get(v1); 
  real x2 = interfaceMesh.X.Get(v2); 
  real x3 = interfaceMesh.X.Get(v3); 

  real y1 = interfaceMesh.Y.Get(v1); 
  real y2 = interfaceMesh.Y.Get(v2); 
  real y3 = interfaceMesh.Y.Get(v3); 

  real z1 = interfaceMesh.Z.Get(v1);
  real z2 = interfaceMesh.Z.Get(v2);
  real z3 = interfaceMesh.Z.Get(v3);

  surfMesh.X.Set( v1,x1 );
  surfMesh.X.Set( v2,x2 );
  surfMesh.X.Set( v3,x3 );

  surfMesh.Y.Set( v1,y1 );
  surfMesh.Y.Set( v2,y2 );
  surfMesh.Y.Set( v3,y3 );

  surfMesh.Z.Set( v1,z1 );
  surfMesh.Z.Set( v2,z2 );
  surfMesh.Z.Set( v3,z3 );

  surfMesh.Marker.Set( v1,0.5 );
  surfMesh.Marker.Set( v2,0.5 );
  surfMesh.Marker.Set( v3,0.5 );
 }

 for( int elem=0;elem<convexMesh.numElems;elem++ )
 {
  int v1 = convexMesh.IEN.Get(elem,0);
  int v2 = convexMesh.IEN.Get(elem,1);
  int v3 = convexMesh.IEN.Get(elem,2);

  real x1 = convexMesh.X.Get(v1); 
  real x2 = convexMesh.X.Get(v2); 
  real x3 = convexMesh.X.Get(v3); 

  real y1 = convexMesh.Y.Get(v1); 
  real y2 = convexMesh.Y.Get(v2); 
  real y3 = convexMesh.Y.Get(v3); 

  real z1 = convexMesh.Z.Get(v1);
  real z2 = convexMesh.Z.Get(v2);
  real z3 = convexMesh.Z.Get(v3);

  surfMesh.X.Set( v1,x1 );
  surfMesh.X.Set( v2,x2 );
  surfMesh.X.Set( v3,x3 );
                       
  surfMesh.Y.Set( v1,y1 );
  surfMesh.Y.Set( v2,y2 );
  surfMesh.Y.Set( v3,y3 );
                       
  surfMesh.Z.Set( v1,z1 );
  surfMesh.Z.Set( v2,z2 );
  surfMesh.Z.Set( v3,z3 );

  surfMesh.Marker.Set( v1,0.0 );
  surfMesh.Marker.Set( v2,0.0 );
  surfMesh.Marker.Set( v3,0.0 );
 }
}

/*
 * este metodo organiza uma estrutura do tipo IEN ordenando a partir de
 * um numero todos os outros vertices sem deixar 'buracos' na estrutura.
 * Isso quer dizer que podemos recriar uma malha comecando a numeracao
 * dos vertices a partir de um numero qualquer.
 * Para isso eh necessario criar um vetor de listas do tipo setNeighbour()
 * e entao fazer um mapeamento de cada vertice para o numero
 * correspondente.
 * !!!metodo ainda nao testado 100%!!!
 * */
SurfaceMesh Model3D::arrangeMesh(SurfaceMesh _mesh,int _nVerts,int _begin)
{
 int nElems = _mesh.numElems;
 int end = _mesh.IEN.Max();

 SurfaceMesh meshReturn;
 meshReturn.IEN.Dim(nElems,3);
 meshReturn.X.Dim(_nVerts);
 meshReturn.Y.Dim(_nVerts);
 meshReturn.Z.Dim(_nVerts);
 meshReturn.Marker.Dim(_nVerts);

 // lista de elementos que contem o vertice indicado pela linha do array
 vector< list<int> > test;  
 test.resize (0);
 test.resize (end+1);
 for( int i=0;i<nElems;i++ )
  for( int j= 0;j<3;j++ )
   test.at( (int) _mesh.IEN.Get(i,j) ).push_back(i);

 int count = 0;
 for( int i=0;i<end+1;i++ )
 {
  // este teste serve para pular todos os 'buracos' da triangulacao
  if( test.at(i).size() > 0 ) 
  {
   //--------------------------------------------------
   //    cout << count << " " << i << " ";
   //    std::ostream_iterator< int > output( cout, " " );
   //    std::copy( test.at(i).begin(),test.at(i).end(), output );
   //    cout << endl;
   //-------------------------------------------------- 

   list<int>::iterator it;
   for( it=test.at(i).begin();it!=test.at(i).end();++it )
   {
	int v1 = _mesh.IEN.Get(*it,0);
	int v2 = _mesh.IEN.Get(*it,1);
	int v3 = _mesh.IEN.Get(*it,2);
	if (v1 == i )
	{
	 meshReturn.IEN.Set(*it,0,count+_begin);
	 meshReturn.X.Set(count,X.Get(v1));
	 meshReturn.Y.Set(count,Y.Get(v1));
	 meshReturn.Z.Set(count,Z.Get(v1));
	 meshReturn.Marker.Set(count,cc.Get(v1));
	}
	else if (v2 == i )
	{
	 meshReturn.IEN.Set(*it,1,count+_begin);
	 meshReturn.X.Set(count,X.Get(v2));
	 meshReturn.Y.Set(count,Y.Get(v2));
	 meshReturn.Z.Set(count,Z.Get(v2));
	 meshReturn.Marker.Set(count,cc.Get(v2));
	}
	else 
	{
	 meshReturn.IEN.Set(*it,2,count+_begin);
	 meshReturn.X.Set(count,X.Get(v3));
	 meshReturn.Y.Set(count,Y.Get(v3));
	 meshReturn.Z.Set(count,Z.Get(v3));
	 meshReturn.Marker.Set(count,cc.Get(v3));
	}
   }
   count++;
  }
 }
 return meshReturn;
}

// este metodo cria duas listas com os vertices do convex hull (outVert)
// e todos os outros menos os vertices do convex hull (inVert). Este
// metodo eh especialmente interessante para aplicacao das condicoes de
// contorno para geometrias complexas.
void Model3D::setInOutVert()
{
 inVert.resize (0);
 outVert.resize (0);

 for(int i=0;i<freeFace.DimI();i++ )
 {
  int v1 = freeFace.Get(i,2);
  int v2 = freeFace.Get(i,3);
  int v3 = freeFace.Get(i,4);
  outVert.push_back(v1);
  outVert.push_back(v2);
  outVert.push_back(v3);
 }
 outVert.sort();
 outVert.unique();

 // retira de inVert todos os vertices presents em outVert.
 list<int>::iterator it;
 for( int vert=0;vert<numVerts;vert++ )
  inVert.push_back(vert);
 for( it=outVert.begin();it!=outVert.end();++it )
  inVert.remove(*it);

//--------------------------------------------------
//  cout << "outVert contains:";
//  for (it=outVert.begin(); it!=outVert.end(); ++it)
//   cout << " " << *it;
//  cout << endl;
//-------------------------------------------------- 
}

/*
 * cc inside = 1.0  |  cc interface = 0.5 | cc outside = 0.0
 * if sum of cc (element) > 2.0, the element is inside of bubble
 * if sum of cc (element) < 2.0, the element is outside
 * if sum of cc (element) = 2.0, by convention, the element is inside
*/ 
void Model3D::setInOutElem()
{
 inElem.resize (0);
 outElem.resize (0);
 for(int i=0;i<IEN.DimI();i++ )
 {
  int v1 = IEN.Get(i,0);
  int v2 = IEN.Get(i,1);
  int v3 = IEN.Get(i,2);
  int v4 = IEN.Get(i,3);
  if( cc.Get(v1) + cc.Get(v2) + cc.Get(v3) + cc.Get(v4) < 2.0 )
   outElem.push_back(i);
  else
   inElem.push_back(i);
 }
}
	   
// cria matrizes de mapeamentos de vizinhos opostos ao vertice em
// questao. Este metodo funciona com Semi-Lagrangian pois permite a
// procura de elementos considerando o vertice em questao. Uma descricao
// completa das matrizes eh definida dentro do proprio metodo.
void Model3D::setOFace()
{
 clMatrix mapViz(numElems,numElems);
 clMatrix faceFace(1000,numGLEP+1);
 //clMatrix freeFace(1000,numGLEP+1);
 freeFace.Dim(1000,numGLEP+1);
 clMatrix faceFaceAux;
 clMatrix freeFaceAux;
 clMatrix mapVizAux(numElems*numGLEP,numGLEP);
 clMatrix mapVizAux2(numElems*numGLEP,numGLEP);
 clMatrix comb(numGLEP,numGLEP-1);  // triangulo = 2
 clVector verts(numGLEP);           // tetraedro = 3
 clVector face(numGLEP-1);
 clVector index(numElems*numGLEP);
 clVector index2(numElems*numGLEP);
 clVector index3(numElems*numGLEP);
 clVector idcol(numElems*numGLEP);
 clVector idcol2(numElems*numGLEP);
 clVector idcol3(numElems*numGLEP);
 clVector vect1(numGLEP+1);
 clVector vert(numGLEP-1);

 //        - nome: comb
 //        - definicao:  matriz de aresta/face para um elemento de referencia
 //
 //            +---+---+---+
 //            | a | a | a |   a = identificacao dos vertices
 //  ---   +---+---+---+---+
 //   |    | b |   |   |   |   b = identificacao das faces   
 //   |    +---+---+---+---+
 //   |    | b |   |   |   |   c = quantidade de faces    
 //   c    +---+---+---+---+
 //   |    | b |   |   |   |   d = quantidade de vertices para formar face
 //   |    +---+---+---+---+
 //   |    | b |   |   |   |
 //  ---   |---+---+---+---+
 //
 //        |______ d ______|
 //        |               |
 //

 comb.Set(0,0,0);
 comb.Set(0,1,1);
 comb.Set(0,2,2);
 comb.Set(1,0,0);
 comb.Set(1,1,1);
 comb.Set(1,2,3);
 comb.Set(2,0,1);
 comb.Set(2,1,2);
 comb.Set(2,2,3);
 comb.Set(3,0,0);
 comb.Set(3,1,2);
 comb.Set(3,2,3);

 int k = 0;
 neighbourElem.resize (0);
 neighbourElem.resize (numVerts);
 for( int i=0;i<numElems;i++ )
 {
  for( int j=0;j<numGLEP;j++ )
  {
   neighbourElem.at( (int)IEN.Get(i,j) ).push_back(i);
  }

  verts.CopyRow(i,IEN);
  for( int j=0;j<numGLEP;j++ )
  {
   face.Set(0,comb.Get(j,0));
   face.Set(1,comb.Get(j,1));
   face.Set(2,comb.Get(j,2));
   vert.Set(0, (int) verts.Get( (int) face.Get(0) ) );
   vert.Set(1, (int) verts.Get( (int) face.Get(1) ) );
   vert.Set(2, (int) verts.Get( (int) face.Get(2) ) );
   vert.Sort();
   mapVizAux.Set(k,0,vert.Get(0));
   mapVizAux.Set(k,1,vert.Get(1));
   mapVizAux.Set(k,2,vert.Get(2));
   mapVizAux.Set(k,3,i);
   k++;
  }
 }

//--------------------------------------------------
//  // ------------------ HELP
//  int v1,v2,v3,v4;
//  list<int> plist5;
//  list<int>::iterator mele5;
//  for( int ii=0;ii<numVerts;ii++ )
//  {
//   if( ii== 675 )
//   {
//    plist5 = neighbourElem.at(ii);
//    cout << "surfaceNode " << ii << endl;
//    for( mele5=plist5.begin(); mele5 != plist5.end();++mele5 )
//    {
// 	v1 = (int) IEN.Get(*mele5,0);
// 	v2 = (int) IEN.Get(*mele5,1);
// 	v3 = (int) IEN.Get(*mele5,2);
// 	v4 = (int) IEN.Get(*mele5,3);
// 	cout << v1 << " " << cc.Get(v1) << " " 
// 	     << v2 << " " << cc.Get(v2) << " "
// 	     << v3 << " " << cc.Get(v3) << " "
// 	     << v4 << " " <<  cc.Get(v4) << " " << *mele5 << endl;
//    }
//    cout << endl;
//   }
//  }
//-------------------------------------------------- 
//--------------------------------------------------
//  for( int i=0;i<numVerts;i++ )
//  {
//   int surfaceNode = i;
//   if( surfaceNode == 659 )
//   {
//    cout << "---------" << surfaceNode << "------------" << endl;
//    std::ostream_iterator< int > output( cout, " " );
//    std::copy( neighbourElem.at(surfaceNode).begin(), neighbourElem.at(surfaceNode).end(), output );
//    cout << endl;
//   }
//  }
//  // ------------------ HELP
//-------------------------------------------------- 
 




 //cout << "mapVizAux: " << endl;
 //mapVizAux.Display();
 //cout << " -------------- " << endl;

 clVector faceaux(3);
 IFACE3D *faces = NULL;
 faces = new IFACE3D[4*(int)mapVizAux.DimI()];
 for( int i=0;i<mapVizAux.DimI();i++ )
 {
  int v1 = (int) mapVizAux.Get(i,0);
  int v2 = (int) mapVizAux.Get(i,1);
  int v3 = (int) mapVizAux.Get(i,2);
  int v4 = (int) mapVizAux.Get(i,3);

  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Set(2,v3);
  faceaux.Sort(); // ordena a linha
  faces[i].p1 = (int)faceaux(0);
  faces[i].p2 = (int)faceaux(1);
  faces[i].p3 = (int)faceaux(2);
  faces[i].p4 = v4;
 }

 // rotina em C que compara vertices de 2 elementos, ordenando-os
 qsort(faces,(int)mapVizAux.DimI(),sizeof(IFACE3D),IFACE3DCompare);

 for( int i=0;i<(int) mapVizAux.DimI();i++ )
 {
  mapVizAux2.Set(i,0,faces[i].p1);
  mapVizAux2.Set(i,1,faces[i].p2);
  mapVizAux2.Set(i,2,faces[i].p3);
  mapVizAux2.Set(i,3,faces[i].p4);
 }

 //        - nome: mapViz
 //        - definicao: matrix com mapeamento de vizinhos de cada
 //                     elemento
 //  					Cada linha da matriz so podera apresentar no 
 //  					maximo o numero de vizinhos de cada elemento, no
 //  					caso do tetraedro os valores diferentes de zero por
 //  					linha nao ultrapassa de 3. A estrutura desta matriz
 //  					deve ser esparsa.
 //
 //              +---+---+---+ ... +---+
 //              | a | a | a |     | a |
 //    ---   +---+---+---+---+ ... +---+
 //     |    | b | e | e | e |     | e |  a = identificacao do elemento
 //     |    +---+---+---+---+ ... +---+
 //     |    .   .   .   .   . ... .   .  b = identificacao do elemento viz.
 //     |	 .   .   .   .   . ... .   .
 //     c    .   .   .   .   . ... .   .  c = numero total de elementos
 //     |    +---+---+---+---+ ... +---+                               
 //     |    | b | e | e | e |     | e |  d = numero total de elementos
 //     |    +---+---+---+---+ ... +---+
 //     |    | b | e | e | e |     | e |  e = identificacao da aresta
 //    ---   +---+---+---+---+ ... +---+
 //
 //          |___________ d ___________|
 //   	     |                         |
 //
 
 //          - nome: faceFace
 //          - definicao: matrix com mapeamento de faces vizinhas de cada
 //   	                elemento. As arestas dobradas sao retiradas da
 //   					matriz mapViz. Se apresenta aresta dobrada quer 
 //   					dizer que 2 elementos possuem a mesma aresta, 
 //   					entao sao elementos vizinhos
 //
 //              +---+---+---+---+---+
 //              | a | b | c | d | e |  a = identificacao do 1o. elemento
 //    ---   +---+---+---+---+---+---+  b = identificacao do 2o. elemento
 //     |    | f |   |   |   |   |   |  c = 1o. vertice da face em comum
 //     |    +---+---+---+---+---+---+  d = 2o. vertice da face em comum
 //     |    .   .   .   .   .   .   .  e = 3o. vertice da face comum
 //     |	 .   .   .   .   .   .   .  f = ident da face dobrada
 //     g    .   .   .   .   .   .   .
 //     |    +---+---+---+---+---+---+
 //     |    | f |   |   |   |   |   |
 //     |    +---+---+---+---+---+---+  g = numero de arestas dobradas
 //     |    | f |   |   |   |   |   |  h = 2 elem + 3 vertices de aresta = 5
 //    ---   +---+---+---+---+---+---+
 //
 //          |________ h ________|
 //          |                   |
 //
 
 //        - nome: freeFace
 //          - definicao: matrix com mapeamento de faces de fronteira 
 //   	                de cada elemento, pois nao possuem arestas 
 //   					dobradas. Como a identificacao segue o 
 //   					padram da matriz faceFace, a 1a coluna eh igual a 
 //   					zero pois em condicao de contorno o elemento nao
 //   					tem vizinho
 //
 //              +---+---+---+---+---+
 //              | a | b | c | d | e |  a = identificacao do 1o. elemento = 0
 //    ---   +---+---+---+---+---+---+  b = identificacao do 2o. elemento
 //     |    | f |   |   |   |   |   |  c = 1o. vertice da face em comum
 //     |    +---+---+---+---+---+---+  d = 2o. vertice da face em comum
 //     |    .   .   .   .   .   .   .  e = 3o. vertice da face comum
 //     |	 .   .   .   .   .   .   .  f = ident da face dobrada
 //     g    .   .   .   .   .   .   .
 //     |    +---+---+---+---+---+---+
 //     |    | f |   |   |   |   |   |
 //     |    +---+---+---+---+---+---+  g = numero de arestas dobradas
 //     |    | f |   |   |   |   |   |  h = 2 elem + 3 vertices de aresta = 5
 //    ---   +---+---+---+---+---+---+
 //
 //          |________ h ________|
 //          |                   |
 //
 
 // procura de elementos que apresentam 3 vertices iguais -> faceFace
 // procura de elementos que nao possuem vizinhos -> freeFace
 mapVizAux = mapVizAux2;
 int iFace = 0;
 int iFree = 0;
 for ( int ii=0;ii<mapVizAux.DimI()-1;ii++ )
 {
  if( (mapVizAux.Get(ii,0)==mapVizAux.Get(ii+1,0)) &&
	  (mapVizAux.Get(ii,1)==mapVizAux.Get(ii+1,1)) &&  
	  (mapVizAux.Get(ii,2)==mapVizAux.Get(ii+1,2)) )  
  {
   // a matriz faceFace eh pre-alocada com 1000 linhas, caso ultrapasse
   // esse valor realoca-se para uma dimensao maior
   if( iFace == faceFace.DimI() ) 
   {
	faceFaceAux = faceFace;
	faceFace.Dim(iFace+faceFaceAux.DimI(),numGLEP+1);
	faceFace.CopyFrom(0,0,faceFaceAux);
   }
   faceFace.Set(iFace,0,mapVizAux.Get(ii+1,3));
   faceFace.Set(iFace,1,mapVizAux.Get(ii,3));
   faceFace.Set(iFace,2,mapVizAux.Get(ii,0));
   faceFace.Set(iFace,3,mapVizAux.Get(ii,1));
   faceFace.Set(iFace,4,mapVizAux.Get(ii,2));
   ii++; // pois existem 2 faces
   iFace++;
  } 
  if((ii==0) || 
	(   (ii!=0) && ( mapVizAux.Get(ii,0)!=mapVizAux.Get(ii-1,0) ||
					 mapVizAux.Get(ii,1)!=mapVizAux.Get(ii-1,1) ||     
					 mapVizAux.Get(ii,2)!=mapVizAux.Get(ii-1,2) ) ))     
  {    
   // a matriz freeFace eh pre-alocada com 1000 linhas, caso ultrapasse
   // esse valor realoca-se para uma dimensao maior
   if( iFree == freeFace.DimI() ) 
   {
	freeFaceAux = freeFace;
	freeFace.Dim(iFree+freeFaceAux.DimI(),numGLEP+1);
	freeFace.CopyFrom(0,0,freeFaceAux);
   }
   freeFace.Set(iFree,0,0);
   freeFace.Set(iFree,1,mapVizAux.Get(ii,3));
   freeFace.Set(iFree,2,mapVizAux.Get(ii,0));
   freeFace.Set(iFree,3,mapVizAux.Get(ii,1));
   freeFace.Set(iFree,4,mapVizAux.Get(ii,2));
   iFree++;
  }
 }
 int ii = mapVizAux.DimI()-1;
 if( ( mapVizAux.Get(ii,0)!=mapVizAux.Get(ii-1,0) ||
	   mapVizAux.Get(ii,1)!=mapVizAux.Get(ii-1,1) ||     
	   mapVizAux.Get(ii,2)!=mapVizAux.Get(ii-1,2) ) )     
 {
   if( iFree == freeFace.DimI() ) 
   {
	freeFaceAux = freeFace;
	freeFace.Dim(iFree+freeFaceAux.DimI(),numGLEP+1);
	freeFace.CopyFrom(0,0,freeFaceAux);
   }
   freeFace.Set(iFree,0,0);
   freeFace.Set(iFree,1,mapVizAux.Get(ii,3));
   freeFace.Set(iFree,2,mapVizAux.Get(ii,0));
   freeFace.Set(iFree,3,mapVizAux.Get(ii,1));
   freeFace.Set(iFree,4,mapVizAux.Get(ii,2));
   iFree++;
 }

 // como as matrizes faceFace e freeFace possuem dimensoes variaveis,
 // precisa-se aloca-las com um numero maior que o necessario. Apos suas
 // atribuicoes, necessita-se dimensiona-las para o tamanho final
 // redimensionalizando...
 faceFaceAux = faceFace;
 freeFaceAux = freeFace;
 faceFace.Dim(iFace,numGLEP+1);
 freeFace.Dim(iFree,numGLEP+1);
 faceFaceAux.CopyTo(0,0,faceFace);
 freeFaceAux.CopyTo(0,0,freeFace);

 mapVizAux.Dim(0,0);
 mapVizAux2.Dim(0,0);

 //        - nome: oFace
 //        - definicao: matrix com mapeamento de elementos opostos ao 
 //   	                vertice em questao. Se encontrar o valor (-1) quer
 //   					dizer que o vertice nao apresenta elemento oposto,
 //   					ou seja, eh um elemento de fronteira
 //
 //              +---+---+---+---+
 //              | a | a | a | a |  a = identificacao dos vertices do elemento
 //    ---   +---+---+---+---+---+      seguindo a ordem local 0-1-2-3
 //     |    | b | e | e | e | e |  
 //     |    +---+---+---+---+---+  b = idendificacao dos elementos
 //     |    .   .   .   .   .   .  
 //     |	 .   .   .   .   .   .  e = identificacao do elemento oposto ao 
 //     c    .   .   .   .   .   .      vertice
 //     |    +---+---+---+---+---+
 //     |    | b | e | e | e | e |
 //     |    +---+---+---+---+---+  c = numero total de elementos
 //     |    | b | e | e | e | e |  d = numero de total de vertices do elemento 
 //    ---   +---+---+---+---+---+
 //
 //          |________ d ________|
 //     	 |                   |
 
 oFace.Dim(numElems,numGLEP);
 for( int i=0;i<numElems;i++ )
  for( int j=0;j<numGLEP;j++ )
   oFace.Set(i,j,-1);

 clVector tetra(IEN.DimJ());

 for( int ii=0;ii<faceFace.DimI();ii++ )
 {
  int elem1 = (int) faceFace.Get(ii,0);
  int elem2 = (int) faceFace.Get(ii,1);
  vert.Set(0, (int) faceFace.Get(ii,2) );
  vert.Set(1, (int) faceFace.Get(ii,3) );
  vert.Set(2, (int) faceFace.Get(ii,4) );

  // armazena o numero da aresta + 1, para manter a matriz esparsa
  mapViz.Set(elem1,elem2,ii+1);
  mapViz.Set(elem2,elem1,ii+1);

  tetra.CopyRow(elem1,IEN);

  for( int kk=0;kk<numGLEP;kk++)
   if( (tetra.Get(kk) != vert.Get(0)) && 
	   (tetra.Get(kk) != vert.Get(1)) && 
	   (tetra.Get(kk) != vert.Get(2)) )
   {
	oFace.Set(elem1,kk,elem2);
	break;
   }

  tetra.CopyRow(elem2,IEN);
  for ( int kk=0;kk<numGLEP;kk++ )
   if( (tetra.Get(kk) != vert.Get(0)) && 
	   (tetra.Get(kk) != vert.Get(1)) && 
	   (tetra.Get(kk) != vert.Get(2)) )
   {
	oFace.Set(elem2,kk,elem1);
	break;
   }
 }

 // correcao da oFace (funcionando para malha step40-20-10)
 int vec[3];
 int v[4];
 int aresta = 0;
 int jAresta = 0;
 int vertOp1 = 0;
 int vertOp2 = 0;
 int kFace = 0;
 int nele;
 int count;
 list<int> plist;
 list<int>::iterator mele; // definicao do iterador
 for( int i=0;i<freeFace.DimI();i++ )
 {
  nele = (int) freeFace.Get(i,1);
  vec[0] = (int) freeFace.Get(i,2); 
  vec[1] = (int) freeFace.Get(i,3); 
  vec[2] = (int) freeFace.Get(i,4); 
  for( int j=0;j<numGLEP-1;j++ )
  {
   ii = (int) freeFace.Get(i,j+2);
   plist = neighbourElem.at(ii);
   for (mele=plist.begin(); mele != plist.end(); ++mele)
   {
	v[0] = (int) IEN.Get(*mele,0);
	v[1] = (int) IEN.Get(*mele,1);
	v[2] = (int) IEN.Get(*mele,2);
	v[3] = (int) IEN.Get(*mele,3);
	count = 0;
	for( int m=0;m<3;m++ )
	 for( int n=0;n<4;n++ )
	  if( (vec[m] == v[n]) && (vec[m] != ii) ){ aresta = vec[m];count++; }
	  
	if( count == 1 ) // aresta localizada!
	{
	 for( int m=0;m<3;m++ )
	  if( vec[m] == aresta ) jAresta = m; 
	  else if( vec[m] != ii ) kFace = m; 
	 
	 for( int m=0;m<3;m++ )
	 {
	  if( (v[m] != ii) && (v[m] != aresta) )
	  {
	   if( testFace(ii,aresta,vec[kFace],v[m]) )
	   {
		for( int n=0;n<4;n++ )
		{
		 if( (IEN.Get(nele,n) != ii) && 
		     (IEN.Get(nele,n) != aresta) && 
			 (IEN.Get(nele,n) != vec[kFace]) )  
		  vertOp1 = n;
		}
		for( int n=0;n<4;n++ )
		{
		 if( (IEN.Get(*mele,n) != ii) && 
		     (IEN.Get(*mele,n) != aresta) && 
			 (IEN.Get(*mele,n) != v[m]) )  
		  vertOp2 = n;
		}
		if( oFace.Get(nele,vertOp1) == -1 )
		 oFace.Set(nele,vertOp1,*mele);

		if( oFace.Get(*mele,vertOp2) == -1 )
		 oFace.Set(*mele,vertOp2,nele);
	   }
	  }
	 }
	}
   }
  }
 }
 delete[] faces;

}

// aplica configuracoes referentes a superficie da modelagem 2 fases.
void Model3D::setSurfaceConfig()
{
 setVertNeighbour(); // neighbourVert
 setInOutVert(); // inVert e outVert
 setInOutElem(); // inElem e outElem
 setSurface(); // surface e nonSurface
 setSurfaceFace(); // elemSurface e neighbourFaceVert
 setSurfaceTri(); // triang superficie - interfaceMesh
 setConvexTri(); // triang parte externa do dominio - convexMesh
 //buildSurfMesh();

 setTriEdge(); 
 setNeighbourSurface(); 
 setTriangleMinEdge(); // minEdge check and config
 saveVTKConvex("./vtk/","conv",0);
 saveVTKSurface("./vtk/","surf",0);
}

bool Model3D::testFace(int v1, int v2, int v3, int v4)
{
 real V,Ax1,Ax2,Ay1,Ay2,Az1,Az2;
 real prodEsc;

 V = getVolume(v1,v2,v3,v4);

  if( fabs(V)<1e-10)
  {
     Az1=0.5*(((X.Get(v2)-X.Get(v1))*(Y.Get(v3)-Y.Get(v1)))
	         -((X.Get(v3)-X.Get(v1))*(Y.Get(v2)-Y.Get(v1))));

     Ay1=0.5*(((X.Get(v2)-X.Get(v1))*(Z.Get(v3)-Z.Get(v1)))
	         -((X.Get(v3)-X.Get(v1))*(Z.Get(v2)-Z.Get(v1))));

     Ax1=0.5*(((Z.Get(v2)-Z.Get(v1))*(Y.Get(v3)-Y.Get(v1)))
	         -((Z.Get(v3)-Z.Get(v1))*(Y.Get(v2)-Y.Get(v1))));

     Az2=0.5*(((X.Get(v2)-X.Get(v1))*(Y.Get(v4)-Y.Get(v1)))
	         -((X.Get(v4)-X.Get(v1))*(Y.Get(v2)-Y.Get(v1))));

     Ay2=0.5*(((X.Get(v2)-X.Get(v1))*(Z.Get(v4)-Z.Get(v1)))
	         -((X.Get(v4)-X.Get(v1))*(Z.Get(v2)-Z.Get(v1))));

     Ax2=0.5*(((Z.Get(v2)-Z.Get(v1))*(Y.Get(v4)-Y.Get(v1)))
	         -((Z.Get(v4)-Z.Get(v1))*(Y.Get(v2)-Y.Get(v1))));

	 prodEsc = Az1*Az2 + Ay1*Ay2 + Ax1*Ax2;
	 return (prodEsc>0);
  }
  else
   return false;
}

// monta um vetor de comprimento de arestas do elemento de menor area e 
// retorna o menor deltaX diferente de 0
// o valor de deltaX retornado nao e necessariamente o menor deltaX da malha,
// porem para malhas regulares o valor tende a ser o menor deltaX da malha
real Model3D::getDeltaXMin()
{
 int v[4];
 int elemIndex=(int) V.MinIndex();
 clVector deltaX(4);

 v[0]=(int) IEN.Get(elemIndex,0);
 v[1]=(int) IEN.Get(elemIndex,1);
 v[2]=(int) IEN.Get(elemIndex,2);
 v[3]=(int) IEN.Get(elemIndex,3);

 deltaX.Set(0,fabs(X.Get(v[0])-X.Get(v[1])));
 deltaX.Set(1,fabs(X.Get(v[1])-X.Get(v[2])));
 deltaX.Set(2,fabs(X.Get(v[2])-X.Get(v[3])));
 deltaX.Set(3,fabs(X.Get(v[3])-X.Get(v[0])));
 deltaX.Sort();

 if( deltaX.Get(0) == 0.0 )
  if( deltaX.Get(1) == 0.0 )
   return deltaX.Get(2);
  else 
   return deltaX.Get(1);
 else
  return deltaX.Get(0);
}

// monta um vetor de comprimento de arestas do elemento de menor area e 
// retorna o menor deltaY diferente de 0
// o valor de deltaY retornado nao e necessariamente o menor deltaY da malha,
// porem para malhas regulares o valor tende a ser o menor deltaY da malha
real Model3D::getDeltaYMin()
{
 int v[4];
 int elemIndex=(int) V.MinIndex();
 clVector deltaY(4);

 v[0] = (int) IEN.Get(elemIndex,0);
 v[1] = (int) IEN.Get(elemIndex,1);
 v[2] = (int) IEN.Get(elemIndex,2);
 v[3] = (int) IEN.Get(elemIndex,3);

 deltaY.Set(0,fabs(Y.Get(v[0])-Y.Get(v[1])));
 deltaY.Set(1,fabs(Y.Get(v[1])-Y.Get(v[2])));
 deltaY.Set(2,fabs(Y.Get(v[2])-Y.Get(v[3])));
 deltaY.Set(3,fabs(Y.Get(v[3])-Y.Get(v[0])));
 deltaY.Sort();

 if( deltaY.Get(0) == 0.0 )
  if( deltaY.Get(1) == 0.0 )
   return deltaY.Get(2);
  else 
   return deltaY.Get(1);
 else
  return deltaY.Get(0);
}

// monta um vetor de comprimento de arestas do elemento de menor area e 
// retorna o menor deltaZ diferente de 0
// o valor de deltaZ retornado nao e necessariamente o menor deltaZ da malha,
// porem para malhas regulares o valor tende a ser o menor deltaZ da malha
real Model3D::getDeltaZMin()
{
 int v[4];
 int elemIndex=(int) V.MinIndex();
 clVector deltaZ(4);

 v[0] = (int) IEN.Get(elemIndex,0);
 v[1] = (int) IEN.Get(elemIndex,1);
 v[2] = (int) IEN.Get(elemIndex,2);
 v[3] = (int) IEN.Get(elemIndex,3);

 deltaZ.Set(0,fabs(Y.Get(v[0])-Y.Get(v[1])));
 deltaZ.Set(1,fabs(Y.Get(v[1])-Y.Get(v[2])));
 deltaZ.Set(2,fabs(Y.Get(v[2])-Y.Get(v[3])));
 deltaZ.Set(3,fabs(Y.Get(v[3])-Y.Get(v[0])));
 deltaZ.Sort();

 if( deltaZ.Get(0) == 0.0 )
  if( deltaZ.Get(1) == 0.0 )
   return deltaZ.Get(2);
  else 
   return deltaZ.Get(1);
 else
  return deltaZ.Get(0);
}

real Model3D::getMaxAbsUC()
{
 clVector aux = uc.Abs();
 real r = aux.Max();
 return r;
}

real Model3D::getMinAbsUC()
{
 clVector aux = uc.Abs();
 real r = aux.Min();
 return r;
}

real Model3D::getMaxAbsVC()
{
 clVector aux = vc.Abs();
 real r = aux.Max();
 return r;
}

real Model3D::getMinAbsVC()
{
 clVector aux = vc.Abs();
 real r = aux.Min();
 return r;
}

real Model3D::getMaxAbsWC()
{
 clVector aux = wc.Abs();
 real r = aux.Max();
 return r;
}

real Model3D::getMinAbsWC()
{
 clVector aux = wc.Abs();
 real r = aux.Min();
 return r;
}

real Model3D::getVolume(int _v1,int _v2,int _v3,int _v4)
{
 return (-1.0/6.0) * (+1*( (X.Get(_v2)*Y.Get(_v3)*Z.Get(_v4)) 
	                      +(Y.Get(_v2)*Z.Get(_v3)*X.Get(_v4)) 
						  +(Z.Get(_v2)*X.Get(_v3)*Y.Get(_v4)) 
						  -(Y.Get(_v2)*X.Get(_v3)*Z.Get(_v4)) 
						  -(X.Get(_v2)*Z.Get(_v3)*Y.Get(_v4)) 
						  -(Z.Get(_v2)*Y.Get(_v3)*X.Get(_v4)) )
	         -X.Get(_v1)*( +Y.Get(_v3)*Z.Get(_v4)
		                   +Y.Get(_v2)*Z.Get(_v3) 
						   +Z.Get(_v2)*Y.Get(_v4)
						   -Y.Get(_v2)*Z.Get(_v4)
						   -Z.Get(_v3)*Y.Get(_v4) 
						   -Z.Get(_v2)*Y.Get(_v3) )
		     +Y.Get(_v1)*( +X.Get(_v3)*Z.Get(_v4)
	            	       +X.Get(_v2)*Z.Get(_v3)
						   +Z.Get(_v2)*X.Get(_v4)
						   -X.Get(_v2)*Z.Get(_v4)
						   -Z.Get(_v3)*X.Get(_v4) 
						   -Z.Get(_v2)*X.Get(_v3) )
		     -Z.Get(_v1)*( +X.Get(_v3)*Y.Get(_v4)
		                   +X.Get(_v2)*Y.Get(_v3) 
						   +Y.Get(_v2)*X.Get(_v4)
						   -X.Get(_v2)*Y.Get(_v4)
						   -Y.Get(_v3)*X.Get(_v4) 
						   -Y.Get(_v2)*X.Get(_v3) ) );
}

real Model3D::getVolume(int _elem)
{
 int v1 = (int)IEN.Get(_elem,0);
 int v2 = (int)IEN.Get(_elem,1);
 int v3 = (int)IEN.Get(_elem,2);
 int v4 = (int)IEN.Get(_elem,3);

 return (-1.0/6.0) * (+1*( (X.Get(v2)*Y.Get(v3)*Z.Get(v4)) 
	                      +(Y.Get(v2)*Z.Get(v3)*X.Get(v4)) 
						  +(Z.Get(v2)*X.Get(v3)*Y.Get(v4)) 
						  -(Y.Get(v2)*X.Get(v3)*Z.Get(v4)) 
						  -(X.Get(v2)*Z.Get(v3)*Y.Get(v4)) 
						  -(Z.Get(v2)*Y.Get(v3)*X.Get(v4)) )
	          -X.Get(v1)*( +Y.Get(v3)*Z.Get(v4)
		                   +Y.Get(v2)*Z.Get(v3) 
						   +Z.Get(v2)*Y.Get(v4)
						   -Y.Get(v2)*Z.Get(v4)
						   -Z.Get(v3)*Y.Get(v4) 
						   -Z.Get(v2)*Y.Get(v3) )
	    	  +Y.Get(v1)*( +X.Get(v3)*Z.Get(v4)
	         	   	       +X.Get(v2)*Z.Get(v3)
						   +Z.Get(v2)*X.Get(v4)
						   -X.Get(v2)*Z.Get(v4)
						   -Z.Get(v3)*X.Get(v4) 
						   -Z.Get(v2)*X.Get(v3) )
			 -Z.Get(v1)*( +X.Get(v3)*Y.Get(v4)
		                   +X.Get(v2)*Y.Get(v3) 
						   +Y.Get(v2)*X.Get(v4)
						   -X.Get(v2)*Y.Get(v4)
						   -Y.Get(v3)*X.Get(v4) 
						   -Y.Get(v2)*X.Get(v3) ) );
}

real Model3D::getArea(int _v1,int _v2,int _v3)
{
 // vectors
 real x1 = X.Get(_v2) - X.Get(_v1);
 real y1 = Y.Get(_v2) - Y.Get(_v1);
 real z1 = Z.Get(_v2) - Z.Get(_v1);

 real x2 = X.Get(_v3) - X.Get(_v1);
 real y2 = Y.Get(_v3) - Y.Get(_v1);
 real z2 = Z.Get(_v3) - Z.Get(_v1);

 real crossX = (y1*z2)-(z1*y2);
 real crossY = -( (x1*z2)-(z1*x2) );
 real crossZ = (x1*y2)-(y1*x2);

 return 0.5*sqrt( crossX*crossX+crossY*crossY+crossZ*crossZ );
}

real Model3D::getArea(int _elem)
{
 int v1=(int) surfMesh.IEN.Get(_elem,0);
 int v2=(int) surfMesh.IEN.Get(_elem,1);
 int v3=(int) surfMesh.IEN.Get(_elem,2);

 
 // vectors
 real x1 = X.Get(v2) - X.Get(v1);
 real y1 = Y.Get(v2) - Y.Get(v1);
 real z1 = Z.Get(v2) - Z.Get(v1);

 real x2 = X.Get(v3) - X.Get(v1);
 real y2 = Y.Get(v3) - Y.Get(v1);
 real z2 = Z.Get(v3) - Z.Get(v1);

 real crossX = (y1*z2)-(z1*y2);
 real crossY = -( (x1*z2)-(z1*x2) );
 real crossZ = (x1*y2)-(y1*x2);

 return 0.5*sqrt( crossX*crossX+crossY*crossY+crossZ*crossZ );
}

real Model3D::getLength(int _v1,int _v2)
{
  real x1 = X.Get(_v1);
  real y1 = Y.Get(_v1);
  real z1 = Z.Get(_v1);

  real x2 = X.Get(_v2);
  real y2 = Y.Get(_v2);
  real z2 = Z.Get(_v2);

  return sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );
}

real Model3D::getAreaHeron(int _elem)
{
 int v1=(int) surfMesh.IEN.Get(_elem,0);
 int v2=(int) surfMesh.IEN.Get(_elem,1);
 int v3=(int) surfMesh.IEN.Get(_elem,2);

 real a = sqrt( (X.Get(v2) - X.Get(v1))*(X.Get(v2) - X.Get(v1)) +
                (Y.Get(v2) - Y.Get(v1))*(Y.Get(v2) - Y.Get(v1)) +
				(Z.Get(v2) - Z.Get(v1))*(Z.Get(v2) - Z.Get(v1)) );

 real b = sqrt( (X.Get(v3) - X.Get(v1))*(X.Get(v3) - X.Get(v1)) +
                (Y.Get(v3) - Y.Get(v1))*(Y.Get(v3) - Y.Get(v1)) +
				(Z.Get(v3) - Z.Get(v1))*(Z.Get(v3) - Z.Get(v1)) );

 real c = sqrt( (X.Get(v3) - X.Get(v2))*(X.Get(v3) - X.Get(v2)) +
                (Y.Get(v3) - Y.Get(v2))*(Y.Get(v3) - Y.Get(v2)) +
				(Z.Get(v3) - Z.Get(v2))*(Z.Get(v3) - Z.Get(v2)) );

 real s = (a+b+c)/2.0;

 return sqrt( s*(s-a)*(s-b)*(s-c) );
}

void Model3D::clearBC()
{
 uc.Dim(numNodes);
 vc.Dim(numNodes);
 wc.Dim(numNodes);
 pc.Dim(numVerts);
 idbcu.Dim(0);
 idbcv.Dim(0);
 idbcw.Dim(0);
 idbcp.Dim(0);
 // nos metodos com Concentracao cc.Dim(numVerts) esta definido, nao
 // precisando defini-lo aqui.
 //cc.Dim(numVerts);
 outflow.Dim(numNodes,1); // usado no metodo Galerkin
}

// re-allocation of vectors and IEN matrix
void Model3D::reAllocStruct()
{
 clVector aux;
 aux = X;
 X.Dim(numNodes);
 X.CopyFrom(0,aux);
 aux = Y;
 Y.Dim(numNodes);
 Y.CopyFrom(0,aux);
 aux = Z;
 Z.Dim(numNodes);
 Z.CopyFrom(0,aux);
 clMatrix IENaux;
 IENaux = IEN;
 IEN.Dim(numElems,numGLEU); // 4 nos por elemento + 1 centroide
 IEN.CopyFrom(0,0,IENaux);
}

// this method checks the orientation of a triangular surface
// _surfaceNode = surface node
// _v1 = 1st vertice of triangle
// _v2 = 2nd vertice of triangle
// _v3 = vertice of tetrahedra inside the bubble that contanis 
// the triangle of surface
bool Model3D::checkNormal(int _surfaceNode,int _v1,int _v2,int _vIn)
{
 real xSurfaceNode = X.Get(_surfaceNode);
 real ySurfaceNode = Y.Get(_surfaceNode);
 real zSurfaceNode = Z.Get(_surfaceNode);
 real xV1 = X.Get(_v1);
 real yV1 = Y.Get(_v1);
 real zV1 = Z.Get(_v1);
 real xV2 = X.Get(_v2);
 real yV2 = Y.Get(_v2);
 real zV2 = Z.Get(_v2);
 real xIn = X.Get(_vIn);
 real yIn = Y.Get(_vIn);
 real zIn = Z.Get(_vIn);

 // vetor surfaceNode -> in
 real x1 = xIn - xSurfaceNode;
 real y1 = yIn - ySurfaceNode;
 real z1 = zIn - zSurfaceNode;
 
 // vetores surfaceNode -> (v1 e v2)
 real xVec1 = xV1 - xSurfaceNode; 
 real yVec1 = yV1 - ySurfaceNode; 
 real zVec1 = zV1 - zSurfaceNode; 
 real xVec2 = xV2 - xSurfaceNode; 
 real yVec2 = yV2 - ySurfaceNode; 
 real zVec2 = zV2 - zSurfaceNode; 

 // normal vector
 // produto vetorial: surfaceNode -> v1 X surfaceNode -> v2
 real x2 = (yVec1*zVec2)-(zVec1*yVec2);
 real y2 = -( (xVec1*zVec2)-(zVec1*xVec2) );
 real z2 = (xVec1*yVec2)-(yVec1*xVec2);

 // scalar product of the normal vector and surfaceNode->in to check if
 // they have the value < 90 degrees
 // produto escalar para saber se os vetores estao no mesmo sentido
 real prodEsc = x1*x2 + y1*y2 + z1*z2;

 if( prodEsc < 0 )
  return true;
 else
  return false;
}

void Model3D::moveXPoints(clVector &_vec,real _dt)
{
 X = X + _vec*_dt;

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  real aux = surfMesh.X.Get(surfaceNode)+(_vec.Get(surfaceNode)*_dt);
  surfMesh.X.Set(surfaceNode,aux);
 }
}

void Model3D::moveYPoints(clVector &_vec,real _dt)
{
 Y = Y + _vec*_dt;

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  real aux = surfMesh.Y.Get(surfaceNode)+(_vec.Get(surfaceNode)*_dt);
  surfMesh.Y.Set(surfaceNode,aux);
 }
}

void Model3D::moveZPoints(clVector &_vec,real _dt)
{
 Z = Z + _vec*_dt;

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  real aux = surfMesh.Z.Get(surfaceNode)+(_vec.Get(surfaceNode)*_dt);
  surfMesh.Z.Set(surfaceNode,aux);
 }
}


SurfaceMesh* Model3D::getSurfMesh(){ return &surfMesh; }
SurfaceMesh* Model3D::getInterfaceMesh(){ return &interfaceMesh; }
SurfaceMesh* Model3D::getConvexMesh(){ return &convexMesh; }
clVector* Model3D::getX(){ return &X; }
real Model3D::getMaxX(){ return X.Max(); }
real Model3D::getMinX(){ return X.Min(); }
void Model3D::setX(clVector _X){ X = _X; }
clVector* Model3D::getY(){ return &Y; }
real Model3D::getMinY(){ return Y.Min(); }
real Model3D::getMaxY(){ return Y.Max(); }
void Model3D::setY(clVector _Y){ Y = _Y; }
real Model3D::getMaxZ(){ return Z.Max(); }
real Model3D::getMinZ(){ return Z.Min(); }
clVector* Model3D::getZ(){ return &Z; }
void Model3D::setZ(clVector _Z){ Z = _Z; }
clVector* Model3D::getUC(){ return &uc; }
clVector* Model3D::getVC(){ return &vc; }
clVector* Model3D::getWC(){ return &wc; }
clVector* Model3D::getPC(){ return &pc; }
clVector* Model3D::getCC(){ return &cc; }
clVector* Model3D::getOutflow(){ return &outflow; }
clVector* Model3D::getIdbcu(){ return &idbcu; }
clVector* Model3D::getIdbcv(){ return &idbcv; }
clVector* Model3D::getIdbcw(){ return &idbcw; }
clVector* Model3D::getIdbcp(){ return &idbcp; }
clVector* Model3D::getIdbcc(){ return &idbcc; }
clMatrix* Model3D::getIEN(){ return &IEN; }
int Model3D::getNumVerts(){ return numVerts; }
int Model3D::getNumNodes(){ return numNodes; }
int Model3D::getNumElems(){ return numElems; }
int Model3D::getNumGLEU(){ return numGLEU; }
int Model3D::getNumGLEP(){ return numGLEP; }
int Model3D::getNumGLEC(){ return numGLEC; }
//clMatrix Model3D::getMapViz(){ return mapViz; }
//clMatrix Model3D::getFaceFace(){ return faceFace; }
clMatrix* Model3D::getOFace(){ return &oFace; }
real Model3D::getXCenter(){ return xCenter; }
real Model3D::getYCenter(){ return yCenter; }
real Model3D::getZCenter(){ return zCenter; }
real Model3D::getBubbleRadius(){ return bubbleRadius; }
clVector* Model3D::getSurface(){ return &surface; }
vector< list<int> >* Model3D::getNeighbourElem(){return &neighbourElem;}
vector< list<int> >* Model3D::getNeighbourVert(){return &neighbourVert;}
vector< list<int> >* Model3D::getNeighbourFace(){return &neighbourFace;}
vector< list<int> >* Model3D::getElemSurface(){return &elemSurface;}
vector< list<int> >* Model3D::getNeighbourFaceVert(){return &neighbourFaceVert;}
vector< list<int> >* Model3D::getSurfaceViz(){return &surfaceViz;}
vector< list<int> >* Model3D::getFaceIEN(){return &faceIEN;}
list<int>* Model3D::getInVert(){return &inVert;}
list<int>* Model3D::getOutVert(){return &outVert;}
list<int>* Model3D::getInElem(){return &inElem;}
list<int>* Model3D::getOutElem(){return &outElem;}

//-------------------------------------------------- 
// Atribui o Model3D do argumento no corrente
//-------------------------------------------------- 
void Model3D::operator=(Model3D &_mRight) 
{
  // ints and floats
  numVerts = _mRight.numVerts;
  numNodes = _mRight.numNodes;
  numGLEU = _mRight.numGLEU;
  numGLEP = _mRight.numGLEP;
  numGLEC = _mRight.numGLEC;
  rMax = _mRight.rMax;
  xCenter = _mRight.xCenter;
  yCenter = _mRight.yCenter;
  zCenter = _mRight.zCenter;
  bubbleRadius = _mRight.bubbleRadius;

  // clVector and clMatrix
  surface = _mRight.surface;
  nonSurface = _mRight.nonSurface;
  uc = _mRight.uc;
  vc = _mRight.vc;
  wc = _mRight.wc;
  pc = _mRight.pc;
  cc = _mRight.cc;
  X = _mRight.X;
  Y = _mRight.Y;
  Z = _mRight.Z;
  outflow = _mRight.outflow;
  idbcu = _mRight.idbcu;
  idbcv = _mRight.idbcv;
  idbcw = _mRight.idbcw;
  idbcp = _mRight.idbcp;
  idbcc = _mRight.idbcc;
  V = _mRight.V;
  IEN = _mRight.IEN;
  faceFace = _mRight.faceFace;
  freeFace = _mRight.freeFace;
  mapViz = _mRight.mapViz;
  oFace = _mRight.oFace;
  surfMesh = _mRight.surfMesh;
  interfaceMesh = _mRight.interfaceMesh;
  convexMesh = _mRight.convexMesh;

  // STL: list and vectors
  neighbourElem = _mRight.neighbourElem; 
  neighbourVert = _mRight.neighbourVert;
  neighbourFace = _mRight.neighbourFace;
  neighbourFaceVert = _mRight.neighbourFaceVert;
  elemSurface = _mRight.elemSurface;
  neighbourSurfaceElem = _mRight.neighbourSurfaceElem;
  neighbourPoint = _mRight.neighbourPoint;
  faceIEN = _mRight.faceIEN;
  surfaceViz = _mRight.surfaceViz;
  outVert = _mRight.outVert;
  inVert = _mRight.inVert;
  outElem = _mRight.outElem;
  inElem = _mRight.inElem;
}

void Model3D::saveVTKConvex( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = (string) _dir + (string) _filename + "TRI" + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename );

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "3D Simulation C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << endl;


 vtkFile << "POINTS " << convexMesh.numVerts << " double" << endl;
 for( int i=0;i<convexMesh.numVerts;i++ )
  vtkFile << X.Get(i) << " " 
          << Y.Get(i) << " " 
		  << Z.Get(i) << endl;

 vtkFile << endl;

 int numTri = convexMesh.numElems;

 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 for( int i=0;i<numTri;i++ )
 {
   vtkFile << "3 " << convexMesh.IEN.Get(i,0) << " "
	               << convexMesh.IEN.Get(i,1) << " "
				   << convexMesh.IEN.Get(i,2) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;

 vtkFile.close();
}

void Model3D::saveVTKSurface( const char* _dir,const char* _filename, int _iter )
{
 stringstream ss;  //convertendo int --> string
 string str;
 ss << _iter;
 ss >> str;

 string file = (string) _dir + (string) _filename + "TRI" + "-" + str + ".vtk";
 const char* filename = file.c_str();

 ofstream vtkFile( filename );

 vtkFile << "# vtk DataFile Version 1.0" << endl;
 vtkFile << "3D Simulation C++" << endl;
 vtkFile << "ASCII" << endl;
 vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;
 vtkFile << endl;


 vtkFile << "POINTS " << interfaceMesh.numVerts<< " double" << endl;
 for( int i=0;i<interfaceMesh.numVerts;i++ )
  vtkFile << interfaceMesh.X.Get(i) << " " 
          << interfaceMesh.Y.Get(i) << " " 
		  << interfaceMesh.Z.Get(i) << endl;

 vtkFile << endl;

 int numTri = interfaceMesh.numElems;

 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 for( int i=0;i<interfaceMesh.numElems;i++ )
 {
   vtkFile << "3 " << interfaceMesh.IEN.Get(i,0) << " "
	               << interfaceMesh.IEN.Get(i,1) << " "
				   << interfaceMesh.IEN.Get(i,2) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;

 vtkFile.close();
}

clVector Model3D::dsearchn(clVector _X,clVector _Y,clVector _Z,
                           clVector &_XI,clVector &_YI,clVector &_ZI)
{
 clVector vmin(_XI.Dim());
 vmin.SetAll(0);
 real dist, dmin;

 // loop sobre todos os pontos a serem procurados
 for( int i=0;i<_XI.Dim();i++ )
 {

  // distancia de um ponto qualquer x(0),y(o) para o ponto informado xi,yi
  dmin = (_XI.Get(i)-_X.Get(0))*(_XI.Get(i)-_X.Get(0))+
         (_YI.Get(i)-_Y.Get(0))*(_YI.Get(i)-_Y.Get(0))+
         (_ZI.Get(i)-_Z.Get(0))*(_ZI.Get(i)-_Z.Get(0));

  // loop sobre todos os pontos da malha (menos o ponto inicial) para
  // saber qual o ponto apresenta a menor distancia com relacao ao ponto
  // informado xi,yi
  for (int j=1;j<_X.Dim();j++)
  {
   //calculo da distancia de outro ponto da malha x(i),y(i)
   dist = (_XI.Get(i)-_X.Get(j))*(_XI.Get(i)-_X.Get(j))+
	      (_YI.Get(i)-_Y.Get(j))*(_YI.Get(i)-_Y.Get(j))+
	      (_ZI.Get(i)-_Z.Get(j))*(_ZI.Get(i)-_Z.Get(j));
   if (dist<dmin)
   {
	dmin=dist;
	vmin.Set(i,j); 
   }
  }
 }

 // retorna o indice do ponto da malha
 return vmin;
}
