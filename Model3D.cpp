// =================================================================== //
// this is file Model3D.cpp, created at 23-Ago-2007                    //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Model3D.h"

using namespace std;

Model3D::Model3D()
{
 numVerts = 0;                
 numElems = 0;                
 numNodes = 0;                
 dVerts = 0;                  
 rMax = 0;                    
 xCenter = 0;
 yCenter = 0;
 zCenter = 0;
 minEdge = 0.10;

 surfMesh.numVerts = 0;
 surfMesh.numElems = 0;
 surfMesh.numNodes = 0;
 surfMesh.numInterfaces = 0;
 surfMesh.numBoundaries = 0;
}

Model3D::Model3D(const Model3D &_mRight)
{
  // ints and floats
  numVerts = _mRight.numVerts;
  numNodes = _mRight.numNodes;
  numElems = _mRight.numElems;
  rMax = _mRight.rMax;
  xCenter = _mRight.xCenter;
  yCenter = _mRight.yCenter;
  zCenter = _mRight.zCenter;
  dVerts = _mRight.dVerts;                  
  averageTriLength = _mRight.averageTriLength;
  averageTriArea = _mRight.averageTriArea;
  averageTetVolume = _mRight.averageTetVolume;
  oper = _mRight.oper;
  opersurf = _mRight.opersurf;

  isp = _mRight.isp;
  ispc = _mRight.ispc;
  rsp = _mRight.rsp;        
  rspn = _mRight.rspn;        
  rspc = _mRight.rspc;
  flip = _mRight.flip;
  spc = _mRight.spc;
  spp = _mRight.spp;
  intet = _mRight.intet;
  maxArea = _mRight.maxArea;
  minArea = _mRight.minArea;
  idMaxArea = _mRight.idMaxArea;
  idMinArea = _mRight.idMinArea;
  minLength = _mRight.minLength;
  maxLength = _mRight.minLength;
  numSurfElems = _mRight.numSurfElems;
  numSurfVerts = _mRight.numSurfVerts;

  ip = _mRight.ip;                    
  ipd = _mRight.ipd;                    
  rp = _mRight.rp;              
  rpi = _mRight.rpi;                   
  rpd = _mRight.rpd;                   
  rpdist = _mRight.rpdist; 
  rph = _mRight.rph;                   
  rpv = _mRight.rpv;                   
  csp = _mRight.csp;                   
  maxVolume = _mRight.maxVolume;
  minVolume = _mRight.minVolume;
  idMaxVolume = _mRight.idMaxVolume;
  idMinVolume = _mRight.idMinVolume;

  // clVector and clMatrix
  surface = _mRight.surface;
  uc = _mRight.uc;
  vc = _mRight.vc;
  wc = _mRight.wc;
  pc = _mRight.pc;
  cc = _mRight.cc;
  heaviside = _mRight.heaviside;
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
  oFace = _mRight.oFace;
  surfMesh = _mRight.surfMesh;
  interfaceMesh = _mRight.interfaceMesh;
  vertIdRegion = _mRight.vertIdRegion;
  elemIdRegion = _mRight.elemIdRegion;
  triEdge = _mRight.triEdge;
  tetVol = _mRight.tetVol;
  edgeSize = _mRight.edgeSize;

  // STL: list and vectors
  initSurfaceVolume = _mRight.initSurfaceVolume;
  surfaceVolume = _mRight.surfaceVolume;
  dVolume = _mRight.dVolume;
  errorVolume = _mRight.errorVolume;
  initSurfaceArea = _mRight.initSurfaceArea;
  surfaceArea = _mRight.surfaceArea;
  dArea = _mRight.dArea;
  errorArea = _mRight.errorArea;
  neighbourElem = _mRight.neighbourElem; 
  neighbourVert = _mRight.neighbourVert;
  neighbourFace = _mRight.neighbourFace;
  neighbourSurfaceElem = _mRight.neighbourSurfaceElem;
  neighbourPoint = _mRight.neighbourPoint;
  faceIEN = _mRight.faceIEN;
  boundaryVert = _mRight.boundaryVert;
  inVert = _mRight.inVert;
  outElem = _mRight.outElem;
  inElem = _mRight.inElem;
}

Model3D::~Model3D(){}

void Model3D::readVTK( const char* filename )
{
 char auxstr[255];
 double coords[3];
 int i,j,k,idnv;

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
	}
   }
  }
 }
 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readVTKHeaviside( const char* filename )
{
 char auxstr[255];
 double fl;
 triEdge.clear();
 heaviside.Dim(numVerts);
 vertIdRegion.Dim(numVerts);
 elemIdRegion.Dim(numElems);

 ifstream vtkFile( filename,ios::in );

 if( !vtkFile )
 {
  cerr << "Esta faltando o arquivo de leitura de heaviside!" << endl;
  exit(1);
 }

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"CHARACTERISTICLENGTH") != 0) )
  vtkFile >> auxstr;

 int nRegions = 0;
 vtkFile >> auxstr; // 1
 vtkFile >> nRegions;
 vtkFile >> auxstr; // float
 triEdge.resize(nRegions);
 for( int nb=0;nb<nRegions;nb++ )
  vtkFile >> triEdge[nb];
 setTriEdge(triEdge);
 
 while( ( !vtkFile.eof())&&(strcmp(auxstr,"elemId") != 0) )
  vtkFile >> auxstr;

 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> auxstr;

 for( int i=0; i < numElems; i++ )
 {
  vtkFile >> fl;
  elemIdRegion.Set(i,fl);
 }

 while( (! vtkFile.eof())&&(strcmp(auxstr,"heaviside") != 0) )
  vtkFile >> auxstr;

 if( !vtkFile.eof() )
 {
  vtkFile >> auxstr;
  vtkFile >> auxstr;
  vtkFile >> auxstr;

  for( int i=0; i < numVerts; i++ )
  {
   vtkFile >> fl;
   heaviside.Set(i,fl);
  }
 }

 while( ( !vtkFile.eof())&&(strcmp(auxstr,"vertId") != 0) )
  vtkFile >> auxstr;

 vtkFile >> auxstr;
 vtkFile >> auxstr;
 vtkFile >> auxstr;

 for( int i=0; i < numVerts; i++ )
 {
  vtkFile >> fl;
  vertIdRegion.Set(i,fl);
 }

 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readVTKSurface( const char* filename )
{
 char auxstr[255];
 double coords[3];
 int i,j,k,idnv;

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
	}
   }
  }
 }
 vtkFile.close();
} // fim do metodo vtkRead

void Model3D::readMSH( const char* filename )
{
 char auxstr[255];
 double coords[3];
 int i,j,k,id;
 int elemNumber,type,numberOfTags;

 ifstream mshFile( filename,ios::in );

 if( !mshFile )
 {
  cerr << "Esta faltando o arquivo de Malha .msh!" << endl;
  exit(1);
 }

 int numberOfPhyNames;
 while( ( !mshFile.eof())&&(strcmp(auxstr,"$PhysicalNames") != 0) )
  mshFile >> auxstr;

 mshFile >> numberOfPhyNames;

 surfMesh.phyNames.resize(numberOfPhyNames);
 if( ( !mshFile.eof())&&(strcmp(auxstr,"$EndPhysicalNames") != 0) )
 {
  for (i=0; i < numberOfPhyNames; i++)
  {
   mshFile >> auxstr;
   mshFile >> id;
   id = id-1;
   mshFile >> auxstr;
   surfMesh.phyNames.at(id)=auxstr;
  }
 }

 /* How many boundaries and interfaces do we have?
  * 
  * numBoundaries = number of walls (ex. wallInflow, wallOutflow etc)
  * numInterfaces = number of interface (ex. bubble, bubble2 etc)
  *
  * */
 for( int i=0;i<numberOfPhyNames;i++ )
 {
  if( surfMesh.phyNames[i].find("bubble") == true )
   surfMesh.numInterfaces++;
  else
   surfMesh.numBoundaries++;
 }
 surfMesh.numInterfaces = numberOfPhyNames - surfMesh.numBoundaries;

 cout << endl;
 cout << "number of boundaries detected: " << surfMesh.numBoundaries 
      << endl;
 cout << "number of phases detected: " << surfMesh.numInterfaces+1
      << endl;
 cout << "number of interfaces: " << surfMesh.numInterfaces
      << endl;
 cout << endl;

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
  surfMesh.elemIdRegion.Dim(surfMesh.numElems);
  surfMesh.idRegion.Dim(surfMesh.numElems);

  for( i=0; i < surfMesh.numElems; i++ )
  {
   mshFile >> elemNumber;
   mshFile >> type; // 2-2D or 3-3D
   mshFile >> numberOfTags;  
   if( numberOfTags == 3 ) // msh file version 2.1
   {
	// surfMesh.elemIdRegion 0 = wall
	// surfMesh.elemIdRegion 1 = surface 1
	// surfMesh.elemIdRegion 2 = surface 2 (if it has more than 1 bubble)
	// surfMesh.elemIdRegion 3 = surface 3 (if it has more than 2 bubbles)
	mshFile >> id;
	id = id-1;
	surfMesh.idRegion.Set(i,id);
	mshFile >> auxstr;
	mshFile >> auxstr;
	if( surfMesh.phyNames.at(id).compare(1,4,"wall") == 0 )
	 surfMesh.elemIdRegion.Set(i,0);
	else 
	{
	 char buffer[10];
	 surfMesh.phyNames.at(id).copy(buffer,2,7);
	 int idBubble = atoi(buffer);
	 surfMesh.elemIdRegion.Set(i,idBubble);
	}
   }
   else // msh file vertion 2.2
   {
	mshFile >> id;
	id = id-1;
	surfMesh.idRegion.Set(i,id);
	mshFile >> auxstr;
	if( surfMesh.phyNames.at(id).compare(1,4,"wall") == 0 )
	 surfMesh.elemIdRegion.Set(i,0);
	else 
	{
	 char buffer[10];
	 surfMesh.phyNames.at(id).copy(buffer,2,7);
	 int idBubble = atoi(buffer);
	 surfMesh.elemIdRegion.Set(i,idBubble);
	}
   }

   for( j=0; j < type+1 ; j++ )
   {
	mshFile >> k;
	k=k-1; // elem .msh comecando com 1
	surfMesh.IEN.Set(i,j,k);
   }
  }
 }
 mshFile.close();

 // filling surfMesh.phyBounds
 surfMesh.phyBounds.clear();
 surfMesh.phyBounds.resize(surfMesh.numVerts);

 for( int i=0; i < surfMesh.numElems; i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  int id = surfMesh.idRegion.Get(i);

  string aux = surfMesh.phyNames.at(id);
  surfMesh.phyBounds.at(v1) = aux;
  surfMesh.phyBounds.at(v2) = aux;
  surfMesh.phyBounds.at(v3) = aux;
 }
} // fim do metodo readMsh

void Model3D::setInterfaceBC()
{
 surfMesh.vertIdRegion.Dim(surfMesh.numVerts);
 surfMesh.Marker.Dim(surfMesh.numVerts);

 // mudar para surfMesh.boundary
 boundaryVert.resize (0);

 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);

  surfMesh.vertIdRegion.Set(v1,surfMesh.elemIdRegion.Get(i));
  surfMesh.vertIdRegion.Set(v2,surfMesh.elemIdRegion.Get(i));
  surfMesh.vertIdRegion.Set(v3,surfMesh.elemIdRegion.Get(i));

  if( surfMesh.elemIdRegion.Get(i) > 0 )
  {
   double aux = 0.5;
   surfMesh.Marker.Set(v1,aux);
   surfMesh.Marker.Set(v2,aux);
   surfMesh.Marker.Set(v3,aux);
  }
  else
  {
   boundaryVert.push_back(v1);
   boundaryVert.push_back(v2);
   boundaryVert.push_back(v3);
  }
 }
 boundaryVert.sort();
 boundaryVert.unique();
//--------------------------------------------------
//  cout << "boundaryVert contains:";
//  for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
//   cout << " " << *it;
//  cout << endl;
//-------------------------------------------------- 
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
 double coords[9];
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
void Model3D::setMeshStep(int nX,int nY,int nZ,const char* _param)
{
 // clean and init tetgen mesh object
 in.initialize();
 out.initialize();

 in.mesh_dim = 3;
 in.numberofpoints = nX*nY*nZ;
 in.pointlist = new REAL[in.numberofpoints * 3];

 double dx = nX/(nX-1.0);
 double dy = nY/(nY-1.0);
 double dz = nZ/(nZ-1.0);
 double aux;
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

 /*
  * Q: Quiet: No terminal output except errors.
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                | meshing 3D points... ";
 tetrahedralize( (char*) _param,&in,&out );
 cout << "finished | " << resetColor() << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
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

 // boundary surface configuration
 surfMesh.numVerts = out.numberofpoints;
 surfMesh.numElems = out.numberoftrifaces;
 surfMesh.IEN.Dim(surfMesh.numElems,3);
 surfMesh.X.Dim(surfMesh.numVerts);
 surfMesh.Y.Dim(surfMesh.numVerts);
 surfMesh.Z.Dim(surfMesh.numVerts);
 for(int i=0; i<out.numberoftrifaces; i++ )
 {
  int v1 = out.trifacelist[3*i+0];
  int v2 = out.trifacelist[3*i+1];
  int v3 = out.trifacelist[3*i+2];
  surfMesh.IEN.Set(i,0,v1);
  surfMesh.IEN.Set(i,1,v3);
  surfMesh.IEN.Set(i,2,v2);
  surfMesh.X.Set(v1,out.pointlist[3*v1+0]);
  surfMesh.Y.Set(v1,out.pointlist[3*v1+1]);
  surfMesh.Z.Set(v1,out.pointlist[3*v1+2]);
  surfMesh.X.Set(v2,out.pointlist[3*v2+0]);
  surfMesh.Y.Set(v2,out.pointlist[3*v2+1]);
  surfMesh.Z.Set(v2,out.pointlist[3*v2+2]);
  surfMesh.X.Set(v3,out.pointlist[3*v3+0]);
  surfMesh.Y.Set(v3,out.pointlist[3*v3+1]);
  surfMesh.Z.Set(v3,out.pointlist[3*v3+2]);
 }
 surfMesh.vertIdRegion.Dim(numVerts,0.0);
 surfMesh.elemIdRegion.Dim(numElems,0.0);
 surfMesh.idRegion.Dim(numElems,0.0);
 surfMesh.Marker.Dim(numVerts,0.0);
 surfMesh.numInterfaces = 0;
 //surfMesh.numBoundaries 3;


 // elemIdRegion and vertIdRegion for single-phase flow is 0, however it
 // is required to be initialized.
 vertIdRegion.Dim(numVerts,0.0);
 elemIdRegion.Dim(numElems,0.0);
 heaviside.Dim(numVerts,0.0);

 in.initialize();
 out.initialize();
}

void Model3D::setMeshStep(int nX,int nY,int nZ)
{
 setMeshStep(nX,nY,nZ,"Q");
}

void Model3D::setStepBC()
{
#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max()) )
	uc.Set(i,1.0);
  }
  if( (Z.Get(i)==Z.Min()) || (Z.Get(i) == Z.Max()) )
  {
   idbcw.AddItem(i);

   wc.Set(i,0.0);
  }
  if( X.Get(i)==X.Max() && Y.Get(i)<Y.Max() && Y.Get(i)>Y.Min() )
   outflow.Set(i,0);
 }
 for( int i=0;i<numVerts;i++ )
 {
  if( X.Get(i)==X.Max() && Y.Get(i)<Y.Max() && Y.Get(i)>Y.Min() )
  {
   idbcp.AddItem(i);

   pc.Set(i,0.0);
  }
 }
}

void Model3D::setStepPBC()
{
#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) || (Z.Get(i)==Z.Min()) || (Z.Get(i)==Z.Max())  )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   wc.Set(i,0.0);
  }
  
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>Y.Min()) && (Y.Get(i)<Y.Max()) )
   {
	idbcu.AddItem(i);
	uc.Set(i,1.0);
   }

 }
 for( int i=0;i<numVerts;i++ )
 {
  if( X.Get(i)==X.Min() && Y.Get(i)>Y.Min() && Y.Get(i)<Y.Max() )
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
  }
 }
}

void Model3D::setWallStepBC()
{
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  if( X.Get(*it)==X.Max() && Y.Get(*it)<Y.Max() && Y.Get(*it)>Y.Min() &&
     (Z.Get(*it)>Z.Min()) && (Z.Get(*it)<Z.Max()) )
  {
   idbcp.AddItem(*it);

   outflow.Set(*it,0);
   pc.Set(*it,0.0);
  }
//--------------------------------------------------
//   else if( (Z.Get(*it)==Z.Min()) || (Z.Get(*it) == Z.Max()) )
//   {
//    idbcw.AddItem(*it);
// 
//    wc.Set(*it,0.0);
//   }
//-------------------------------------------------- 
  else if( (X.Get(*it)==X.Min()) && (Y.Get(*it)>(Y.Max()/2.0)) && 
	       (Y.Get(*it)<Y.Max()) && 
           (Z.Get(*it)>Z.Min()) && (Z.Get(*it)<Z.Max()) )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else  
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double aux = 0.0;
   uc.Set(*it,aux);
   vc.Set(*it,aux);
   wc.Set(*it,aux);
  }
 }
}

void Model3D::setCStepBC()
{
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  if( (X.Get(*it)==X.Min()) || (Y.Get(*it)==Y.Min()) || (Y.Get(*it)==Y.Max()) )
  {
   if( (X.Get(*it)==X.Min()) && (Y.Get(*it)>(Y.Max()/2.0)) && (Y.Get(*it)<Y.Max()) )
   {
	idbcc.AddItem(*it);

	cc.Set(*it,1.0);
   }
  }
 }
}

void Model3D::setStepReservBC()
{
#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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
#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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

void Model3D::setWallNormalVWBC()
{
	for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
	{
		if (surfMesh.phyBounds.at(*it) == "\"wallNormalV\"")
		{
			idbcv.AddItem(*it);
			vc.Set(*it,0.0);
		}
		if (surfMesh.phyBounds.at(*it) == "\"wallNormalW\"")
		{
			idbcw.AddItem(*it);
			wc.Set(*it,0.0);
		}

	}

}

void Model3D::setWallMovingPBC(double _velInf, double _velSup)
{
	for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
	{
		if (surfMesh.phyBounds.at(*it) == "\"wallMovingInf\"")
		{
			idbcu.AddItem(*it);
			idbcv.AddItem(*it);
			idbcw.AddItem(*it);

			uc.Set(*it,_velInf);
			vc.Set(*it,0.0);
			wc.Set(*it,0.0);
		}
		if (surfMesh.phyBounds.at(*it) == "\"wallMovingSup\"")
		{
			idbcu.AddItem(*it);
			idbcv.AddItem(*it);
			idbcw.AddItem(*it);

			uc.Set(*it,_velSup);
			vc.Set(*it,0.0);
			wc.Set(*it,0.0);
		}

	}

}

void Model3D::setCubeVortexBC()
{
#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
	
  if( (X.Get(i)==X.Max()) || (X.Get(i)==X.Min()) )
  {
   	//idbcu.AddItem(i);
  	//uc.Set(i,0.0);
  }
 
  if( (Y.Get(i)==Y.Max()) || (Y.Get(i)==Y.Min()) )
  {
   	idbcv.AddItem(i);
  	vc.Set(i,0.0);
  }

  if( (Z.Get(i)==Z.Max()) || (Z.Get(i)==Z.Min()) )
  {
   	idbcw.AddItem(i);
  	wc.Set(i,0.0);
  }

 }
}

void Model3D::setWallSlipEdgesBC()
{
  for ( int i = 0; i < numVerts; ++i )
  {
	double x = X.Get(i);
	double y = Y.Get(i);
	double z = Z.Get(i);

	double xm = X.Min();
	double ym = Y.Min();
	double zm = Z.Min();

	double xM = X.Max();
	double yM = Y.Max();
	double zM = Z.Max();
  
	if ( ( (x == xm) || (x == xM) ) && ( (z == zm) || (z == zM) ) && ( (y > ym) || (y < yM) ) )
	{
		idbcu.AddItem(i);
		idbcw.AddItem(i);
		uc.Set(i,0.0);
		wc.Set(i,0.0);
	}
	if ( ( (x == xm) || (x == xM) ) && ( (y == ym) || (y == yM) ) && ( (z > zm) || (z < zM) ) )
	{
		idbcu.AddItem(i);
		idbcv.AddItem(i);
		uc.Set(i,0.0);
		vc.Set(i,0.0);
	}

  }
}
void Model3D::setAdimenStep()
{
 double aux;
 double factor = 1.0/(Y.Max()-Y.Min());

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

void Model3D::setMeshDisk(int nLados1Poli,int nCircMax,int nZ,const
  char* _param)
{
 double aux;

 clVector xCirc(1);
 clVector yCirc(1);
 xCirc.Set(0,0.0);
 yCirc.Set(0,0.0);

 double dr = 1;
 double r = dr;
 double dl = ( (2*3.141592)/nLados1Poli)*dr;
 double theta = 0;
 double dTheta = 0;
 double points2d = 1;
 for( int nCirc=1;nCirc<=nCircMax;nCirc++ )
 {
  theta = 0.0;
  dTheta = (dl/nCirc)*dr;
  for( int k=0;k<(nLados1Poli*nCirc);k++ )
  {
   aux = r*cos(theta);
   xCirc.AddItem(aux);
   aux = r*sin(theta);
   yCirc.AddItem(aux);

   theta = theta + dTheta;
   points2d++;

   if( theta >= 2*3.141592 ) break;
  }
  r=r+dr;
 }

 numVerts = points2d*nZ;

 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);

 double z = 0;
 double points3d = 0;
 for( int jz=0;jz<nZ;jz++ )
 {
  double jzNorm = (double) jz/(nZ-1);

  z = jzNorm;                  // linear
  //z = jzNorm*jzNorm;           // quadratic
  //z = jzNorm*jzNorm*jzNorm;    // cubic
  //z = exp(jzNorm)-1;             // exponential
  //z = exp(4*jzNorm)-1; // exponential (6 points conc boundary layer)
   
  for( int jCirc=1;jCirc<=xCirc.Dim();jCirc++ )
  {
   aux = xCirc.Get(jCirc-1);
   X.Set(points3d,aux);
   aux = yCirc.Get(jCirc-1);
   Y.Set(points3d,aux);
   aux = z;
   Z.Set(points3d,aux);

   points3d++;
  }
 }

 // clean and init tetgen mesh object
 in.initialize();
 out.initialize();

 in.mesh_dim = 3;
 in.numberofpoints = numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];

 for( int i=0;i<numVerts;i++ )
 {
  in.pointlist[3*i+0] = X.Get(i);
  in.pointlist[3*i+1] = Y.Get(i);
  in.pointlist[3*i+2] = Z.Get(i);
 }

 /*
  * Q: Quiet: No terminal output except errors.
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                | meshing 3D points... ";
 tetrahedralize( (char*) _param,&in,&out );
 cout << "finished | " << resetColor() << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
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

 // boundary surface configuration
 surfMesh.numVerts = out.numberofpoints;
 surfMesh.numElems = out.numberoftrifaces;
 surfMesh.IEN.Dim(surfMesh.numElems,3);
 surfMesh.X.Dim(surfMesh.numVerts);
 surfMesh.Y.Dim(surfMesh.numVerts);
 surfMesh.Z.Dim(surfMesh.numVerts);
 for(int i=0; i<out.numberoftrifaces; i++ )
 {
  int v1 = out.trifacelist[3*i+0];
  int v2 = out.trifacelist[3*i+1];
  int v3 = out.trifacelist[3*i+2];
  surfMesh.IEN.Set(i,0,v1);
  surfMesh.IEN.Set(i,1,v3);
  surfMesh.IEN.Set(i,2,v2);
  surfMesh.X.Set(v1,out.pointlist[3*v1+0]);
  surfMesh.Y.Set(v1,out.pointlist[3*v1+1]);
  surfMesh.Z.Set(v1,out.pointlist[3*v1+2]);
  surfMesh.X.Set(v2,out.pointlist[3*v2+0]);
  surfMesh.Y.Set(v2,out.pointlist[3*v2+1]);
  surfMesh.Z.Set(v2,out.pointlist[3*v2+2]);
  surfMesh.X.Set(v3,out.pointlist[3*v3+0]);
  surfMesh.Y.Set(v3,out.pointlist[3*v3+1]);
  surfMesh.Z.Set(v3,out.pointlist[3*v3+2]);
 }
 surfMesh.vertIdRegion.Dim(numVerts,0.0);
 surfMesh.elemIdRegion.Dim(numElems,0.0);
 surfMesh.idRegion.Dim(numElems,0.0);
 surfMesh.Marker.Dim(numVerts,0.0);
 surfMesh.numInterfaces = 0;
 //surfMesh.numBoundaries 3;


 // elemIdRegion and vertIdRegion for single-phase flow is 0, however it
 // is required to be initialized.
 vertIdRegion.Dim(numVerts,0.0);
 elemIdRegion.Dim(numElems,0.0);
 heaviside.Dim(numVerts,0.0);
 
 // como nao ha regiao predefinida, este metodo nao funciona aqui!
 //mesh3d = convertTetgenToMesh3d(out);

 in.initialize();
 out.initialize();
}

void Model3D::setMeshDisk(int nLados1Poli,int nCircMax,int nZ)
{
 setMeshDisk(nLados1Poli,nCircMax,nZ,"Q");
}

/* Method to transform the coordinates of the disk problem to the
 * sphere.
 *
 * Obs.: this method keeps the edge equally distanced
 * */
void Model3D::transformDiskToSphere()
{
 // ndr = number of radius intervals or 2nd input in
 // setMeshDisk(...,ndr,...)
 int ndr = 0;
 for( int i=0;i<numVerts;i++ )
  if( Z.Get(i) == Z.Min() && Y.Get(i) == 0 && X.Get(i) > 0 )
   ndr++;

 double aux = 0;
 double maxRadius = (X.Max()+Y.Max())/2.0;
 for( int i=0;i<numNodes;i++ )
 {
  double radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );
  double height = Z.Get(i);
  double dr = maxRadius/ndr;
  double nr = radius/dr;
  double theta = nr*(3.14159265359/2.0)/ndr;

  aux = (sin(theta)*(maxRadius+height))*(X.Get(i)/(radius+1E-10));
  X.Set(i,aux);

  aux = (sin(theta)*(maxRadius+height))*(Y.Get(i)/(radius+1E-10));
  Y.Set(i,aux);

  double newRadius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );
  aux = sqrt( (maxRadius+height)*(maxRadius+height) - newRadius*newRadius );
  Z.Set(i,aux);
 }
}


void Model3D::mesh2Dto3D(const char* _param)
{
 // clean and init tetgen mesh object
 in.initialize();
 out.initialize();

 in.mesh_dim = 3;
 in.numberofpoints = surfMesh.numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 convertModel3DtoTetgen(in);
 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");

 /*
  * Q: Quiet: No terminal output except errors.
  * Y: mesh boundary preserving
  * YY: no mesh boundary points inserted
  * R: mesh coarsening
  * C: Checks the consistency of the final mesh.
  * A: Assigns attributes to identify tetrahedra in certain regions.
  * a: Applies a maximum tetrahedron volume constraint.
  * p:  Tetrahedralizes a picecwise linear complex
  * q: Quality mesh generation. Min radius-edge ratio may be specifyed
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                 | meshing surface to 3D domain... ";
 tetrahedralize( (char*) _param,&in,&out ); // no insertion of points
 cout << "finished | " << resetColor() << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 convertTetgenToModel3D(out);

 mesh3d = convertTetgenToMesh3d(out);

 in.initialize();
 out.initialize();
}

void Model3D::mesh2Dto3D()
{
 mesh2Dto3D("QYYApa");
}

/*
 * mapEdgeTri.Get(i,0) // tamanho da aresta
 * mapEdgeTri.Get(i,1) // numero do 1o. vertice da aresta
 * mapEdgeTri.Get(i,2) // numero do 2o. vertice da areata
 * mapEdgeTri.Get(i,3) // numero do 3o. vertice do 1o. elemento
 * mapEdgeTri.Get(i,4) // numero do 3o. vertice do 2o. elemento
 * mapEdgeTri.Get(i,5) // 1o. elemento
 * mapEdgeTri.Get(i,6) // 2o. elemento 
 *
 * */
void Model3D::setMapEdgeTri()
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
  double x1=surfMesh.X.Get(faces[j].p1);
  double y1=surfMesh.Y.Get(faces[j].p1);
  double z1=surfMesh.Z.Get(faces[j].p1);
  double x2=surfMesh.X.Get(faces[j].p2);
  double y2=surfMesh.Y.Get(faces[j].p2);
  double z2=surfMesh.Z.Get(faces[j].p2);
  double length = vectorLength(x1-x2,y1-y2,z1-z2);

  mapEdgeTri.Set(i,0,length); // tamanho da aresta
  mapEdgeTri.Set(i,1,faces[j].p1 ); // numero do 1o. vertice da aresta
  mapEdgeTri.Set(i,2,faces[j].p2 ); // numero do 2o. vertice da areata
  mapEdgeTri.Set(i,3,faces[j].p3 );   // numero do 3o. vertice do 1o. elemento
  mapEdgeTri.Set(i,4,faces[j+1].p3 ); // numero do 3o. vertice do 2o. elemento
  mapEdgeTri.Set(i,5,faces[j].p4 ); // 1o. elemento
  mapEdgeTri.Set(i,6,faces[j+1].p4 ); // 2o. elemento 

  j=j+2; // pois cada aresta eh dividida com apenas 2 elementos
 }
 delete[] faces;
}

int Model3D::findEdge(int _v1,int _v2)
{
 int aux=-1;
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  if( (mapEdgeTri.Get(i,1) == _v1 && mapEdgeTri.Get(i,2) == _v2) ||
	  (mapEdgeTri.Get(i,1) == _v2 && mapEdgeTri.Get(i,2) == _v1) )
  {
   aux = i;
   break;
  }
 }
 return aux;
}

void Model3D::insertPointsByLength(const char* _interpolation,
                                   double _param)
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int edge=0;edge<mapEdgeTri.DimI();edge++ )
 {
  // edge length
  double edgeLength = mapEdgeTri.Get(edge,0);
  int v1 = mapEdgeTri.Get(edge,1); // v1
  //int v2 = mapEdgeTri.Get(edge,2); // v2
  //int v3elem1 = mapEdgeTri.Get(edge,3); // v3elem1
  //int v3elem2 = mapEdgeTri.Get(edge,4); // v3elem2
  double vertID = surfMesh.vertIdRegion.Get(v1); 

//--------------------------------------------------
//   double curv1 = fabs(surfMesh.curvature.Get(v1));
//   double curv2 = fabs(surfMesh.curvature.Get(v2));
//-------------------------------------------------- 
//--------------------------------------------------
//   double maxCurv = max(curv1,curv2);
//   double erro = maxCurv*edgeLength;
//-------------------------------------------------- 
  // this works, but is not consistent!!! CHANGE IT SOON!
//--------------------------------------------------
//   double curv1 = fabs(surfMesh.curvature.Get(v1));
//   double curv2 = fabs(surfMesh.curvature.Get(v2));
//   double curv3 = fabs(surfMesh.curvature.Get(v3elem1));
//   double curv4 = fabs(surfMesh.curvature.Get(v3elem2));
//-------------------------------------------------- 

  // to avoid contraction at high curvature regions
  //bool curvTest = (curv1 < 40 && curv2 < 40 && curv3 < 40 && curv4 < 40);

  if( vertID > 0 &&
	  //Z.Get(v1) != Z.Min() && Z.Get(v2) != Z.Min() &&
	  //Z.Get(v1) != Z.Max() && Z.Get(v2) != Z.Max() &&
      //erro > 0.7 )
      //curvTest &&
	  edgeLength > _param*triEdge[vertID] ) 
  {
   insertSurfacePoint(edge,_interpolation);

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   isp[vertID]++;
   opersurf[vertID]++;
  }
 }
}

void Model3D::insertPointsByLength(const char* _interpolation)
{
 insertPointsByLength(_interpolation,1.4);
}

/* Method to surface points when the curvature of the points is higher
 * then a specific value. This method should not be used often due to
 * the excessive deletion of points. The example below shows an 
 * initial surface mesh that contains angles of 90 degrees 
 * (points 1 and 2), this method will delete such necessery points.
 *
 *      ------------------------------------------------
 *            1 o ----------o 
 *              |            \     flow
 *              |             )    ---->
 *              |            /
 *            2 o ----------o 
 *      ------------------------------------------------
 *
 *  OBS.: Maximum curvature that a edge can handle. (empiric)
 *
 *  edgeLength = 0.08 ---> curvature < 60
 *  edgeLength = 0.10 ---> curvature < 50
 *
 *  edgeLength*curv < 5.0
 * */
void Model3D::removePointsByCurvature()
{
 for( int surfaceNode=0;surfaceNode<surfMesh.numVerts;surfaceNode++ )
 {
  // Checking the largest neighbour edge of surfaceNode
  // maxEdgeLength is then compared to surfaceNode's curvature
  double maxEdgeLength = 0;
  int listSize = neighbourPoint.at(surfaceNode).size();
  list<int> plist = neighbourPoint.at(surfaceNode);
  list<int>::iterator vert=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   double P0x = surfMesh.X.Get(surfaceNode);
   double P0y = surfMesh.Y.Get(surfaceNode);
   double P0z = surfMesh.Z.Get(surfaceNode);

   int v1 = *vert;++vert;
   double P1x = surfMesh.X.Get(v1);
   double P1y = surfMesh.Y.Get(v1);
   double P1z = surfMesh.Z.Get(v1);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

   if( edgeLength > maxEdgeLength )
	maxEdgeLength = edgeLength;
  }

  // Still don't know about the curv value!!!
  // edge length
  double vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
  double curv = fabs(surfMesh.curvature.Get(surfaceNode));

  // if the mesh is highly distorted (Re=160) and there are no
  // tetrahedron elements enough inside the sphere, then (I don't know
  // why) the surface mesh is sucked in, thus forming strange vertex
  // with curvature smaller then -30, thus this can be avoided using the
  // IF below.
//--------------------------------------------------
//   if( vertID > 0 && 
// 	  (maxEdgeLength*curv > 4.0 || surfMesh.curvature.Get(surfaceNode) < -30) )
//-------------------------------------------------- 
//--------------------------------------------------
//   if( vertID > 0 && maxEdgeLength*curv > 4.0 )
//-------------------------------------------------- 
  if( vertID > 0 && curv > 80 )
  {
   cout << "----------------- " << color(none,red,black) 
	    << "removing vertex with curvature (" 
		<< resetColor() << curv
		<< color(none,red,black) 
		<< ") at (" 
		<< resetColor()
		<< surfMesh.vertIdRegion.Get(surfaceNode)
		<< color(none,red,black) 
		<< "): "
		<< resetColor() << surfaceNode << endl;
   //saveVTKSurface("./vtk/","deleteBefore",surfaceNode);

   // delete surfaceNode from surface, xSurface, ySurface, zSurface vectors
   // surface is not used to add/remove/flip elements before the remeshing
   // but it's used on removePointsByInterfaceDistance
   // to be implemented 

   removeSurfacePoint(surfaceNode);

   // updating curvature value
   setNormalAndKappa();

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   rspc[vertID]++;
   opersurf[vertID]++;
  }
 }
}

//--------------------------------------------------
// void Model3D::removePointsByCurvature()
// {
//  // number of removed surface points by Curvature
//  fill(rspc.begin(),rspc.end(),0);
// 
//  for( int i=0;i<surfMesh.numVerts;i++ )
//  {
//   // edge length
//   int surfaceNode = i;
//   double vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
//   double curv = fabs(surfMesh.curvature.Get(surfaceNode));
//   double edgeLength = triEdge[vertID];
//   // Still don't know about the curv value!!!
//   if( vertID > 0 && edgeLength*curv > 4.5 )
//   {
//    cout << "----------------- " << color(none,red,black) 
// 	    << "removing vertex with curvature (" 
// 		<< resetColor() << curv
// 		<< color(none,red,black) 
// 		<< ") at (" 
// 		<< resetColor()
// 		<< surfMesh.vertIdRegion.Get(surfaceNode)
// 		<< color(none,red,black) 
// 		<< "): "
// 		<< resetColor() << surfaceNode << endl;
//    //saveVTKSurface("./vtk/","deleteBefore",surfaceNode);
// 
//    // delete surfaceNode from surface, xSurface, ySurface, zSurface vectors
//    // surface is not used to add/remove/flip elements before the remeshing
//    // but it's used on removePointsByInterfaceDistance
//    // to be implemented 
// 
//    // marking the desired elements for deletion
//    list<int> plist = neighbourSurfaceElem.at(surfaceNode);
//    for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
// 	markSurfElemForDeletion(*mele);
// 
//    // deleting elements
//    deleteSurfaceElements();
// 
//    // after the deletion process it's mandatory to create new elements
//    // to fill the space left by the deleting process
//    //surfaceTriangulator(surfaceNode);
//    surfaceTriangulatorEarClipping(surfaceNode,
//                                   neighbourPoint.at(surfaceNode),
//                                   "no");
// 
//    // deleting X,Y and Z coordinate; deleting the point maker funcition
//    deleteSurfacePoint(surfaceNode);
//    
//    // update surface, edge matrix, surface neigh elems and points
//    restoreMappingArrays();
//    
//    // updating curvature value
//    setNormalAndKappa();
// 
//    saveVTKSurface("./vtk/","remCurv",rspc[vertID]);
//    rspc[vertID]++;
//   }
//  }
// }
//-------------------------------------------------- 

void Model3D::insertPointsByCurvature(const char* _interpolation)
{
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge length
  //double edgeLength = mapEdgeTri.Get(i,0);
  int v1 = mapEdgeTri.Get(i,1);
  int v2 = mapEdgeTri.Get(i,2);
  double curv1 = fabs(surfMesh.curvature.Get(v1));
  double curv2 = fabs(surfMesh.curvature.Get(v2));
  double vertID = surfMesh.vertIdRegion.Get(v1); 
  double length = mapEdgeTri.Get(i,0);
  //double erro = max(curv1,curv2);

  if( curv1*length > 2.5 || curv2*length > 2.5  )
  {
   cout << "----------- Inserting by curvature..." << endl;
   insertSurfacePoint(i,_interpolation);

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   ispc[vertID]++;
   opersurf[vertID]++;
  }
 }
}

void Model3D::insertPointsByInterfaceDistance(const char* _interpolation)
{
 double Ymax1=100;
 double Ymin1=-100;
 double Ymax2=-100;
 double Ymin2=100;

 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  // bubble 1 (Y<0)
  if( Y.Get(i) < 0 && heaviside.Get(i)==0.5 )
  {
   if(Y.Get(i)>Ymin1) Ymin1=Y.Get(i);
   if(Y.Get(i)<Ymax1) Ymax1=Y.Get(i);
  }
  // bubble 2 (Y>0)
  if( Y.Get(i) > 0 && heaviside.Get(i)==0.5 )
  {
   if(Y.Get(i)<Ymin2) Ymin2=Y.Get(i);
   if(Y.Get(i)>Ymax2) Ymax2=Y.Get(i);
  }
 }

 int ny = 4;
 double distMin = Ymin2-Ymin1;
 //double distMax = Ymax2-Ymax1;
 double dy = (Ymin2-Ymin1)/(ny-1);

 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge length
  double edgeLength = mapEdgeTri.Get(i,0);
  double aux = (distMin/2.0)+2*dy;
  //double aux = (distMin/2.0)+(distMax-distMin)/3.0;
  if( surfMesh.Marker.Get(mapEdgeTri.Get(i,1)) == 0.5 && 
	  surfMesh.Y.Get(mapEdgeTri.Get(i,1)) < 1.0*aux && 
	  surfMesh.Y.Get(mapEdgeTri.Get(i,1)) > -1.0*aux &&
	  edgeLength > 4*dy )
  {
   insertSurfacePoint(i,_interpolation);
  }
 }
}

/* triangulate the interface using triangles with sorted vertices.
 * This method picks the first vertex of the neighbourPoint list and
 * connect it to each pair of vertex to mesh locally the empty space.
 * Ex. neighbourPoint = 0 1 2 3 4 5 0
 *     tri1 = 0 1 2
 *     tri2 = 0 2 3
 *     tri3 = 0 3 4
 *     tri4 = 0 4 5
 *
 * _v is the removed vertex located in the center of the polyhedron
 * neighbourPoint
 *
 * OBS: this methods does not work for conncave polyhedrons
 * */
void Model3D::surfaceTriangulator(int _v)
{
 int listSize = neighbourPoint.at(_v).size();
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 int vert0 = *point;++point;

 for( int i=0;i<listSize-3;i++ )
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

/*
 * This Method mesh an empty space (where _v is the removed vertex)
 * using the ear clipping technique which uses every 3 consective
 * vertices, thus forming the triangulation.
 * Ex.: neighbourPoint = 0 1 2 3 4 5 0
 *      tri1 = 0 1 2
 *      tri2 = 2 3 4
 *      tri3 = 4 5 0
 *      tri4 = 0 2 4
 *
 * */
void Model3D::surfaceTriangulatorEarClipping(int _v,
                                             list<int> _list,
											 const char* _mode)
{
 int vert1,vert2,vert3;
 list<int> myList = _list;
 int listSize = myList.size();

 // listSize = number of points of polyhedron + 1 (to close the
 // polyhedron)
 // ex.: 0 1 2 3 4 5 0 (closed polyhedron)
 while( listSize>4 )
 {
  if( strcmp( _mode,"quality" ) == 0 )
  {
   clVector bestTri = triangleQuality(_v);
   vert1 = bestTri.Get(0);
   vert2 = bestTri.Get(1);
   vert3 = bestTri.Get(2);
  }
  else
  {
   list<int>::iterator point=myList.begin();
   vert1 = *point;++point;
   vert2 = *point;++point;
   vert3 = *point;
  }

  myList.remove(vert2); // removing 2nd vertice
  listSize--;

  // adding new element
  surfMesh.IEN.AddRow();
  int elem = surfMesh.IEN.DimI()-1;
  surfMesh.IEN.Set(elem,0,vert1);
  surfMesh.IEN.Set(elem,1,vert2);
  surfMesh.IEN.Set(elem,2,vert3);

  // add new elemIdRegion
  list<int> plist2 = neighbourSurfaceElem.at(_v);
  list<int>::iterator mele=plist2.begin();
  surfMesh.elemIdRegion.AddItem(surfMesh.elemIdRegion.Get(*mele));
  surfMesh.idRegion.AddItem(surfMesh.idRegion.Get(*mele));
  mele = plist2.end();

  surfMesh.numElems++;
 }
 // adding last remaning element
 list<int>::iterator point=myList.begin();
 vert1 = *point;++point;
 vert2 = *point;++point;
 vert3 = *point;
 // adding new element
 surfMesh.IEN.AddRow();
 int elem = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem,0,vert1);
 surfMesh.IEN.Set(elem,1,vert2);
 surfMesh.IEN.Set(elem,2,vert3);

  // add new elemIdRegion
 list<int> plist2 = neighbourSurfaceElem.at(_v);
 list<int>::iterator mele=plist2.begin();
 surfMesh.elemIdRegion.AddItem(surfMesh.elemIdRegion.Get(*mele));
 surfMesh.idRegion.AddItem(surfMesh.idRegion.Get(*mele));
 mele = plist2.end();

 surfMesh.numElems++;
}

/* input: vertex to be deleted
 *
 *
 *
 * output: vector with 3 nodes to create 1 element
 * description: this method check the quality of all ear triangles that
 * can be created from the resulting polyhedron after deletion of vertex
 * _v.
 * 
 * Problem to be solved:
 *
 *          0 *
 *            |\
 *			  | \
 *		    1 *  \
 *			  |   \
 *			  |    \
 *		    2 *-----* 3
 *
 * 0-1-2-3 are coplanar and 0-1-2 are aligned. This method will chose
 * the triangulation as following:
 * 
 * 1st. best triangle: 0-3-2
 * 2nd. and final triangle: 2-1-0 (invalid triangle!)
 *
 * This method should triangulate like:
 * 
 * 1st. 0-3-1
 * 2nd. 2-1-3
 * 
 * HOW CAN WE DO IT?
 */
clVector Model3D::triangleQuality(int _v)
{
 int listSize = neighbourPoint.at(_v).size();
 list<int> plist = neighbourPoint.at(_v);
 list<int>::iterator point=plist.begin();
 int vert1,vert2,vert3;
 double qMax = 0; 
 clVector bestTri(3);

 for( int i=0;i<listSize-1;i++ )
 {
  vert1 = *point;++point;
  vert2 = *point;++point;
  vert3 = *point;--point;

  double P1x = surfMesh.X.Get(vert1);
  double P1y = surfMesh.Y.Get(vert1);
  double P1z = surfMesh.Z.Get(vert1);

  double P2x = surfMesh.X.Get(vert2);
  double P2y = surfMesh.Y.Get(vert2);
  double P2z = surfMesh.Z.Get(vert2);

  double P3x = surfMesh.X.Get(vert3);
  double P3y = surfMesh.Y.Get(vert3);
  double P3z = surfMesh.Z.Get(vert3);

  double q = triangleQualityMeasure(P1x,P1y,P1z,
	                              P2x,P2y,P2z,
			    				  P3x,P3y,P3z);

  if( qMax<q )
  {
   qMax = q; 
   bestTri.Set(0,vert1);
   bestTri.Set(1,vert2);
   bestTri.Set(2,vert3);
  }
 }

 return bestTri;
}

void Model3D::deleteSurfacePoint(int _v)
{
 X.Delete(_v);
 Y.Delete(_v);
 Z.Delete(_v);
 heaviside.Delete(_v);
 curvature.Delete(_v);
 edgeSize.Delete(_v);
 numVerts--;

 surfMesh.X.Delete(_v);
 surfMesh.Y.Delete(_v);
 surfMesh.Z.Delete(_v);
 surfMesh.xNormal.Delete(_v);
 surfMesh.yNormal.Delete(_v);
 surfMesh.zNormal.Delete(_v);
 surfMesh.vertIdRegion.Delete(_v);
 surfMesh.phyBounds.erase(surfMesh.phyBounds.begin()+_v);

 // Eh preciso atualizar os valores de boundaryVert que sao maios que _v
 // caso _v < boundaryVert.Dimension
 // pois boundaryVert nao eh uma lista continua de valores.
//--------------------------------------------------
//  if( (unsigned int) _v < boundaryVert.max_size() )
//   setInterfaceBC();
//-------------------------------------------------- 

 surfMesh.Marker.Delete(_v);
 surfMesh.numVerts--;

 // updating surfMesh.IEN
 for( int i=0;i<surfMesh.IEN.DimI();i++ )
  for( int j=0;j<surfMesh.IEN.DimJ();j++ )
   if( surfMesh.IEN.Get(i,j)>_v )
	surfMesh.IEN.Set(i,j,surfMesh.IEN.Get(i,j)-1);
}

void Model3D::markSurfElemForDeletion(int _elem)
{
 surfMesh.IEN.Set(_elem,0,-1);
 surfMesh.IEN.Set(_elem,1,-1);
 surfMesh.IEN.Set(_elem,2,-1);
}

void Model3D::deleteSurfaceElements()
{
 // deleting elements and elemIdRegion
 for( int i=0;i<surfMesh.IEN.DimI();i++ )
  if( surfMesh.IEN.Get(i,0) == -1 && 
	  surfMesh.IEN.Get(i,1) == -1 && 
	  surfMesh.IEN.Get(i,2) == -1 )
  {
   surfMesh.IEN.DelLine(i);
   surfMesh.elemIdRegion.Delete(i);
   surfMesh.idRegion.Delete(i);
   //surfMesh.xNormalElem.Delete(i);
   //surfMesh.yNormalElem.Delete(i);
   //surfMesh.zNormalElem.Delete(i);
   surfMesh.numElems--;
   i--; // should go back to verify the next element as well
  }
}

/* esta rotina nao esta otimizada; criamos uma matriz clMatrix para
 * organizar e ordenar os vertices do poliedro resultante da retirada
 * de um ponto da malha de superficie. Para isso configuramos a tal
 * matriz com as linhas representando cada aresta do poliedro sendo que
 * as 2 colunas representam o 1o. vertice e o 2o. vertice da aresta.
 * Esta aresta por sinal nao esta ordenada em algum sentido (horario ou
 * anti-horario) e para isso eh preciso organiza-la e ordenar as linhas
 * da matriz. Isto eh feito buscando o 2 vertice da 1a linha e
 * verificando nas linhas remanescentes da matriz qual delas contem o
 * mesmo vertices. Se for a 1a coluna basta substituir a 2a. linha da
 * matriz pela linha em que o elemento se encontra. Caso esteja na 2a.
 * coluna fazemos o mesmo passo acima porem invertemos tambem o vertice
 * da 1a. coluna com o da 2o coluna. Fazemos este loop ate a ultima
 * linha e suprimimos os elementos repetidos. No final temos uma matriz
 * organizada por arestas e nos.
 * Exemplo:
 *
 * test = [ 0  1 ]   -->   [ 0  1 ]
 *        [ 2  0 ]   -->   [ 1  3 ]
 *        [ 1  3 ]   -->   [ 3  2 ]
 *        [ 3  2 ]   -->   [ 2  0 ]                               
 */
list<int> Model3D::setPolyhedron(list<int> _myList)
{
 list<int> myList = _myList;

  int listSize = myList.size();

  clMatrix test(listSize/2,2);
  list<int>::iterator mele=myList.begin();
  for( int i=0;i<listSize/2;i++ )
  {
   test.Set(i,0,*mele);++mele;
   test.Set(i,1,*mele);++mele;
  }

//--------------------------------------------------
//   cout << "-------- " << endl;
//   test.Display();
//   cout << endl;
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
	 cerr << endl;
	 cerr << endl;
	 cerr << color(blink,red,black);
	 cerr << "              *----------------------------------------*" << endl;
	 cerr << "                       The orientation is wrong!        " << endl;
	 cerr << "              *----------------------------------------*" << endl;
	 cerr << resetColor();
	 cerr << endl;
	 cerr << endl;

	 test.Set(k+1,0,vSwap2);
	 test.Set(k+1,1,vSwap1);
	 test.Set(z,0,v3);
	 test.Set(z,1,v2);
	 break;
	}
   }
  }
//--------------------------------------------------
//   cout << endl;
//   cout << endl;
//   test.Display();
//   cout << "-------- " << endl;
//-------------------------------------------------- 

  /* 
   * old neighbourPoint = [ 0 1 3 1 2 0 3 2 ] 
   * new neighbourPoint = [ 0  1  3  2  0 ]
   */
  myList.clear();
  for( int i=0;i<test.DimI();i++ )
   for( int j=0;j<test.DimJ();j++ )
	myList.push_back(test.Get(i,j));

  // removing all but the first element from every consecutive group of
  // equal elements in the list container
  myList.unique();

  //--------------------------------------------------
  //   cout << "---------" << _v << "------------" << endl;
  //   std::ostream_iterator< int > output( cout, " " );
  //   std::copy( neighbourPoint.at(_v).begin(), 
  //              neighbourPoint.at(_v).end(), output );
  //   cout << endl;
  //-------------------------------------------------- 

  return myList;
}

void Model3D::setNeighbourSurfaceElem()
{
 // list of element neighbours
 neighbourSurfaceElem.resize (0);
 neighbourSurfaceElem.resize (surfMesh.numVerts);
 for( int elem=0;elem<surfMesh.numElems;elem++ )
  for( int vert=0;vert<3;vert++ )
   neighbourSurfaceElem.at( (int)surfMesh.IEN.Get(elem,vert) ).push_back(elem);
}

/*
 *
 * Orientation is always: v1 - v2 - v3
 *
 *
 * */
list<int> Model3D::getNeighbourSurfacePoint(int _node)
{
 int node = _node;
 list<int> myList;

 list<int> plist = neighbourSurfaceElem.at(node);
 for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
 {
  int v1 = (int) surfMesh.IEN.Get(*mele,0);
  int v2 = (int) surfMesh.IEN.Get(*mele,1);
  int v3 = (int) surfMesh.IEN.Get(*mele,2);

  if( v1 == node )
  {
   myList.push_back(v2);
   myList.push_back(v3);
  }
  else if( v2 == node )
  {
   myList.push_back(v3);
   myList.push_back(v1);
  }
  else
  {
   myList.push_back(v1);
   myList.push_back(v2);
  }
 }
 
 // to sort the surface neighbour list (see setPolyhedron)
 myList = setPolyhedron(myList);
 
 return myList;
}

void Model3D::setNeighbourSurfacePoint()
{
 // list of point neighbours
 neighbourPoint.resize (0);
 neighbourPoint.resize (surfMesh.numVerts);

 for( int i=0;i<surfMesh.numVerts;i++ )
  neighbourPoint.at(i) = getNeighbourSurfacePoint(i);
}

 /* Method to flip edge (v1 and v2).
  *
  *                 v1                                  v1
  *                  o                                   o 
  *                 /|\                                 / \
  *                / | \                               /   \
  *               /  |  \            flip             /     \
  *              /   |   \         -------->         /       \
  *             /    |    \                         /    1    \
  *            /     |     \                       /           \
  *           /      |      \                     /             \
  *  v3elem1 o   1   |   2   o v3elem2   v3elem1 o ------------- o v3elem2        *           \      |      /                     \             /              
  *            \     |     /                       \           /      
  *             \    |    /                         \    2    /               
  *              \   |   /                           \       /               
  *               \  |  /                             \     /        
  *                \ | /                               \   /                
  *                 \|/                                 \ /         
  *                  o                                   o                 
  *                 v2                                  v2 
  *
  * */
void Model3D::flipTriangleEdges()
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
  double edge = i;
  int v1 = mapEdgeTri.Get(edge,1); // v1
  int v2 = mapEdgeTri.Get(edge,2); // v2
  int v3elem1 = mapEdgeTri.Get(edge,3); // v3elem1
  int v3elem2 = mapEdgeTri.Get(edge,4); // v3elem2
  int elem1 = mapEdgeTri.Get(edge,5); // elem1
  int elem2 = mapEdgeTri.Get(edge,6); // elem2
  int elemID = surfMesh.elemIdRegion.Get(elem1);

  /* --- Checking the orientation of the mapEdgeTri Matrix ---- */
  clVector normalSurfMesh = getNormalElem(elem1);
  clVector normalToCompare = getNormalElem(v1,v3elem1,v2);
  if( dotProd( normalSurfMesh.Get(0),
	           normalSurfMesh.Get(1),
			   normalSurfMesh.Get(2),
               normalToCompare.Get(0),
			   normalToCompare.Get(1),
			   normalToCompare.Get(2) ) < 0 )
  {
   swap(v1,v2);
  }
  /* -------------------- END of CHECKING --------------------- */

  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);

  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);

  double P3elem1x = surfMesh.X.Get(v3elem1);
  double P3elem1y = surfMesh.Y.Get(v3elem1);
  double P3elem1z = surfMesh.Z.Get(v3elem1);

  double P3elem2x = surfMesh.X.Get(v3elem2);
  double P3elem2y = surfMesh.Y.Get(v3elem2);
  double P3elem2z = surfMesh.Z.Get(v3elem2);

  // elem1
  double area1 = getArea(P1x,P1y,P1z,
	                   P3elem1x,P3elem1y,P3elem1z,
					   P2x,P2y,P2z);
  double q1 = triangleQualityMeasure(P1x,P1y,P1z,
	                               P3elem1x,P3elem1y,P3elem1z,
								   P2x,P2y,P2z);
  double c1 = getCircumRadius(P1x,P1y,P1z,
	                        P3elem1x,P3elem1y,P3elem1z,
							P2x,P2y,P2z);
  clVector normalElem1 = getNormalElem(elem1);

  // elem2
  double area2 = getArea(P2x,P2y,P2z,
	                   P3elem2x,P3elem2y,P3elem2z,
					   P1x,P1y,P1z);
  double q2 = triangleQualityMeasure(P2x,P2y,P2z,
	                               P3elem2x,P3elem2y,P3elem2z,
								   P1x,P1y,P1z);
  double c2 = getCircumRadius(P2x,P2y,P2z,
	                        P3elem2x,P3elem2y,P3elem2z,
							P1x,P1y,P1z);
  clVector normalElem2 = getNormalElem(elem2);

  // elem3 (new elem1)
  double area3 = getArea(P1x,P1y,P1z,
	                   P3elem1x,P3elem1y,P3elem1z,
					   P3elem2x,P3elem2y,P3elem2z);
  double q3 = triangleQualityMeasure(P1x,P1y,P1z,
	                               P3elem1x,P3elem1y,P3elem1z,
	                               P3elem2x,P3elem2y,P3elem2z);
  double c3 = getCircumRadius(P1x,P1y,P1z,
	                        P3elem1x,P3elem1y,P3elem1z,
	                        P3elem2x,P3elem2y,P3elem2z);
  clVector normalElem3 = getNormalElem(v1,v3elem1,v3elem2);

  // elem4 (new elem2)
  double area4 = getArea(P2x,P2y,P2z,
	                   P3elem2x,P3elem2y,P3elem2z,
					   P3elem1x,P3elem1y,P3elem1z);
  double q4 = triangleQualityMeasure(P2x,P2y,P2z,
	                               P3elem2x,P3elem2y,P3elem2z,
	                               P3elem1x,P3elem1y,P3elem1z);
  double c4 = getCircumRadius(P2x,P2y,P2z,
	                        P3elem2x,P3elem2y,P3elem2z,
	                        P3elem1x,P3elem1y,P3elem1z);
  clVector normalElem4 = getNormalElem(v2,v3elem2,v3elem1);

//--------------------------------------------------
//   double angleOld = angle3D(normalElem1.Get(0),normalElem1.Get(1),
// 	                      normalElem1.Get(2),
// 						  normalElem2.Get(0),normalElem2.Get(1),
// 						  normalElem2.Get(2) );
// 
//   double angleNew = angle3D(normalElem3.Get(0),normalElem3.Get(1),
// 	                      normalElem3.Get(2),
// 						  normalElem4.Get(0),normalElem4.Get(1),
// 						  normalElem4.Get(2) );
//-------------------------------------------------- 

  /*                   v1              
   *                 . o            
   *                "  |  " .    y  
   *     v3elem2   "   |     " . 
   *             o"  x |         " .
   *                   |             " .
   *                   o --------------- o  v2
   *                 v3elem1
   *                 
   *                 x: normal to the future plane
   *                 y: vector v1-v2
   *
   *  OBS.: need to check with the angle between x and y allows the
   *  flipping
   *
   *
   * */
//--------------------------------------------------
//   double v12x = P2x-P1x;
//   double v12y = P2y-P1y;
//   double v12z = P2z-P1z;
// 
//   double v13elem1x = P3elem1x-P1x;
//   double v13elem1y = P3elem1y-P1y;
//   double v13elem1z = P3elem1z-P1z;
// 
//   double v13elem2x = P3elem2x-P1x;
//   double v13elem2y = P3elem2y-P1y;
//   double v13elem2z = P3elem2z-P1z;
// 
//   double angle1 = angle3D(v12x,v12y,v12z,
// 	                    v13elem1x,v13elem1y,v13elem1z);
//   double angle2 = angle3D(v12x,v12y,v12z,
// 	                    v13elem2x,v13elem2y,v13elem2z);
// 
//   double v21x = P1x-P2x;
//   double v21y = P1y-P2y;
//   double v21z = P1z-P2z;
// 
//   double v23elem1x = P3elem1x-P2x;
//   double v23elem1y = P3elem1y-P2y;
//   double v23elem1z = P3elem1z-P2z;
// 
//   double v23elem2x = P3elem2x-P2x;
//   double v23elem2y = P3elem2y-P2y;
//   double v23elem2z = P3elem2z-P2z;
// 
//   double angle3 = dotProd(v21x,v21y,v21z,
// 	                    v23elem1x,v23elem1y,v23elem1z);
//   double angle4 = dotProd(v21x,v21y,v21z,
// 	                    v23elem2x,v23elem2y,v23elem2z);
//-------------------------------------------------- 

  // this works, but is not consistent!!! CHANGE IT SOON!
  //double curv1 = fabs(surfMesh.curvature.Get(v1));
  //double curv2 = fabs(surfMesh.curvature.Get(v2));
  //double curv3 = fabs(surfMesh.curvature.Get(v3elem1));
  //double curv4 = fabs(surfMesh.curvature.Get(v3elem2));

  // to avoid contraction at high curvature regions
  //bool curvTest = (curv1 < 20 && curv2 < 20 && curv3 < 20 && curv4 < 20);

  // to avoid 3 elem neighbors vertex
  bool neighTest = (neighbourSurfaceElem.at( v1 ).size() > 4 &&  
                    neighbourSurfaceElem.at( v2 ).size() > 4);   

  /* FLIPPING requirements:
   * - sum of quality of old triangles < sum of quality of new triangles
   * - curvature of all the 4 vertices < 40
   * - sum of area of old triangles > sum of area of new triangles
   * - angle between plane of old triangles > 90 degrees
   * - sum of circumcenters of old triangles > sum of circumcenters of
   *   new triangles
   * */
  if( surfMesh.Marker.Get(v1) == 0.5 && 
	  q1+q2 < q3+q4 &&              // quality sum
	  area1+area2  > area3+area4 && // area sum
	  c1+c2 > c3+c4 &&              // circum radius
	  //angleNew < angleOld &&
	  //curvTest &&
	  neighTest ) 
  {
   cout << "------------- " << color(none,green,black) 
	    << "flipping edge at (" << resetColor()
		<< surfMesh.elemIdRegion.Get(elem1)
		<< color(none,green,black) 
		<< "): " << resetColor() 
		<< v1 << " " << v2 
		<< color(none,green,black) 
		<< " --> " << resetColor()
		<< v3elem1 << " " << v3elem2 << endl;

   // new elem1
   surfMesh.IEN.Set(elem1,0,v1);
   surfMesh.IEN.Set(elem1,1,v3elem1);
   surfMesh.IEN.Set(elem1,2,v3elem2);

   // new elem2
   surfMesh.IEN.Set(elem2,0,v2);
   surfMesh.IEN.Set(elem2,1,v3elem2);
   surfMesh.IEN.Set(elem2,2,v3elem1);

   // update surface, edge matrix, surface neigh elems and points
   restoreMappingArrays();

   saveVTKSurface("./vtk/","surface",opersurf[elemID]);
   flip[elemID]++;
   opersurf[elemID]++;

   // removing low quality elements
   removePointByNeighbourCheck(v1);
   removePointByNeighbourCheck(v2);
   removePointByNeighbourCheck(v3elem1);
   removePointByNeighbourCheck(v3elem2);
  }
 }
}

 /* Method to insert a new point (vAdd) between 2 vertices (v1 and v2).
  *
  *           v3elem1                          v3elem1          
  *              o                                o                   
  *             / \                              /|\
  *            /   \                            / | \
  *           /     \        Add vertex        /  |  \
  *          /       \       --------->       /   |   \
  *         /    1    \                      /    |    \
  *        /           \                    /  1  |  3  \
  *       /             \                  /      | vAdd \
  *   v1 o ------------- o v2          v1 o ----- o ----- o v2               
  *       \             /                  \      |      /              
  *        \           /                    \  2  |  4  /              
  *         \    2    /                      \    |    /                 
  *          \       /                        \   |   /                 
  *           \     /                          \  |  /                 
  *            \   /                            \ | /                   
  *             \ /                              \|/                   
  *              o                                o                     
  *           v3elem2                          v3elem2                    
  *
  * This method also considers the curvature of v1 and v2 to pull up or
  * down (h value) the vAdd.    
  *
  *                                            vAdd
  *    N1                         N2             - x -           <-
  *     ^                         ^            -       -           |
  *      \                       /           -           -         | h
  *       \                     /           /             \        |
  *        \                   /           /               \       |
  *         o ------ x ------ o           o --------------- o    <-
  *        v1       vAdd      v2         v1                 v2
  *       
  * */
void Model3D::insertSurfacePoint(int _edge,const char* _interpolation)
{
 int vAdd = surfMesh.numVerts; // aditional vertice

 cout << "------------- " << color(none,yellow,black) 
      << "inserting vertex at (" << resetColor()
	  << surfMesh.elemIdRegion.Get(mapEdgeTri.Get(_edge,5)) 
	  << color(none,yellow,black) 
	  << "): "
      << resetColor() << vAdd << endl;
 //cout << mapEdgeTri.Get(_edge,0) << endl;
 //saveVTKSurface("./vtk/","insertBefore",vAdd);

 // edge vertices
 int v1 = mapEdgeTri.Get(_edge,1);
 int v2 = mapEdgeTri.Get(_edge,2);
 int v3elem1 = mapEdgeTri.Get(_edge,3);
 int v3elem2 = mapEdgeTri.Get(_edge,4);

 // elements
 int elem1 = mapEdgeTri.Get(_edge,5);
 int elem2 = mapEdgeTri.Get(_edge,6);

 /* --- Checking the orientation of the mapEdgeTri Matrix ---- */
 clVector normalSurfMesh = getNormalElem(elem1);
 clVector normalToCompare = getNormalElem(v1,v2,v3elem1);
 if( dotProd( normalSurfMesh.Get(0),
	           normalSurfMesh.Get(1),
			   normalSurfMesh.Get(2),
               normalToCompare.Get(0),
			   normalToCompare.Get(1),
			   normalToCompare.Get(2) ) < 0 )
  {
   swap(v1,v2);
  }
  /* -------------------- END of CHECKING --------------------- */

 double XvAdd = 0.0;
 double YvAdd = 0.0;
 double ZvAdd = 0.0;

 if( strcmp( _interpolation,"bi-curvature") == 0 ) 
 {
  clVector coordAdd1 = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                              X.Get(v2),Y.Get(v2),Z.Get(v2),
								  surfMesh.xNormal.Get(v1)+
								  surfMesh.xNormal.Get(v3elem1),
								  surfMesh.yNormal.Get(v1)+
								  surfMesh.yNormal.Get(v3elem1),
								  surfMesh.zNormal.Get(v1)+
								  surfMesh.zNormal.Get(v3elem1),
								  surfMesh.xNormal.Get(v2)+
								  surfMesh.xNormal.Get(v3elem2),
								  surfMesh.yNormal.Get(v2)+
								  surfMesh.yNormal.Get(v3elem2),
								  surfMesh.zNormal.Get(v2)+
								  surfMesh.zNormal.Get(v3elem2) );

  XvAdd = coordAdd1.Get(0);
  YvAdd = coordAdd1.Get(1);
  ZvAdd = coordAdd1.Get(2);
 }
 else if( strcmp( _interpolation,"curvature") == 0 ) 
 {
  clVector coordAdd = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                             X.Get(v2),Y.Get(v2),Z.Get(v2),
								 surfMesh.xNormal.Get(v1),
								 surfMesh.yNormal.Get(v1),
								 surfMesh.zNormal.Get(v1),
								 surfMesh.xNormal.Get(v2),
								 surfMesh.yNormal.Get(v2),
								 surfMesh.zNormal.Get(v2) );

  XvAdd = coordAdd.Get(0);
  YvAdd = coordAdd.Get(1);
  ZvAdd = coordAdd.Get(2);
 }
 else // flat
 {
  // add point in the middle of a edge (not consider curvature)
  XvAdd = ( surfMesh.X.Get(v1)+ surfMesh.X.Get(v2) )*0.5;
  YvAdd = ( surfMesh.Y.Get(v1)+ surfMesh.Y.Get(v2) )*0.5;
  ZvAdd = ( surfMesh.Z.Get(v1)+ surfMesh.Z.Get(v2) )*0.5;
 }


//--------------------------------------------------
//   cout << "Flat: " << endl;
//   cout << "x: " << ( surfMesh.X.Get(v1)+ surfMesh.X.Get(v2) )*0.5 << endl;
//   cout << "y: " << ( surfMesh.Y.Get(v1)+ surfMesh.Y.Get(v2) )*0.5 << endl;
//   cout << "z: " << ( surfMesh.Z.Get(v1)+ surfMesh.Z.Get(v2) )*0.5 << endl;
//   cout << endl;
//   cout << "curvature: " << endl;
//   cout << "x: " << XvAdd << endl;
//   cout << "y: " << YvAdd << endl;
//   cout << "z: " << ZvAdd << endl;
//   cout << " ----------------- " << endl;
//-------------------------------------------------- 

//--------------------------------------------------
//  cout << "v1: " << v1 << " " << "v2: " << v2 << endl;
//  cout << "1: " << Xcurv1 << " " << Ycurv1 << " " << Zcurv1 << endl;
//  cout << "2: " << Xcurv2 << " " << Ycurv2 << " " << Zcurv2 << endl;
//  cout << "add: " << XvAdd << " " << YvAdd << " " << ZvAdd << endl;
//  cout << "radius: " << r2D << endl;
//  cout << "xMid2D-Xc: " << xMid2D-Xc << endl;
//  cout << "center: " << Xc << " " << Yc << endl;
//  cout << "midNew: " << yMidNew1 << " " << yMidNew2 << endl;
//  cout << "curvatures v1 - v2: " << surfMesh.curvature.Get(v1) << " " 
//       << surfMesh.curvature.Get(v2) << endl;
//  cout << "normal: " << normalXUnit << " " 
//       << normalYUnit << " " << normalZUnit << endl;
//-------------------------------------------------- 

 // insert aditional vertice coordinate
 X.AddItem(vAdd,XvAdd);
 Y.AddItem(vAdd,YvAdd);
 Z.AddItem(vAdd,ZvAdd);
 heaviside.AddItem(vAdd,heaviside.Get(v1));
 curvature.AddItem(vAdd,curvature.Get(v1));
 edgeSize.AddItem(vAdd,edgeSize.Get(v1));

 surfMesh.X.AddItem(XvAdd);
 surfMesh.Y.AddItem(YvAdd);
 surfMesh.Z.AddItem(ZvAdd);
 surfMesh.Marker.AddItem(surfMesh.Marker.Get(v1)); 
 surfMesh.vertIdRegion.AddItem(surfMesh.vertIdRegion.Get(v1));
 surfMesh.phyBounds.push_back(surfMesh.phyBounds.at(v1));

 // incremeting the number of points
 surfMesh.numVerts++;
 numVerts++;


 /* by adding 1 point on the edge it is necessary to divide the
  * original element and also the oposite element by 2, becoming 4
  * elements in total. */
 /* 
  *           v3elem1                          v3elem1          
  *              o                                o                   
  *             / \                              /|\
  *            /   \                            / | \
  *           /     \        Add vertex        /  |  \
  *          /       \       --------->       /   |   \
  *         /    1    \                      /    |    \
  *        /           \                    /  1  |  3  \
  *       /             \                  /      | vAdd \
  *   v1 o ------------- o v2          v1 o ----- o ----- o v2               
  *       \             /                  \      |      /              
  *        \           /                    \  2  |  4  /              
  *         \    2    /                      \    |    /                 
  *          \       /                        \   |   /                 
  *           \     /                          \  |  /                 
  *            \   /                            \ | /                   
  *             \ /                              \|/                   
  *              o                                o                     
  *           v3elem2                          v3elem2                    
  *
  * */
 
 // 1st. new element (v1 - vAdd - v3elem1) 
 // on the same position of the OLD 1st. element (v1 - v2 - v3elem1)
 // OLD ELEM1 //
 surfMesh.IEN.Set(elem1,0,v1);
 surfMesh.IEN.Set(elem1,1,vAdd);
 surfMesh.IEN.Set(elem1,2,v3elem1);
 // add new elemIdRegion
 surfMesh.elemIdRegion.Set(elem1,surfMesh.elemIdRegion.Get(elem1));
 surfMesh.idRegion.Set(elem1,surfMesh.idRegion.Get(elem1));

 // 2nd. new element (v1 - vAdd - v3elem2) 
 // on the same position of the OLD 2nd. element (v1 - v2 - v3elem2)
 // OLD ELEM2 //
 surfMesh.IEN.Set(elem2,0,v1);
 surfMesh.IEN.Set(elem2,1,v3elem2);
 surfMesh.IEN.Set(elem2,2,vAdd);
 // add new elemIdRegion
 surfMesh.elemIdRegion.Set(elem2,surfMesh.elemIdRegion.Get(elem1));
 surfMesh.idRegion.Set(elem2,surfMesh.idRegion.Get(elem1));

 // 3rd. new element (v2 - vAdd - v3elem1) on the last row
 // OLD ELEM1 //
 surfMesh.IEN.AddRow();
 int elem3 = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem3,0,v2);
 surfMesh.IEN.Set(elem3,1,v3elem1);
 surfMesh.IEN.Set(elem3,2,vAdd);
 // add new elemIdRegion
 surfMesh.elemIdRegion.AddItem(surfMesh.elemIdRegion.Get(elem1));
 surfMesh.idRegion.AddItem(surfMesh.idRegion.Get(elem1));
 surfMesh.numElems++;

 // 4th. new element (v2 - vAdd - v3elem2) on the last row
 // OLD ELEM2 //
 surfMesh.IEN.AddRow();
 int elem4 = surfMesh.IEN.DimI()-1;
 surfMesh.IEN.Set(elem4,0,v2);
 surfMesh.IEN.Set(elem4,1,vAdd);
 surfMesh.IEN.Set(elem4,2,v3elem2);
 // add new elemIdRegion
 surfMesh.elemIdRegion.AddItem(surfMesh.elemIdRegion.Get(elem1));
 surfMesh.idRegion.AddItem(surfMesh.idRegion.Get(elem1));
 surfMesh.numElems++;
 
 // update surface, edge matrix, surface neigh elems and points
 restoreMappingArrays();

 // curvature is approx. the average between vertices
//--------------------------------------------------
//double curv = (surfMesh.curvature.Get(v1)+surfMesh.curvature.Get(v2))*0.5;
//  surfMesh.curvature.AddItem(curv);
//  curvature.AddItem(vAdd,curv);
//-------------------------------------------------- 

 // computing curvature
 clVector myVec = getNormalAndKappa(vAdd,
                    getNeighbourSurfacePoint(vAdd));
 curvature.AddItem(vAdd,myVec.Get(0));
 surfMesh.curvature.AddItem(vAdd,myVec.Get(0));
 surfMesh.xNormal.AddItem(vAdd,myVec.Get(1));
 surfMesh.yNormal.AddItem(vAdd,myVec.Get(2));
 surfMesh.zNormal.AddItem(vAdd,myVec.Get(3));

//--------------------------------------------------
// cout << "curv(v1):      " << surfMesh.curvature.Get(v1) << endl;
// cout << "curv(v2):      " << surfMesh.curvature.Get(v2) << endl;
// cout << "curv(v3elem1): " << surfMesh.curvature.Get(v3elem1) << endl;
// cout << "curv(v3elem2): " << surfMesh.curvature.Get(v3elem2) << endl;
// cout << "new curv:        " << curv << endl;
// cout << "calculated curv: " << myVec.Get(0) << endl;
//-------------------------------------------------- 
}

void Model3D::removeSurfacePoint(int _node)
{
 //saveVTKSurface("./vtk/","deleteBefore",v1);

 // delete v1 from surface, xSurface, ySurface, zSurface vectors
 // surface is not used to add/remove/flip elements before the remeshing
 // but it's used on removePointsByInterfaceDistance
 // to be implemented 

 // marking the desired elements for deletion
 list<int> plist = neighbourSurfaceElem.at(_node);
 for( list<int>::iterator mele=plist.begin(); mele != plist.end();++mele )
  markSurfElemForDeletion(*mele);

 // deleting elements
 deleteSurfaceElements();

 // after the deletion process it's mandatory to create new elements
 // to fill the space left by the deleting process
 //surfaceTriangulator(_node);
 surfaceTriangulatorEarClipping(_node,
                                neighbourPoint.at(_node),
								"no");

 // deleting X,Y and Z coordinate; deleting the point maker funcition
 deleteSurfacePoint(_node);

 // update surface, edge matrix, surface neigh elems and points
 restoreMappingArrays();
}


/* Method to contract edge (v1 and v2).
 *
 *                v3elem1                          v3elem1          
 *                   o                                o                   
 *                  / \                               |
 *                 /   \                              |
 *                /     \                             |
 *       o       /       \       o                o   |   o
 *        \     /    1    \     /                  \  |  /
 *         \   /           \   /                    \ | /
 *          \ /             \ /                      \| 
 *        v1 o ------------- o v2                     o v1
 *          / \             / \                      /|\
 *         /   \           /   \                    / | \
 *        /     \    2    /     \                  /  |  \
 *       o       \       /       o                o   |   o
 *                \     /                             |
 *                 \   /                              |
 *                  \ /                               |                   
 *                   o                                o                     
 *                v3elem2                          v3elem2                    
 *    
 * */
void Model3D::contractEdgesByLength(const char* _interpolation,
                                    double _param)
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int edge=0;edge<mapEdgeTri.DimI();edge++ )
 {
  // int length = mapEdgeTri.Get(edge,0); // length
  int v1 = mapEdgeTri.Get(edge,1); 
  int v2 = mapEdgeTri.Get(edge,2); 
  int v3elem1 = mapEdgeTri.Get(edge,3);
  int v3elem2 = mapEdgeTri.Get(edge,4);

  double curv1 = fabs(surfMesh.curvature.Get(v1));
  double curv2 = fabs(surfMesh.curvature.Get(v2));
  double curv3 = fabs(surfMesh.curvature.Get(v3elem1));
  double curv4 = fabs(surfMesh.curvature.Get(v3elem2));
  int elem1 = mapEdgeTri.Get(edge,5);
  int elem2 = mapEdgeTri.Get(edge,6);
  double edgeLength = mapEdgeTri.Get(edge,0);
  clVector normalElem1 = getNormalElem(elem1);
  clVector normalElem2 = getNormalElem(elem2);

//--------------------------------------------------
//   // radii
//   double angle = angle3D(normalElem1.Get(0),normalElem1.Get(1),
// 	                   normalElem1.Get(2),
// 					   normalElem2.Get(0),normalElem2.Get(1),
// 					   normalElem2.Get(2) );
//-------------------------------------------------- 

  // verifying the length of each surface edge
  int elemID = surfMesh.elemIdRegion.Get(mapEdgeTri.Get(edge,5));

  //double minCurv = min(curv1,curv2);
  //double erro = minCurv*edgeLength;
  
  // angle test
  // bool angleTest = angle > 0.0;

  // elemID out of boudary
  bool elemIDTest = elemID > 0;

  // to avoid contraction at high curvature regions
  bool curvTest = (curv1 < 40 && curv2 < 40 && curv3 < 40 && curv4 < 40);

  // to avoid 3 elem neighbors vertex
  bool neighTest = (neighbourSurfaceElem.at( v3elem1 ).size() > 4 &&  
                    neighbourSurfaceElem.at( v3elem2 ).size() > 4);   

  //if( elemID > 0 && erro < 0.5*erroS )//&&
  if( edgeLength < _param*triEdge[elemID] &&
      //erro < 0.03 &&
	  elemIDTest && 
	  curvTest && 
	  neighTest //&& 
	  //angleTest 
	)
  {
//--------------------------------------------------
//    cout << " ----------------- " << endl;
//    cout << "v1: " << v1 << endl;
//    cout << "v2: " << v2 << endl;
//    cout << "v3elem1: " << v3elem1 << endl;
//    cout << "v3elem2: " << v3elem2 << endl;
//    cout << "elem1: " << elem1 << endl;
//    cout << "elem2: " << elem2 << endl;
//    cout << " ----------------- " << endl;
//-------------------------------------------------- 

   if( strcmp( _interpolation,"curvature") == 0 ) 
   {
	// using curvature
	clVector coordAdd = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                               X.Get(v2),Y.Get(v2),Z.Get(v2),
								   surfMesh.xNormal.Get(v1),
								   surfMesh.yNormal.Get(v1),
								   surfMesh.zNormal.Get(v1),
								   surfMesh.xNormal.Get(v2),
								   surfMesh.yNormal.Get(v2),
								   surfMesh.zNormal.Get(v2) );

	double XvAdd = coordAdd.Get(0);
	double YvAdd = coordAdd.Get(1);
	double ZvAdd = coordAdd.Get(2);
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else if( strcmp( _interpolation,"bi-curvature") == 0 ) 
   {
	// using bi-curvature
	clVector coordAdd1 = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                                X.Get(v2),Y.Get(v2),Z.Get(v2),
									surfMesh.xNormal.Get(v1),
									surfMesh.yNormal.Get(v1),
									surfMesh.zNormal.Get(v1),
									surfMesh.xNormal.Get(v2),
									surfMesh.yNormal.Get(v2),
									surfMesh.zNormal.Get(v2) );

	clVector coordAdd2 = fitEllipse( X.Get(v3elem1),
	                                Y.Get(v3elem1),
									Z.Get(v3elem1),
									X.Get(v3elem2),
									Y.Get(v3elem2),
									Z.Get(v3elem2),
									surfMesh.xNormal.Get(v3elem1),
									surfMesh.yNormal.Get(v3elem1),
									surfMesh.zNormal.Get(v3elem1),
									surfMesh.xNormal.Get(v3elem2),
									surfMesh.yNormal.Get(v3elem2),
									surfMesh.zNormal.Get(v3elem2) );

	double XvAdd = (coordAdd1.Get(0)+coordAdd2.Get(0))*0.5;
	double YvAdd = (coordAdd1.Get(1)+coordAdd2.Get(1))*0.5;
	double ZvAdd = (coordAdd1.Get(2)+coordAdd2.Get(2))*0.5;
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else // flat
   {
	//--------------------------------------------------
	// double XvNew = ( surfMesh.X.Get(v1)+ surfMesh.X.Get(v2) )*0.5;
	// double YvNew = ( surfMesh.Y.Get(v1)+ surfMesh.Y.Get(v2) )*0.5;
	// double ZvNew = ( surfMesh.Z.Get(v1)+ surfMesh.Z.Get(v2) )*0.5;
	//-------------------------------------------------- 
	double XvNew = surfMesh.X.Get(v1);
	double YvNew = surfMesh.Y.Get(v1);
	double ZvNew = surfMesh.Z.Get(v1);
	surfMesh.X.Set(v1, XvNew );
	surfMesh.Y.Set(v1, YvNew );
	surfMesh.Z.Set(v1, ZvNew );
	X.Set(v1, XvNew );
	Y.Set(v1, YvNew );
	Z.Set(v1, ZvNew );
   }

   // changing surfMesh.IEN from v2 to v1
   for( int i=0;i<surfMesh.IEN.DimI();i++ )
	for( int j=0;j<surfMesh.IEN.DimJ();j++ )
	 if( surfMesh.IEN.Get(i,j)==v2 )
	  surfMesh.IEN.Set(i,j,v1);

   // deleting element 1 and 2 (see comments above)
   markSurfElemForDeletion(elem1);
   markSurfElemForDeletion(elem2);
   deleteSurfaceElements();

   // delete v2 because v1 is always lower then v2
   deleteSurfacePoint(v2);

   // update surface, edge matrix, surface neigh elems and points
   restoreMappingArrays();

   // computing curvature
   clVector myVec = getNormalAndKappa(v1,getNeighbourSurfacePoint(v1));
   curvature.Set(v1,myVec.Get(0));
   surfMesh.curvature.Set(v1,myVec.Get(0));
   surfMesh.xNormal.Set(v1,myVec.Get(1));
   surfMesh.yNormal.Set(v1,myVec.Get(2));
   surfMesh.zNormal.Set(v1,myVec.Get(3));

   // removing low quality elements
   if( v3elem1 > v2 )
	v3elem1 = v3elem1-1;
   if( v3elem2 > v2 )
	v3elem2 = v3elem2-1;
   removePointByNeighbourCheck(v1);
   removePointByNeighbourCheck(v3elem1);
   removePointByNeighbourCheck(v3elem2);

   cout << "------------- " << color(none,blue,black) 
	    << "contracting edge at (" << resetColor()
		<< surfMesh.elemIdRegion.Get(elem1)
		<< color(none,blue,black) 
		<< "): " << resetColor() 
		<< v2 << color(none,blue,black) 
		<< " --> " << resetColor()
		<< v1 << endl;
   saveVTKSurface("./vtk/","surface",opersurf[elemID]);
   csp[elemID]++;
   opersurf[elemID]++;

//--------------------------------------------------
//    cout << "curv1: " << curv1  
//         << "  curv2: " << curv2 
//         << "  new: " <<  myVec.Get(0) << endl;
//-------------------------------------------------- 
  }
 }
}

void Model3D::contractEdgesByLength2(const char* _interpolation,
                                    double _param)
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int edge=0;edge<mapEdgeTri.DimI();edge++ )
 {
  // int length = mapEdgeTri.Get(edge,0); // length
  int v1 = mapEdgeTri.Get(edge,1); 
  int v2 = mapEdgeTri.Get(edge,2); 
  int v3elem1 = mapEdgeTri.Get(edge,3);
  int v3elem2 = mapEdgeTri.Get(edge,4);

  int elem1 = mapEdgeTri.Get(edge,5);
  int elem2 = mapEdgeTri.Get(edge,6);
  double edgeLength = mapEdgeTri.Get(edge,0);
  clVector normalElem1 = getNormalElem(elem1);
  clVector normalElem2 = getNormalElem(elem2);

//--------------------------------------------------
//   // radii
//   double angle = angle3D(normalElem1.Get(0),normalElem1.Get(1),
// 	                   normalElem1.Get(2),
// 					   normalElem2.Get(0),normalElem2.Get(1),
// 					   normalElem2.Get(2) );
//-------------------------------------------------- 

  // verifying the length of each surface edge
  int elemID = surfMesh.elemIdRegion.Get(mapEdgeTri.Get(edge,5));

  //double minCurv = min(curv1,curv2);
  //double erro = minCurv*edgeLength;
  
  // angle test
  // bool angleTest = angle > 0.0;

  // elemID out of boudary
  bool elemIDTest = elemID > 0;

  // to avoid 3 elem neighbors vertex
  bool neighTest = (neighbourSurfaceElem.at( v3elem1 ).size() > 4 &&  
                    neighbourSurfaceElem.at( v3elem2 ).size() > 4);   

  //if( elemID > 0 && erro < 0.5*erroS )//&&
  if( edgeLength < _param*triEdge[elemID] &&
      //erro < 0.03 &&
	  elemIDTest && 
	  neighTest //&& 
	  //angleTest 
	)
  {
//--------------------------------------------------
//    cout << " ----------------- " << endl;
//    cout << "v1: " << v1 << endl;
//    cout << "v2: " << v2 << endl;
//    cout << "v3elem1: " << v3elem1 << endl;
//    cout << "v3elem2: " << v3elem2 << endl;
//    cout << "elem1: " << elem1 << endl;
//    cout << "elem2: " << elem2 << endl;
//    cout << " ----------------- " << endl;
//-------------------------------------------------- 

   if( strcmp( _interpolation,"curvature") == 0 ) 
   {
	// using curvature
	clVector coordAdd = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                               X.Get(v2),Y.Get(v2),Z.Get(v2),
								   surfMesh.xNormal.Get(v1),
								   surfMesh.yNormal.Get(v1),
								   surfMesh.zNormal.Get(v1),
								   surfMesh.xNormal.Get(v2),
								   surfMesh.yNormal.Get(v2),
								   surfMesh.zNormal.Get(v2) );

	double XvAdd = coordAdd.Get(0);
	double YvAdd = coordAdd.Get(1);
	double ZvAdd = coordAdd.Get(2);
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else if( strcmp( _interpolation,"bi-curvature") == 0 ) 
   {
	// using bi-curvature
	clVector coordAdd1 = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                                X.Get(v2),Y.Get(v2),Z.Get(v2),
									surfMesh.xNormal.Get(v1),
									surfMesh.yNormal.Get(v1),
									surfMesh.zNormal.Get(v1),
									surfMesh.xNormal.Get(v2),
									surfMesh.yNormal.Get(v2),
									surfMesh.zNormal.Get(v2) );

	clVector coordAdd2 = fitEllipse( X.Get(v3elem1),
	                                Y.Get(v3elem1),
									Z.Get(v3elem1),
									X.Get(v3elem2),
									Y.Get(v3elem2),
									Z.Get(v3elem2),
									surfMesh.xNormal.Get(v3elem1),
									surfMesh.yNormal.Get(v3elem1),
									surfMesh.zNormal.Get(v3elem1),
									surfMesh.xNormal.Get(v3elem2),
									surfMesh.yNormal.Get(v3elem2),
									surfMesh.zNormal.Get(v3elem2) );

	double XvAdd = (coordAdd1.Get(0)+coordAdd2.Get(0))*0.5;
	double YvAdd = (coordAdd1.Get(1)+coordAdd2.Get(1))*0.5;
	double ZvAdd = (coordAdd1.Get(2)+coordAdd2.Get(2))*0.5;
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else // flat
   {
	//--------------------------------------------------
	// double XvNew = ( surfMesh.X.Get(v1)+ surfMesh.X.Get(v2) )*0.5;
	// double YvNew = ( surfMesh.Y.Get(v1)+ surfMesh.Y.Get(v2) )*0.5;
	// double ZvNew = ( surfMesh.Z.Get(v1)+ surfMesh.Z.Get(v2) )*0.5;
	//-------------------------------------------------- 
	double XvNew = surfMesh.X.Get(v1);
	double YvNew = surfMesh.Y.Get(v1);
	double ZvNew = surfMesh.Z.Get(v1);
	surfMesh.X.Set(v1, XvNew );
	surfMesh.Y.Set(v1, YvNew );
	surfMesh.Z.Set(v1, ZvNew );
	X.Set(v1, XvNew );
	Y.Set(v1, YvNew );
	Z.Set(v1, ZvNew );
   }

   // changing surfMesh.IEN from v2 to v1
   for( int i=0;i<surfMesh.IEN.DimI();i++ )
	for( int j=0;j<surfMesh.IEN.DimJ();j++ )
	 if( surfMesh.IEN.Get(i,j)==v2 )
	  surfMesh.IEN.Set(i,j,v1);

   // deleting element 1 and 2 (see comments above)
   markSurfElemForDeletion(elem1);
   markSurfElemForDeletion(elem2);
   deleteSurfaceElements();

   // delete v2 because v1 is always lower then v2
   deleteSurfacePoint(v2);

   // update surface, edge matrix, surface neigh elems and points
   restoreMappingArrays();

   // computing curvature
   clVector myVec = getNormalAndKappa(v1,getNeighbourSurfacePoint(v1));
   curvature.Set(v1,myVec.Get(0));
   surfMesh.curvature.Set(v1,myVec.Get(0));
   surfMesh.xNormal.Set(v1,myVec.Get(1));
   surfMesh.yNormal.Set(v1,myVec.Get(2));
   surfMesh.zNormal.Set(v1,myVec.Get(3));

   // removing low quality elements
   if( v3elem1 > v2 )
	v3elem1 = v3elem1-1;
   if( v3elem2 > v2 )
	v3elem2 = v3elem2-1;
   removePointByNeighbourCheck(v1);
   removePointByNeighbourCheck(v3elem1);
   removePointByNeighbourCheck(v3elem2);

   cout << "------------- " << color(none,blue,black) 
	    << "contracting edge at (" << resetColor()
		<< surfMesh.elemIdRegion.Get(elem1)
		<< color(none,blue,black) 
		<< "): " << resetColor() 
		<< v2 << color(none,blue,black) 
		<< " --> " << resetColor()
		<< v1 << endl;
   saveVTKSurface("./vtk/","surface",opersurf[elemID]);
   csp[elemID]++;
   opersurf[elemID]++;

//--------------------------------------------------
//    cout << "curv1: " << curv1  
//         << "  curv2: " << curv2 
//         << "  new: " <<  myVec.Get(0) << endl;
//-------------------------------------------------- 
  }
 }
}

/* If only the _interpolation input is given, use 0.6 as contraction
 * condition
 * */
void Model3D::contractEdgesByLength(const char* _interpolation)
{
 contractEdgesByLength(_interpolation,0.6);
}

void Model3D::removePointsByLength(double _param)
{
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
 // edge vertices
 double edgeLength = mapEdgeTri.Get(i,0);
 int v1 = mapEdgeTri.Get(i,1);
 int v2 = mapEdgeTri.Get(i,2);
 //int v3elem1 = mapEdgeTri.Get(i,3);
 //int v3elem2 = mapEdgeTri.Get(i,4);
 int elemID = surfMesh.elemIdRegion.Get(mapEdgeTri.Get(i,5));
 int vertID = surfMesh.vertIdRegion.Get(v1);
 //double curv1 = fabs(surfMesh.curvature.Get(v1));
 //double curv2 = fabs(surfMesh.curvature.Get(v2));
 //double curv3 = fabs(surfMesh.curvature.Get(v3elem1));
 //double curv4 = fabs(surfMesh.curvature.Get(v3elem2));
 
 //double minCurv = min(curv1,curv2);
 //double erro = minCurv*edgeLength;
 
 // angle test
 // bool angleTest = angle > 0.0;
 
 // elemID out of boudary
 bool elemIDTest = elemID > 0;
 
 // to avoid contraction at high curvature regions
 //bool curvTest = (curv1 < 40 && curv2 < 40 && curv3 < 40 && curv4 < 40);
 
 
 if( edgeLength < _param*triEdge[elemID] &&
     //erro < 0.03 &&
     elemIDTest //&& 
     //curvTest 
     //angleTest 
   )
  {
    // sum of all neighbour edge length of the 1st. point
    double sumLength1=0;
    int listSize1 = neighbourPoint.at(v1).size();
    list<int> plist1 = neighbourPoint.at(v1);
    list<int>::iterator vert1=plist1.begin();
    for( int i=0;i<listSize1-1;i++ )
    {
     double P0x = surfMesh.X.Get(v1);
     double P0y = surfMesh.Y.Get(v1);
     double P0z = surfMesh.Z.Get(v1);
     
     int v = *vert1;++vert1;
     double P1x = surfMesh.X.Get(v);
     double P1y = surfMesh.Y.Get(v);
     double P1z = surfMesh.Z.Get(v);
     
     sumLength1 += distance(P0x,P0y,P0z,P1x,P1y,P1z);
    }

   // sum of all neighbour edge length of the 1st. point
   double sumLength2=0;
   int listSize2 = neighbourPoint.at(v2).size();
   list<int> plist2 = neighbourPoint.at(v2);
   list<int>::iterator vert2=plist2.begin();
   for( int i=0;i<listSize2-1;i++ )
   {
    double P0x = surfMesh.X.Get(v2);
    double P0y = surfMesh.Y.Get(v2);
    double P0z = surfMesh.Z.Get(v2);

    int v = *vert2;++vert2;
    double P2x = surfMesh.X.Get(v);
    double P2y = surfMesh.Y.Get(v);
    double P2z = surfMesh.Z.Get(v);

    sumLength2 += distance(P0x,P0y,P0z,P2x,P2y,P2z);
   }

   // check which node has the smallest length sum and proceed
   if( sumLength1 < sumLength2 )
   {
    cout << "------------- " << color(none,red,black) 
         << "removing vertex at (" 
         << resetColor()
         << surfMesh.vertIdRegion.Get(v1)
         << color(none,red,black) 
         << "): "
      << resetColor() << v1 << endl;

    removeSurfacePoint(v1);

    saveVTKSurface("./vtk/","surface",opersurf[vertID]);
    rsp[vertID]++;
    opersurf[vertID]++;
   }
   else // if the 2nd. node has the smallest edge length sum
   {
    cout << "------------ " << color(none,red,black) 
         << "removing vertex at(" 
         << resetColor()
         << surfMesh.vertIdRegion.Get(v2)
         << color(none,red,black) 
         << "): "
      << resetColor() << v2 << endl;

    removeSurfacePoint(v2);

    saveVTKSurface("./vtk/","surface",opersurf[vertID]);
    rsp[vertID]++;
    opersurf[vertID]++;
   }
  }
 }
}

void Model3D::removePointsByLength()
{
 removePointsByLength(0.6);
}

/*
 * Remove point according to a fixed distance between the interface and
 * the point itself. The fixed distance is the variable "d" and d>0
 * means that it cannot remove points lied on the surface. 
 *
 * OBS.: this method should be used carefuly due to an inconsitency with
 * the TETGEN mesh generator. If this method removes too many points,
 * the TETGEN creates tetrahedral elements with all the 4 nodes in the
 * interface, thus the element ratio quality is low and its volume is
 * too small.
 *
 * */
void Model3D::removePointsByInterfaceDistance()
{
 /*     
  *                l*sqrt(6)
  *   height = h = ---------- = l*0.8164
  *                    3
  * */

 double triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));

 double h = triEdgeMin*0.8864; 
 for( int i=0;i<numVerts;i++ )
 {
  int vertID = vertIdRegion.Get(i);
  double d = interfaceDistance.Get(i);

  if( d>0 && d<h*1.0 ) // hiRe
  {
//--------------------------------------------------
//    cout << "--- " << color(none,red,black) << "removing vertex by distance: "
// 	    << resetColor() << i << endl;
//-------------------------------------------------- 

   mark3DPointForDeletion(i);
   rpi[vertID]++;
  }
 }
 //cout << "  removed by interface Distance: " << rpi[vertID] << endl;
}

/* 
 * Remove point according to its distance compared to its neighbours.
 * The idea of such a method was to avoid cluster of points in one
 * specific area. But with the implementation of the
 * remove3dMeshPointsByDiffusion, this method is obsolete.
 *
 * */
void Model3D::remove3dMeshPointsByDistance()
{
 for( int i=surfMesh.numVerts;i<numVerts;i++ )
 {
  int vertID = vertIdRegion.Get(i);

  for( int j=surfMesh.numVerts;j<numVerts;j++ )
  {
   if( heaviside.Get(i) != 0.5 && heaviside.Get(j) != 0.5 )
   {

	double d = distance( X.Get(i),Y.Get(i),Z.Get(i),
	                   X.Get(j),Y.Get(j),Z.Get(j) );
	//--------------------------------------------------
	// if( interfaceDistance.Get(i) > 3.0 &&
	//     interfaceDistance.Get(j) > 3.0 &&
	// 	d>0 && d<3.0*triEdge[2] )
	//-------------------------------------------------- 
	//if( d>0 && d<1.0*triEdge[1] )
	if( d>0 && d<2.0*triEdge[1] )
	{
	//--------------------------------------------------
	//  cout << "- " << color(none,blue,black) 
	//       << "removing dense vertex cluster: "
	//       << resetColor() << i << " " << heaviside.Get(i) << endl;
	//-------------------------------------------------- 
	 mark3DPointForDeletion(i);
	 rpdist[vertID]++;
	}
   }
  }
//--------------------------------------------------
//   cout << "  removed by distance: " << rpdist[vertID] << endl;
//-------------------------------------------------- 
 }
}

/*
 * Insert point(s) according to the solution of the diffusion equation
 * given by the class Helmholtz3D.
 *
 * */
void Model3D::insert3dMeshPointsByDiffusion()
{
 /*
  * mapEdge.Set(edge,0,numVerts+edge); // numero da aresta
  * mapEdge.Set(edge,1,xMid ); // coordenada X do centro da aresta
  * mapEdge.Set(edge,2,yMid ); // coordenada Y do centro da aresta
  * mapEdge.Set(edge,3,zMid ); // coordenada Y do centro da aresta
  * mapEdge.Set(edge,4,faces[i].p1 ); // 1o noh
  * mapEdge.Set(edge,5,faces[i].p2 ); // 2o noh
  * */
 for( int e=0;e<mapEdge.DimI();e++ )
 {
  double XvAdd = mapEdge.Get(e,1);
  double YvAdd = mapEdge.Get(e,2);
  double ZvAdd = mapEdge.Get(e,3);

  int v1 = mapEdge.Get(e,4);
  int v2 = mapEdge.Get(e,5);

  int vertID = vertIdRegion.Get(v1);

  double x1=X.Get(v1);
  double y1=Y.Get(v1);
  double z1=Z.Get(v1);
  double x2=X.Get(v2);
  double y2=Y.Get(v2);
  double z2=Z.Get(v2);
  double length = distance(x1,y1,z1,x2,y2,z2);

  // minVert should be bigger then surfMesh.numVerts because we are
  // treating only 3D vertices after the surface mesh vertices.
  //int minVert = min(v1,v2);
  int maxVert = max(v1,v2);
  double maxEdge = max(edgeSize.Get(v1),edgeSize.Get(v2));
  //double hSum = heaviside.Get(v1) + heaviside.Get(v2);

  // edgeSize is the result of \nabla^2 edge = 0
  if( length > 2.5*maxEdge && 
	//--------------------------------------------------
	//   interfaceDistance.Get(v1) > 2*triEdge[1] &&
	//   interfaceDistance.Get(v2) > 2*triEdge[1] &&
	//-------------------------------------------------- 
	  //ipd[vertID] < 200 &&
	  maxVert > surfMesh.numVerts )
  {
   cout << v1 << " (" << edgeSize.Get(v1) << ") " 
	    << v2 << " (" << edgeSize.Get(v2) << ") " << endl;
   cout << x1 << " " << y1 << " " << z1 << endl;
   cout << x2 << " " << y2 << " " << z2 << endl;
   cout << X.Get(v2) << " " << Y.Get(v2) << " " << Z.Get(v2) << endl;
   cout << e << " " << length << " " << edgeSize.Get(v1) << endl;
   int vAdd = numVerts; // aditional vertice

   cout << "- " << color(none,blue,black) 
	            << "inserting vertex: "
				<< resetColor() << vAdd << " " 
				<< heaviside.Get(maxVert) << endl;

   X.AddItem(vAdd,XvAdd);
   Y.AddItem(vAdd,YvAdd);
   Z.AddItem(vAdd,ZvAdd);
   heaviside.AddItem(vAdd,heaviside.Get(maxVert));
   vertIdRegion.AddItem(vAdd,vertIdRegion.Get(maxVert));
   edgeSize.AddItem(vAdd,edgeSize.Get(maxVert));

   numVerts++;
   dVerts++;
   ipd[vertID]++;
  }
 }
 //cout << "  inserted by diffusion: " << ipd[vertID] << endl;
}

/*
 * Remove point(s) according to the solution of the diffusion equation
 * given by the class Helmholtz3D.
 *
 * */
void Model3D::remove3dMeshPointsByDiffusion()
{
 int vertID = 0;
 for( int e=0;e<mapEdge.DimI();e++ )
 {
  int v1 = mapEdge.Get(e,4);
  int v2 = mapEdge.Get(e,5);

  vertID = vertIdRegion.Get(v1);

  double x1=X.Get(v1);
  double y1=Y.Get(v1);
  double z1=Z.Get(v1);
  double x2=X.Get(v2);
  double y2=Y.Get(v2);
  double z2=Z.Get(v2);
  double length = distance(x1,y1,z1,x2,y2,z2);

  //int maxVert = max(v1,v2);
  int minVert = min(v1,v2);

  //double size = (edgeSize.Get(v1)+edgeSize.Get(v2))/2.0;
  double size = min(edgeSize.Get(v1),edgeSize.Get(v2));

  //cout << e << " " << length << " " << edgeSize.Get(v1) << endl;
  // edgeSize is the result of \nabla^2 edge = f
  if( length < 0.3*size && 
  //if( length < 0.7*size && 
	  minVert > surfMesh.numVerts )
  {
   cout << "  removed vert: " << minVert << endl;
   mark3DPointForDeletion(minVert);
   rpd[vertID]++;
  }
 }
 cout << "  removed by diffusion: " << rpd[vertID] << endl;
}

/* 
 * Implementation of the break-up model of 2 bubbles. Still not
 * finished.
 * */
void Model3D::breakup()
{
 for( int i=0;i<IEN.DimI();i++ )
 {
  int v1 = IEN.Get(i,0);
  int v2 = IEN.Get(i,1);
  int v3 = IEN.Get(i,2);
  int v4 = IEN.Get(i,3);
  if( heaviside.Get(v1) == 0.5 && heaviside.Get(v2) == 0.5 &&
	  heaviside.Get(v3) == 0.5 && heaviside.Get(v4) == 0.5 )
  {
   heaviside.Set(v1,0.0);
   heaviside.Set(v2,0.0);
   heaviside.Set(v3,0.0);
   heaviside.Set(v4,0.0);
  }
 }
}

/*
 * Insert point where the area of the surface triangle is bigger than a
 * value given by the "test" variable.
 *
 * */
void Model3D::insertPointsByArea()
{
 int lastRow;
 double test = 0.008;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  // P1
  int v1 = surfMesh.IEN.Get(i,0);
  double p1x = surfMesh.X.Get(v1);
  double p1y = surfMesh.Y.Get(v1);
  double p1z = surfMesh.Z.Get(v1);

  // P2
  int v2 = surfMesh.IEN.Get(i,1);
  double p2x = surfMesh.X.Get(v2);
  double p2y = surfMesh.Y.Get(v2);
  double p2z = surfMesh.Z.Get(v2);

  // P3
  int v3 = surfMesh.IEN.Get(i,2);
  double p3x = surfMesh.X.Get(v3);
  double p3y = surfMesh.Y.Get(v3);
  double p3z = surfMesh.Z.Get(v3);

  if(  getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z) > test )
  {
   int v4 = surfMesh.numVerts;

   // coords of triangle centroid
   double centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3) )/3.0;
   double centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3) )/3.0;
   double centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3) )/3.0;

   X.AddItem(v4,centroidX);
   Y.AddItem(v4,centroidY);
   Z.AddItem(v4,centroidZ);
   heaviside.AddItem(v4,heaviside.Get(v1)); 
   edgeSize.AddItem(v4,edgeSize.Get(v1));

   surfMesh.X.AddItem(v4,centroidX);
   surfMesh.Y.AddItem(v4,centroidY);
   surfMesh.Z.AddItem(v4,centroidZ);
   surfMesh.Marker.AddItem(v4,surfMesh.Marker.Get(v1));
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

/* 
 * strategy to ADD points - to be validated
 * This method works where two bubbles interact.
 *
 */
void Model3D::insertPointsBetweenBubblesByPosition()
{
 int ny = 4; // number of points between interfaces
 int nPoints = 20;

 double Xmax1=0.0;  double Ymax1=0.0; 
 double Xmin1=0.0;  double Ymin1=Y.Min();     
 double Zmax1=0.0;  double Ymax2=0.0; 
 double Zmin1=0.0;  double Ymin2=Y.Max();     

 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  // bubble 1 (Y<0)
  if( surfMesh.Y.Get(i) < 0 && surfMesh.Marker.Get(i)==0.5 )
  {
   if(X.Get(i)>Xmax1) Xmax1=X.Get(i);
   if(X.Get(i)<Xmin1) Xmin1=X.Get(i);

   if(Y.Get(i)>Ymin1) Ymin1=Y.Get(i);
   if(Y.Get(i)<Ymax1) Ymax1=Y.Get(i);

   if(Z.Get(i)>Zmax1) Zmax1=Z.Get(i);
   if(Z.Get(i)<Zmin1) Zmin1=Z.Get(i);
  }
  // bubble 2 (Y>0)
  if( surfMesh.Y.Get(i) > 0 && surfMesh.Marker.Get(i)==0.5 )
  {
   if(Y.Get(i)<Ymin2) Ymin2=Y.Get(i);
   if(Y.Get(i)>Ymax2) Ymax2=Y.Get(i);
  }
 }

 // initial position
 double xi = Xmin1;
 double yi = Ymin1;
 double zi = Zmin1;

 // distance between points
 double dx = (Xmax1-Xmin1)/(nPoints-1);
 double dy = (Ymin2-Ymin1)/(ny+1);
 double dz = (Zmax1-Zmin1)/(nPoints-1);

 // counter to numberize added points
 int count = surfMesh.numVerts;

 for( int i=0;i<nPoints;i++ )
 {
  for( int j=1;j<=ny;j++ )
  {
   for( int k=0;k<nPoints;k++ )
   {
	in.pointlist[3*count+0] = xi + dx*i;
	in.pointlist[3*count+1] = yi + dy*j;
	in.pointlist[3*count+2] = zi + dz*k;
	in.pointmarkerlist[count] = 11;

	count++;
   }
  }
 }
}
 
/* 
 * This method remesh completly the domain preserving only the points
 * located at surface and convex-hull. To do so, surfMesh.numVerts and
 * surfMesh.IEN need to be set on the beginning of the running program,
 * usually when the mesh is created from the .MSH file 
 *
 * */
void Model3D::mesh2Dto3DOriginal(const char* _param)
{
 // clean and init tetgen mesh object
 in.initialize();
 out.initialize();

 in.mesh_dim = 3;
 //in.numberofpoints = surfMesh.numVerts + 1600; // num of add points
 in.numberofpoints = surfMesh.numVerts; // num of add points
 numVerts = in.numberofpoints;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 convertModel3DtoTetgen(in);
 
 // add points between bubbles
 //insertPointsBetweenBubblesByPosition();

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");

 /*
  * Q: Quiet: No terminal output except errors.
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                | complete re-meshing the domain... ";
 tetrahedralize( (char*) _param,&in,&out );
 cout << "finished | " << resetColor() << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 convertTetgenToModel3D(out);
 mesh3d = convertTetgenToMesh3d(out);

 cout << endl;
 while( checkMeshQuality(out) == true ) 
 {
  in.initialize();
  out.initialize();

  in.mesh_dim = 3;
  in.numberofpoints = numVerts;
  in.pointlist = new REAL[in.numberofpoints * 3];
  in.pointmarkerlist = new int[in.numberofpoints];

  convertModel3DtoTetgen(in);

  cout << "----> fixing 3D mesh points... ";
  tetrahedralize( (char*) _param,&in,&out );
  cout << "finished <---- " << endl;;

  convertTetgenToModel3D(out);
  mesh3d = convertTetgenToMesh3d(out);
 }
 cout << endl;

 in.initialize();
 out.initialize();
}

void Model3D::mesh2Dto3DOriginal()
{
 mesh2Dto3DOriginal("QYYApa");
}

/* 
 * Convert Model3D data structure to TETGEN.
 * 
 * */
void Model3D::convertModel3DtoTetgen(tetgenio &_tetmesh)
{
 /* ------------ pontos da malha separados em 2 loops ------------ */
 // adiciona na estrutura tetgen as coordenadas dos pontos da 
 // superficie e do convex-hull
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  in.pointlist[3*i+0] = surfMesh.X.Get(i); 
  in.pointlist[3*i+1] = surfMesh.Y.Get(i); 
  in.pointlist[3*i+2] = surfMesh.Z.Get(i); 
  if( surfMesh.Marker.Get(i) == 0.0 ) // convex-hull
   in.pointmarkerlist[i] = 11;
  if( surfMesh.Marker.Get(i) == 0.5 ) // interface
   in.pointmarkerlist[i] = 22; // interface
 }

 // adicionando pontos que nao sao da interface e do convex-hull
 for( int i=surfMesh.numVerts;i<_tetmesh.numberofpoints;i++ )
 {
  in.pointlist[3*i+0] = X.Get(i);
  in.pointlist[3*i+1] = Y.Get(i);
  in.pointlist[3*i+2] = Z.Get(i);
  if( heaviside.Get(i) == 0.0 ) // fora da bolha
   in.pointmarkerlist[i] = 11;
  if( heaviside.Get(i) == 1.0 ) // dentro da bolha
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
  *  */
 /* 
  * This procedure defines specific regions in the volumetric mesh. To
  * do so it must define each region by a number starting to 1. 
  * Ex.: 1 --> wall
  *      2 --> bubble1
  *      3 --> bubble2, etc.
  * Note that the TETGEN indexes for attributes like this one does NOT
  * start at 0, thus we should sum 1 for each attribule. Thus, it is
  * possible to localize a point inside a specific region after the
  * meshing process by TETGEN.
  *
  * */

 // fluido interior + fluido exterior + superficie
 in.numberofregions = surfMesh.elemIdRegion.Max()+1; 
 in.regionlist = new REAL[in.numberofregions*5];

 setNeighbourSurfaceElem();

 // fora e dentro das bolhas
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
 {
  int node;
  //double curv = fabs(surfMesh.curvature.Get(node));
  // find the first vertex with region == nb
  //for( int i=0;i<surfMesh.numVerts;i++ )
  for( int i=surfMesh.numVerts-1;i>=0;i-- )
  {
   clVector myVec = getNormalAndKappa(i,getNeighbourSurfacePoint(i));
   double curv = fabs(myVec.Get(0));
  //for( int i=0;i<surfMesh.numVerts;i++ )
   if( surfMesh.vertIdRegion.Get(i) == nb && curv < 20 )
   //if( surfMesh.vertIdRegion.Get(i) == nb )
   {
	node = i;
	break;
   }
  }

  clVector myVec = getNormalAndKappa(node,getNeighbourSurfacePoint(node));
  double xIn = surfMesh.X.Get(node)-0.1*triEdge[nb]*myVec.Get(1);
  double yIn = surfMesh.Y.Get(node)-0.1*triEdge[nb]*myVec.Get(2);
  double zIn = surfMesh.Z.Get(node)-0.1*triEdge[nb]*myVec.Get(3);

  double edge = triEdge[nb];
  if( edgeSize.Dim() > 0 )
   edge = edgeSize.Max();

  in.regionlist[5*nb+0] = xIn;
  in.regionlist[5*nb+1] = yIn;
  in.regionlist[5*nb+2] = zIn;
  in.regionlist[5*nb+3] = nb+1;
  in.regionlist[5*nb+4] = 10*edge*edge*edge*1.4142/12.0;
  //in.regionlist[5*nb+4] = tetVol[nb];
//--------------------------------------------------
//   cout << " ------ " << node << " ------" << endl;
//   cout << "triEdge: " << triEdge[nb] << endl;
//   cout << "xIn: " << xIn << endl; 
//   cout << "yIn: " << yIn << endl; 
//   cout << "zIn: " << zIn << endl; 
//   cout << surfMesh.X.Get(node) << "      " << "xNormal: " << myVec.Get(1) << endl;
//   cout << surfMesh.Y.Get(node) << "      " << "yNormal: " << myVec.Get(2) << endl;
//   cout << surfMesh.Z.Get(node) << "      " << "zNormal: " << myVec.Get(3) << endl;
//-------------------------------------------------- 
 }

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
  double mSum = surfMesh.Marker.Get(v1) + 
	          surfMesh.Marker.Get(v2) + 
			  surfMesh.Marker.Get(v3);
  // melhorar esta configuracao de facet para bolha e convex hull
  if( mSum == 0.5 ) // bubble/drop
   in.facetmarkerlist[i] = 10;
  else if( mSum == 0.0 ) // out
   in.facetmarkerlist[i] = 20;
  else // in
   in.facetmarkerlist[i] = 30;
 }

//--------------------------------------------------
//  // hole in 3D mesh
//  in.numberofholes = 1; 
//  in.holelist = new REAL[3];
//  in.holelist[0] = 5.0;
//  in.holelist[1] = 5.0;
//  in.holelist[2] = 5.0;
//-------------------------------------------------- 
}

/* 
 * Convert Model3D.SurfaceMesh data structure to TETGEN.
 * 
 * */
tetgenio Model3D::convertSurfaceMeshToTetGen(SurfaceMesh _mesh,
                                             tetgenio &_tetmesh)
{
 // add to tetgen struct the point coords
 for( int i=0;i<_mesh.numVerts;i++ )
 {
  _tetmesh.pointlist[3*i+0] = _mesh.X.Get(i);
  _tetmesh.pointlist[3*i+1] = _mesh.Y.Get(i);
  _tetmesh.pointlist[3*i+2] = _mesh.Z.Get(i);
  if( _mesh.Marker.Get(i) == 0.0 )
   _tetmesh.pointmarkerlist[i] = 11;
  if( _mesh.Marker.Get(i) == 0.5 )
   _tetmesh.pointmarkerlist[i] = 22; // same id of facetmarker
 }

 // definindo regiao fora da bolha 
 _tetmesh.regionlist[0] = surfMesh.X.Min()+0.01;
 _tetmesh.regionlist[1] = surfMesh.Y.Min()+0.01;
 _tetmesh.regionlist[2] = surfMesh.Z.Min()+0.01;
 _tetmesh.regionlist[3] = 1;

 tetgenio::facet *f;   // Define a pointer of facet. 
 tetgenio::polygon *p; // Define a pointer of polygon.

 // defining the interface and convex-hull
 // definindo a superficie da bolha e convex-hull
 for( int i=0;i<_mesh.numElems;i++ )
 {
  int v1 = _mesh.IEN.Get(i,0);
  int v2 = _mesh.IEN.Get(i,1);
  int v3 = _mesh.IEN.Get(i,2);
  f = &_tetmesh.facetlist[i];
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
  double mSum = surfMesh.Marker.Get(v1) + 
	          surfMesh.Marker.Get(v2) + 
			  surfMesh.Marker.Get(v3);
  // melhorar esta configuracao de facet para bolha e convex hull
  if( mSum == 0.5 ) // bubble/drop
   in.facetmarkerlist[i] = 10;
  else if( mSum == 0.0 ) // out
   in.facetmarkerlist[i] = 20;
  else // in
   in.facetmarkerlist[i] = 30;
 }
 return _tetmesh;
}

/* 
 * Convert TETGEN data structure to Model3D.Mesh3D
 * 
 * */
Mesh3D Model3D::convertTetgenToMesh3d(tetgenio &_tetmesh)
{
 Mesh3D mesh;
 mesh.numElems = _tetmesh.numberoftetrahedra;
 mesh.numNodes = _tetmesh.numberofpoints+_tetmesh.numberoftetrahedra;
 mesh.numVerts = _tetmesh.numberofpoints;
 mesh.IEN.Dim(mesh.numElems,4);
 mesh.Marker.Dim(mesh.numVerts);
 mesh.Marker.SetAll(0.0);

 // varre lista de elementos e passa para estrutura IEN
 for( int i=0;i<_tetmesh.numberoftetrahedra;i++ )
 {
  // setting de heaviside = 0 para fora da bolha e heaviside = 0.5 para interface
  if( _tetmesh.tetrahedronattributelist[i] == 1 )
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = _tetmesh.tetrahedronlist[i*4+j];
	mesh.IEN.Set(i,j,vertice);
	mesh.Marker.Set(vertice,0.0);
   }
  }
  // setting de heaviside = 1 para dentro da bolha e heaviside = 0.5 para interface
  else 
  {
   for( int j=0;j<4;j++ )
   {
	int vertice = _tetmesh.tetrahedronlist[i*4+j];
	mesh.IEN.Set(i,j,vertice);
	mesh.Marker.Set(vertice,1.0);
   }
  }
 }

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 mesh.X.Dim(numNodes);
 mesh.Y.Dim(numNodes);
 mesh.Z.Dim(numNodes);
 for( int i=0;i<mesh.numVerts;i++ )
 {
  mesh.X.Set(i,_tetmesh.pointlist[3*i+0]);
  mesh.Y.Set(i,_tetmesh.pointlist[3*i+1]);
  mesh.Z.Set(i,_tetmesh.pointlist[3*i+2]);
  if( _tetmesh.pointmarkerlist[i] == 10 ||
	  _tetmesh.pointmarkerlist[i] == 22 )
   mesh.Marker.Set(i,0.5);
 }

 return mesh;
}

/* 
 * Convert TETGEN data structure to Model3D
 * 
 * */
void Model3D::convertTetgenToModel3D(tetgenio &_tetmesh)
{
 int numVertsOld = numVerts;
 numElems = _tetmesh.numberoftetrahedra;
 numNodes = _tetmesh.numberofpoints+_tetmesh.numberoftetrahedra;
 numVerts = _tetmesh.numberofpoints;
 
 // varre lista de elementos e passa para estrutura IEN
 IEN.Dim(numElems,4);
 heaviside.Dim(numVerts);
 heaviside.SetAll(0.0);
 vertIdRegion.Dim(numVerts);
 elemIdRegion.Dim(numElems);
 elemIdRegion.SetAll(0.0);
 //edgeSize.Dim(numVerts);
 for( int i=0;i<_tetmesh.numberoftetrahedra;i++ )
 {
  elemIdRegion.Set(i,_tetmesh.tetrahedronattributelist[i]-1);

  // set:
  // heaviside = 0 outside bubble
  // elemIdRegion = out.tetrahedronattributelist[i]
  for( int j=0;j<4;j++ )
  {
   int vertice = _tetmesh.tetrahedronlist[i*4+j];
   IEN.Set(i,j,vertice);
   vertIdRegion.Set(vertice,_tetmesh.tetrahedronattributelist[i]-1);
   //edgeSize.Set(vertice,triEdge[_tetmesh.tetrahedronattributelist[i]-1]);

   // setting heaviside
   if( _tetmesh.tetrahedronattributelist[i] == 1 ) // wall and out mesh
	heaviside.Set(vertice,0.0);
   else
	heaviside.Set(vertice,1.0);
  }
 }

 // atualizando valores de X,Y,Z,uc,vc,wc e pc
 X.Dim(numNodes);
 Y.Dim(numNodes);
 Z.Dim(numNodes);
 for( int i=0;i<numVerts;i++ )
 {
  X.Set(i,_tetmesh.pointlist[3*i+0]);
  Y.Set(i,_tetmesh.pointlist[3*i+1]);
  Z.Set(i,_tetmesh.pointlist[3*i+2]);
  if( _tetmesh.pointmarkerlist[i] == 10 ||
	  _tetmesh.pointmarkerlist[i] == 22 )
  {
   heaviside.Set(i,0.5);
   vertIdRegion.Set(i,0.5);
  }
 }

 // updating edgeSize vector due to increasing of mesh points
 clVector zeros(numVerts-numVertsOld);
 edgeSize.Append(zeros);
}

/*
 * Test tetrahedron mesh quality.
 *
 * */
bool Model3D::checkMeshQuality(tetgenio &_tetmesh)
{
 bool flag = false;

 for( int i=0;i<_tetmesh.numberoftetrahedra;i++ )
 {
  int v1 = _tetmesh.tetrahedronlist[i*4+0];
  int v2 = _tetmesh.tetrahedronlist[i*4+1];
  int v3 = _tetmesh.tetrahedronlist[i*4+2];
  int v4 = _tetmesh.tetrahedronlist[i*4+3];
  int region = _tetmesh.tetrahedronattributelist[i]-1;

  int elemID = elemIdRegion.Get(i);

  if( _tetmesh.pointmarkerlist[v1] == 22 &&
	  _tetmesh.pointmarkerlist[v2] == 22 &&
	  _tetmesh.pointmarkerlist[v3] == 22 &&
	  _tetmesh.pointmarkerlist[v4] == 22 )
  {
   double xMid = ( _tetmesh.pointlist[3*v1+0]+
                 _tetmesh.pointlist[3*v2+0]+
				 _tetmesh.pointlist[3*v3+0]+
				 _tetmesh.pointlist[3*v4+0] )/4.0;

   double yMid = ( _tetmesh.pointlist[3*v1+1]+
                 _tetmesh.pointlist[3*v2+1]+
				 _tetmesh.pointlist[3*v3+1]+
				 _tetmesh.pointlist[3*v4+1] )/4.0;

   double zMid = ( _tetmesh.pointlist[3*v1+2]+
                 _tetmesh.pointlist[3*v2+2]+
				 _tetmesh.pointlist[3*v3+2]+
				 _tetmesh.pointlist[3*v4+2] )/4.0;

   X.AddItem(numVerts,xMid);
   Y.AddItem(numVerts,yMid);
   Z.AddItem(numVerts,zMid);
   heaviside.AddItem(fabs(region-1));
   edgeSize.AddItem(numVerts,edgeSize.Get(v1));
   numVerts++;

   badtet[elemID]++;
   flag = true;
  }
 }
 return flag;
}

/*
 * Remove point(s) according to the tetrahedron volume. If such a volume
 * is smaller then the variable "tetVol", one of the four nodes are
 * deleted.
 *
 * */
void Model3D::remove3dMeshPointsByVolume()
{
 double vSum;
 double vertSum;
 int v[NUMGLE];
 int elemID;
 int vert=0;

 tetVol.clear();
 tetVol.resize(triEdge.size()); // number of surfaces + 1
 double triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));

 // set tetVol ---> wall,bubble1, bubble2 etc.
 for( int v=0;v<(int) triEdge.size();v++ )
  tetVol[v] = triEdgeMin*triEdgeMin*triEdgeMin*sqrt(2.0)/12.0;

 for( int elem=0;elem<numElems;elem++ )
 {
  v[0] = IEN.Get(elem,0);
  v[1] = IEN.Get(elem,1);
  v[2] = IEN.Get(elem,2);
  v[3] = IEN.Get(elem,3);

  double hSum = heaviside.Get(v[0])+heaviside.Get(v[1])+
              heaviside.Get(v[2])+heaviside.Get(v[3]);

  elemID = elemIdRegion.Get(elem);

  double maxVol = max(1.0E-05,0.01*tetVol[elemIdRegion.Get(elem)]);

  int count=0;
  if( hSum != 2.0 && 
	  fabs(getVolume(elem)) < maxVol ) 
  {
   // add to checkVert only non surface vertex
   list<int> checkVert;
   for( int j=0;j<NUMGLE;j++ )
	if( heaviside.Get(v[j]) != 0.5 )
	  checkVert.push_back( v[j] );

   vertSum = 1.0E+17; // initial value

   // check sum of volumes of each non surface vertex neighbours
   list<int> plist = checkVert;
   for(list<int>::iterator mvert=plist.begin(); mvert!= plist.end();++mvert )
   {
	vSum = 0;
	list<int> plist2 = neighbourElem.at(*mvert);
	for(list<int>::iterator mele=plist2.begin(); mele != plist2.end();++mele )
	 vSum += fabs(getVolume(*mele));

	// The problem is: if all the vertices are lower then
	// surfMesh.numVerts, 
	if( vSum < vertSum && *mvert > surfMesh.numVerts )
	{
	 vertSum = vSum;
	 vert = *mvert;
	 count++;
	}
   }
   // mark points to delete
   if( count > 0 )
   {
	mark3DPointForDeletion(vert);
	rpv[elemID]++;
   }
  }
 }
 cout << "  removed by volume: " << rpv[elemID] << endl;
}

/*
 * Mark point of the 3D mesh structure to be deleted by setting the
 * heaviside vector to -1
 *
 * */
void Model3D::mark3DPointForDeletion(int _vert)
{
 heaviside.Set(_vert,-1);
}

/* 
 * Perform 3D mesh point deletion according to the vector heaviside (set
 * -1)
 *
 * */
void Model3D::delete3DPoints()
{
 for( int dp=0;dp<heaviside.Dim();dp++ )
 {
  if( heaviside.Get(dp) == -1 )
  {
   X.Delete(dp);
   Y.Delete(dp);
   Z.Delete(dp);
   heaviside.Delete(dp);
   //interfaceDistance.Delete(dp);
   edgeSize.Delete(dp);
   numVerts--;
   dVerts--;
   dp--;
  }
 }
}

/* 
 * TETGEN output mesh report.
 *
 * */
void Model3D::triMeshStats()
{
 // fora e dentro das bolhas
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 vector<double> sumLength;sumLength.clear();
 sumLength.resize(surfMesh.elemIdRegion.Max()+1);
 vector<double> sumArea;sumArea.clear();
 sumArea.resize(surfMesh.elemIdRegion.Max()+1);
 for( int e=0;e<surfMesh.numElems;e++ )
 {
  int v1 = surfMesh.IEN.Get(e,0);
  double p1x = surfMesh.X.Get(v1); 
  double p1y = surfMesh.Y.Get(v1); 
  double p1z = surfMesh.Z.Get(v1); 

  int v2=surfMesh.IEN.Get(e,1);
  double p2x = surfMesh.X.Get(v2); 
  double p2y = surfMesh.Y.Get(v2); 
  double p2z = surfMesh.Z.Get(v2); 

  int v3=surfMesh.IEN.Get(e,2);
  double p3x = surfMesh.X.Get(v3); 
  double p3y = surfMesh.Y.Get(v3); 
  double p3z = surfMesh.Z.Get(v3); 

  int elemID = surfMesh.elemIdRegion.Get(e);

  double length12 = distance(p1x,p1y,p1z,p2x,p2y,p2z);
  double length13 = distance(p1x,p1y,p1z,p3x,p3y,p3z);
  double length23 = distance(p2x,p2y,p2z,p3x,p3y,p3z);

  sumLength[elemID] += ( length12+length13+length23 )/3.0;

  double auxMinLength = min( length12,length13 );
  auxMinLength = min( auxMinLength,length23 );

  double auxMaxLength = max( length12,length13 );
  auxMaxLength = max( auxMaxLength,length23 );

  double area = fabs(getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z));
  sumArea[elemID] += area;

  // areas
  if( area < minArea[elemID] ) 
  {
   minArea[elemID] = area;
   idMinArea[elemID] = e;
  }
  if( area > maxArea[elemID] ) 
  {
   maxArea[elemID] = area;
   idMaxArea[elemID] = e;
  }

  // lengths
  if( auxMinLength < minLength[elemID] ) 
   minLength[elemID] = auxMinLength;

  if( auxMaxLength > maxLength[elemID] ) 
   maxLength[elemID] = auxMaxLength;

  numSurfElems[elemID]++;
 }
 
 // average length
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
 {
  averageTriLength[nb] = sumLength[nb]/numSurfElems[nb];
  averageTriArea[nb] = sumArea[nb]/numSurfElems[nb];
 }
}

void Model3D::tetMeshStats()
{
 /* ******************************************** *
  *                                              *
  * regular    a^3 sqrt(2)    triEdge^3 sqrt(2)  *
  *   tet    = ----------- = ------------------  *
  * volume         12                12          *
  *                                              *
  * ******************************************** */

 vector<double> sumVolume;sumVolume.clear();
 sumVolume.resize(elemIdRegion.Max()+1);
 vector<int> count;count.clear();
 count.resize(elemIdRegion.Max()+1);
 for( int e=0;e<numElems;e++ )
 {
  int v1 = IEN.Get(e,0);
  int v2 = IEN.Get(e,1);
  int v3 = IEN.Get(e,2);
  int v4 = IEN.Get(e,3);

  int elemID = elemIdRegion.Get(e);
  double volume = fabs(getVolume(e));

  if( heaviside.Get(v1) == 0.5 &&
	  heaviside.Get(v2) == 0.5 &&
	  heaviside.Get(v3) == 0.5 &&
	  heaviside.Get(v4) == 0.5 )
   intet[elemID]++;

  sumVolume[elemID] += volume;

  if( volume < minVolume[elemID] ) 
  {
   minVolume[elemID] = volume;
   idMinVolume[elemID] = e;
  }
  if( volume > maxVolume[elemID] ) 
  {
   maxVolume[elemID] = volume;
   idMaxVolume[elemID] = e;
  }
  count[elemID]++;
 }
 for( int nb=0;nb<=elemIdRegion.Max();nb++ )
  averageTetVolume[nb] = sumVolume[nb]/count[nb];
}

/*
 * Mesh a set of points using TETGEN.
 *
 * */
void Model3D::mesh3DPoints(const char* _param)
{
 // init tetgen mesh object
 in.initialize();
 out.initialize();

 in.mesh_dim = 3;
 in.numberofpoints = numVerts;
 in.pointlist = new REAL[in.numberofpoints * 3];
 in.pointmarkerlist = new int[in.numberofpoints];

 convertModel3DtoTetgen(in);

 //in.save_poly("bubble");
 //in.save_nodes("bubble");
 //in.save_elements("in");
 /*
  * Q: Quiet: No terminal output except errors.
  * Y: mesh boundary preserving
  * YY: no mesh boundary points inserted
  * R: mesh coarsening
  * C: Checks the consistency of the final mesh.
  * A: Assigns attributes to identify tetrahedra in certain regions.
  * a: Applies a maximum tetrahedron volume constraint.
  * p:  Tetrahedralizes a picecwise linear complex
  * q: Quality mesh generation. Min radius-edge ratio may be specifyed
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                     | re-meshing 3D points... ";
 tetrahedralize( (char*) _param,&in,&out ); 
 cout << "finished | " << resetColor() << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 //out.save_elements("out");
 //out.save_nodes("out");
 //out.save_poly("out");
 //out.save_faces("out");

 dVerts = out.numberofpoints - numVerts;

 convertTetgenToModel3D(out);
 mesh3d = convertTetgenToMesh3d(out);

 //breakup();
 
 in.initialize();
 out.initialize();
}

void Model3D::mesh3DPoints()
{
 //tetrahedralize( (char*) "QYYRCApq1.414q10a",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApqq10a",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApqq10",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYAp",&in,&out ); // no insertion of points
 mesh3DPoints("QYYApa");
}

void Model3D::setDiskCouetteBC()
{
 double aux;
 rMax = Y.Max(); 

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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

void Model3D::setInfiniteDiskBC(double _F,double _G, double _H)
{
 double omega,aux,radius;
 rMax = Y.Max();

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

  if( Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   ////idbcp.AddItem(i); // caso com c.c. livre em w

   uc.Set(i,radius*_F); // Z=10
   vc.Set(i,radius*_G); // Z=10
   wc.Set(i,(-1.0)*_H); // Z=10
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
   outflow.Set(i,aux);

   if( i < numVerts )
   {
	idbcp.AddItem(i);
	aux = 0.0;
	pc.Set(i,aux);
   }
  }
 }
}

void Model3D::setInfiniteSphereBC(double _F,double _G, double _H)
{
 //double radius;
 double omega,aux;
 rMax = Y.Max();

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  //radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

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

  if( (Z.Get(i)<Z.Max() && 
	   Z.Get(i)>Z.Min() && 
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) > 
	   (rMax*rMax - 0.001) ) ) ||
	  (Z.Get(i) == Z.Max()))
  {
   outflow.Set(i,aux);

   if( i < numVerts )
   {
	idbcp.AddItem(i);
	aux = 0.0;
	pc.Set(i,aux);
   }
  }
 }
}

void Model3D::setFiniteDiskBC()
{
 double omega,aux;
 rMax = Y.Max();

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  if( Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0); 
   vc.Set(i,0.0); 
   wc.Set(i,0.0); 
  }

  if( Z.Get(i) == Z.Min() && 
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)<rMax/4.0*rMax/4.0) )
  {
   //radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   omega=1.0;

   //aux = (-1.0)*omega*radius;
   aux = (-1.0)*Y.Get(i)*omega;
   uc.Set(i,aux);
   //aux = omega*radius;
   aux = X.Get(i)*omega;
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);
  }

  // free-surface boundary condition
//--------------------------------------------------
//   if( Z.Get(i) == Z.Min() && 
// 	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)<rMax/4.0*rMax/4.0) )
//   {
//    idbcw.AddItem(i);
// 
//    wc.Set(i,0.0); 
//   }
//-------------------------------------------------- 

  // casca do cilindro
  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   uc.Set(i,0.0); 
   vc.Set(i,0.0); 
   wc.Set(i,0.0); 
  }
 }

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Min() && 
	 (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>=rMax/4.0*rMax/4.0) )
  {
   idbcp.AddItem(i);
   aux = 0.0;
   pc.Set(i,aux);
  }
 }
}

void Model3D::setCDiskBC()
{
 double aux;
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


void Model3D::setGenericBC()
{    
 /* This IF selects the priority boundary conditions to be set on the
  * phyBounds vector. Note that only these names will be written on top
  * of the others. For instance if a corner point has 2 types of
  * boundary condition, the phyNames below will be written.
  * 
  * \remark PBC has priority over all the conditions already declared, 
  * namely:
  * 
  * - slip walls    : 3rd priority
  * - inflow walls  : 2nd priority
  * - no-slip walls : 1st priority
  *
  * For this reason, any corner point will be set with periodic condition
  * and the new ordering is:
  *
  * - periodic walls : 4th priority
  * - slip walls     : 3rd priority
  * - inflow walls   : 2nd priority
  * - no-slip walls  : 1st priority
  *
  * */ 

  /* This algorithm sets up physical groups of periodic boundary
  * conditions for a unique pair of master/slave boundaries defined as
  * 
  * "wallLeft"  and   "wallRight" 
  *
  * at ".geo" file. 
  *
  */

 // auxiliary x,y,z pbc coordinates
 vector<double> xpbcMaster(0);
 vector<double> ypbcMaster(0);
 vector<double> zpbcMaster(0);
 vector<double> xpbcSlave(0);
 vector<double> ypbcSlave(0);
 vector<double> zpbcSlave(0);
 
 for( int i = 0; i < surfMesh.numElems; i++ )
 {		
  	int v1 = surfMesh.IEN.Get(i,0);
  	int v2 = surfMesh.IEN.Get(i,1);
	int v3 = surfMesh.IEN.Get(i,2);
  	int id = surfMesh.idRegion.Get(i);
	
	// left indices
  	if ( surfMesh.phyNames.at(id).compare(5,4,"Left") == 0 )
  	{
		pbcIndicesLeft.push_back(v1);
		pbcIndicesLeft.push_back(v2);
		pbcIndicesLeft.push_back(v3);
	}

	// right indices
  	if ( surfMesh.phyNames.at(id).compare(5,5,"Right") == 0 )
	{
		pbcIndicesRight.push_back(v1);
		pbcIndicesRight.push_back(v2);
		pbcIndicesRight.push_back(v3);

	}
	
 }

 // removing duplicatas
 set<int> setOne( pbcIndicesLeft.begin(), pbcIndicesLeft.end() );
 set<int> setTwo( pbcIndicesRight.begin(), pbcIndicesRight.end() );
 set<int>::iterator itsetOne;
 set<int>::iterator itsetTwo;

 // resizing to reallocate
 pbcIndicesLeft.resize(0);
 pbcIndicesRight.resize(0);
			
 for (itsetOne = setOne.begin(); itsetOne != setOne.end(); ++itsetOne)
 {
		cout << "Index: " << *itsetOne << endl;
		pbcIndicesLeft.push_back(*itsetOne);
				 
		// master coordinates
		double xL = surfMesh.X.Get(*itsetOne);
		double yL = surfMesh.Y.Get(*itsetOne);
		double zL = surfMesh.Z.Get(*itsetOne);
		xpbcMaster.push_back(xL);
		ypbcMaster.push_back(yL);
		zpbcMaster.push_back(zL);
 }
			 
	cout << "Master nodes stored." << endl;
			
 for (itsetTwo = setTwo.begin(); itsetTwo != setTwo.end(); ++itsetTwo)
 {
		cout << "Index: " << *itsetTwo << endl;
		pbcIndicesRight.push_back(*itsetTwo);

		// slave coordinates
		double xR = surfMesh.X.Get(*itsetTwo);
		double yR = surfMesh.Y.Get(*itsetTwo);
		double zR = surfMesh.Z.Get(*itsetTwo);
		xpbcSlave.push_back(xR);
		ypbcSlave.push_back(yR);
		zpbcSlave.push_back(zR);
 }

	cout << "Slave nodes stored." << endl;

	int sizeM = pbcIndicesLeft.size();
	int sizeS = pbcIndicesRight.size();
			
	// spatial correspondence checking
	if ( sizeM != sizeS )
	{
		string warn = "Master/Slave nodes differing in quantity! PBC not applicable.";
		cerr << warn << endl;
	}
	// pairing checking
	else
	{
		for (int is = 0; is < sizeM; ++is)
		{
			double yM = ypbcMaster.at(is);
			double zM = zpbcMaster.at(is);
			double yS = ypbcSlave.at(is);
			double zS = zpbcSlave.at(is);
			
			// simple extrusion assumed
			if ( ( fabs( yS - yM ) > EPS ) && 
			     ( fabs( zS - zM ) > EPS ) )
			{
				cout << "Entry not periodic: " << is << endl;	
				cout << "ibL: "  << pbcIndicesLeft.at(is) << "\t yM = " << yM << "\t zM = " << zM << endl;
				cout << "ibR: " << pbcIndicesRight.at(is) << "\t yS = " << yS << "\t zS = " << zS << endl;
			}
			
	 	}
		
	 }

 for( int i=0; i < surfMesh.numElems; i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  int id = surfMesh.idRegion.Get(i);

  // 3nd. priority
  if( surfMesh.phyNames.at(id).compare(5,7,"NormalU") == 0 || 
      surfMesh.phyNames.at(id).compare(5,7,"NormalV") == 0 ||
      surfMesh.phyNames.at(id).compare(5,7,"NormalW") == 0 )
  {
   string aux = surfMesh.phyNames.at(id);
   surfMesh.phyBounds.at(v1) = aux;
   surfMesh.phyBounds.at(v2) = aux;
   surfMesh.phyBounds.at(v3) = aux;
  }
 }
 
 for( int i=0; i < surfMesh.numElems; i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  int id = surfMesh.idRegion.Get(i);

  // 2nd. priority
  if( surfMesh.phyNames.at(id).compare(5,7,"InflowU") == 0 || 
      surfMesh.phyNames.at(id).compare(5,7,"InflowV") == 0 || 
      surfMesh.phyNames.at(id).compare(5,7,"InflowW") == 0 || 
      surfMesh.phyNames.at(id).compare(5,16,"InflowUParabolic") == 0 || 
      surfMesh.phyNames.at(id).compare(5,16,"InflowVParabolic") == 0 || 
      surfMesh.phyNames.at(id).compare(5,16,"InflowWParabolic") == 0 )
  {
   string aux = surfMesh.phyNames.at(id);
   surfMesh.phyBounds.at(v1) = aux;
   surfMesh.phyBounds.at(v2) = aux;
   surfMesh.phyBounds.at(v3) = aux;
  }
 }

 for( int i=0; i < surfMesh.numElems; i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  int id = surfMesh.idRegion.Get(i);
  
  // 1st. priority
  if( surfMesh.phyNames.at(id).compare(5,6,"NoSlip") == 0 || 
      surfMesh.phyNames.at(id).compare(5,19,"NoSlipConcentration") == 0 || 
      surfMesh.phyNames.at(id).compare(5,14,"NoSlipPressure") == 0 || 
      surfMesh.phyNames.at(id).compare(5,4,"InvU") == 0 || 
      surfMesh.phyNames.at(id).compare(5,4,"InvV") == 0 || 
      surfMesh.phyNames.at(id).compare(5,4,"InvW") == 0 ||
      surfMesh.phyNames.at(id).compare(5,14,"Inflow2Bubbles") == 0 ||
      surfMesh.phyNames.at(id).compare(5,17,"Inflow2AxiBubbles") == 0 )
  {
   string aux = surfMesh.phyNames.at(id);
   surfMesh.phyBounds.at(v1) = aux;
   surfMesh.phyBounds.at(v2) = aux;
   surfMesh.phyBounds.at(v3) = aux;
  }
 }

 // calculating channel's diameter.
 double diameterXY = ( dist(X.Min(),X.Max()) + 
                     dist(Y.Min(),Y.Max()) ) / 2.0;
 double diameterXZ = ( dist(X.Min(),X.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
 double diameterYZ = ( dist(Y.Min(),Y.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
//--------------------------------------------------
//  double diameterYZ = distance(Y.Min(),Z.Min(),Y.Max(),Z.Max());
//-------------------------------------------------- 

 int count = 0;
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  // outflow condition
  if( surfMesh.phyBounds.at(*it) == "\"wallOutflow\"" )
  {
   idbcp.AddItem(*it);

   pc.Set(*it,0.0);
  }

  // inflow condition U
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowUParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( Y.Get(*it)*Y.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   double Umax = 1.0;
   double aux = 1.0*Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
	                                         (diameterYZ/2.0)) );

   //aux=1.0;
   uc.Set(*it,aux);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }

  // inflow condition V
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowVParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( X.Get(*it)*X.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   double Vmax = 2.0;
   double aux = Vmax*( 1.0-radius*radius/((diameterXZ/2.0)*
	                                    (diameterXZ/2.0)) );

   uc.Set(*it,0.0);
   vc.Set(*it,aux);
   wc.Set(*it,0.0);
  }

  // inflow condition W
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowWParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( X.Get(*it)*X.Get(*it) + Y.Get(*it)*Y.Get(*it) );

   // Parabolic profile
   double Wmax = 2.0;
   double aux = Wmax*( 1.0-radius*radius/((diameterXY/2.0)*
	                                    (diameterXY/2.0)) );

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,aux);
  }

  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowU\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  
  // inflow condition V
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowV\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,1.0);
   wc.Set(*it,0.0);
  }

  // inflow condition W
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,1.0);
  }

  // 2 bubbles inflow condition
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflow2Bubbles\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double aux = X.Get(*it);
   uc.Set(*it,aux);
   aux = (-1.0)*Y.Get(*it);
   vc.Set(*it,aux);
   aux = Z.Get(*it);
   wc.Set(*it,aux);
  }

  // 2 Axi bubbles inflow condition
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflow2AxiBubbles\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double aux = X.Get(*it);
   uc.Set(*it,aux);
   aux = (-1.0)*Y.Get(*it);
   vc.Set(*it,aux);
   aux = 0.0;
   wc.Set(*it,aux);
  }

  // moving boundary condition as inflow set to Zero
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowZeroU\"" || 
           surfMesh.phyBounds.at(*it) == "\"wallInflowZeroV\"" ||
		   surfMesh.phyBounds.at(*it) == "\"wallInflowZeroW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }

  // moving boundary U
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvU\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,-1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }

  // moving boundary V
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvV\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,-1.0);
   wc.Set(*it,0.0);
  }

  // moving boundary W
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,-1.0);
  }
  
  // NoSlip with Concentration b.c.
  else if( surfMesh.phyBounds.at(*it) == "\"wallNoSlipConcentration\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
   idbcc.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
   cc.Set(*it,1.0);
  }

  // moving boundary U
  else if( surfMesh.phyBounds.at(*it) == "\"wallNoSlipPressure\"" )
  {
   if( X.Get(*it) == X.Max() && 
	   Z.Get(*it) == Z.Min() && 
	   count < 1 )
   {
	idbcp.AddItem(*it);

	pc.Set(*it,0.0);
   }
   else
   {
	idbcu.AddItem(*it);
	idbcv.AddItem(*it);
	idbcw.AddItem(*it);

	uc.Set(*it,-1.0);
	vc.Set(*it,0.0);
	wc.Set(*it,0.0);
   }
   if( Z.Get(*it) == Z.Min() )
   {
	idbcc.AddItem(*it);

	cc.Set(*it,1.0);
   }
  }

  // symmetry boundary U
  else if( surfMesh.phyBounds.at(*it) == "\"wallNormalU\"" )
  {
   idbcu.AddItem(*it);
   uc.Set(*it,0.0);
  }

  // symmetry boundary V
  else if( surfMesh.phyBounds.at(*it) == "\"wallNormalV\"" )
  {
   idbcv.AddItem(*it);
   vc.Set(*it,0.0);
  }

  // symmetry boundary W
  else if( surfMesh.phyBounds.at(*it) == "\"wallNormalW\"" )
  {
   idbcw.AddItem(*it);
   wc.Set(*it,0.0);
  }
  // no slip condition if any other is imposed
  //if( surfMesh.phyBounds.at(*it) == "\"wallNoSlip\"" )
  else
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
 }
 //integralParabolic();
}

void Model3D::setBubbleArrayPeriodicBC()
{

	for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
    {
        if ( surfMesh.phyBounds.at(*it) == "\"wallMoving\"" )
        {
	        idbcu.AddItem(*it);
	        idbcv.AddItem(*it);
	        idbcw.AddItem(*it);
	        idbcp.AddItem(*it);
	 
	        uc.Set(*it,-1.0);
	        vc.Set(*it,0.0);
	        wc.Set(*it,0.0);
	        pc.Set(*it,0.0);
	     }
	 }
}

void Model3D::setOnePointPressureBC()
{
	for (int i = 0; i < numVerts; ++i)
	{
	    // Neumman for pressure ("one-point")
		if ( ( fabs( X.Get(i) - 0.5*X.Max() ) < 0.01 ) &&
		     ( fabs( Y.Get(i) - Y.Min()     ) < 0.01 ) &&
		     ( fabs( Z.Get(i) - 0.5*Z.Max() ) < 0.01 ) )
		{
			idbcp.AddItem(i);
			pc.Set(i,0.0);
			cout << "Pressure index set: "<< i << endl;
			break;
		}
		   
	}
}


void Model3D::setGenericBC(double _vel)
{    
 clearBC();

 // calculating channel's diameter.
 double diameterXY = ( dist(X.Min(),X.Max()) + 
                     dist(Y.Min(),Y.Max()) ) / 2.0;
 double diameterXZ = ( dist(X.Min(),X.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
 double diameterYZ = ( dist(Y.Min(),Y.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
//--------------------------------------------------
//  double diameterYZ = distance(Y.Min(),Z.Min(),Y.Max(),Z.Max());
//-------------------------------------------------- 

 int count = 0;
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  // outflow condition
  if( surfMesh.phyBounds.at(*it) == "\"wallOutflow\"" )
  {
   idbcp.AddItem(*it);
  
   pc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowU\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowUParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( Y.Get(*it)*Y.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   double Umax = 1.0;
   double aux = Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
	                                    (diameterYZ/2.0)) );
   //aux=1.0;

   uc.Set(*it,aux-1.0-_vel);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }

  // inflow condition V
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowVParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( X.Get(*it)*X.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   double Vmax = 2.0;
   double aux = Vmax*( 1.0-radius*radius/((diameterXZ/2.0)*
	                                    (diameterXZ/2.0)) );

   uc.Set(*it,0.0);
   vc.Set(*it,aux-_vel);
   wc.Set(*it,0.0);
  }

  // inflow condition W
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowWParabolic\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double radius = sqrt( X.Get(*it)*X.Get(*it) + Y.Get(*it)*Y.Get(*it) );

   // Parabolic profile
   double Wmax = 2.0;
   double aux = Wmax*( 1.0-radius*radius/((diameterXY/2.0)*
	                                    (diameterXY/2.0)) );

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,aux-_vel);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowV\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,1.0);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,1.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowZeroU\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0-_vel);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowZeroV\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0-_vel);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInflowZeroW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0-_vel);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvU\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,-1.0-_vel);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvV\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,-1.0-_vel);
   wc.Set(*it,0.0);
  }
  else if( surfMesh.phyBounds.at(*it) == "\"wallInvW\"" )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,-1.0-_vel);
  }

  // moving boundary U with 1 node pressure
  else if( surfMesh.phyBounds.at(*it) == "\"wallNoSlipPressure\"" )
  {
   if( X.Get(*it) == X.Max() && 
	   Z.Get(*it) == Z.Min() &&
	   count < 1 )
   {
	idbcp.AddItem(*it);

	pc.Set(*it,0.0);
	count++;
   }
   else
   {
	idbcu.AddItem(*it);
	idbcv.AddItem(*it);
	idbcw.AddItem(*it);

	uc.Set(*it,0.0-1.0-_vel);
	vc.Set(*it,0.0);
	wc.Set(*it,0.0);
   }
   if( Z.Get(*it) == Z.Min() )
   {
	idbcc.AddItem(*it);

	cc.Set(*it,1.0);
   }
  }

  //if( surfMesh.phyBounds.at(*it) == "\"wallNoSlip\"" )
  else
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);
  
   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
 }
}

void Model3D::setWallBC()
{    
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  idbcu.AddItem(*it);
  idbcv.AddItem(*it);
  idbcw.AddItem(*it);

  double aux = 0.0;
  uc.Set(*it,aux);
  vc.Set(*it,aux);
  wc.Set(*it,aux);
 }
}

void Model3D::setWallBC(double _vel)
{    
 clearBC();

 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  idbcu.AddItem(*it);
  idbcv.AddItem(*it);
  idbcw.AddItem(*it);

  double aux = 0.0;
  uc.Set(*it,aux-_vel);
  vc.Set(*it,aux);
  wc.Set(*it,aux);
 }
}

void Model3D::setMovingWallBC()
{    
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
//--------------------------------------------------
//   if( X.Get(*it) == X.Min() &&
// 	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
// 	  Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
//-------------------------------------------------- 
  if( X.Get(*it)==X.Min() && 
    ( Z.Get(*it)*Z.Get(*it)+Y.Get(*it)*Y.Get(*it) < (0.5*0.5 - 0.01) ) )
  {
   idbcp.AddItem(*it);

   pc.Set(*it,0.0);
  }
//--------------------------------------------------
//   else if( X.Get(*it) == X.Max() &&
// 	       Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
// 		   Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
//-------------------------------------------------- 
  else if( X.Get(*it)==X.Max() && 
	     ( Z.Get(*it)*Z.Get(*it)+Y.Get(*it)*Y.Get(*it) < (0.5*0.5 - 0.01) ) )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   double aux = 0.0;
   uc.Set(*it,aux);
   vc.Set(*it,aux);
   wc.Set(*it,aux);
  }
  else
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,-1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
 }
}

void Model3D::setMovingWallBC(double _vel)
{    
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
//--------------------------------------------------
//   if( X.Get(*it) == X.Min() &&
// 	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
// 	  Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
//-------------------------------------------------- 
  if( X.Get(*it)==X.Min() && 
    ( Z.Get(*it)*Z.Get(*it)+Y.Get(*it)*Y.Get(*it) < (0.5*0.5 - 0.01) ) )
  {
   pc.Set(*it,0.0);
  }
//--------------------------------------------------
//   else if( X.Get(*it) == X.Max() &&
// 	       Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
// 		   Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
//-------------------------------------------------- 
  else if( X.Get(*it)==X.Max() && 
	     ( Z.Get(*it)*Z.Get(*it)+Y.Get(*it)*Y.Get(*it) < (0.5*0.5 - 0.01) ) )
  {
   uc.Set(*it,0.0-_vel);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else
  {
   uc.Set(*it,-1.0-_vel);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
 }
}

void Model3D::setMicroWallBC()
{    
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  if( X.Get(*it) == X.Min() &&
	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
	  Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  if( (Y.Get(*it) == Y.Min()) || (Y.Get(*it) == Y.Max()) || 
	  (Z.Get(*it) == Z.Min()) || (Z.Get(*it) == Z.Max()) )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,-1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  if( X.Get(*it) == X.Max() &&
	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() && 
	  Z.Get(*it) > Z.Min() && Z.Get(*it) < Z.Max() )
  {
   idbcp.AddItem(*it);

   pc.Set(*it,0.0);
  }
 }
}

void Model3D::setCircularWallBC()
{    
 double radius = 0.5;
 for (list<int>::iterator it=boundaryVert.begin(); 
                          it!=boundaryVert.end(); 
						  ++it)
 {
  if( X.Get(*it) == X.Min() &&
	 (Y.Get(*it)*Y.Get(*it)+Z.Get(*it)*Z.Get(*it) < 
	  (radius*radius - 0.001) ) )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,1.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  else if( X.Get(*it) == X.Max() &&
	      (Y.Get(*it)*Y.Get(*it)+Z.Get(*it)*Z.Get(*it) < 
		   (radius*radius - 0.001) ) )
  {
   idbcp.AddItem(*it);

   pc.Set(*it,0.0);
  }
  else
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
 }
}

void Model3D::setWallCouetteBC()
{    
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  if( Z.Get(*it) == Z.Min() &&
	  X.Get(*it) > X.Min() && X.Get(*it) < X.Max() && 
	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,1.0);
  }
  if( Z.Get(*it) > Z.Min() && 
    ( X.Get(*it) == X.Min() || X.Get(*it) == X.Max() || 
      Y.Get(*it) == Y.Min() || Y.Get(*it) == Y.Max() ) )
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   uc.Set(*it,0.0);
   vc.Set(*it,0.0);
   wc.Set(*it,0.0);
  }
  if( Z.Get(*it) == Z.Max() &&
	  X.Get(*it) > X.Min() && X.Get(*it) < X.Max() && 
	  Y.Get(*it) > Y.Min() && Y.Get(*it) < Y.Max() )
  {
   idbcp.AddItem(*it);

   pc.Set(*it,0.0);
  }
 }
}

void Model3D::set2AxiBubblesBC()
{
 double aux = 0;
 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  // condicao de velocidade
  if( (Y.Get(*it)==Y.Min()) || (Y.Get(*it)==Y.Max()) )  
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   aux = X.Get(*it);
   uc.Set(*it,aux);
   aux = (-1.0)*Y.Get(*it);
   vc.Set(*it,aux);
   aux = 0.0;
   wc.Set(*it,aux);
  }
  else if( (Z.Get(*it)==Z.Min()) || Z.Get(*it)==Z.Max() )  
  {
   idbcw.AddItem(*it);

   aux = 0.0;
   wc.Set(*it,aux);
  }
  // condicao de outflow
  else
  {
   idbcp.AddItem(*it);

   aux = 0.0;
   pc.Set(*it,aux);
  }
 }
}

void Model3D::set2BubblesBC()
{
 double aux;

 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  if( (Y.Get(*it)==Y.Min()) || (Y.Get(*it)==Y.Max()) )  
  {
   idbcu.AddItem(*it);
   idbcv.AddItem(*it);
   idbcw.AddItem(*it);

   aux = X.Get(*it);
   uc.Set(*it,aux);
   aux = (-1.0)*Y.Get(*it);
   vc.Set(*it,aux);
   aux = Z.Get(*it);
   wc.Set(*it,aux);
  }
  // condicao de outflow
  else
  {
   idbcp.AddItem(*it);

   aux = 0.0;
   pc.Set(*it,aux);
  }
 }
}

// FSBC - Free Surface Boundary Condition
void Model3D::setDiskFSBC()
{
 double aux;
 heaviside.Dim(numVerts);
 rMax = Y.Max();

#if NUMGLEU == 5
 double numBCPoints = numVerts;
#else
 double numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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

   // heaviside = 0.5 -> noh da superficie
   heaviside.Set(i,0.5); // para funcionamento do ALE

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   pc.Set(i,aux);
   outflow.Set(i,aux);
  }
 }
}

void Model3D::setDiskCFSBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) == Z.Min() )
  {
   idbcc.AddItem(i);

   double aux = 1.0;
   cc.Set(i,aux);
  }
 }
}

void Model3D::setAdimenDisk()
{
 double aux;
 double Red = 100;
 double factorz = 1.0/(Z.Max()-Z.Min());
 double xrMax = X.Max();
 double yrMax = Y.Max();

 for( int i=0;i<numVerts;i++ )
 {
  aux = (X.Get(i)/xrMax)*Red;
  X.Set(i,aux);
  aux = (Y.Get(i)/yrMax)*Red;
  Y.Set(i,aux);
  //aux = Z.Get(i)*factorz*4;
  //aux = Z.Get(i)*factorz*6;
  aux = Z.Get(i)*factorz*10;
  Z.Set(i,aux);
 }
}

void Model3D::setAdimenDiskCouette()
{
 double aux;
 double Red = 100;
 double factorz = 1.0/(Z.Max()-Z.Min());

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
 double aux;
 double Red = 1;
 double factorz = 1.0/(Z.Max()-Z.Min());

 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i) != Z.Min() &&
	  (X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < (rMax*rMax - 0.1) ) )
  {
   aux = Z.Get(i)*( 1.0 + (X.Get(i)/Red)*factorz*0.3 );
   Z.Set(i,aux);
  }
 }
}

void Model3D::setPerturbSurf2()
{
 double aux;
 double Red = 1;
 double factorz = 1.0/(Z.Max()-Z.Min());
 double xMid = X.Min()+(X.Max()-X.Min())/2.0;

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
 double aux;
 double Red = 1;
 double factorz = 1.0/(Z.Max()-Z.Min());
 double rMid = rMax/10.0;

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

/*
 * This method perturbs the radius of a sphere to became an 
 * axisymmetric ellipsoid where x=y and differs from z.
 *
 * _factor = perturbed z diameter
 * 
 * dx*dy*dz = 1 (volume constant)
 * dx = dy, 
 * dz = 1.01
 * dx = sqrt(1/dz)
 *
 * */
void Model3D::setSphereToEllipsoid(double _factor)
{
 double aux;
 double dz = _factor;
 double dx = sqrt(1.0/dz);

 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  aux = surfMesh.X.Get(i)*dx;
  surfMesh.X.Set(i,aux);

  aux = surfMesh.Y.Get(i)*dx;
  surfMesh.Y.Set(i,aux);

  aux = surfMesh.Z.Get(i)*dz;
  surfMesh.Z.Set(i,aux);
 }

 for( int i=0;i<numVerts;i++ )
 {
  aux = X.Get(i)*dx;
  X.Set(i,aux);

  aux = Y.Get(i)*dx;
  Y.Set(i,aux);

  aux = Z.Get(i)*dz;
  Z.Set(i,aux);
 }
}

/* 
 * This method enlarges the sphere by a factor (_factor) * radius
 * keeping the same shape
 *
 * dx*dy*dz = 1 (volume constant)
 * dx = dy = dz 
 * Ex. _factor = 1 --> V = (4/3)*PI*0.5^3
 *     _factor = 2 --> V = (4/3)*PI*1.0^3
 *     _factor = 4 --> V = (4/3)*PI*2.0^3
 *
 * */
void Model3D::setBiggerSphere(double _factor)
{
 double aux;
 double ds = _factor;

 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  aux = surfMesh.X.Get(i)*ds;
  surfMesh.X.Set(i,aux);

  aux = surfMesh.Y.Get(i)*ds;
  surfMesh.Y.Set(i,aux);

  aux = surfMesh.Z.Get(i)*ds;
  surfMesh.Z.Set(i,aux);
 }

 for( int i=0;i<numVerts;i++ )
 {
  aux = X.Get(i)*ds;
  X.Set(i,aux);

  aux = Y.Get(i)*ds;
  Y.Set(i,aux);

  aux = Z.Get(i)*ds;
  Z.Set(i,aux);
 }
}

void Model3D::setMiniElement()
{
 numElems = numElems;
 numVerts = numVerts;
 numNodes = numVerts + numElems;

 clearBC();
 reAllocStruct();
 checkTetrahedronOrientation();

 setCentroid();
}

void Model3D::setCentroid()
{
 for( int i=0;i<numElems;i++ )
 {
  int v1 = (int) IEN.Get(i,0);
  int v2 = (int) IEN.Get(i,1);
  int v3 = (int) IEN.Get(i,2);
  int v4 = (int) IEN.Get(i,3);

  int vAdd = numVerts+i;

  int pos = IEN.DimJ()-1;
  IEN.Set(i,pos,vAdd);
  double centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
  double centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
  double centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;
  X.Set(vAdd,centroidX);
  Y.Set(vAdd,centroidY);
  Z.Set(vAdd,centroidZ);
 }
}

void Model3D::centroidPositionCorrection()
{
 for( int i=0;i<numElems;i++ )
 {
  int v1 = (int)IEN.Get(i,0);
  int v2 = (int)IEN.Get(i,1);
  int v3 = (int)IEN.Get(i,2);
  int v4 = (int)IEN.Get(i,3);
  int v5 = (int)IEN.Get(i,4);

  double centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
  double centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
  double centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;

  X.Set(v5,centroidX);
  Y.Set(v5,centroidY);
  Z.Set(v5,centroidZ);
 }
}

void Model3D::edgeMidPointPositionCorrection()
{
 double xMid,yMid,zMid;

 for( int i=0;i<numElems;i++ )
 {
  int v1  = (int)IEN.Get(i,0);
  int v2  = (int)IEN.Get(i,1);
  int v3  = (int)IEN.Get(i,2);
  int v4  = (int)IEN.Get(i,3);
  int v5  = (int)IEN.Get(i,4);
  int v6  = (int)IEN.Get(i,5);
  int v7  = (int)IEN.Get(i,6);
  int v8  = (int)IEN.Get(i,7);
  int v9  = (int)IEN.Get(i,8);
  int v10 = (int)IEN.Get(i,9);

  // v5 position correction
  xMid = ( X.Get(v1)+X.Get(v2) )/2.0;
  yMid = ( Y.Get(v1)+Y.Get(v2) )/2.0;
  zMid = ( Z.Get(v1)+Z.Get(v2) )/2.0;
  X.Set(v5,xMid);
  Y.Set(v5,yMid);
  Z.Set(v5,zMid);

  // v6 position correction
  xMid = ( X.Get(v1)+X.Get(v3) )/2.0;
  yMid = ( Y.Get(v1)+Y.Get(v3) )/2.0;
  zMid = ( Z.Get(v1)+Z.Get(v3) )/2.0;
  X.Set(v6,xMid);
  Y.Set(v6,yMid);
  Z.Set(v6,zMid);

  // v7 position correction
  xMid = ( X.Get(v1)+X.Get(v4) )/2.0;
  yMid = ( Y.Get(v1)+Y.Get(v4) )/2.0;
  zMid = ( Z.Get(v1)+Z.Get(v4) )/2.0;
  X.Set(v7,xMid);
  Y.Set(v7,yMid);
  Z.Set(v7,zMid);

  // v8 position correction
  xMid = ( X.Get(v2)+X.Get(v4) )/2.0;
  yMid = ( Y.Get(v2)+Y.Get(v4) )/2.0;
  zMid = ( Z.Get(v2)+Z.Get(v4) )/2.0;
  X.Set(v8,xMid);
  Y.Set(v8,yMid);
  Z.Set(v8,zMid);

  // v9 position correction
  xMid = ( X.Get(v3)+X.Get(v4) )/2.0;
  yMid = ( Y.Get(v3)+Y.Get(v4) )/2.0;
  zMid = ( Z.Get(v3)+Z.Get(v4) )/2.0;
  X.Set(v9,xMid);
  Y.Set(v9,yMid);
  Z.Set(v9,zMid);

  // v10 position correction
  xMid = ( X.Get(v2)+X.Get(v3) )/2.0;
  yMid = ( Y.Get(v2)+Y.Get(v3) )/2.0;
  zMid = ( Z.Get(v2)+Z.Get(v3) )/2.0;
  X.Set(v10,xMid);
  Y.Set(v10,yMid);
  Z.Set(v10,zMid);
 }
}

void Model3D::checkTetrahedronOrientation()
{
 V.Dim(numElems);

 for( int i=0;i<numElems;i++ )
 {
  double volume = getVolume(i);
  V.Set(i,volume);

  if( fabs(volume)<1.0E-14)
   cerr << " -------- " << "tetrahedon singular (" << i << ")" 
        << " -------- " << endl;

  if( volume<0.0 )
  {
   int v2 = (int)IEN.Get(i,1);
   int v3 = (int)IEN.Get(i,2);

   IEN.Set(i,1,v3);
   IEN.Set(i,2,v2);
  }
 }
}

/*
 * mapEdge.Set(edge,0,numVerts+edge); // numero da aresta
 * mapEdge.Set(edge,1,xMid ); // coordenada X do centro da aresta
 * mapEdge.Set(edge,2,yMid ); // coordenada Y do centro da aresta
 * mapEdge.Set(edge,3,zMid ); // coordenada Y do centro da aresta
 * mapEdge.Set(edge,4,faces[i].p1 ); // 1o noh
 * mapEdge.Set(edge,5,faces[i].p2 ); // 2o noh
 *
 * */
void Model3D::setMapping()
{
 setNeighbour();

 int numFaces = 4; // teraedro tem 6 arestas
 clVector faceaux(3);
 IFACE3D *faces = NULL;
 faces = new IFACE3D[numFaces*numElems];

 int numEdges = 6; // teraedro tem 6 arestas
 clVector edgeaux(2);
 IFACE2D *edges = NULL;
 edges = new IFACE2D[numEdges*numElems];

 checkTetrahedronOrientation();

 for( int mele=0;mele<numElems;mele++ )
 {
  int v1 = (int) IEN.Get(mele,0);
  int v2 = (int) IEN.Get(mele,1);
  int v3 = (int) IEN.Get(mele,2);
  int v4 = (int) IEN.Get(mele,3);

  // 1st. face 0-1-2
  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Set(2,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFaces*mele+0].p1 = (int) faceaux.Get(0);
  faces[numFaces*mele+0].p2 = (int) faceaux.Get(1);
  faces[numFaces*mele+0].p3 = (int) faceaux.Get(2);
  faces[numFaces*mele+0].p4 = mele;
  faces[numFaces*mele+0].p5 = 3; // ID of v4

  // 2nd. face 0-1-3
  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Set(2,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFaces*mele+1].p1 = (int) faceaux.Get(0);
  faces[numFaces*mele+1].p2 = (int) faceaux.Get(1);
  faces[numFaces*mele+1].p3 = (int) faceaux.Get(2);
  faces[numFaces*mele+1].p4 = mele;
  faces[numFaces*mele+1].p5 = 2; // ID of v3

  // 3rd. face 0-2-3
  faceaux.Set(0,v1);
  faceaux.Set(1,v3);
  faceaux.Set(2,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFaces*mele+2].p1 = (int) faceaux.Get(0);
  faces[numFaces*mele+2].p2 = (int) faceaux.Get(1);
  faces[numFaces*mele+2].p3 = (int) faceaux.Get(2);
  faces[numFaces*mele+2].p4 = mele;
  faces[numFaces*mele+2].p5 = 1; // ID of v2

  // 4th. face 1-2-3
  faceaux.Set(0,v2);
  faceaux.Set(1,v3);
  faceaux.Set(2,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFaces*mele+3].p1 = (int) faceaux.Get(0);
  faces[numFaces*mele+3].p2 = (int) faceaux.Get(1);
  faces[numFaces*mele+3].p3 = (int) faceaux.Get(2);
  faces[numFaces*mele+3].p4 = mele;
  faces[numFaces*mele+3].p5 = 0; // ID of v1

  // 1st. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+0].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+0].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+0].p3 = mele;

  // 2nd. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+1].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+1].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+1].p3 = mele;

  // 3rd. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+2].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+2].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+2].p3 = mele;

  // 4th. edge
  faceaux.Set(0,v2);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+3].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+3].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+3].p3 = mele;

  // 5th. edge
  faceaux.Set(0,v2);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+4].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+4].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+4].p3 = mele;

  // 6th. edge
  faceaux.Set(0,v3);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  edges[numEdges*mele+5].p1 = (int) edgeaux.Get(0);
  edges[numEdges*mele+5].p2 = (int) edgeaux.Get(1);
  edges[numEdges*mele+5].p3 = mele;
 }

 // ordena uma estrutura (matriz) em ordem crescente na linha e coluna
 // as faces continuam repetidas neste ponto.
 qsort(faces,numFaces*numElems,sizeof(IFACE3D),IFACE3DCompare);
 qsort(edges,numEdges*numEdges,sizeof(IFACE2D),IFACE2DCompare);

//--------------------------------------------------
//  for( int i=0;i<numFaces*numElems;i++ )
//   cout << faces[i].p1 << " "
//        << faces[i].p2 << " "
//        << faces[i].p3 << endl;
//-------------------------------------------------- 

 /*        - nome: oFace
           - definicao: matrix com mapeamento de elementos opostos ao 
		                vertice em questao. Se encontrar o valor (-1) quer
						dizer que o vertice nao apresenta elemento oposto,
						ou seja, eh um elemento de fronteira

               +---+---+---+
               | a | a | a |  a = identificacao dos vertices do elemento
     ---   +---+---+---+---+      seguindo a ordem local 0-1-2
      |    | b | e | e | e |  
      |    +---+---+---+---+  b = idendificacao dos elementos
	  |    .   .   .   .   .  
	  |	   .   .   .   .   .  e = identificacao do elemento oposto ao 
	  c    .   .   .   .   .      vertice
      |    +---+---+---+---+
      |    | b | e | e | e |
      |    +---+---+---+---+  c = numero total de elementos
      |    | b | e | e | e |  d = numero de total de vertices do elemento 
     ---   +---+---+---+---+

	       |________ d ________|
		   |                   |
 */

 // vertices belonging to the boundary
 boundaryVert.clear();
 oFace.Dim(numElems,NUMGLE);
 for( int i=0;i<numElems;i++ )
  for( int j=0;j<NUMGLE;j++ )
   oFace.Set(i,j,-1);

 /* p1 = v1
  * p2 = v2
  * p3 = v3
  * p4 = elem
  * p5 = v postition (0,1,2 or 3)
  * */
 int i = 0;
 while( i < numFaces*numElems-1 )
 {
  if( faces[i].p1 == faces[i+1].p1 && 
	  faces[i].p2 == faces[i+1].p2 &&
	  faces[i].p3 == faces[i+1].p3 )
  {
   oFace.Set(faces[i].p4,faces[i].p5,faces[i+1].p4);
   oFace.Set(faces[i+1].p4,faces[i+1].p5,faces[i].p4);
   i += 2;
  }
  else
  {
   boundaryVert.push_back(faces[i].p1);
   boundaryVert.push_back(faces[i].p2);
   boundaryVert.push_back(faces[i].p3);
   i += 1;
  }
 }
 boundaryVert.sort();
 boundaryVert.unique();

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

 int edge=0;
 minEdge = 1000000;

 /*
  *  neighbourEdge is a vector list related to the elements ID that
  *  shares the same edge. In 2D, one edge is shared only by 2 elements,
  *  but in 3D there is no fixed number of neighbours elements.
  * */
 neighbourEdge.clear();
 neighbourEdge.resize (numEdges*numElems);
 mapEdge.Dim(numEdges*numElems,6);

 // numeracao de arestas a partir de numVerts e associacao das arestas
 // aos elementos para inclusao na IEN OBS.: the for loop is from 0 to
 // numEdges*numElems-1. This can be done because the last line of edges
 // (edges[numEdges*numElems-1]) is always a repetion of the previous one
 // (edges[numEdges*numElems-2]).
 for( int i=0;i<numEdges*numElems-1;i++ )
 {
  double x1=X.Get(edges[i].p1);
  double y1=Y.Get(edges[i].p1);
  double z1=Z.Get(edges[i].p1);
  double x2=X.Get(edges[i].p2);
  double y2=Y.Get(edges[i].p2);
  double z2=Z.Get(edges[i].p2);

  double xMid = (x1+x2)*0.5;
  double yMid = (y1+y2)*0.5;
  double zMid = (z1+z2)*0.5;
  double length = vectorLength(x1-x2,y1-y2,z1-z2);

  // checking the minimum edge length
  if( length < minEdge )
   minEdge = length;

  mapEdge.Set(edge,0,numVerts+edge); // numero da aresta
  mapEdge.Set(edge,1,xMid ); // coordenada X do centro da aresta
  mapEdge.Set(edge,2,yMid ); // coordenada Y do centro da aresta
  mapEdge.Set(edge,3,zMid ); // coordenada Y do centro da aresta
  mapEdge.Set(edge,4,edges[i].p1 ); // 1o noh
  mapEdge.Set(edge,5,edges[i].p2 ); // 2o noh

  neighbourEdge.at(edge).push_back(edges[i].p3);

  // edges[numEdges*numElems]
  while( (i<numEdges*numElems-1 ) &&
         (edges[i].p1 == edges[i+1].p1) &&
         (edges[i].p2 == edges[i+1].p2) )
		 
  {
   neighbourEdge.at(edge).push_back(edges[i+1].p3);
   i++; // pula 1 linha
  }
  edge++; // numero total de arestas
 }
 clMatrix resizedMapEdge(edge,6);
 mapEdge.CopyTo(0,0,resizedMapEdge);
 mapEdge = resizedMapEdge;
 neighbourEdge.resize (edge); // trim vector para numero real de itens

//--------------------------------------------------
//  for( int i=0;i<edge;i++ )
//  {
//   cout << "---------" << i << "------------" << endl;
//   std::ostream_iterator< int > output( cout, " " );
//   std::copy( neighbourEdge.at(i).begin(), 
// 	         neighbourEdge.at(i).end(), output );
//   cout << endl;
//  }
//-------------------------------------------------- 

 delete[] faces;
 delete[] edges;
}

void Model3D::setQuadElement()
{
 // atualizado vetores com numero total de nos
 numElems = numElems;
 numVerts = numVerts;
 numNodes = numVerts+mapEdge.DimI(); // atualizando numNodes

 clearBC();
 reAllocStruct();
 checkTetrahedronOrientation();


 /*
  *  Tetrahedron divided by 4 facets:
  *
  *      facet 1        facet 2        facet 3        facet 4
  *
  *         v3             v3             v4             v3       
  *         o              o              o              o       
  *        / \            / \            / \            / \
  *    v6 o   o v8    v6 o   o v10   v7 o   o v9    v8 o   o v10   
  *      /     \        /     \        /     \        /     \     
  *     o - o - o      o - o - o      o - o - o      o - o - o   
  *    v1   v5   v2   v1   v7   v4   v1   v5   v2   v2   v9   v4  
  *
  * */

 for( int i=0;i<mapEdge.DimI();i++ )
 {
  int edge   = mapEdge.Get(i,0);
  double xc    = mapEdge.Get(i,1);
  double yc    = mapEdge.Get(i,2);
  double zc    = mapEdge.Get(i,3);
  int vEdge1 = mapEdge.Get(i,4); 
  int vEdge2 = mapEdge.Get(i,5);

  list<int> plist = neighbourEdge.at (i);
  for(list<int>::iterator elem=plist.begin();elem!=plist.end();++elem )
  {
   // add at coords to X and Y vectors
   X.Set(edge,xc); 
   Y.Set(edge,yc);
   Z.Set(edge,zc);

   int v1 = IEN.Get(*elem,0);
   int v2 = IEN.Get(*elem,1);
   int v3 = IEN.Get(*elem,2);
   int v4 = IEN.Get(*elem,3);

   // vertex v5 between v1 e v2
   if( (vEdge1 == v1 && vEdge2 == v2) || 
	   (vEdge2 == v1 && vEdge1 == v2) )
   {
	int v5 = edge;
	IEN.Set(*elem,4,v5);
   }

   // vertex v6 between v2 e v3
   if( (vEdge1 == v2 && vEdge2 == v3) || 
	   (vEdge2 == v2 && vEdge1 == v3) )
   {
	int v6 = edge;
	IEN.Set(*elem,5,v6);
   }

   // vertex v7 between v3 e v1
   if( (vEdge1 == v3 && vEdge2 == v1) || 
	   (vEdge2 == v3 && vEdge1 == v1) )
   {
	int v7 = edge;
	IEN.Set(*elem,6,v7);
   }

   // vertex v8 between v1 e v4
   if( (vEdge1 == v1 && vEdge2 == v4) || 
	   (vEdge2 == v1 && vEdge1 == v4) )
   {
	int v8 = edge;
	IEN.Set(*elem,7,v8);
   }

   // vertex v9 between v2 e v4
   if( (vEdge1 == v2 && vEdge2 == v4) || 
	   (vEdge2 == v2 && vEdge1 == v4) )
   {
	int v9 = edge;
	IEN.Set(*elem,8,v9);
   }

   // vertex v10 between v3 e v4
   if( (vEdge1 == v3 && vEdge2 == v4) || 
	   (vEdge2 == v3 && vEdge1 == v4) )
   {
	int v10 = edge;
	IEN.Set(*elem,9,v10);
   }
  }
 }
}

void Model3D::setNeighbour()
{
 neighbourElem.resize (0);
 neighbourElem.resize (numVerts);

 for( int i=0;i<numElems;i++ )
  for( int j= 0;j<NUMGLEU;j++ )
   neighbourElem.at( (int)IEN.Get(i,j) ).push_back(i);
}

/* create neighbourVert vector list, which is the mapping of neighbour
 * vertices of each 3D mesh vertex (numVerts size). In this way, the
 * vector is numVerts long and each entry of this vector is a list that
 * can change the number of elements because the mesh is unstructured.
 *
 * input: neighbourElem
 * output: neighbourVert
 * */
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

/* finds the surface and non-surface nodes and the corresponding x,y,z
 * coordinates (only for surface nodes).
 * */
void Model3D::setSurface()
{
 list<int> plist;
 list<int>::iterator vert;

 // procurando vertices da bolha
 clVector surfaceAux = surfMesh.Marker==0.5;
 surface = surfaceAux.Find();
 
 // {x,y,z}Surface representam os valores de {X,Y,Z} dos nos da interface
 xSurface.Dim( surface.Dim() );
 ySurface.Dim( surface.Dim() );
 zSurface.Dim( surface.Dim() );
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  xSurface.Set(i,surfMesh.X.Get( surfaceNode ));
  ySurface.Set(i,surfMesh.Y.Get( surfaceNode ));
  zSurface.Set(i,surfMesh.Z.Get( surfaceNode ));
 }
} // fecha metodo setSurface

/* cria matriz IEN para os elementos da superficie que no caso 3D sao
 * triagulos. Este metodo utiliza neighbourPoint 
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

  int listSize = neighbourPoint.at(surfaceNode).size();
  list<int> plist = neighbourPoint.at(surfaceNode);
  list<int>::iterator mele=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   int v1 = *mele;++mele;
   int v2 = *mele;

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
 * este metodo organiza uma estrutura do tipo IEN ordenando a partir de
 * um numero todos os outros vertices sem deixar 'buracos' na estrutura.
 * Isso quer dizer que podemos recriar uma malha comecando a numeracao
 * dos vertices a partir de um numero qualquer.
 * Para isso eh necessario criar um vetor de listas do tipo setNeighbour()
 * e entao fazer um mapeamento de cada vertice para o numero
 * correspondente.
 * !!!metodo ainda nao testado 100%!!!
 * */
SurfaceMesh Model3D::arrangeMesh(SurfaceMesh _tetmesh,int _nVerts,int _begin)
{
 int nElems = _tetmesh.numElems;
 int end = _tetmesh.IEN.Max();

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
   test.at( (int) _tetmesh.IEN.Get(i,j) ).push_back(i);

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
	int v1 = _tetmesh.IEN.Get(*it,0);
	int v2 = _tetmesh.IEN.Get(*it,1);
	int v3 = _tetmesh.IEN.Get(*it,2);
	if (v1 == i )
	{
	 meshReturn.IEN.Set(*it,0,count+_begin);
	 meshReturn.X.Set(count,X.Get(v1));
	 meshReturn.Y.Set(count,Y.Get(v1));
	 meshReturn.Z.Set(count,Z.Get(v1));
	 meshReturn.Marker.Set(count,heaviside.Get(v1));
	}
	else if (v2 == i )
	{
	 meshReturn.IEN.Set(*it,1,count+_begin);
	 meshReturn.X.Set(count,X.Get(v2));
	 meshReturn.Y.Set(count,Y.Get(v2));
	 meshReturn.Z.Set(count,Z.Get(v2));
	 meshReturn.Marker.Set(count,heaviside.Get(v2));
	}
	else 
	{
	 meshReturn.IEN.Set(*it,2,count+_begin);
	 meshReturn.X.Set(count,X.Get(v3));
	 meshReturn.Y.Set(count,Y.Get(v3));
	 meshReturn.Z.Set(count,Z.Get(v3));
	 meshReturn.Marker.Set(count,heaviside.Get(v3));
	}
   }
   count++;
  }
 }
 return meshReturn;
}

// este metodo cria duas listas com os vertices do convex hull (boundaryVert)
// e todos os outros menos os vertices do convex hull (inVert). Este
// metodo eh especialmente interessante para aplicacao das condicoes de
// contorno para geometrias complexas.
void Model3D::setInOutVert()
{
 inVert.resize (0);

 // retira de inVert todos os vertices presents em boundaryVert.
 list<int>::iterator it;
 for( int vert=0;vert<numVerts;vert++ )
  inVert.push_back(vert);
 for( it=boundaryVert.begin();it!=boundaryVert.end();++it )
  inVert.remove(*it);

//--------------------------------------------------
//  cout << "boundaryVert contains:";
//  for (it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
//   cout << " " << *it;
//  cout << endl;
//-------------------------------------------------- 
}

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

  double hsum = heaviside.Get(v1)+heaviside.Get(v2)+
              heaviside.Get(v3)+heaviside.Get(v4);

  if( hsum > 1.5 )
   inElem.push_back(i);
  else
   outElem.push_back(i);
 }
}

// aplica configuracoes referentes a superficie da modelagem 2 fases.
void Model3D::setSurfaceConfig()
{
 setVertNeighbour(); // neighbourVert (3D mesh)
 setInOutVert(); // inVert e boundaryVert
 setInOutElem(); // inElem e outElem

 // update surface, edge matrix, surface neigh elems and points
 restoreMappingArrays();

 //setSurfaceTri(); // triang superficie - interfaceMesh

 setInterfaceDistance();
 setNormalAndKappa();
 setKappaSurface();

 setSurfaceVolume();
 setSurfaceArea();
 triMeshStats();
 tetMeshStats();
}

/* restores all the matrix and vectors mappings when a point is inserted
 * or deleted, and an edge is contracted.
 * This method encapsulate the updating routines.
 *
 * */
void Model3D::restoreMappingArrays()
{
 // update surface
 setSurface();

 // updating edge matrix
 setMapEdgeTri();

 // updating surface neighbour elems
 setNeighbourSurfaceElem();

 // updating surface neighbour points
 setNeighbourSurfacePoint();
}


bool Model3D::testFace(int v1, int v2, int v3, int v4)
{
 double V,Ax1,Ax2,Ay1,Ay2,Az1,Az2;
 double prodEsc;

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

double Model3D::getVolume(int _v1,int _v2,int _v3,int _v4)
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

double Model3D::getVolume(int _elem)
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

void Model3D::clearBC()
{
 uc.Dim(numNodes);
 vc.Dim(numNodes);
 wc.Dim(numNodes);
 pc.Dim(numVerts);
 cc.Dim(numVerts);
 idbcu.Dim(0);
 idbcv.Dim(0);
 idbcw.Dim(0);
 idbcp.Dim(0);
 idbcc.Dim(0);
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
 IEN.Dim(numElems,NUMGLEU); // 4 nos por elemento + 1 centroide
 IEN.CopyFrom(0,0,IENaux);
}

void Model3D::moveXPoints(clVector &_vec,double _dt)
{
 //X = X + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  double aux = X.Get(i)+(_vec.Get(i)*_dt);
  X.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  double aux = surfMesh.X.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.X.Set(i,aux);
 }
}

void Model3D::moveYPoints(clVector &_vec,double _dt)
{
 //Y = Y + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  double aux = Y.Get(i)+(_vec.Get(i)*_dt);
  Y.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  double aux = surfMesh.Y.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.Y.Set(i,aux);
 }
}

void Model3D::moveZPoints(clVector &_vec,double _dt)
{
 //Z = Z + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  double aux = Z.Get(i)+(_vec.Get(i)*_dt);
  Z.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  double aux = surfMesh.Z.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.Z.Set(i,aux);
 }
}


SurfaceMesh* Model3D::getSurfMesh(){ return &surfMesh; }
SurfaceMesh* Model3D::getInterfaceMesh(){ return &interfaceMesh; }
Mesh3D* Model3D::getMesh3d(){ return &mesh3d; }
clVector* Model3D::getX(){ return &X; }
double Model3D::getMaxX(){ return X.Max(); }
double Model3D::getMinX(){ return X.Min(); }
void Model3D::setX(clVector _X){ X = _X; }
clVector* Model3D::getY(){ return &Y; }
double Model3D::getMinY(){ return Y.Min(); }
double Model3D::getMaxY(){ return Y.Max(); }
void Model3D::setY(clVector _Y){ Y = _Y; }
double Model3D::getMaxZ(){ return Z.Max(); }
double Model3D::getMinZ(){ return Z.Min(); }
clVector* Model3D::getZ(){ return &Z; }
void Model3D::setZ(clVector _Z){ Z = _Z; }
clVector* Model3D::getUC(){ return &uc; }
clVector* Model3D::getVC(){ return &vc; }
clVector* Model3D::getWC(){ return &wc; }
clVector* Model3D::getPC(){ return &pc; }
clVector* Model3D::getCC(){ return &cc; }
clVector* Model3D::getHeaviside(){ return &heaviside; }
clVector* Model3D::getOutflow(){ return &outflow; }
clVector* Model3D::getIdbcu(){ return &idbcu; }
clVector* Model3D::getIdbcv(){ return &idbcv; }
clVector* Model3D::getIdbcw(){ return &idbcw; }
clVector* Model3D::getIdbcp(){ return &idbcp; }
clVector* Model3D::getIdbcc(){ return &idbcc; }
clMatrix* Model3D::getIEN(){ return &IEN; }
clMatrix* Model3D::getMapEdge(){ return &mapEdge; }
clMatrix* Model3D::getMapEdgeTri(){ return &mapEdgeTri; }
clVector* Model3D::getInterfaceDistance(){ return &interfaceDistance; }
clVector* Model3D::getElemIdRegion(){ return &elemIdRegion; }
clVector* Model3D::getVertIdRegion(){ return &vertIdRegion; }
clDMatrix* Model3D::getCurvature(){ return &curvature; }
int Model3D::getNumVerts(){ return numVerts; }
int Model3D::getNumNodes(){ return numNodes; }
int Model3D::getNumElems(){ return numElems; }
clMatrix* Model3D::getOFace(){ return &oFace; }
double Model3D::getXCenter(){ return xCenter; }
double Model3D::getYCenter(){ return yCenter; }
double Model3D::getZCenter(){ return zCenter; }
clVector* Model3D::getSurface(){ return &surface; }
vector< list<int> >* Model3D::getNeighbourElem(){return &neighbourElem;}
vector< list<int> >* Model3D::getNeighbourVert(){return &neighbourVert;}
vector< list<int> >* Model3D::getNeighbourFace(){return &neighbourFace;}
vector< list<int> >* Model3D::getNeighbourPoint(){return &neighbourPoint;}
vector<int>* Model3D::getPbcIndicesLeft(){return &pbcIndicesLeft;} // PBC
vector<int>* Model3D::getPbcIndicesRight(){return &pbcIndicesRight;} // PBC
vector< vector<int> >* Model3D::getPbcIndicesMaster(){return &pbcIndicesMaster;} // PBC
vector< vector<int> >* Model3D::getPbcIndicesSlave(){return &pbcIndicesSlave;} // PBC
vector< list<int> >* Model3D::getFaceIEN(){return &faceIEN;}
list<int>* Model3D::getInVert(){return &inVert;}
list<int>* Model3D::getBoundaryVert(){return &boundaryVert;}
list<int>* Model3D::getInElem(){return &inElem;}
list<int>* Model3D::getOutElem(){return &outElem;}
double Model3D::getMinEdge(){return minEdge;}
vector<double> Model3D::getTriEdge(){return triEdge;}
void Model3D::setTetVol(vector< double > _tetVol){tetVol= _tetVol;}
vector<double> Model3D::getTetVol(){return tetVol;}
clVector* Model3D::getEdgeSize(){ return &edgeSize; }
void Model3D::setEdgeSize(clVector _edgeSize){ edgeSize = _edgeSize; }
vector<double> Model3D::getAverageTriLength(){ return averageTriLength; }
vector<double> Model3D::getAverageTriArea(){ return averageTriArea; }
vector<double> Model3D::getAverageTetVolume(){ return averageTetVolume; }

vector<double> Model3D::getInitSurfaceArea(){return initSurfaceArea;}
vector<double> Model3D::getInitSurfaceVolume(){return initSurfaceVolume;}
vector<double> Model3D::getDArea(){return dArea;}
vector<double> Model3D::getErrorArea(){return errorArea;}
vector<double> Model3D::getSurfaceArea(){return surfaceArea;}
vector<double> Model3D::getDVolume(){return dVolume;}
vector<double> Model3D::getErrorVolume(){return errorVolume;}
vector<double> Model3D::getSurfaceVolume(){return surfaceVolume;}
clVector* Model3D::getCloser(){ return &closer; }

// Mesh indexes:
vector<int> Model3D::getISP(){return isp;}
vector<int> Model3D::getISPC(){return ispc;}
vector<int> Model3D::getRSP(){return rsp;}
vector<int> Model3D::getRSPN(){return rspn;}
vector<int> Model3D::getRSPC(){return rspc;}
vector<int> Model3D::getFLIP(){return flip;}
vector<int> Model3D::getSPC(){return spc;}
vector<int> Model3D::getSPP(){return spp;}
vector<int> Model3D::getINTET(){return intet;}
vector<double> Model3D::getMinArea(){return minArea;}
vector<double> Model3D::getMaxArea(){return maxArea;}
vector<double> Model3D::getMinLength(){return minLength;}
vector<double> Model3D::getMaxLength(){return maxLength;}
vector<int> Model3D::getIdMinArea(){return idMinArea;}
vector<int> Model3D::getIdMaxArea(){return idMaxArea;}
vector<int> Model3D::getNumSurfElems(){return numSurfElems;}
vector<int> Model3D::getNumSurfVerts(){return numSurfVerts;}
vector<int> Model3D::getOPER(){return oper;}
vector<int> Model3D::getOPERSURF(){return opersurf;}
vector<int> Model3D::getIP(){return ip;}
vector<int> Model3D::getIPD(){return ipd;}
vector<int> Model3D::getRP(){return rp;}
vector<int> Model3D::getRPI(){return rpi;}
vector<int> Model3D::getRPD(){return rpd;}
vector<int> Model3D::getRPDist(){return rpdist;}
vector<int> Model3D::getRPH(){return rph;}
vector<int> Model3D::getRPV(){return rpv;}
vector<int> Model3D::getCSP(){return csp;}
vector<double> Model3D::getMinVolume(){return minVolume;}
vector<double> Model3D::getMaxVolume(){return maxVolume;}
vector<int> Model3D::getIdMinVolume(){return idMinVolume;}
vector<int> Model3D::getIdMaxVolume(){return idMaxVolume;}

//-------------------------------------------------- 
// Atribui o Model3D do argumento no corrente
//-------------------------------------------------- 
void Model3D::operator=(Model3D &_mRight) 
{
  // ints and floats
  numVerts = _mRight.numVerts;
  numNodes = _mRight.numNodes;
  numElems = _mRight.numElems;
  rMax = _mRight.rMax;
  xCenter = _mRight.xCenter;
  yCenter = _mRight.yCenter;
  zCenter = _mRight.zCenter;
  dVerts = _mRight.dVerts;                  
  minEdge = _mRight.minEdge;
  averageTriLength = _mRight.averageTriLength;
  averageTriArea = _mRight.averageTriArea;
  averageTetVolume = _mRight.averageTetVolume;
  isp = _mRight.isp;
  ispc = _mRight.ispc;
  rsp = _mRight.rsp;        
  rspn = _mRight.rspn;        
  rspc = _mRight.rspc;
  flip = _mRight.flip;
  spc = _mRight.spc;
  spp = _mRight.spp;
  intet = _mRight.intet;
  maxArea = _mRight.maxArea;
  minArea = _mRight.minArea;
  maxLength = _mRight.maxLength;
  minLength = _mRight.minLength;
  idMaxArea = _mRight.idMaxArea;
  numSurfElems = _mRight.numSurfElems;
  numSurfVerts = _mRight.numSurfVerts;
  idMinArea = _mRight.idMinArea;

  oper = _mRight.oper;                    
  opersurf = _mRight.opersurf;                    

  ip = _mRight.ip;                    
  ipd = _mRight.ipd;                    
  rp = _mRight.rp;              
  rpi = _mRight.rpi;                   
  rpd = _mRight.rpd;                   
  rpdist = _mRight.rpdist;                   
  rph = _mRight.rph;                   
  rpv = _mRight.rpv;                   
  csp = _mRight.csp;                   
  maxVolume = _mRight.maxVolume;
  minVolume = _mRight.minVolume;
  idMaxVolume = _mRight.idMaxVolume;
  idMinVolume = _mRight.idMinVolume;

  // clVector and clMatrix
  surface = _mRight.surface;
  uc = _mRight.uc;
  vc = _mRight.vc;
  wc = _mRight.wc;
  pc = _mRight.pc;
  cc = _mRight.cc;
  heaviside = _mRight.heaviside;
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
  oFace = _mRight.oFace;
  surfMesh = _mRight.surfMesh;
  interfaceMesh = _mRight.interfaceMesh;
  vertIdRegion = _mRight.vertIdRegion;
  elemIdRegion = _mRight.elemIdRegion;
  triEdge = _mRight.triEdge;
  tetVol = _mRight.tetVol;
  edgeSize = _mRight.edgeSize;

  // STL: list and vectors
  initSurfaceVolume = _mRight.initSurfaceVolume;
  surfaceVolume = _mRight.surfaceVolume;
  initSurfaceArea = _mRight.initSurfaceArea;
  surfaceArea = _mRight.surfaceArea;
  neighbourElem = _mRight.neighbourElem; 
  neighbourVert = _mRight.neighbourVert;
  neighbourFace = _mRight.neighbourFace;
  neighbourSurfaceElem = _mRight.neighbourSurfaceElem;
  neighbourPoint = _mRight.neighbourPoint;
  faceIEN = _mRight.faceIEN;
  boundaryVert = _mRight.boundaryVert;
  inVert = _mRight.inVert;
  outElem = _mRight.outElem;
  inElem = _mRight.inElem;
  pbcIndicesLeft = _mRight.pbcIndicesLeft; // PBC
  pbcIndicesRight = _mRight.pbcIndicesRight; // PBC
  pbcIndicesMaster = _mRight.pbcIndicesMaster; // PBC
  pbcIndicesSlave = _mRight.pbcIndicesSlave; // PBC
}

void Model3D::saveVTKSurface( const char* _dir,
                              const char* _filename, 
							  int _iter )
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

 vtkFile << "POINTS " << surfMesh.numVerts<< " double" << endl;
 for( int i=0;i<surfMesh.numVerts;i++ )
  vtkFile << surfMesh.X.Get(i) << " " 
          << surfMesh.Y.Get(i) << " " 
		  << surfMesh.Z.Get(i) << endl;

 vtkFile << endl;

 int numTri = 0;
 for( int i=0;i<surfMesh.numElems;i++ )
  if( surfMesh.elemIdRegion.Get(i) > 0 )
   numTri++;

 vtkFile << "CELLS " << numTri << " " << 4*numTri << endl;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  if( surfMesh.elemIdRegion.Get(i) > 0 )
   vtkFile << "3 " << surfMesh.IEN.Get(i,0) << " "
	               << surfMesh.IEN.Get(i,1) << " "
				   << surfMesh.IEN.Get(i,2) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numTri << endl;
 for( int i=0;i<numTri;i++ )
  vtkFile << "5 ";

 vtkFile << endl;

 vtkFile.close();
}

void Model3D::saveVTK( const char* _dir,const char* _filename, int _iter )
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
  vtkFile << X.Get(i) << " " 
          << Y.Get(i) << " " 
		  << Z.Get(i) << endl;

 vtkFile << endl;

 vtkFile << "CELLS " << numElems << " " << 4*numElems << endl;
 for( int i=0;i<numElems;i++ )
 {
   vtkFile << "4 " << IEN.Get(i,0) << " "
	               << IEN.Get(i,1) << " "
	               << IEN.Get(i,2) << " "
				   << IEN.Get(i,3) << endl;
 }
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
  vtkFile << "10 ";

 vtkFile << endl;

 vtkFile.close();
}

/* This method computes the curvature of a surface point using tangent 
 * vectors of triangles belonging to the star of such a point and
 * summing for all surface elements, then dividing to the sum of 1/3 the
 * total area.
 * Description:
 *  1_ using neighbourPoint we compute the unit vectors of each element
 *  inside the star of the surface node.
 *  2_ [normal vector] the cross product of both unit vectors 
 *  (P1 - surfaceNode and P2 - surfaceNode) gives a unit normal vector 
 *  of the respective triangle.
 *  3_ [curvature calculation] the normal unit vector of the opposite
 *  edge (lied on the surface element plane) is integrated in the
 *  half length of the specific edge, thus it is computed the curvature
 *  of such a plane.
 *  4_ summing all cross vectors and all curvature vectors with respect
 *  to the elements that contains the surface node, it is possible to
 *  calculate with precision the normal unit vector and curvature of
 *  each surface node.
 *
 *  NORMAL: OUTWARD direction
 *
 *  requirement: neighbourSurfaceElem,neighbourPoint
 *  input: node, list of neighbour elements of the node
 *  output: vector(4) = curvature,
 *                      x normal,
 *                      Y normal,
 *                      Z normal
 * */
clVector Model3D::getNormalAndKappa(int _node,list<int> _myList)
{
 double P0x = surfMesh.X.Get(_node);
 double P0y = surfMesh.Y.Get(_node);
 double P0z = surfMesh.Z.Get(_node);

 //int c1 = 0;
 double fx = 0;
 double fy = 0;
 double fz = 0;
 double sumArea = 0;
 double sumXCrossUnit = 0;
 double sumYCrossUnit = 0;
 double sumZCrossUnit = 0;

 int listSize = _myList.size();
 list<int>::iterator mele=_myList.begin();
 for( int i=0;i<listSize-1;i++ )
 {
  int v1 = *mele;++mele;
  int v2 = *mele;

  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);
  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);

  // distance do ponto 0 ate metade do segmento 01
  double a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance do ponto 0 ate metade do segmento 02
  double b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // distance da metade do segmento 01 ate metade do segmento 02
  double c = distance(P1x,P1y,P1z,P2x,P2y,P2z);

  // vetores 
  double x1 = P1x-P0x;
  double y1 = P1y-P0y;
  double z1 = P1z-P0z;

  double x2 = P2x-P0x;
  double y2 = P2y-P0y;
  double z2 = P2z-P0z;

  double xReta = P2x-P1x;
  double yReta = P2y-P1y;
  double zReta = P2z-P1z;

  // vetores unitarios deslocados para origem do sistema (0,0,0)
  double x1Unit = x1/a;
  double y1Unit = y1/a;
  double z1Unit = z1/a;

  double x2Unit = x2/b;
  double y2Unit = y2/b;
  double z2Unit = z2/b;

  // normal to each triangular face
  clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

  // somatorio NAO ponderado pela area dos vetores unitarios normais 
  // aos triangulos encontrados na estrela do vertice
  sumXCrossUnit += cross.Get(0);
  sumYCrossUnit += cross.Get(1);
  sumZCrossUnit += cross.Get(2);

  // soma dos vetores 1Unit + 2Unit = resultante
  double xPlaneRes = x1Unit+x2Unit;
  double yPlaneRes = y1Unit+y2Unit;
  double zPlaneRes = z1Unit+z2Unit;

  clVector proj = projection(xPlaneRes,yPlaneRes,zPlaneRes,
	xReta,yReta,zReta);
  double xPlaneTang = proj.Get(0);
  double yPlaneTang = proj.Get(1);
  double zPlaneTang = proj.Get(2);

  /* subtraindo vetor tangente do vetor unitario para encontrar as
   * coordenadas do vetor normal SITUADO NA SUPERFICE (face do
   * tetrahedro = triangulo)
   * 
   *                 P0
   *                  ^
   *                 / \
   *                /   \
   *               /     \
   *            P1 ------- P2
   *                  ----> PlaneTang
   *                  |\
   *                  | \  PlaneRes
   *                  |  \
   *                PlaneNormal
   */
  double xPlaneNormal = xPlaneRes - xPlaneTang;
  double yPlaneNormal = yPlaneRes - yPlaneTang;
  double zPlaneNormal = zPlaneRes - zPlaneTang;

  double len = vectorLength(xPlaneNormal,yPlaneNormal,zPlaneNormal);

  // Unitario do vetor resultante situado no plano do triangulo
  // combinacao linear dos vetores unitarios das arestas do triangulo
  double xPlaneNormalUnit = xPlaneNormal/len;
  double yPlaneNormalUnit = yPlaneNormal/len;
  double zPlaneNormalUnit = zPlaneNormal/len;

  // normal ao plano integrada na distancia (MOD) dos 2 vertices medianos
  // force = resultante das componentes * tamanho da aresta que sera
  // usada como referencia no calculo da area do triangulo 
  double base = c/2.0;

  fx += xPlaneNormalUnit*base;
  fy += yPlaneNormalUnit*base;
  fz += zPlaneNormalUnit*base;

  // 1/3 of area P0-Pm01-Pm02
  double area = getArea(P0x,P0y,P0z,P1x,P1y,P1z,P2x,P2y,P2z);

  sumArea += (1.0/3.0)*area;

  // norb's correction
  //double fact = (a*b*c*c)/ (4*area*area);
  //sumArea += (1.0/4.0)*area*fact;
 }
 mele=_myList.end();

 double len = vectorLength(sumXCrossUnit,sumYCrossUnit,sumZCrossUnit);
 double xNormalUnit = sumXCrossUnit/len;
 double yNormalUnit = sumYCrossUnit/len;
 double zNormalUnit = sumZCrossUnit/len;

 // intensidade da forca resultante
 double force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) )/sumArea;

//--------------------------------------------------
//  if( _node == 55 )
//  {
//   cout << "force: " << fx << "," << fy << "," << fz << endl;
//   cout << "normal: " << xNormalUnit << "," << yNormalUnit << "," <<
//    zNormalUnit << endl;
//  cout << "dot: " << dotProd(fx,fy,fz,xNormalUnit,yNormalUnit,zNormalUnit) << endl;
//  }
//-------------------------------------------------- 

 /* This if statement test weather the curvature (force) will be set
  * to negative or positive following the three cases based on the
  * normal componentes ({x,y,z}NormalUnit) and the computed force
  * (f{x,y,z}):
  *
  *       \              |             /
  *        \     N       |    N       /    N
  *   <---  ) --->       | --->      (  --->
  *   f    /             |            \ --->
  *       /              |  f=0        \   f
  *
  *   dotProd < 0     dotProd = 0   dotProd > 0
  *   force > 0       force = 0     force < 0
  *
  * */
 if( dotProd(fx,fy,fz,xNormalUnit,yNormalUnit,zNormalUnit) > 0.0 )
  force = -force;

 clVector vec(4);
 vec.Set(0,force);       // curvature
 vec.Set(1,xNormalUnit); // x normal
 vec.Set(2,yNormalUnit); // y normal
 vec.Set(3,zNormalUnit); // z normal

 return vec;
} // fecha metodo getNormalAndKappa

/*
 * To compute normalAndKappa at all surfMesh vertices
 *
 * */
void Model3D::setNormalAndKappa()
{
 setSurface();
 setNeighbourSurfaceElem();
 setNeighbourSurfacePoint();

 surfMesh.xNormal.Dim(surfMesh.numVerts);
 surfMesh.yNormal.Dim(surfMesh.numVerts);
 surfMesh.zNormal.Dim(surfMesh.numVerts);
 surfMesh.curvature.Dim(surfMesh.numVerts);

 // loop sobre todos os nos da superficie 
 for( int i=0;i<surface.Dim();i++ )
 {
  int node = surface.Get(i);

//--------------------------------------------------
//   clVector vec = getNormalAndKappa(node,
// 	              getNeighbourSurfacePoint(node));
//-------------------------------------------------- 

  clVector vec = getNormalAndKappa(node,
	              getNeighbourSurfacePoint(node));

  double force       = vec.Get(0);
  double xNormalUnit = vec.Get(1);
  double yNormalUnit = vec.Get(2);
  double zNormalUnit = vec.Get(3);

  surfMesh.xNormal.Set(node,xNormalUnit);
  surfMesh.yNormal.Set(node,yNormalUnit);
  surfMesh.zNormal.Set(node,zNormalUnit);
  surfMesh.curvature.Set(node,force);

//--------------------------------------------------
//   if( node == 2 )
//   {
//    cout << " ---> " << node << " <---" << endl;
//    cout << "x: " << surfMesh.X.Get(node) << endl;
//    cout << "y: " << surfMesh.Y.Get(node) << endl;
//    cout << "z: " << surfMesh.Z.Get(node) << endl;
//    cout << "xNormal: " << xNormalUnit << endl;
//    cout << "yNormal: " << yNormalUnit << endl;
//    cout << "zNormal: " << zNormalUnit << endl;
//    cout << "xN: " << surfMesh.X.Get(node) + xNormalUnit << endl;
//    cout << "yN: " << surfMesh.Y.Get(node) + yNormalUnit << endl;
//    cout << "zN: " << surfMesh.Z.Get(node) + zNormalUnit << endl;
//    cout << "force: " << force << endl;
//    list<int> neighPoint = getNeighbourSurfacePoint(node);
//    for (list<int>::iterator it=neighPoint.begin(); it!=neighPoint.end(); ++it)
// 	cout << *it << " ";
//    cout <<endl;
//   }
//-------------------------------------------------- 
 }
 //setNormalAndKappa2D();
} // fecha metodo setNormalAndKappa

/* input: setSurface()
 * output: closer,xCloser,yCloser,zCloser
 *
 * */
void Model3D::setCloser()
{
 // closer=surface(dsearchn(X(surface),Y(surface),X,Y));
 // esta funcao retorna o noh da interface (surface) mais 
 // proximo de cada noh da malha (vertices)
 closer = dsearchn(xSurface,ySurface,zSurface,X,Y,Z);
 xCloser.Dim( closer.Dim() );
 yCloser.Dim( closer.Dim() );
 zCloser.Dim( closer.Dim() );
 for( int i=0;i<closer.Dim();i++ )
 {
  int aux = closer.Get(i);
  closer.Set(i,surface.Get(aux)); // alterando os valores de closer(i)
  xCloser.Set(i,X.Get( closer.Get(i) ));
  yCloser.Set(i,Y.Get( closer.Get(i) ));
  zCloser.Set(i,Z.Get( closer.Get(i) ));
 }
}

/* Set the distance of each 3D mesh node to the closest surface node
 *
 * input: setCloser (xCloser,yCloser,zCloser)
 * output: interfaceDistance
 *
 * */
void Model3D::setInterfaceDistance()
{
 setCloser();

 interfaceDistance.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  double aux = distance( X.Get(i),Y.Get(i),Z.Get(i),
	                   xCloser.Get(i),yCloser.Get(i),zCloser.Get(i) );
  interfaceDistance.Set(i,aux);
 }
}

/* espalhando kappa calculado na superfice para todos os pontos da
 * malha, com isso garantimos uma forca distribuida
 * */
void Model3D::setKappaSurface()
{
 setCloser();
 curvature.Dim(3*numNodes);
 mesh3d.curvature.Dim(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  int aux = closer.Get(i);
  curvature.Set( i,surfMesh.curvature.Get(aux) );
  curvature.Set( i+numNodes,surfMesh.curvature.Get(aux) );
  curvature.Set( i+2*numNodes,surfMesh.curvature.Get(aux) );
  mesh3d.curvature.Set( i,surfMesh.curvature.Get(aux) );
  mesh3d.curvature.Set( i+numNodes,surfMesh.curvature.Get(aux) );
  mesh3d.curvature.Set( i+2*numNodes,surfMesh.curvature.Get(aux) );
 }
}

// espalhando kappa calculado na superfice para todos os pontos da
// malha, com isso garantimos uma forca distribuida
void Model3D::setKappaSurface(clVector &_kappa)
{
 setCloser();
 curvature.Dim(3*numNodes);
 mesh3d.curvature.Dim(3*numNodes);
 for( int i=0;i<numNodes;i++ )
 {
  int aux = closer.Get(i);
  curvature.Set( i,_kappa.Get(aux) );
  curvature.Set( i+numNodes,_kappa.Get(aux) );
  curvature.Set( i+2*numNodes,_kappa.Get(aux) );
  mesh3d.curvature.Set( i,_kappa.Get(aux) );
  mesh3d.curvature.Set( i+numNodes,_kappa.Get(aux) );
  mesh3d.curvature.Set( i+2*numNodes,_kappa.Get(aux) );
 }
}

void Model3D::setInitSurfaceVolume()
{
 initSurfaceVolume.clear();
 initSurfaceVolume.resize((int) surfMesh.elemIdRegion.Max()+1);

 // surfMesh.elemIdRegion == 0 --> none
 // surfMesh.elemIdRegion == 1 --> wall
 // surfMesh.elemIdRegion == 2 --> bubble 1
 // surfMesh.elemIdRegion == 3 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
  initSurfaceVolume[nb] = getSurfaceVolume(nb);
}

void Model3D::setSurfaceVolume()
{
 surfaceVolume.clear();
 surfaceVolume.resize((int) surfMesh.elemIdRegion.Max()+1);

 // surfMesh.elemIdRegion == 0 --> none
 // surfMesh.elemIdRegion == 1 --> wall
 // surfMesh.elemIdRegion == 2 --> bubble 1
 // surfMesh.elemIdRegion == 3 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
  surfaceVolume[nb] = getSurfaceVolume(nb);
}

void Model3D::setInitSurfaceArea()
{
 initSurfaceArea.clear();
 initSurfaceArea.resize((int) surfMesh.elemIdRegion.Max()+1);

 // surfMesh.elemIdRegion == 0 --> none
 // surfMesh.elemIdRegion == 1 --> wall
 // surfMesh.elemIdRegion == 2 --> bubble 1
 // surfMesh.elemIdRegion == 3 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
  initSurfaceArea[nb] = getSurfaceArea(nb);
}

void Model3D::setSurfaceArea()
{
 surfaceArea.clear();
 surfaceArea.resize((int) surfMesh.elemIdRegion.Max()+1);

 // surfMesh.elemIdRegion == 0 --> none
 // surfMesh.elemIdRegion == 1 --> wall
 // surfMesh.elemIdRegion == 2 --> bubble 1
 // surfMesh.elemIdRegion == 3 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
  surfaceArea[nb] = getSurfaceArea(nb);
}

/* 
 * Computing the volume of a surface V by the divergence theorem:
 *
 * V = \int_V dV = (1/3) * \oint_S x \cdot n dS
 *
 * The centroid X is given by:
 *
 * S = (1/V) * \int_V x dV = (1/2V) * \oint_S (x \cdot x) n dS 
 *
 * The centroid velocity W is:
 *
 * W = (1/V) \int w dV = (1/V) * \oint (zu) \cdot n dS
 *
 * OBS: The mesh needs to be oriented, otherwise the method doesn't work!
 * */
double Model3D::getSurfaceVolume(int _region)
{
 double sumVolume = 0;
 //double sumCentroidX = 0;
 //double sumCentroidY = 0;
 //double sumCentroidZ = 0;
 for( int mele=0;mele<surfMesh.numElems;mele++ )
 {
  if( _region == surfMesh.elemIdRegion.Get(mele) )
  {
   // P1
   int v1 = surfMesh.IEN.Get(mele,0);
   double p1x = surfMesh.X.Get(v1);
   double p1y = surfMesh.Y.Get(v1);
   double p1z = surfMesh.Z.Get(v1);

   // P2
   int v2 = surfMesh.IEN.Get(mele,1);
   double p2x = surfMesh.X.Get(v2);
   double p2y = surfMesh.Y.Get(v2);
   double p2z = surfMesh.Z.Get(v2);

   // P3
   int v3 = surfMesh.IEN.Get(mele,2);
   double p3x = surfMesh.X.Get(v3);
   double p3y = surfMesh.Y.Get(v3);
   double p3z = surfMesh.Z.Get(v3);

   // element centroid
   double xCentroid = (p1x+p2x+p3x)/3.0;
   double yCentroid = (p1y+p2y+p3y)/3.0;
   double zCentroid = (p1z+p2z+p3z)/3.0;

   // distance from point 1 to 2
   double a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

   // distance from point 2 to 3
   double b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

   // unit vectors
   double x1Unit = (p2x-p1x)/a;
   double y1Unit = (p2y-p1y)/a;
   double z1Unit = (p2z-p1z)/a;

   double x2Unit = (p3x-p2x)/b;
   double y2Unit = (p3y-p2y)/b;
   double z2Unit = (p3z-p2z)/b;

   // calculando o produto vetorial de cada elemento triangular da superficie
   clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

   // somatorio ponderado pela area dos vetores unitarios normais 
   // aos triangulos encontrados na estrela do vertice
   double xNormalElem = cross.Get(0);
   double yNormalElem = cross.Get(1);
   double zNormalElem = cross.Get(2);

   double len = vectorLength(xNormalElem,yNormalElem,zNormalElem);

   double xNormalElemUnit = xNormalElem/len;
   double yNormalElemUnit = yNormalElem/len;
   double zNormalElemUnit = zNormalElem/len;

   double area = getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z);

   sumVolume += ( xCentroid*xNormalElemUnit + 
	              yCentroid*yNormalElemUnit +
	              zCentroid*zNormalElemUnit ) * area;

//--------------------------------------------------
//    sumCentroidX += ( xCentroid*xCentroid + 
// 	                 yCentroid*yCentroid +
// 	                 zCentroid*zCentroid ) * xNormalElemUnit * area;
//    sumCentroidY += ( xCentroid*xCentroid + 
// 	                 yCentroid*yCentroid +
// 	                 zCentroid*zCentroid ) * yNormalElemUnit * area;
//    sumCentroidZ += ( xCentroid*xCentroid + 
// 	                 yCentroid*yCentroid +
// 	                 zCentroid*zCentroid ) * zNormalElemUnit * area;
//-------------------------------------------------- 
  }
 }
 double vol = (1.0/3.0)*sumVolume;
 //double xc = sumCentroidX/(2*vol);
 //double yc = sumCentroidY/(2*vol);
 //double zc = sumCentroidZ/(2*vol);

 return vol;
}

/*
 * Calculates the volume of a surface identified by (int _region). This
 * method computes the sum of all inner tetrahdron elements for a 
 * specific region.
 *
 * input: region defined by elemIdRegion.
 *        1 none
 *        2 wall
 *        3 bubble1
 *        4 bubble2, etc.
 *
 * return: volume of a surface
 *
 * */
double Model3D::getSurfaceVolumeTET(int _region)
{
 double sumVolume=0;
 for( int i=0;i<numElems;i++ )
  if( _region == elemIdRegion.Get(i) )
   sumVolume += getVolume(i);
 return sumVolume;
}

/*
 * Calculates the area of a surface identified by (int _region). This
 * method computes the sum of all surface triangle elements for a 
 * specific region.
 *
 * input: region defined by elemIdRegion.
 *        1 none
 *        2 wall
 *        3 bubble1
 *        4 bubble2, etc.
 *
 * return: area of a surface
 *
 * */
double Model3D::getSurfaceArea(int _region)
{
 double sumArea = 0;
 for(int tri=0;tri<surfMesh.numElems;tri++ )
  if( surfMesh.elemIdRegion.Get(tri) == _region )
  {
   // P1
   int v1 = surfMesh.IEN.Get(tri,0);
   double p1x = surfMesh.X.Get(v1);
   double p1y = surfMesh.Y.Get(v1);
   double p1z = surfMesh.Z.Get(v1);

   // P2
   int v2 = surfMesh.IEN.Get(tri,1);
   double p2x = surfMesh.X.Get(v2);
   double p2y = surfMesh.Y.Get(v2);
   double p2z = surfMesh.Z.Get(v2);

   // P3
   int v3 = surfMesh.IEN.Get(tri,2);
   double p3x = surfMesh.X.Get(v3);
   double p3y = surfMesh.Y.Get(v3);
   double p3z = surfMesh.Z.Get(v3);

   sumArea += getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z);
  }
 return sumArea;
}

void Model3D::setSingleElement()
{
 numVerts = 4;
 numElems = 1;
 numNodes = 4;

 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 IEN.Dim(numElems);

 // point 0
 X.Set(0,0.0);
 Y.Set(0,0.0);
 Z.Set(0,0.0);

 // point 1
 X.Set(1,1.0);
 Y.Set(1,0.0);
 Z.Set(1,0.0);

 // point 2
 X.Set(2,0.0);
 Y.Set(2,1.0);
 Z.Set(2,0.0);

 // point 3
 X.Set(3,0.0);
 Y.Set(3,0.0);
 Z.Set(3,1.0);

 // elem 1
 IEN.Set(0,0,0);
 IEN.Set(0,1,1);
 IEN.Set(0,2,2);
 IEN.Set(0,3,3);
}

void Model3D::setTwoElements()
{
 numVerts = 5;
 numElems = 2;
 numNodes = 5;

 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 IEN.Dim(numElems);

 // point 0
 X.Set(0,0.0);
 Y.Set(0,0.0);
 Z.Set(0,0.0);

 // point 1
 X.Set(1,1.0);
 Y.Set(1,0.0);
 Z.Set(1,0.0);

 // point 2
 X.Set(2,0.0);
 Y.Set(2,1.0);
 Z.Set(2,0.0);

 // point 3
 X.Set(3,0.0);
 Y.Set(3,0.0);
 Z.Set(3,1.0);

 // point 4
 X.Set(4,0.0);
 Y.Set(4,-1.0);
 Z.Set(4,0.0);

 // elem 1
 IEN.Set(0,0,0);
 IEN.Set(0,1,1);
 IEN.Set(0,2,2);
 IEN.Set(0,3,3);

 // elem 2
 IEN.Set(1,0,0);
 IEN.Set(1,1,1);
 IEN.Set(1,2,3);
 IEN.Set(1,3,4);
}

void Model3D::setThreeElements()
{
 numVerts = 6;
 numElems = 3;
 numNodes = 6;

 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 IEN.Dim(numElems);

 // point 0
 X.Set(0,0.0);
 Y.Set(0,0.0);
 Z.Set(0,0.0);

 // point 1
 X.Set(1,1.0);
 Y.Set(1,0.0);
 Z.Set(1,0.0);

 // point 2
 X.Set(2,0.0);
 Y.Set(2,1.0);
 Z.Set(2,0.0);

 // point 3
 X.Set(3,0.0);
 Y.Set(3,0.0);
 Z.Set(3,1.0);

 // point 4
 X.Set(4,0.0);
 Y.Set(4,-1.0);
 Z.Set(4,0.0);

 // point 5
 X.Set(4,-1.0);
 Y.Set(4,0.0);
 Z.Set(4,0.0);

 // elem 0
 IEN.Set(0,0,0);
 IEN.Set(0,1,1);
 IEN.Set(0,2,2);
 IEN.Set(0,3,3);

 // elem 1
 IEN.Set(1,0,0);
 IEN.Set(1,1,1);
 IEN.Set(1,2,3);
 IEN.Set(1,3,4);

 // elem 2
 IEN.Set(2,0,0);
 IEN.Set(2,1,2);
 IEN.Set(2,2,3);
 IEN.Set(2,3,5);
}

void Model3D::setFourElements()
{
 numVerts = 6;
 numElems = 4;
 numNodes = 6;

 X.Dim(numVerts);
 Y.Dim(numVerts);
 Z.Dim(numVerts);
 IEN.Dim(numElems);

 // point 0     // point 1    // point 2    // point 3
 X.Set(0,0.0);  X.Set(1,1.0); X.Set(2,0.0); X.Set(3,0.0);
 Y.Set(0,0.0);  Y.Set(1,0.0); Y.Set(2,1.0); Y.Set(3,0.0);
 Z.Set(0,0.0);  Z.Set(1,0.0); Z.Set(2,0.0); Z.Set(3,1.0);

 // point 4     // point 5
 X.Set(4,0.0);  X.Set(4,-1.0);
 Y.Set(4,-1.0); Y.Set(4,0.0);
 Z.Set(4,0.0);  Z.Set(4,0.0);

 // elem 0       // elem 1       // elem 2       // elem 3
 IEN.Set(0,0,0); IEN.Set(1,0,0); IEN.Set(2,0,0); IEN.Set(3,0,0);
 IEN.Set(0,1,1); IEN.Set(1,1,1); IEN.Set(2,1,2); IEN.Set(3,1,3);
 IEN.Set(0,2,2); IEN.Set(1,2,3); IEN.Set(2,2,3); IEN.Set(3,2,4);
 IEN.Set(0,3,3); IEN.Set(1,3,4); IEN.Set(2,3,5); IEN.Set(3,3,5);
}


/*
 * Check the orientation of each surface triagle (convex-hull + bubbles)
 * by comparing (dot product) the orientation of the element's normal vector 
 * and the direction of the element to the centroid relative to the
 * geometrical shape. For example, if the simulation has 1 bubble, then
 * this method calculates the centroid of the convex-hull and then
 * change the orientation of each wrong oriented triangle, finally it
 * computes the centroid of the bubble and change the orientation of
 * each wrong oriented triangle that is part of the bubble.
 *
 *       o ------------- o         
 *       |               |         
 *       |               |            -      
 *       |               |          /   \
 *       |       x       |    +    |  x  |   
 *       |               |          \   /             
 *       |               |            -               
 *       |               |         
 *       o ------------- o                      x centroid
 *
 *
 * OBS: works only for convex shapes

 * ALTERNATIVE: GMsh can verify the orientation of a triangle surface
 * mesh. To do so, select TOOLS->OPTIONS->MESH-VISIBILITY and set a
 * non-zero (10 or 100) value on the NORMAL BOX. It will be able to
 * check the normal directions of each PHYSICAL SURFACE, with something
 * is wrong, you can edit the .geo file and change the sign (+ or -) of the
 * specific PHYSICAL SURFACE.
 * 
 * If GMsh delivers an oriented surface mesh, this method becomes
 * obsolete.
 * */
void Model3D::checkTriangleOrientation()
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int nb=0;nb<=surfMesh.elemIdRegion.Max();nb++ )
 {
  clVector centroid = computeConvexRegionCentroid(nb);
  double xc = centroid.Get(0);
  double yc = centroid.Get(1);
  double zc = centroid.Get(2);

  for( int elem=0;elem<surfMesh.numElems;elem++ )
  {
   if( surfMesh.elemIdRegion.Get(elem) == nb )
   {
	int v1 = surfMesh.IEN.Get(elem,0);
	double p1x = surfMesh.X.Get(v1);
	double p1y = surfMesh.Y.Get(v1);
	double p1z = surfMesh.Z.Get(v1);

	int v2 = surfMesh.IEN.Get(elem,1);
	double p2x = surfMesh.X.Get(v2);
	double p2y = surfMesh.Y.Get(v2);
	double p2z = surfMesh.Z.Get(v2);

	int v3 = surfMesh.IEN.Get(elem,2);
	double p3x = surfMesh.X.Get(v3);
	double p3y = surfMesh.Y.Get(v3);
	double p3z = surfMesh.Z.Get(v3);

	/*               
	 *               v3
	 *               o 
	 *              / \
	 *             /   \
	 *            /     \
	 *           o ----- o 
	 *         v1         v2
	 *
	 * */

	double vx = p1x-xc;
	double vy = p1y-yc;
	double vz = p1z-zc;

	double v1x = p2x-p1x;
	double v1y = p2y-p1y;
	double v1z = p2z-p1z;

	double v2x = p3x-p2x;
	double v2y = p3y-p2y;
	double v2z = p3z-p2z;

	// normal to each triangular face
	clVector cross = crossProd(v1x,v1y,v1z,v2x,v2y,v2z);

	if( dotProd(vx,vy,vz,cross.Get(0),cross.Get(1),cross.Get(2)) < 0.0 )
	{
	 surfMesh.IEN.Set(elem,0,v2);
	 surfMesh.IEN.Set(elem,1,v1);
	}
   }
  }
 }
}

clVector Model3D::computeConvexRegionCentroid(int _region)
{
 double sumX = 0;
 double sumY = 0;
 double sumZ = 0;
 int count = 0;
 for( int elem=0;elem<surfMesh.numElems;elem++ )
 {
  if( _region == surfMesh.elemIdRegion.Get(elem) )
  {
   int v1 = surfMesh.IEN.Get(elem,0);
   int v2 = surfMesh.IEN.Get(elem,1);
   int v3 = surfMesh.IEN.Get(elem,2);

   sumX += ( surfMesh.X.Get(v1)+surfMesh.X.Get(v2)+surfMesh.X.Get(v3) )/3.0;
   sumY += ( surfMesh.Y.Get(v1)+surfMesh.Y.Get(v2)+surfMesh.Y.Get(v3) )/3.0;
   sumZ += ( surfMesh.Z.Get(v1)+surfMesh.Z.Get(v2)+surfMesh.Z.Get(v3) )/3.0;
   count++;
  }
 }
 double xc = sumX/count;
 double yc = sumY/count;
 double zc = sumZ/count;

 clVector centroid(3);
 centroid.Set(0,xc);
 centroid.Set(1,yc);
 centroid.Set(2,zc);

 return centroid;
}

//--------------------------------------------------
// clVector Model3D::getNormalAndKappaByGauss(int _node,list<int> _myList)
// {
//  double P0x = surfMesh.X.Get(_node);
//  double P0y = surfMesh.Y.Get(_node);
//  double P0z = surfMesh.Z.Get(_node);
// 
//  // Normal Vector
//  double xNormal = 0;
//  double yNormal = 0;
//  double zNormal = 0;
// 
//  // Fundamental Form Coefficients
//  double E = 0;
//  double F = 0;
//  double G = 0;
//  double L = 0;
//  double M = 0;
//  double N = 0;
// 
//  int listSize = _myList.size();
//  list<int>::iterator mele=_myList.begin();
//  for( int i=0;i<listSize-1;i++ )
//  {
//   int v1 = *mele;++mele;
//   int v2 = *mele;
// 
//   double P1x = surfMesh.X.Get(v1);
//   double P1y = surfMesh.Y.Get(v1);
//   double P1z = surfMesh.Z.Get(v1);
//   double P2x = surfMesh.X.Get(v2);
//   double P2y = surfMesh.Y.Get(v2);
//   double P2z = surfMesh.Z.Get(v2);
// 
//   // distance do ponto 0 ate metade do segmento 01
//   double a = distance(P0x,P0y,P0z,P1x,P1y,P1z);
// 
//   // distance do ponto 0 ate metade do segmento 02
//   double b = distance(P0x,P0y,P0z,P2x,P2y,P2z);
// 
//   // vetores 
//   double x1 = P1x-P0x;
//   double y1 = P1y-P0y;
//   double z1 = P1z-P0z;
// 
//   double x2 = P2x-P0x;
//   double y2 = P2y-P0y;
//   double z2 = P2z-P0z;
// 
//   // vetores unitarios deslocados para origem do sistema (0,0,0)
//   double x1Unit = x1/a;
//   double y1Unit = y1/a;
//   double z1Unit = z1/a;
// 
//   double x2Unit = x2/b;
//   double y2Unit = y2/b;
//   double z2Unit = z2/b;
// 
//   // normal to each triangular face
//   clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);
// 
//   // somatorio NAO ponderado pela area dos vetores unitarios normais 
//   // aos triangulos encontrados na estrela do vertice
//   xNormal += cross.Get(0);
//   yNormal += cross.Get(1);
//   zNormal += cross.Get(2);
// 
//   // vector Ru
//   double ux = x1Unit - xNormal;
//   double uy = y1Unit - yNormal;
//   double uz = z1Unit - zNormal;
// 
//   // vector Rv
//   double vx = x2Unit - xNormal;
//   double vy = y2Unit - yNormal;
//   double vz = z2Unit - zNormal;
// 
//   // First Fundamental Form
//   E = ux*ux + uy*uy + uz*uz;
//   F = ux*vx + uy*vy + uz*vz;
//   G = vx*vx + vy*vy + vz*vz;
// 
//   // vector Ruu
//   double uux = 0;
//   double uuy = 0;
//   double uuz = 0;
// 
//   // vector Ruv
//   double uvx = 0;
//   double uvy = 0;
//   double uvz = 0;
// 
//   // vector Rvv
//   double vvx = 0;
//   double vvy = 0;
//   double vvz = 0;
//   
//   // Second Fundamental Form
//   L = uux*xNormal + uuy*yNormal + uuz*zNormal;
//   M = uvx*xNormal + uvy*yNormal + uvz*zNormal;
//   N = vvx*xNormal + vvy*yNormal + vvz*zNormal;
// 
//  }
//  mele=_myList.end();
// 
//  double len = vectorLength(xNormal,yNormal,zNormal);
//  double xNormalUnit = xNormal/len;
//  double yNormalUnit = yNormal/len;
//  double zNormalUnit = zNormal/len;
// 
//  double meanCurvature = (E*N + G*L - 2*F*M)/2*(E*G-F*F);
// 
//  clVector vec(4);
//  vec.Set(0,meanCurvature);
//  vec.Set(1,xNormalUnit);
//  vec.Set(2,yNormalUnit);
//  vec.Set(3,zNormalUnit);
// 
//  return vec;
// } // fecha metodo getNormalAndKappaByGauss
// 
//-------------------------------------------------- 
clVector Model3D::getNormalAndKappaByDesbrun(int _node,list<int> _myList)
{
 double P0x = surfMesh.X.Get(_node);
 double P0y = surfMesh.Y.Get(_node);
 double P0z = surfMesh.Z.Get(_node);

 //int c1 = 0;
 double fx = 0;
 double fy = 0;
 double fz = 0;
 double sumMixedArea = 0;
 double sumXCrossUnit = 0;
 double sumYCrossUnit = 0;
 double sumZCrossUnit = 0;

 int listSize = _myList.size();
 list<int>::iterator mele=_myList.begin();

 // adding 2nd vertex to the end of the list
 // old: 0 1 2 3 4 5 0 
 // new: 0 1 2 3 4 5 0 1
 list<int> li = _myList;
 list<int>::iterator mele2=li.begin();
 ++mele2;
 _myList.push_back(*mele2);

 for( int i=0;i<listSize-1;i++ )
 {
  int v1 = *mele;++mele;
  int v2 = *mele;++mele;
  int v3 = *mele;--mele;

  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);
  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);
  double P3x = surfMesh.X.Get(v3);
  double P3y = surfMesh.Y.Get(v3);
  double P3z = surfMesh.Z.Get(v3);

  // distance from point 0 to 1  
  double a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance from point 0 to 2  
  double b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // distance from point 0 to 3  
  double c = distance(P0x,P0y,P0z,P3x,P3y,P3z);

  // distance from point 1 to 2  
  double d = distance(P1x,P1y,P1z,P2x,P2y,P2z);

  // distance from point 2 to 3  
  double e = distance(P2x,P2y,P2z,P3x,P3y,P3z);

  // vetores 
  double x1 = P1x-P0x;
  double y1 = P1y-P0y;
  double z1 = P1z-P0z;

  double x2 = P2x-P0x;
  double y2 = P2y-P0y;
  double z2 = P2z-P0z;

  // vetores unitarios deslocados para origem do sistema (0,0,0)
  double x1Unit = x1/a;
  double y1Unit = y1/a;
  double z1Unit = z1/a;

  double x2Unit = x2/b;
  double y2Unit = y2/b;
  double z2Unit = z2/b;

  // normal to each triangular face
  clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

  // somatorio NAO ponderado pela area dos vetores unitarios normais 
  // aos triangulos encontrados na estrela do vertice
  sumXCrossUnit += cross.Get(0);
  sumYCrossUnit += cross.Get(1);
  sumZCrossUnit += cross.Get(2);

  // angles (law of cosine)
  double alpha = acos( -(b*b-d*d-a*a)/(2*d*a) );
  double beta = acos( -(b*b-c*c-e*e)/(2*c*e) );

  double a1 = acos( -(d*d-a*a-b*b)/(2*a*b) );
  double a2 = alpha;
  double a3 = acos( -(a*a-b*b-d*d)/(2*b*d) );

//--------------------------------------------------
//   cout << a1 << " " << a2 << " " << " " << a3 << " " << a1+a2+a3 <<endl;
//   cout << alpha << " " << beta << endl;
//   cout << "----"<< endl;
//-------------------------------------------------- 

  fx += ( 1.0/tan(alpha)+1.0/tan(beta) )*(x2);
  fy += ( 1.0/tan(alpha)+1.0/tan(beta) )*(y2);
  fz += ( 1.0/tan(alpha)+1.0/tan(beta) )*(z2);

  // area P0-P1-P2
  double area = getArea(P0x,P0y,P0z,P1x,P1y,P1z,P2x,P2y,P2z);
  double voronoiArea = (1.0/8.0)*( (a*a/(tan(a3))) + (b*b/(tan(a2))));

  // voronoi area
  if( a1 < 3.1415/2.0 && a2 < 3.1415/2.0 && a3 < 3.1415/2.0 )
   sumMixedArea += voronoiArea;
  else
   if( a1 > 3.1415/2.0 )
	sumMixedArea += area/2.0;
   else
	sumMixedArea += area/4.0;
 }
 mele=_myList.end();

 fx = fx/(2*sumMixedArea);
 fy = fy/(2*sumMixedArea);
 fz = fz/(2*sumMixedArea);

 double len = vectorLength(sumXCrossUnit,sumYCrossUnit,sumZCrossUnit);
 double xNormalUnit = sumXCrossUnit/len;
 double yNormalUnit = sumYCrossUnit/len;
 double zNormalUnit = sumZCrossUnit/len;

 // resulting force (already divided by 2*area)
 double force = sqrt( (fx*fx)+(fy*fy)+(fz*fz) );

 /* This if statement test weather the curvature (force) will be set
  * to negative or positive following the three cases based on the
  * normal componentes ({x,y,z}NormalUnit) and the computed force
  * (f{x,y,z}):
  *
  *       \              |             /
  *        \     N       |    N       /    N
  *   <---  ) --->       | --->      (  --->
  *   f    /             |            \ --->
  *       /              |  f=0        \   f
  *
  *   dotProd < 0     dotProd = 0   dotProd > 0
  *   force > 0       force = 0     force < 0
  *
  * */
 if( dotProd(fx,fy,fz,xNormalUnit,yNormalUnit,zNormalUnit) > 0.0 )
  force = -force;

 clVector vec(4);
 vec.Set(0,force);       // curvature
 vec.Set(1,xNormalUnit); // x normal
 vec.Set(2,yNormalUnit); // y normal
 vec.Set(3,zNormalUnit); // z normal

 return vec;
} // fecha metodo getNormalAndKappaByDesbrun

clVector Model3D::getNormalElem(int _elem)
{
 int v1 = surfMesh.IEN.Get(_elem,0);
 double p1x = surfMesh.X.Get(v1);
 double p1y = surfMesh.Y.Get(v1);
 double p1z = surfMesh.Z.Get(v1);

 int v2 = surfMesh.IEN.Get(_elem,1);
 double p2x = surfMesh.X.Get(v2);
 double p2y = surfMesh.Y.Get(v2);
 double p2z = surfMesh.Z.Get(v2);

 int v3 = surfMesh.IEN.Get(_elem,2);
 double p3x = surfMesh.X.Get(v3);
 double p3y = surfMesh.Y.Get(v3);
 double p3z = surfMesh.Z.Get(v3);
 
 // distance from point 1 to 2
 double a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

 // distance from point 2 to 3
 double b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

 /*               
  *               v3
  *               o 
  *              / \
  *             /   \ b
  *            /     \
  *           o ----- o 
  *         v1    a    v2
  *
  * */
 double v1x = (p2x-p1x)/a;
 double v1y = (p2y-p1y)/a;
 double v1z = (p2z-p1z)/a;

 double v2x = (p3x-p2x)/b;
 double v2y = (p3y-p2y)/b;
 double v2z = (p3z-p2z)/b;

 clVector normal = crossProd(v1x,v1y,v1z,v2x,v2y,v2z);
 double length = vectorLength(normal.Get(0),normal.Get(1),normal.Get(2));

 // unit normal to each triangular face
 return normal/length;
} // fecha metodo getNormalElem

clVector Model3D::getNormalElem(int _v1,int _v2,int _v3)
{
 double p1x = surfMesh.X.Get(_v1);
 double p1y = surfMesh.Y.Get(_v1);
 double p1z = surfMesh.Z.Get(_v1);

 double p2x = surfMesh.X.Get(_v2);
 double p2y = surfMesh.Y.Get(_v2);
 double p2z = surfMesh.Z.Get(_v2);

 double p3x = surfMesh.X.Get(_v3);
 double p3y = surfMesh.Y.Get(_v3);
 double p3z = surfMesh.Z.Get(_v3);
 
 // distance from point 1 to 2
 double a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

 // distance from point 2 to 3
 double b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

 /*               
  *               v3
  *               o 
  *              / \
  *             /   \ b
  *            /     \
  *           o ----- o 
  *         v1    a    v2
  *
  * */
 double v1x = (p2x-p1x)/a;
 double v1y = (p2y-p1y)/a;
 double v1z = (p2z-p1z)/a;

 double v2x = (p3x-p2x)/b;
 double v2y = (p3y-p2y)/b;
 double v2z = (p3z-p2z)/b;

 clVector normal = crossProd(v1x,v1y,v1z,v2x,v2y,v2z);
 double length = vectorLength(normal.Get(0),normal.Get(1),normal.Get(2));

 // unit normal to each triangular face
 return normal/length;
} // fecha metodo getNormalElem

/*
 * Loop over all the surface edges, checking the angle between normal
 * vectors of both planes sharing the edge.
 *         
 *         ^
 *         |
 *         | N1
 *         |
 *    x ------- o
 *    v1        |   N2
 *              | ----->     theta = 90
 *              |
 *              |
 *              x v2
 *
 * theta --> 0, the planes are NOT collapsing
 * theta --> 180, the planes are about to collapse
 *
 * */
void Model3D::checkAngleBetweenPlanes()
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int edge=0;edge<mapEdgeTri.DimI();edge++ )
 {
  double v1 = mapEdgeTri.Get(edge,1);
  double v2 = mapEdgeTri.Get(edge,2);
  double v3elem1 = mapEdgeTri.Get(edge,3);
  double v3elem2 = mapEdgeTri.Get(edge,4);
  int elem1 = mapEdgeTri.Get(edge,5);
  int elem2 = mapEdgeTri.Get(edge,6);
  int vertID = surfMesh.vertIdRegion.Get(v1);

  // elem1
  clVector normalElem1 = getNormalElem(elem1);
  normalElem1 = normalElem1/20;
  // centroid 
  clVector centroidTRIElem1 = centroidTRI3D(surfMesh.X.Get(v1),
                                            surfMesh.Y.Get(v1),
                                            surfMesh.Z.Get(v1),
                                            surfMesh.X.Get(v2),
                                            surfMesh.Y.Get(v2),
                                            surfMesh.Z.Get(v2),
                                            surfMesh.X.Get(v3elem1),
                                            surfMesh.Y.Get(v3elem1),
                                            surfMesh.Z.Get(v3elem1) );

  // elem2
  clVector normalElem2 = getNormalElem(elem2);
  normalElem2 = normalElem2/20;
  // centroid 
  clVector centroidTRIElem2 = centroidTRI3D(surfMesh.X.Get(v1),
                                            surfMesh.Y.Get(v1),
                                            surfMesh.Z.Get(v1),
                                            surfMesh.X.Get(v2),
                                            surfMesh.Y.Get(v2),
                                            surfMesh.Z.Get(v2),
                                            surfMesh.X.Get(v3elem2),
                                            surfMesh.Y.Get(v3elem2),
                                            surfMesh.Z.Get(v3elem2) );



  double theta = angle3D( normalElem1.Get(0),
	                    normalElem1.Get(1),
						normalElem1.Get(2),
	                    normalElem2.Get(0),
						normalElem2.Get(1),
						normalElem2.Get(2) );

  if( (180*theta/3.1415) > 120  && vertID > 0)
  {
//--------------------------------------------------
//    cout << "v1: " << mapEdgeTri.Get(edge,1) << endl;
//    cout << "v2: " << mapEdgeTri.Get(edge,2) << endl;
//    cout << "elem1:        " << elem1 << endl; 
//    cout << "    centroid: " << centroidTRIElem1.Get(0) << " " 
// 	                        << centroidTRIElem1.Get(1) << " "
// 	                        << centroidTRIElem1.Get(2) << endl;
//    cout << "    normal:   " << centroidTRIElem1.Get(0)+normalElem1.Get(0) << " " 
// 	                        << centroidTRIElem1.Get(1)+normalElem1.Get(1) << " "
// 	                        << centroidTRIElem1.Get(2)+normalElem1.Get(2) << endl;
//    cout << "elem2:        " << elem2 << endl; 
//    cout << "    centroid: " << centroidTRIElem2.Get(0) << " " 
// 	                        << centroidTRIElem2.Get(1) << " "
// 						    << centroidTRIElem2.Get(2) << endl;
//    cout << "    normal:   " << centroidTRIElem2.Get(0)+normalElem2.Get(0) << " " 
// 	                        << centroidTRIElem2.Get(1)+normalElem2.Get(1) << " "
// 						    << centroidTRIElem2.Get(2)+normalElem2.Get(2) << endl;
//    cout << "  theta: " << 180*theta/3.14159 << endl;
//    cout << " --------------- " << endl;
//-------------------------------------------------- 

   cout << "------------- " << color(none,magenta,black) 
	    << "smooth vertex high angle plane (" 
		<< resetColor() << 180*theta/3.14159
		<< color(none,magenta,black) << ") at (" 
		<< resetColor() << surfMesh.vertIdRegion.Get(v3elem1)
		<< color(none,magenta,black) << "): " 
		<< resetColor() << v3elem1 
		<< color(none,magenta,black) << " and " 
		<< resetColor() << v3elem2 
		<< resetColor() << endl;

   if( surfMesh.curvature.Get(v3elem1) > 0 )
	smoothPoint(v3elem1);
   if( surfMesh.curvature.Get(v3elem2) > 0 )
	smoothPoint(v3elem2);

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   spp[vertID]++;
   opersurf[vertID]++;
  }
 }
}

/* 
 * This if checks whether the point on the surface has only 3 surface
 * elements. This is caused by the re-meshing process that sometimes
 * creates these problematic local mesh. Such a point should be
 * removed to avoid mesh problems.
 * */
void Model3D::removePointsByNeighbourCheck()
{
 for( int i=0;i<surfMesh.numVerts;i++ )
  removePointByNeighbourCheck(i);
}

void Model3D::removePointByNeighbourCheck(int _node)
{
 int vertID = surfMesh.vertIdRegion.Get(_node);
 /* 
  * This if checks whether the point on the surface has only 3 surface
  * elements. This is caused by the re-meshing process that sometimes
  * creates these problematic local mesh. Such a point should be
  * removed to avoid mesh problems.
  * */
 int elemListSize = neighbourSurfaceElem.at( _node ).size();

 if( elemListSize < 4 )
 {
  removeSurfacePoint(_node);

  if( elemListSize < 3 )
  {
   cout << "------------- " << color(none,red,black) 
	<< "removing fake triangle: " << resetColor() 
	<< _node << endl;
   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   rspn[vertID]++;
   opersurf[vertID]++;
  }

  /*
   * This if removes a surface point when the number of its neighbours
   * are equal. Usually these points demage the mesh
   * quality and breaks the simulation flow. 
   * */
  if( elemListSize == 3 )
  {
   cout << "------------- " << color(none,red,black) 
	<< "removing low-quality point cluster: " << resetColor() 
	<< _node << endl;

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   rspn[vertID]++;
   opersurf[vertID]++;
  }
 }
}

void Model3D::smoothPointsByCurvature()
{
 for( int surfaceNode=0;surfaceNode<surfMesh.numVerts;surfaceNode++ )
 {
  double vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
  double curv = fabs(surfMesh.curvature.Get(surfaceNode));

  double maxEdgeLength = 0;
  int listSize = neighbourPoint.at(surfaceNode).size();
  list<int> plist = neighbourPoint.at(surfaceNode);
  list<int>::iterator vert=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   double P0x = surfMesh.X.Get(surfaceNode);
   double P0y = surfMesh.Y.Get(surfaceNode);
   double P0z = surfMesh.Z.Get(surfaceNode);

   int v1 = *vert;++vert;
   double P1x = surfMesh.X.Get(v1);
   double P1y = surfMesh.Y.Get(v1);
   double P1z = surfMesh.Z.Get(v1);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

   if( edgeLength > maxEdgeLength )
	maxEdgeLength = edgeLength;
  }

  //if( vertID > 0 && maxEdgeLength*curv > 3.5 )
  if( vertID > 0 && curv > 65 )
  {
   cout << "------------- " << color(none,magenta,black) 
	<< "smoothing vertex with curvature (" 
	<< resetColor() << curv
	<< color(none,magenta,black) 
	<< ") at (" 
	<< resetColor()
	<< surfMesh.vertIdRegion.Get(surfaceNode)
	<< color(none,magenta,black) 
	<< "): "
	 << resetColor() << surfaceNode << endl;

   smoothPoint(surfaceNode);

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   spc[vertID]++;
   opersurf[vertID]++;
  }
 }
}

void Model3D::smoothPoint(int _node)
{
 double aux;
 double xSum = 0.0;
 double ySum = 0.0;
 double zSum = 0.0;
 double distSum = 0.0;

 int listSize = neighbourPoint.at(_node).size();
 list<int> plist = neighbourPoint.at(_node);
 list<int>::iterator vert=plist.begin();
 for( int i=0;i<listSize-1;i++ )
 {
  double P0x = surfMesh.X.Get(_node);
  double P0y = surfMesh.Y.Get(_node);
  double P0z = surfMesh.Z.Get(_node);

  int v1 = *vert;++vert;
  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);

  double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  distSum += edgeLength;   
  xSum += ( P1x-P0x )*edgeLength;
  ySum += ( P1y-P0y )*edgeLength;
  zSum += ( P1z-P0z )*edgeLength;
 }
 aux = X.Get(_node) + (1.0/distSum)*xSum; 
 X.Set(_node,aux);
 surfMesh.X.Set(_node,aux);

 aux = Y.Get(_node) + (1.0/distSum)*ySum; 
 Y.Set(_node,aux);
 surfMesh.Y.Set(_node,aux);

 aux = Z.Get(_node) + (1.0/distSum)*zSum; 
 Z.Set(_node,aux);
 surfMesh.Z.Set(_node,aux);
}

/* This method sets triEdge value, which is used many times during the
 * simulation: insert points, delete points, volume correction, define
 * in and out regions in two-phase flows etc. Currently, triEdge is
 * being set to the averageTriLength, which is the average length value
 * of each region (elemIdRegion).
 * surfMesh.elemIdRegion 0 = wall
 * surfMesh.elemIdRegion 1 = surface 1
 * surfMesh.elemIdRegion 2 = surface 2 (if it has more than 1 bubble)
 * surfMesh.elemIdRegion 3 = surface 3 (if it has more than 2 bubbles)
 * */
void Model3D::setTriEdge()
{
 initMeshParameters();
 triMeshStats();
 triEdge = averageTriLength;
}

/* This method sets the value of triEdge which will be used on many
 * methods during the simulation.
 * */
void Model3D::setTriEdge(vector< double > _triEdge)
{
 initMeshParameters();
 triEdge = _triEdge;
}

/* initialize index variables. These variables are set according to the
 * region in the domain (surface, 3dmesh, inside bubble, outside
 * bubble). They are not required on the calculations, however they keep
 * important information about the remeshing process.
 *
 * input: index variables
 * output: initialized index variables
 *
 * */
void Model3D::initMeshParameters()
{
 // number of surfaces 
 int numSurface = surfMesh.numInterfaces+1;
 oper.resize(numSurface);    // oper: num of operations of time step
 opersurf.resize(numSurface); // oper: num of operations on the surf
                                      
 isp.resize(numSurface);  // isp: num of inserted surface points by length
 ispc.resize(numSurface); // ispc: num of inserted surface points by curv
 rsp.resize(numSurface);  // rsp: num of removed surface points by length
 rspn.resize(numSurface); // rspn: num of removed surface points by neigh check
 rspc.resize(numSurface); // rspc: num of removed surface points by curv
 csp.resize(numSurface);  // csp: num of contracted surface points
 flip.resize(numSurface); // flip: flipping operations
 spc.resize(numSurface);  // spc: smoothing operations                           
 spp.resize(numSurface);  // spp: smoothing operations
                                                                                  
 ip.resize(numSurface);     // ip: num of inserted 3d mesh points
 ipd.resize(numSurface);    // ipd: by diffusion 
 rp.resize(numSurface);     // rp: num of removed 3d mesh points
 rpi.resize(numSurface);    // rpi: by interface distance
 rpv.resize(numSurface);    // rpv: by volume 
 rpd.resize(numSurface);    // rpd: by diffusion 
 rpdist.resize(numSurface); // rpd: by distance 
 rph.resize(numSurface);    // rph: by height                                      
 badtet.resize(numSurface); // num of shit tetrahedrons
                                                                                            
 // set surface lengths                
 triEdge.resize(numSurface); // surface triangle length by region
 averageTriLength.resize(numSurface); // average surface triangle length
 averageTriArea.resize(numSurface);   // average surface triangle area
 averageTetVolume.resize(numSurface); // average tetrahedron volume
 tetVol.resize(numSurface);  // recommended tet volume for each region
 //initSurfaceVolume.resize(numSurface);  // init surface volume
 //surfaceVolume.resize(numSurface);  // surface volume
 idMinVolume.resize(numSurface);  // ID of min tet volume
 idMaxVolume.resize(numSurface);  // ID of max tet volume
 minVolume.resize(numSurface);    // min tet volume
 maxVolume.resize(numSurface);    // max tet volume     
 dVolume.resize(numSurface);      // delta volume     
 errorVolume.resize(numSurface);  // error volume     
 //initSurfaceArea.resize(numSurface);  // init surface area 
 //surfaceArea.resize(numSurface);  // surface area 
 idMaxArea.resize(numSurface);    // ID of min tri area   
 idMinArea.resize(numSurface);    // ID of max tri area                            
 maxArea.resize(numSurface);      // min triangle area    
 minArea.resize(numSurface);      // max triangle area                                     
 dArea.resize(numSurface);        // delta area 
 errorArea.resize(numSurface);    // error area 
 minLength.resize(numSurface);    // min triangle length
 maxLength.resize(numSurface);    // max triangle length
 numSurfElems.resize(numSurface); // number of surface elements
 numSurfVerts.resize(numSurface); // number of surface points    
 intet.resize(numSurface);        // number of tets with 4 surface nodes

 fill(oper.begin(),oper.end(),0);
 fill(opersurf.begin(),opersurf.end(),0);
 fill(isp.begin(),isp.end(),0);
 fill(ispc.begin(),ispc.end(),0);
 fill(rsp.begin(),rsp.end(),0);
 fill(rspc.begin(),rspc.end(),0);
 fill(rspn.begin(),rspn.end(),0);
 fill(csp.begin(),csp.end(),0);
 fill(flip.begin(),flip.end(),0);
 fill(spc.begin(),spc.end(),0);
 fill(spp.begin(),spp.end(),0);
 fill(ip.begin(),ip.end(),0);
 fill(ipd.begin(),ipd.end(),0);
 fill(rp.begin(),rp.end(),0);
 fill(rpi.begin(),rpi.end(),0);
 fill(rpv.begin(),rpv.end(),0);
 fill(rpd.begin(),rpd.end(),0);
 fill(rpdist.begin(),rpdist.end(),0);
 fill(rph.begin(),rph.end(),0);

 fill(badtet.begin(),badtet.end(),0);
 fill(averageTriLength.begin(),averageTriLength.end(),0);
 fill(averageTriArea.begin(),averageTriArea.end(),0);
 fill(averageTetVolume.begin(),averageTetVolume.end(),0);
 fill(tetVol.begin(),tetVol.end(),0);
 //fill(initSurfaceVolume.begin(),initSurfaceVolume.end(),0);
 //fill(surfaceVolume.begin(),surfaceVolume.end(),0);
 fill(maxVolume.begin(),maxVolume.end(),1.0E-20);
 fill(minVolume.begin(),minVolume.end(),1.0E+20);
 fill(idMaxVolume.begin(),idMaxVolume.end(),0);
 fill(idMinVolume.begin(),idMinVolume.end(),0);
 fill(dVolume.begin(),dVolume.end(),0);
 fill(errorVolume.begin(),errorVolume.end(),0);
 //fill(initSurfaceArea.begin(),initSurfaceArea.end(),0);
 //fill(surfaceArea.begin(),surfaceArea.end(),0);
 fill(maxArea.begin(),maxArea.end(),1.0E-20);
 fill(minArea.begin(),minArea.end(),1.0E+20);
 fill(idMaxArea.begin(),idMaxArea.end(),0);
 fill(idMinArea.begin(),idMinArea.end(),0);
 fill(dArea.begin(),dArea.end(),0);
 fill(errorArea.begin(),errorArea.end(),0);
 fill(maxLength.begin(),maxLength.end(),1.0E-20);
 fill(minLength.begin(),minLength.end(),1.0E+20);
 fill(numSurfElems.begin(),numSurfElems.end(),0);
 fill(numSurfVerts.begin(),numSurfVerts.end(),0);
 fill(intet.begin(),intet.end(),0);
}

/* check and remove vertices that are too close to the surface mesh
 * structure (boundary and interface vertices)
 * In this method, besides the 3 vertices of a triangle surface mesh, 
 * the midEdges and centroid vertices (total 6
 * coordinates) are used to compare the distance to the wall/interface.
 * If it is lower than a predefinided value (if statement in the end of
 * the method), the 3D mesh vertex (vertID3D) is removed.
 * 
 * input: {X,Y,Z} and surfMesh{X,Y,Z}
 * output: mark for deletion vertID3D
 * requirements: neighbourVert
 *
 * */
void Model3D::remove3dMeshPointsByHeight()
{
 // loop at all surface mesh vertices
 for( int elem=0;elem<surfMesh.numElems;elem++ )
 {
  int v1 = surfMesh.IEN.Get(elem,0);
  int v2 = surfMesh.IEN.Get(elem,1);
  int v3 = surfMesh.IEN.Get(elem,2);

  // vertices of the triangular surface mesh
  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);
  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);
  double P3x = surfMesh.X.Get(v3);
  double P3y = surfMesh.Y.Get(v3);
  double P3z = surfMesh.Z.Get(v3);

  // centroid
  double xMid = (P1x+P2x+P3x)/3.0;
  double yMid = (P1y+P2y+P3y)/3.0;
  double zMid = (P1z+P2z+P3z)/3.0;

  // mid 12
  double xMid12 = (P1x+P2x)/2.0;
  double yMid12 = (P1y+P2y)/2.0;
  double zMid12 = (P1z+P2z)/2.0;

  // mid 13
  double xMid13 = (P1x+P3x)/2.0;
  double yMid13 = (P1y+P3y)/2.0;
  double zMid13 = (P1z+P3z)/2.0;

  // mid 23
  double xMid23 = (P2x+P3x)/2.0;
  double yMid23 = (P2y+P3y)/2.0;
  double zMid23 = (P2z+P3z)/2.0;

  // list of neighbouring points
  list<int> plist = neighbourVert.at(v1);
  for(list<int>::iterator vert=plist.begin(); vert != plist.end();++vert )
  {
   // *vert cannot be a surface mesh vertex
   if( *vert > surfMesh.numVerts )
   {
	int vertID = surfMesh.vertIdRegion.Get(v1);
	int vertID3d = vertIdRegion.Get(*vert);

	double Pxvert = X.Get(*vert);
	double Pyvert = Y.Get(*vert);
	double Pzvert = Z.Get(*vert);

	double height1 = distance(P1x,P1y,P1z,Pxvert,Pyvert,Pzvert);
	double height2 = distance(P2x,P2y,P2z,Pxvert,Pyvert,Pzvert);
	double height3 = distance(P3x,P3y,P3z,Pxvert,Pyvert,Pzvert);
	// centroid
	double height4 = distance(xMid,yMid,zMid,Pxvert,Pyvert,Pzvert);
	// xMid12 
	double height5 = distance(xMid12,yMid12,zMid12,Pxvert,Pyvert,Pzvert);
	// xMid13 
	double height6 = distance(xMid13,yMid13,zMid13,Pxvert,Pyvert,Pzvert);
	// xMid23 
	double height7 = distance(xMid23,yMid23,zMid23,Pxvert,Pyvert,Pzvert);

	double minHeight = min(height1,height2);
	minHeight = min(minHeight,height3);
	minHeight = min(minHeight,height4);
	minHeight = min(minHeight,height5);
	minHeight = min(minHeight,height6);
	minHeight = min(minHeight,height7);

	//--------------------------------------------------
	// bool pressureTest=(surfMesh.phyBounds.at(v1) == "\"wallOutflow\"" ) ||
	//                   (surfMesh.phyBounds.at(v2) == "\"wallOutflow\"" ) ||
	//                   (surfMesh.phyBounds.at(v3) == "\"wallOutflow\"" );
	//-------------------------------------------------- 

	if( minHeight < 0.7*triEdge[vertID] //&& 
	//if( minHeight < 0.4*edgeSize.Get(v1)  && 
	    //pressureTest 
	  )
	  // vertID > 0)
	{
	 mark3DPointForDeletion(*vert);
	 cout << "-----> deleting (" << *vert <<  ")" 
	      << "   " << minHeight << " < " << 0.4*triEdge[vertID] <<  endl;
	 rph[vertID3d]++;
	}
   }
  }
 }
}

/* set normal and kappa vectors using 2D plane when the interface is
 * part of the boundary. This calculation is similar to the 2D code and
 * it is applied, for instance, on the two-phase annular flow, where the
 * interface is also part of the boundary domain.
 *
 *
 *
 * */
void Model3D::setNormalAndKappa2D()
{
 vector< list<int> > neighbourLinePoint;  // 
 neighbourLinePoint.clear();

 //for( int vert=0;vert<surfMesh.numVerts;vert++ )
 for (list<int>::iterator vert=boundaryVert.begin(); 
                          vert!=boundaryVert.end(); 
						  ++vert)
 {
  // identifying which *vert is part of the interface
  list<int> auxList; 
  //if( surfMesh.Marker.Get(*vert) == 0.5 )
  if( heaviside.Get(*vert) == 0.5 )
  {
   auxList.push_back(*vert);

   // neighPoint = all neighbour vertices of *vert
   list<int> neighPoint = getNeighbourSurfacePoint(*vert);
   for (list<int>::iterator it=neighPoint.begin(); it!=neighPoint.end(); ++it)
	if( heaviside.Get(*it) == 0.5 && 
	    (Z.Get(*it) == Z.Max() || Z.Get(*it) == Z.Min()) )
	 auxList.push_back(*it);

   neighbourLinePoint.push_back( auxList );
  }
 }

 /*       v1           v0
  *         o -------- o
  *                     \
  *                      \
  *                       \
  *                         o v2
  * */
 int vectorSize = neighbourLinePoint.size();
 for( int i=0;i<vectorSize;i++ )
 {
  list<int> plist = neighbourLinePoint.at(i);
  list<int>::iterator lineVert=plist.begin(); 

  int v0 = *lineVert;++lineVert;
  int v1 = *lineVert;++lineVert; 
  int v2 = *lineVert;

//--------------------------------------------------
//   cout << v0 << " " << v1 << " " << v2 << " "
//        << surfMesh.curvature.Get(v0) << endl;
//-------------------------------------------------- 

  double P0x = surfMesh.X.Get(v0);
  double P0y = surfMesh.Y.Get(v0);
  double P0z = surfMesh.Z.Get(v0);

  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);

  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);

  // distance do ponto 0 ate ponto 1
  double a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance do ponto 0 ate ponto 2
  double b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // vetors
  double x1Unit = (P1x-P0x)/a;
  double y1Unit = (P1y-P0y)/a;

  double x2Unit = (P2x-P0x)/b;
  double y2Unit = (P2y-P0y)/b;

  double fx = x1Unit+x2Unit;
  double fy = y1Unit+y2Unit;

  // 1/2 of length P0-P1 and P0-P2
  double sumLength = (a+b)/2.0;

  /* 2D rotation of z = 90 degrees
   * x' = x*cos(z) - y*sin(z)
   * y' = x*sin(z) + y*cos(z)
   * */

  double xNormalUnit = +y1Unit*1 - y2Unit*1;
  double yNormalUnit = -x1Unit*1 + x2Unit*1;

  double len = vectorLength(xNormalUnit,yNormalUnit);

  xNormalUnit = xNormalUnit/len;
  yNormalUnit = yNormalUnit/len;

  // intensidade da forca resultante
  double force = sqrt( (fx*fx)+(fy*fy) )/sumLength;

 /* This if statement test weather the curvature (force) will be set
  * to negative or positive following the three cases based on the
  * normal componentes ({x,y,z}NormalUnit) and the computed force
  * (f{x,y,z}):
  *
  *       \              |             /
  *        \     N       |    N       /    N
  *   <---  ) --->       | --->      (  --->
  *   f    /             |            \ --->
  *       /              |  f=0        \   f
  *
  *   dotProd < 0     dotProd = 0   dotProd > 0
  *   force > 0       force = 0     force < 0
  *
  * */
 if( dotProd(fx,fy,xNormalUnit,yNormalUnit) > 0.0 )
  force = -force;

  surfMesh.curvature.Set(v0,force);
  surfMesh.xNormal.Set(v0,xNormalUnit);
  surfMesh.yNormal.Set(v0,yNormalUnit);
  //surfMesh.zNormal.Set(v0,surfMesh.zNormal.Get(v0));
  surfMesh.zNormal.Set(v0,0.0);
 }
} // fecha metodo getNormalAndKappa2D

void Model3D::setWallAnnularBC()
{    
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  double rMax = 1.0;
  if( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) < rMax*rMax - 0.001 )
   //--------------------------------------------------
   //    if( X.Get(i) < 0.5 && X.Get(i) > -0.5 &&  
   //        Y.Get(i) < 0.5 && Y.Get(i) > -0.5 ) 
   //-------------------------------------------------- 
  {
   if( surfMesh.Z.Get(i) == surfMesh.Z.Min() ) 
   {
   //idbcu.AddItem(i);
   //idbcv.AddItem(i);
   //idbcw.AddItem(i);

   //double aux = 0.0;
   //uc.Set(i,aux);
   //vc.Set(i,aux);

   wc.Set(i,1.0);
   }
   //else if( heaviside.Get(i) != 0.5 )
   else 
   {
	idbcp.AddItem(i);

	pc.Set(i,0.0);
   }
  }
  else 
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);

   double aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
 }
}

void Model3D::setWallInterfaceBC()
{    
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  double rMax = 0.5;

  if( surfMesh.Z.Get(i)==surfMesh.Z.Max() || surfMesh.Z.Get(i)==surfMesh.Z.Min() ) 
  {
   // adding vapor phase in the wall boundary which is not set by
   // default from the GMsh program
   if( surfMesh.X.Get(i)*surfMesh.X.Get(i)+
	   surfMesh.Y.Get(i)*surfMesh.Y.Get(i) < rMax*rMax - 0.001 ) 
//--------------------------------------------------
//    if( surfMesh.X.Get(i) < 0.5 && surfMesh.X.Get(i) > -0.5 &&  
//        surfMesh.Y.Get(i) < 0.5 && surfMesh.Y.Get(i) > -0.5 ) 
//-------------------------------------------------- 
   {
	//surfMesh.vertIdRegion.Set(i,1.0);
	//surfMesh.Marker.Set(i,1.0);
   }
  }
 }
}

clVector Model3D::computeConvexRegionCentroid2D(double _zPlane)
{
 double sumX = 0;
 double sumY = 0;
 //double sumZ = 0;
 int count = 0;
 //for( int vert=0;vert<surfMesh.numVerts;vert++ )
 for (list<int>::iterator vert=boundaryVert.begin(); 
                          vert!=boundaryVert.end(); 
						  ++vert)
 {
  if( heaviside.Get(*vert) == 0.5 && Z.Get(*vert) == _zPlane )
  {
   sumX += surfMesh.X.Get(*vert);
   sumY += surfMesh.Y.Get(*vert);
   //sumZ += surfMesh.Z.Get(*vert);
   count++;
  }
 }

 double xc = sumX/count;
 double yc = sumY/count;
 //double zc = sumZ/count;

 clVector centroid(2);
 centroid.Set(0,xc);
 centroid.Set(1,yc);

 return centroid;
}

void Model3D::checkLineOrientation()
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for (list<int>::iterator vert=boundaryVert.begin(); 
                          vert!=boundaryVert.end(); 
                          ++vert)
 {
  clVector centroid = computeConvexRegionCentroid2D(Z.Min());
  double xc = centroid.Get(0);
  double yc = centroid.Get(1);

  for( int elem=0;elem<surfMesh.numElems;elem++ )
  {
   int v1 = surfMesh.IEN.Get(elem,0);
   double p1x = surfMesh.X.Get(v1);
   double p1y = surfMesh.Y.Get(v1);

   int v2 = surfMesh.IEN.Get(elem,1);
   double p2x = surfMesh.X.Get(v2);
   double p2y = surfMesh.Y.Get(v2);


   /*               
	*           o ----- o 
	*         v1         v2
	* */

   // reference vector
   double vx = p1x-xc;
   double vy = p1y-yc;

   // line edge vector
   double x1 = p1x-p2x;
   double y1 = p1y-p2y;

   // Normal
   /* 2D rotation of z = 90 degrees
	* x' = x*cos(z) - y*sin(z)
	* y' = x*sin(z) + y*cos(z)
	* */
   double v1x = -y1*1;
   double v1y = +x1*1;

   if( dotProd(vx,vy,v1x,v1y) > 0.0 )
   {
	surfMesh.IEN.Set(elem,0,v2);
	surfMesh.IEN.Set(elem,1,v1);
   }
  }
 }
}

/* 
 * Bubble volume correction, print screen and save in file with
 * iteratios
 * */
void Model3D::applyBubbleVolumeCorrection()
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int nb=1;nb<=surfMesh.elemIdRegion.Max();nb++ )
 {
  stringstream ss;  //convertendo int --> string
  string str;
  ss << nb;
  ss >> str;

  string fileAux = "dat/updateVolume" + str + ".dat";
  const char* filename = fileAux.c_str();
  ifstream testFile( filename );
  ofstream file( filename,ios::app );
  if( testFile )
  {
   testFile.close();
   cout << "appending on file updateVolume" << nb << ".dat" << endl;
  }
  else
  {
   cout << "Creating file updateVolume" << nb << ".dat" << endl;
   file << "#count" << setw(19) << "volume" 
                    << setw(18) << "area" 
                    << setw(18) << "errorv" 
                    << setw(18) << "errora" 
                    << setw(18) << "dv" 
                    << setw(18) << "da" 
     			    << endl;
  }

  double aux = 0;

  dArea[nb] = (initSurfaceArea[nb] - surfaceArea[nb]);
  dVolume[nb] = (initSurfaceVolume[nb] - surfaceVolume[nb]);

  errorArea[nb] = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
  errorVolume[nb] = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);

  double TOL = initSurfaceVolume[nb]*0.00000001;

  int count = 0;
  //while( fabs(dVolume[nb]) > 1E-06 )
  //while( fabs(dArea[nb]) > 1E-06 )
  //while( fabs(dVolume[nb]) > TOL && count < 30 )
  while( fabs(errorVolume[nb]) > TOL && count < 30 )
  {
  file << setprecision(20) << scientific; 
  file << setw(10) << count << " " 
       << setw(17) << surfaceVolume[nb] << " " 
       << setw(17) << surfaceArea[nb] << " " 
       << setw(17) << errorVolume[nb] << " " 
       << setw(17) << errorArea[nb] << " " 
       << setw(17) << dVolume[nb] << " " 
       << setw(17) << dArea[nb] << " " 
       << setw(5) << setprecision(0) << fixed  
       << endl;

   for( int i=0;i<surface.Dim();i++ )
   {
	int surfaceNode = surface.Get(i);

	if( surfMesh.vertIdRegion.Get(surfaceNode) == nb )
	{
	 aux = surfMesh.X.Get(surfaceNode) + 
	       surfMesh.xNormal.Get(surfaceNode)*triEdge[nb]*errorVolume[nb];
	       //surfMesh.xNormal.Get(surfaceNode)*1.1*(dVolume[nb]/fabs(dArea[nb]));
	 X.Set(surfaceNode,aux);
	 surfMesh.X.Set(surfaceNode,aux);

	 aux = surfMesh.Y.Get(surfaceNode) + 
	       surfMesh.yNormal.Get(surfaceNode)*triEdge[nb]*errorVolume[nb];
	       //surfMesh.yNormal.Get(surfaceNode)*1.1*(dVolume[nb]/fabs(dArea[nb]));
	 Y.Set(surfaceNode,aux);
	 surfMesh.Y.Set(surfaceNode,aux);

	 aux = surfMesh.Z.Get(surfaceNode) + 
	       surfMesh.zNormal.Get(surfaceNode)*triEdge[nb]*errorVolume[nb];
	       //surfMesh.zNormal.Get(surfaceNode)*1.1*(dVolume[nb]/fabs(dArea[nb]));
	 Z.Set(surfaceNode,aux);
	 surfMesh.Z.Set(surfaceNode,aux);
	}
   }
   surfaceVolume[nb] = getSurfaceVolume(nb);
   surfaceArea[nb] = getSurfaceArea(nb);
   dArea[nb] = (initSurfaceArea[nb] - surfaceArea[nb]);
   dVolume[nb] = (initSurfaceVolume[nb] - surfaceVolume[nb]);
   errorArea[nb] = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
   errorVolume[nb] = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);

   count++;
   //cout << nb << " " << dVolume[nb] << " " << initSurfaceVolume[nb] << " " << surfaceVolume[nb] << endl;
   
  }
  file.close();

  cout << endl;
  cout << setw(20) << color(none,red,black) 
                   << "|--------- VOLUME CORRECTION ---------|" << endl;
  cout << setw(33) << color(none,white,black) << "|initial: " 
                   << initSurfaceVolume[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|final: " << surfaceVolume[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|dv: " << dVolume[nb] << endl;
  cout << setw(26) << color(none,white,black) 
                   << "volume |error: " << fabs(errorVolume[nb]) << endl;
  cout << setw(21) << color(none,red,black) 
                   << "     ---------------------------- " << endl;
  cout << setw(33) << color(none,white,black) << "|initial: " 
                   << initSurfaceArea[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|final: " << surfaceArea[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|da: " << dArea[nb] << endl;
  cout << setw(26) << color(none,white,black) 
                   << "  area |error: " << fabs(errorArea[nb]) << endl;
  cout << setw(21) << color(none,red,black) 
                   << "     ---------------------------- " << endl;
  cout << setw(28) << color(none,white,black) 
                   << "number of iterations: " << count << endl;
  cout << setw(20) << color(none,red,black) 
                   << "|-------------------------------------|" << endl;
  cout << resetColor() << endl;
 }
}

void Model3D::integralParabolic()
{
 double sumUArea = 0;
 double sumArea = 0;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);

  double P1x = surfMesh.X.Get(v1);
  double P1y = surfMesh.Y.Get(v1);
  double P1z = surfMesh.Z.Get(v1);

  double P2x = surfMesh.X.Get(v2);
  double P2y = surfMesh.Y.Get(v2);
  double P2z = surfMesh.Z.Get(v2);

  double P3x = surfMesh.X.Get(v3); 
  double P3y = surfMesh.Y.Get(v3); 
  double P3z = surfMesh.Z.Get(v3); 

  if( P1x == X.Max() &&
      P2x == X.Max() && 
      P3x == X.Max() )
  {
   double area = getArea(P1x,P1y,P1z,P2x,P2y,P2z,P3x,P3y,P3z);
   sumArea += area;
   sumUArea += area*(uc.Get(v1)+uc.Get(v2)+uc.Get(v3))/3.0;
  }
 }
 cout << sumUArea/sumArea << endl;
}

/*
 * Insert point(s) according to the solution of the diffusion equation
 * given by the class Helmholtz3D.
 *
 * */
void Model3D::insert3dMeshPointsByVolume()
{
 int vertID = 0;

 for( int elem=0;elem<IEN.DimI();elem++ )
 {
  int v1 = (int) IEN.Get(elem,0);
  int v2 = (int) IEN.Get(elem,1);
  int v3 = (int) IEN.Get(elem,2);
  int v4 = (int) IEN.Get(elem,3);

  double vol = fabs(getVolume(elem));

  int maxVert = max(v1,v2);
  maxVert = max(maxVert,v3);
  maxVert = max(maxVert,v4);

  double edgeMean = ( edgeSize.Get(v1) +
                    edgeSize.Get(v2) +
					edgeSize.Get(v3) +
					edgeSize.Get(v4) )/4.0;

  double tet = edgeMean*edgeMean*edgeMean*sqrt(2.0)/12.0;

  if( vol > 4.0*tet )
  {
   int vAdd = numVerts; // aditional vertice
   double XvAdd = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
   double YvAdd = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
   double ZvAdd = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;

   cout << "- " << color(none,blue,black) 
	            << "inserting vertex: "
				<< resetColor() << vAdd << endl;

   X.AddItem(vAdd,XvAdd);
   Y.AddItem(vAdd,YvAdd);
   Z.AddItem(vAdd,ZvAdd);
   heaviside.AddItem(vAdd,0);
   vertIdRegion.AddItem(vAdd,vertIdRegion.Get(v1));
   edgeSize.AddItem(vAdd,edgeSize.Get(v1));

   numVerts++;
   dVerts++;
   ipd[vertID]++;
  }
 }
 cout << "  inserted by diffusion: " << ipd[vertID] << endl;
}
