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
  nonSurface = _mRight.nonSurface;
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
  faceFace = _mRight.faceFace;
  freeFace = _mRight.freeFace;
  mapViz = _mRight.mapViz;
  oFace = _mRight.oFace;
  surfMesh = _mRight.surfMesh;
  interfaceMesh = _mRight.interfaceMesh;
  convexMesh = _mRight.convexMesh;
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
}

Model3D::~Model3D(){}

void Model3D::readVTK( const char* filename )
{
 char auxstr[255];
 real coords[3];
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
 real fl;
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
 real coords[3];
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
 real coords[3];
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
 // clean and init tetgen mesh object
 in.initialize();
 out.initialize();

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

 /*
  * Q: Quiet: No terminal output except errors.
  * */
 cout << endl;
 cout << "            " 
      << "|-----------------------------------------------------|" << endl;
 cout << color(blink,blue,black) 
      << "                | meshing 3D points... ";
 tetrahedralize( (char*) "Q",&in,&out );
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

 in.initialize();
 out.initialize();
}

void Model3D::setStepBC()
{
#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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

   real aux = 0.0;
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
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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
 real aux;

 clVector xCirc(1);
 clVector yCirc(1);
 xCirc.Set(0,0.0);
 yCirc.Set(0,0.0);

 real dr = 1;
 real r = dr;
 real dl = ( (2*3.141592)/nLados1Poli)*dr;
 real theta = 0;
 real dTheta = 0;
 real points2d = 1;
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

 real z = 0;
 real dz = 0;
 real points3d = 0;
 for( int jz=0;jz<nZ;jz++ )
 {
  real jzNorm = (real) jz/(nZ-1);

  //dz = jzNorm;                  // linear
  //dz = jzNorm*jzNorm;           // quadratic
  //dz = jzNorm*jzNorm*jzNorm;    // cubic
  dz = exp(jzNorm);               // exponential
   
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
  z=z+dz;
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
 tetrahedralize( (char*) "Q",&in,&out );
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
 
 // como nao ha regiao predefinida, este metodo nao funciona aqui!
 //mesh3d = convertTetgenToMesh3d(out);

 in.initialize();
 out.initialize();
}


void Model3D::mesh2Dto3D()
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
 //tetrahedralize( (char*) "QYYRCApq1.414q10a",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApqq",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApa0.1",&in,&out ); 
 //tetrahedralize( (char*) "QYYRCApa",&in,&out );
 tetrahedralize( (char*) "QYYApa",&in,&out ); // no insertion of points
 //tetrahedralize( (char*) "QYYAp",&in,&out ); // no insertion of points
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
  real x1=surfMesh.X.Get(faces[j].p1);
  real y1=surfMesh.Y.Get(faces[j].p1);
  real z1=surfMesh.Z.Get(faces[j].p1);
  real x2=surfMesh.X.Get(faces[j].p2);
  real y2=surfMesh.Y.Get(faces[j].p2);
  real z2=surfMesh.Z.Get(faces[j].p2);
  real length = vectorLength(x1-x2,y1-y2,z1-z2);

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

void Model3D::insertPointsByLength()
{
 // surfMesh.elemIdRegion == 0 --> wall
 // surfMesh.elemIdRegion == 1 --> bubble 1
 // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
 for( int edge=0;edge<mapEdgeTri.DimI();edge++ )
 {
  // edge length
  real edgeLength = mapEdgeTri.Get(edge,0);
  real v1 = mapEdgeTri.Get(edge,1);
  //real v2 = mapEdgeTri.Get(edge,2);
  real vertID = surfMesh.vertIdRegion.Get(v1); 

//--------------------------------------------------
//   real curv1 = fabs(surfMesh.curvature.Get(v1));
//   real curv2 = fabs(surfMesh.curvature.Get(v2));
//   real maxCurv = max(curv1,curv2);
//   real erro = maxCurv*edgeLength;
//-------------------------------------------------- 

  if( vertID > 0 &&
	  //Z.Get(v1) != Z.Min() && Z.Get(v2) != Z.Min() &&
	  //Z.Get(v1) != Z.Max() && Z.Get(v2) != Z.Max() &&
      //erro > 0.7 )
	  edgeLength > 1.4*triEdge[vertID] ) 
  {
   //insertSurfacePoint(edge,"flat");
   insertSurfacePoint(edge,"curvature");
   //insertSurfacePoint(edge,"bi-curvature");

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   isp[vertID]++;
   opersurf[vertID]++;
  }
 }
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
  real maxEdgeLength = 0;
  int listSize = neighbourPoint.at(surfaceNode).size();
  list<int> plist = neighbourPoint.at(surfaceNode);
  list<int>::iterator vert=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   real P0x = surfMesh.X.Get(surfaceNode);
   real P0y = surfMesh.Y.Get(surfaceNode);
   real P0z = surfMesh.Z.Get(surfaceNode);

   int v1 = *vert;++vert;
   real P1x = surfMesh.X.Get(v1);
   real P1y = surfMesh.Y.Get(v1);
   real P1z = surfMesh.Z.Get(v1);

   real edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

   if( edgeLength > maxEdgeLength )
	maxEdgeLength = edgeLength;
  }

  // Still don't know about the curv value!!!
  // edge length
  real vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
  real curv = fabs(surfMesh.curvature.Get(surfaceNode));

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
//   real vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
//   real curv = fabs(surfMesh.curvature.Get(surfaceNode));
//   real edgeLength = triEdge[vertID];
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

void Model3D::insertPointsByCurvature()
{
 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge length
  //real edgeLength = mapEdgeTri.Get(i,0);
  int v1 = mapEdgeTri.Get(i,1);
  int v2 = mapEdgeTri.Get(i,2);
  real curv1 = fabs(surfMesh.curvature.Get(v1));
  real curv2 = fabs(surfMesh.curvature.Get(v2));
  real vertID = surfMesh.vertIdRegion.Get(v1); 
  real length = mapEdgeTri.Get(i,0);
  //real erro = max(curv1,curv2);

  if( curv1*length > 2.5 || curv2*length > 2.5  )
  {
   cout << "----------- Inserting by curvature..." << endl;
   insertSurfacePoint(i,"flat");
   //insertSurfacePoint(i,"curvature");
   //insertSurfacePoint(i,"bi-curvature");

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   ispc[vertID]++;
   opersurf[vertID]++;
  }
 }
}

void Model3D::insertPointsByInterfaceDistance()
{
 real Ymax1=100;
 real Ymin1=-100;
 real Ymax2=-100;
 real Ymin2=100;

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
 real distMin = Ymin2-Ymin1;
 //real distMax = Ymax2-Ymax1;
 real dy = (Ymin2-Ymin1)/(ny-1);

 for( int i=0;i<mapEdgeTri.DimI();i++ )
 {
  // edge length
  real edgeLength = mapEdgeTri.Get(i,0);
  real aux = (distMin/2.0)+2*dy;
  //real aux = (distMin/2.0)+(distMax-distMin)/3.0;
  if( surfMesh.Marker.Get(mapEdgeTri.Get(i,1)) == 0.5 && 
	  surfMesh.Y.Get(mapEdgeTri.Get(i,1)) < 1.0*aux && 
	  surfMesh.Y.Get(mapEdgeTri.Get(i,1)) > -1.0*aux &&
	  edgeLength > 4*dy )
  {
   insertSurfacePoint(i,"flat");
   //insertSurfacePoint(i,"curvature");
   //insertSurfacePoint(i,"bi-curvature");
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
 real qMax = 0; 
 clVector bestTri(3);

 for( int i=0;i<listSize-1;i++ )
 {
  vert1 = *point;++point;
  vert2 = *point;++point;
  vert3 = *point;--point;

  real q = triangleQualityMeasure(vert1,vert2,vert3);

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
void Model3D::flipTriangleEdge()
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
  real edge = i;
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

  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);

  real P2x = surfMesh.X.Get(v2);
  real P2y = surfMesh.Y.Get(v2);
  real P2z = surfMesh.Z.Get(v2);

  real P3elem1x = surfMesh.X.Get(v3elem1);
  real P3elem1y = surfMesh.Y.Get(v3elem1);
  real P3elem1z = surfMesh.Z.Get(v3elem1);

  real P3elem2x = surfMesh.X.Get(v3elem2);
  real P3elem2y = surfMesh.Y.Get(v3elem2);
  real P3elem2z = surfMesh.Z.Get(v3elem2);

  // elem1
  real area1 = getAreaVert(v1,v3elem1,v2);
  real q1 = triangleQualityMeasure(v1,v3elem1,v2);
  real c1 = getCircumRadius(P1x,P1y,P1z,
	                        P3elem1x,P3elem1y,P3elem1z,
							P2x,P2y,P2z);
  clVector normalElem1 = getNormalElem(elem1);

  // elem2
  real area2 = getAreaVert(v2,v3elem2,v1);
  real q2 = triangleQualityMeasure(v2,v3elem2,v1);
  real c2 = getCircumRadius(P2x,P2y,P2z,
	                        P3elem2x,P3elem2y,P3elem2z,
							P1x,P1y,P1z);
  clVector normalElem2 = getNormalElem(elem2);

  // elem3 (new elem1)
  real area3 = getAreaVert(v1,v3elem1,v3elem2);
  real q3 = triangleQualityMeasure(v1,v3elem1,v3elem2);
  real c3 = getCircumRadius(P1x,P1y,P1z,
	                        P3elem1x,P3elem1y,P3elem1z,
	                        P3elem2x,P3elem2y,P3elem2z);
  clVector normalElem3 = getNormalElem(v1,v3elem1,v3elem2);

  // elem4 (new elem2)
  real area4 = getAreaVert(v2,v3elem2,v3elem1);
  real q4 = triangleQualityMeasure(v2,v3elem2,v3elem1);
  real c4 = getCircumRadius(P2x,P2y,P2z,
	                        P3elem2x,P3elem2y,P3elem2z,
	                        P3elem1x,P3elem1y,P3elem1z);
  clVector normalElem4 = getNormalElem(v2,v3elem2,v3elem1);

//--------------------------------------------------
//   real angleOld = angle3D(normalElem1.Get(0),normalElem1.Get(1),
// 	                      normalElem1.Get(2),
// 						  normalElem2.Get(0),normalElem2.Get(1),
// 						  normalElem2.Get(2) );
// 
//   real angleNew = angle3D(normalElem3.Get(0),normalElem3.Get(1),
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
//   real v12x = P2x-P1x;
//   real v12y = P2y-P1y;
//   real v12z = P2z-P1z;
// 
//   real v13elem1x = P3elem1x-P1x;
//   real v13elem1y = P3elem1y-P1y;
//   real v13elem1z = P3elem1z-P1z;
// 
//   real v13elem2x = P3elem2x-P1x;
//   real v13elem2y = P3elem2y-P1y;
//   real v13elem2z = P3elem2z-P1z;
// 
//   real angle1 = angle3D(v12x,v12y,v12z,
// 	                    v13elem1x,v13elem1y,v13elem1z);
//   real angle2 = angle3D(v12x,v12y,v12z,
// 	                    v13elem2x,v13elem2y,v13elem2z);
// 
//   real v21x = P1x-P2x;
//   real v21y = P1y-P2y;
//   real v21z = P1z-P2z;
// 
//   real v23elem1x = P3elem1x-P2x;
//   real v23elem1y = P3elem1y-P2y;
//   real v23elem1z = P3elem1z-P2z;
// 
//   real v23elem2x = P3elem2x-P2x;
//   real v23elem2y = P3elem2y-P2y;
//   real v23elem2z = P3elem2z-P2z;
// 
//   real angle3 = dotProd(v21x,v21y,v21z,
// 	                    v23elem1x,v23elem1y,v23elem1z);
//   real angle4 = dotProd(v21x,v21y,v21z,
// 	                    v23elem2x,v23elem2y,v23elem2z);
//-------------------------------------------------- 

  // this works, but is not consistent!!! CHANGE IT SOON!
  //real curv1 = fabs(surfMesh.curvature.Get(v1));
  //real curv2 = fabs(surfMesh.curvature.Get(v2));
  //real curv3_1 = fabs(surfMesh.curvature.Get(v3elem1));
  //real curv3_2 = fabs(surfMesh.curvature.Get(v3elem2));

  /* FLIPPING requirements:
   * - sum of quality of old triangles < sum of quality of new triangles
   * - curvature of all the 4 vertices < 40
   * - sum of area of old triangles > sum of area of new triangles
   * - angle between plane of old triangles > 90 degrees
   * - sum of circumcenters of old triangles > sum of circumcenters of
   *   new triangles
   * */
  if( surfMesh.Marker.Get(v1) == 0.5 && 
	q1+q2 < q3+q4 && // quality sum
	//--------------------------------------------------
	// (curv1 < 40 && curv2 < 40) && // curvature
	// (curv3_1 < 40 && curv3_2 < 40) && // curvature
	//-------------------------------------------------- 
	area1+area2  > area3+area4 && // area sum
	//angleNew < angleOld &&
    neighbourSurfaceElem.at( v1 ).size() > 4 &&  
    neighbourSurfaceElem.at( v2 ).size() > 4 &&  
	c1+c2 > c3+c4 ) // circum radius
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
void Model3D::insertSurfacePoint(int _edge,const char* _mode)
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

 real XvAdd = 0.0;
 real YvAdd = 0.0;
 real ZvAdd = 0.0;

 if( strcmp( _mode,"bi-curvature") == 0 ) 
 {
  clVector coordAdd1 = fitEllipse( X.Get(v1),Y.Get(v1),Z.Get(v1),
	                              X.Get(v2),Y.Get(v2),Z.Get(v2),
								  surfMesh.xNormal.Get(v1),
								  surfMesh.yNormal.Get(v1),
								  surfMesh.zNormal.Get(v1),
								  surfMesh.xNormal.Get(v2),
								  surfMesh.yNormal.Get(v2),
								  surfMesh.zNormal.Get(v2) );

  clVector coordAdd2 = fitEllipse( X.Get(v3elem1),Y.Get(v3elem1),Z.Get(v3elem1),
	                              X.Get(v3elem2),Y.Get(v3elem2),Z.Get(v3elem2),
								  surfMesh.xNormal.Get(v3elem1),
								  surfMesh.yNormal.Get(v3elem1),
								  surfMesh.zNormal.Get(v3elem1),
								  surfMesh.xNormal.Get(v3elem2),
								  surfMesh.yNormal.Get(v3elem2),
								  surfMesh.zNormal.Get(v3elem2) );

  XvAdd = (coordAdd1.Get(0)+coordAdd2.Get(0))*0.5;
  YvAdd = (coordAdd1.Get(1)+coordAdd2.Get(1))*0.5;
  ZvAdd = (coordAdd1.Get(2)+coordAdd2.Get(2))*0.5;
 }
 else if( strcmp( _mode,"curvature") == 0 ) 
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
//real curv = (surfMesh.curvature.Get(v1)+surfMesh.curvature.Get(v2))*0.5;
//  surfMesh.curvature.AddItem(curv);
//  curvature.AddItem(vAdd,curv);
//-------------------------------------------------- 

 clVector myVec = getNormalAndKappaByDesbrun(vAdd,
                    getNeighbourSurfacePoint(vAdd));
 surfMesh.curvature.AddItem(vAdd,myVec.Get(0));
 curvature.AddItem(vAdd,myVec.Get(0));

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
void Model3D::contractEdgeByLength()
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

  real curv1 = fabs(surfMesh.curvature.Get(v1));
  real curv2 = fabs(surfMesh.curvature.Get(v2));
  real curv3 = fabs(surfMesh.curvature.Get(v3elem1));
  real curv4 = fabs(surfMesh.curvature.Get(v3elem2));
  int elem1 = mapEdgeTri.Get(edge,5);
  int elem2 = mapEdgeTri.Get(edge,6);
  real edgeLength = mapEdgeTri.Get(edge,0);
  clVector normalElem1 = getNormalElem(elem1);
  clVector normalElem2 = getNormalElem(elem2);

//--------------------------------------------------
//   // radii
//   real angle = angle3D(normalElem1.Get(0),normalElem1.Get(1),
// 	                   normalElem1.Get(2),
// 					   normalElem2.Get(0),normalElem2.Get(1),
// 					   normalElem2.Get(2) );
//-------------------------------------------------- 

  // verifying the length of each surface edge
  int elemID = surfMesh.elemIdRegion.Get(mapEdgeTri.Get(edge,5));

  //real minCurv = min(curv1,curv2);
  //real erro = minCurv*edgeLength;
  
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
  if( edgeLength < 0.5*triEdge[elemID] &&
      //erro < 0.03 &&
	  elemIDTest && curvTest && neighTest //&& angleTest 
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

   const char* _mode = "flat";

   if( strcmp( _mode,"curvature") == 0 ) 
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

	real XvAdd = coordAdd.Get(0);
	real YvAdd = coordAdd.Get(1);
	real ZvAdd = coordAdd.Get(2);
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else if( strcmp( _mode,"bi-curvature") == 0 ) 
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

	real XvAdd = (coordAdd1.Get(0)+coordAdd2.Get(0))*0.5;
	real YvAdd = (coordAdd1.Get(1)+coordAdd2.Get(1))*0.5;
	real ZvAdd = (coordAdd1.Get(2)+coordAdd2.Get(2))*0.5;
	surfMesh.X.Set(v1, XvAdd );
	surfMesh.Y.Set(v1, YvAdd );
	surfMesh.Z.Set(v1, ZvAdd );
	X.Set(v1, XvAdd );
	Y.Set(v1, YvAdd );
	Z.Set(v1, ZvAdd );
   }
   else // flat
   {
	// flat
	real XvNew = ( surfMesh.X.Get(v1)+ surfMesh.X.Get(v2) )*0.5;
	real YvNew = ( surfMesh.Y.Get(v1)+ surfMesh.Y.Get(v2) )*0.5;
	real ZvNew = ( surfMesh.Z.Get(v1)+ surfMesh.Z.Get(v2) )*0.5;
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
   clVector myVec = getNormalAndKappaByDesbrun(v1,getNeighbourSurfacePoint(v1));
   surfMesh.curvature.Set(v1,myVec.Get(0));
   curvature.Set(v1,myVec.Get(0));

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
  }
 }
}

void Model3D::removePointsByLength()
{
  for( int i=0;i<mapEdgeTri.DimI();i++ )
  {
   // edge vertices
   int v1 = mapEdgeTri.Get(i,1);
   int v2 = mapEdgeTri.Get(i,2);
   int vertID = surfMesh.vertIdRegion.Get(v1);

   // verifying the length of each surface edge
   if( mapEdgeTri.Get(i,0) < 0.2*triEdge[vertID] ) //&&
   {
	// sum of all neighbour edge length of the 1st. point
	real sumLength1=0;
	int listSize1 = neighbourPoint.at(v1).size();
	list<int> plist1 = neighbourPoint.at(v1);
	list<int>::iterator vert1=plist1.begin();
	for( int i=0;i<listSize1-1;i++ )
	{
	 real P0x = surfMesh.X.Get(v1);
	 real P0y = surfMesh.Y.Get(v1);
	 real P0z = surfMesh.Z.Get(v1);

	 int v = *vert1;++vert1;
	 real P1x = surfMesh.X.Get(v);
	 real P1y = surfMesh.Y.Get(v);
	 real P1z = surfMesh.Z.Get(v);

     sumLength1 += distance(P0x,P0y,P0z,P1x,P1y,P1z);
	}

	// sum of all neighbour edge length of the 1st. point
	real sumLength2=0;
	int listSize2 = neighbourPoint.at(v2).size();
	list<int> plist2 = neighbourPoint.at(v2);
	list<int>::iterator vert2=plist2.begin();
	for( int i=0;i<listSize2-1;i++ )
	{
	 real P0x = surfMesh.X.Get(v2);
	 real P0y = surfMesh.Y.Get(v2);
	 real P0z = surfMesh.Z.Get(v2);

	 int v = *vert2;++vert2;
	 real P2x = surfMesh.X.Get(v);
	 real P2y = surfMesh.Y.Get(v);
	 real P2z = surfMesh.Z.Get(v);

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

 real triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));

 real h = triEdgeMin*0.8864; 
 for( int i=0;i<numVerts;i++ )
 {
  int vertID = vertIdRegion.Get(i);
  real d = interfaceDistance.Get(i);

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

	real d = distance( X.Get(i),Y.Get(i),Z.Get(i),
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
  real XvAdd = mapEdge.Get(e,1);
  real YvAdd = mapEdge.Get(e,2);
  real ZvAdd = mapEdge.Get(e,3);

  int v1 = mapEdge.Get(e,4);
  int v2 = mapEdge.Get(e,5);

  int vertID = vertIdRegion.Get(v1);

  real x1=X.Get(v1);
  real y1=Y.Get(v1);
  real z1=Z.Get(v1);
  real x2=X.Get(v2);
  real y2=Y.Get(v2);
  real z2=Z.Get(v2);
  real length = distance(x1,y1,z1,x2,y2,z2);

  // minVert should be bigger then surfMesh.numVerts because we are
  // treating only 3D vertices after the surface mesh vertices.
  //int minVert = min(v1,v2);
  int maxVert = max(v1,v2);
  real maxEdge = max(edgeSize.Get(v1),edgeSize.Get(v2));
  //real hSum = heaviside.Get(v1) + heaviside.Get(v2);

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

  real x1=X.Get(v1);
  real y1=Y.Get(v1);
  real z1=Z.Get(v1);
  real x2=X.Get(v2);
  real y2=Y.Get(v2);
  real z2=Z.Get(v2);
  real length = distance(x1,y1,z1,x2,y2,z2);

  //int maxVert = max(v1,v2);
  int minVert = min(v1,v2);

  //real size = (edgeSize.Get(v1)+edgeSize.Get(v2))/2.0;
  real size = min(edgeSize.Get(v1),edgeSize.Get(v2));

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
 real test = 0.008;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);
  if( getAreaElem(i) > test )
  {
   int v4 = surfMesh.numVerts;

   // coords of triangle centroid
   real centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3) )/3.0;
   real centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3) )/3.0;
   real centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3) )/3.0;

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

 real Xmax1=0.0;  real Ymax1=0.0; 
 real Xmin1=0.0;  real Ymin1=Y.Min();     
 real Zmax1=0.0;  real Ymax2=0.0; 
 real Zmin1=0.0;  real Ymin2=Y.Max();     

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
 real xi = Xmin1;
 real yi = Ymin1;
 real zi = Zmin1;

 // distance between points
 real dx = (Xmax1-Xmin1)/(nPoints-1);
 real dy = (Ymin2-Ymin1)/(ny+1);
 real dz = (Zmax1-Zmin1)/(nPoints-1);

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
void Model3D::mesh2Dto3DOriginal()
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
 tetrahedralize( (char*) "QYYRCApq1.414q10a",&in,&out );
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
  tetrahedralize( (char*) "QYYRCApq1.414q10a",&in,&out );
  cout << "finished <---- " << endl;;

  convertTetgenToModel3D(out);
  mesh3d = convertTetgenToMesh3d(out);
 }
 cout << endl;

 in.initialize();
 out.initialize();
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
  //real curv = fabs(surfMesh.curvature.Get(node));
  // find the first vertex with region == nb
  //for( int i=0;i<surfMesh.numVerts;i++ )
  for( int i=surfMesh.numVerts-1;i>=0;i-- )
  {
   clVector myVec = getNormalAndKappaByDesbrun(i,getNeighbourSurfacePoint(i));
   real curv = fabs(myVec.Get(0));
  //for( int i=0;i<surfMesh.numVerts;i++ )
   if( surfMesh.vertIdRegion.Get(i) == nb && curv < 20 )
   //if( surfMesh.vertIdRegion.Get(i) == nb )
   {
	node = i;
	break;
   }
  }

  clVector myVec = getNormalAndKappaByDesbrun(node,getNeighbourSurfacePoint(node));
  real xIn = surfMesh.X.Get(node)-0.1*triEdge[nb]*myVec.Get(1);
  real yIn = surfMesh.Y.Get(node)-0.1*triEdge[nb]*myVec.Get(2);
  real zIn = surfMesh.Z.Get(node)-0.1*triEdge[nb]*myVec.Get(3);

  real edge = triEdge[nb];
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
  real mSum = surfMesh.Marker.Get(v1) + 
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
  real mSum = surfMesh.Marker.Get(v1) + 
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
   real xMid = ( _tetmesh.pointlist[3*v1+0]+
                 _tetmesh.pointlist[3*v2+0]+
				 _tetmesh.pointlist[3*v3+0]+
				 _tetmesh.pointlist[3*v4+0] )/4.0;

   real yMid = ( _tetmesh.pointlist[3*v1+1]+
                 _tetmesh.pointlist[3*v2+1]+
				 _tetmesh.pointlist[3*v3+1]+
				 _tetmesh.pointlist[3*v4+1] )/4.0;

   real zMid = ( _tetmesh.pointlist[3*v1+2]+
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
void Model3D::removePointByVolume()
{
 real vSum;
 real vertSum;
 int v[NUMGLE];
 int vert=0;

 tetVol.clear();
 tetVol.resize(triEdge.size()); // number of surfaces + 1
 real triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));

 // set tetVol ---> wall,bubble1, bubble2 etc.
 for( int v=0;v<(int) triEdge.size();v++ )
  tetVol[v] = triEdgeMin*triEdgeMin*triEdgeMin*sqrt(2.0)/12.0;

 for( int elem=0;elem<numElems;elem++ )
 {
  v[0] = IEN.Get(elem,0);
  v[1] = IEN.Get(elem,1);
  v[2] = IEN.Get(elem,2);
  v[3] = IEN.Get(elem,3);

  real hSum = heaviside.Get(v[0])+heaviside.Get(v[1])+
              heaviside.Get(v[2])+heaviside.Get(v[3]);

  int elemID = elemIdRegion.Get(elem);

  real maxVol = max(1.0E-05,0.01*tetVol[elemIdRegion.Get(elem)]);

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
 //cout << "  removed by volume: " << rpv[elemID] << endl;
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
 vector<real> sumLength;sumLength.clear();
 sumLength.resize(surfMesh.elemIdRegion.Max()+1);
 vector<real> sumArea;sumArea.clear();
 sumArea.resize(surfMesh.elemIdRegion.Max()+1);
 for( int e=0;e<surfMesh.numElems;e++ )
 {
  int v1 = surfMesh.IEN.Get(e,0);
  real p1x = surfMesh.X.Get(v1); 
  real p1y = surfMesh.Y.Get(v1); 
  real p1z = surfMesh.Z.Get(v1); 

  int v2=surfMesh.IEN.Get(e,1);
  real p2x = surfMesh.X.Get(v2); 
  real p2y = surfMesh.Y.Get(v2); 
  real p2z = surfMesh.Z.Get(v2); 

  int v3=surfMesh.IEN.Get(e,2);
  real p3x = surfMesh.X.Get(v3); 
  real p3y = surfMesh.Y.Get(v3); 
  real p3z = surfMesh.Z.Get(v3); 

  int elemID = surfMesh.elemIdRegion.Get(e);

  real length12 = distance(p1x,p1y,p1z,p2x,p2y,p2z);
  real length13 = distance(p1x,p1y,p1z,p3x,p3y,p3z);
  real length23 = distance(p2x,p2y,p2z,p3x,p3y,p3z);

  sumLength[elemID] += ( length12+length13+length23 )/3.0;

  real auxMinLength = min( length12,length13 );
  auxMinLength = min( auxMinLength,length23 );

  real auxMaxLength = max( length12,length13 );
  auxMaxLength = max( auxMaxLength,length23 );

  real area = fabs(getAreaElem(e));
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

 vector<real> sumVolume;sumVolume.clear();
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
  real volume = fabs(getVolume(e));

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
void Model3D::mesh3DPoints()
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
 //tetrahedralize( (char*) "QYYRCApq1.414q10a",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApqq10a",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYRCApa",&in,&out );
 //tetrahedralize( (char*) "QYYCApa0.5",&in,&out ); 
 tetrahedralize( (char*) "QYYApa",&in,&out ); 
 //tetrahedralize( (char*) "QYYRCApqq10",&in,&out ); // quality
 //tetrahedralize( (char*) "QYYApaq",&in,&out ); // 
 //tetrahedralize( (char*) "QYYAp",&in,&out ); // no insertion of points
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

void Model3D::setDiskCouetteBC()
{
 real aux;
 rMax = Y.Max(); 

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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

void Model3D::setNuCteDiskBC()
{
 real omega,aux;
 real radius;
 rMax = Y.Max();

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
 {
  radius = sqrt( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i) );

  if( Z.Get(i) == Z.Max() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   //idbcp.AddItem(i); // caso com c.c. livre em w
   idbcw.AddItem(i);

   //uc.Set(i,radius*2.5666593e-02); // Z=4
   //vc.Set(i,radius*3.4943977e-02); // Z=4
   //wc.Set(i,-8.2505646e-01); // Z=4

   //uc.Set(i,radius*4.5487756e-03); // Z=6
   //vc.Set(i,radius*5.9598499e-03); // Z=6
   //wc.Set(i,-8.7414071e-01); // Z=6

   uc.Set(i,radius*1.3326987e-04); // Z=10
   vc.Set(i,radius*1.7327920e-04); // Z=10
   wc.Set(i,-8.8416563E-01); // Z=10
   
   //uc.Set(i,0.0); // Z=12
   //vc.Set(i,0.0); // Z=12
   //pc.Set(i,-0.391141); // caso com c.c. livre em w
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

 }

 for( int i=0;i<numVerts;i++ )
 {
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

void Model3D::setNuCDiskBC()
{
 real omega,aux,radius;
 rMax = Y.Max();

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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
   outflow.Set(i,aux);

 }
 for( int i=0;i<numVerts;i++ )
 {
  if( Z.Get(i)<Z.Max() && Z.Get(i)>Z.Min() && 
	( X.Get(i)*X.Get(i)+Y.Get(i)*Y.Get(i)>rMax*rMax - 0.001) )
  {
   idbcp.AddItem(i);
   aux = 0.0;
   pc.Set(i,aux);
  }
 }
}

void Model3D::setFiniteDiskBC()
{
 real omega,aux;
 rMax = Y.Max();

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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

void Model3D::setCubeBC()
{
#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
#endif

 for( int i=0;i<numBCPoints;i++ )
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

void Model3D::setGenericBC()
{    
 /* This IF selects the priority boundary conditions to be set on the
  * phyBounds vector. Note that only these names will be written on top
  * of the others. For instance if a corner point has 2 types of
  * boundary condition, the phyNames below will be written.
  * */ 
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
 real diameterXY = ( dist(X.Min(),X.Max()) + 
                     dist(Y.Min(),Y.Max()) ) / 2.0;
 real diameterXZ = ( dist(X.Min(),X.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
 real diameterYZ = ( dist(Y.Min(),Y.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
//--------------------------------------------------
//  real diameterYZ = distance(Y.Min(),Z.Min(),Y.Max(),Z.Max());
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

   real radius = sqrt( Y.Get(*it)*Y.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   real Umax = 1.0;
   real aux = 1.0*Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
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

   real radius = sqrt( X.Get(*it)*X.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   real Vmax = 2.0;
   real aux = Vmax*( 1.0-radius*radius/((diameterXZ/2.0)*
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

   real radius = sqrt( X.Get(*it)*X.Get(*it) + Y.Get(*it)*Y.Get(*it) );

   // Parabolic profile
   real Wmax = 2.0;
   real aux = Wmax*( 1.0-radius*radius/((diameterXY/2.0)*
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

   real aux = X.Get(*it);
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

   real aux = X.Get(*it);
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

void Model3D::setGenericBC(real _vel)
{    
 clearBC();

 // calculating channel's diameter.
 real diameterXY = ( dist(X.Min(),X.Max()) + 
                     dist(Y.Min(),Y.Max()) ) / 2.0;
 real diameterXZ = ( dist(X.Min(),X.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
 real diameterYZ = ( dist(Y.Min(),Y.Max()) + 
                     dist(Z.Min(),Z.Max()) ) / 2.0;
//--------------------------------------------------
//  real diameterYZ = distance(Y.Min(),Z.Min(),Y.Max(),Z.Max());
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

   real radius = sqrt( Y.Get(*it)*Y.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   real Umax = 1.0;
   real aux = Umax*( 1.0-radius*radius/((diameterYZ/2.0)*
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

   real radius = sqrt( X.Get(*it)*X.Get(*it) + Z.Get(*it)*Z.Get(*it) );

   // Parabolic profile
   real Vmax = 2.0;
   real aux = Vmax*( 1.0-radius*radius/((diameterXZ/2.0)*
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

   real radius = sqrt( X.Get(*it)*X.Get(*it) + Y.Get(*it)*Y.Get(*it) );

   // Parabolic profile
   real Wmax = 2.0;
   real aux = Wmax*( 1.0-radius*radius/((diameterXY/2.0)*
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

  real aux = 0.0;
  uc.Set(*it,aux);
  vc.Set(*it,aux);
  wc.Set(*it,aux);
 }
}

void Model3D::setWallBC(real _vel)
{    
 clearBC();

 for (list<int>::iterator it=boundaryVert.begin(); it!=boundaryVert.end(); ++it)
 {
  idbcu.AddItem(*it);
  idbcv.AddItem(*it);
  idbcw.AddItem(*it);

  real aux = 0.0;
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

   real aux = 0.0;
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

void Model3D::setMovingWallBC(real _vel)
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
 real radius = 0.5;
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
 real aux = 0;
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
 real aux;

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

void Model3D::setCube(real _lim1,real _lim2,real _eps)
{
 real eps2 = _eps;
 real lim1 = _lim1;
 real lim2 = _lim2;

 heaviside.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  // na interface
  if( (X.Get(i)<lim2+eps2) && (X.Get(i)>lim1-eps2) && 
      (Y.Get(i)<lim2+eps2) && (Y.Get(i)>lim1-eps2) && 
      (Z.Get(i)<lim2+eps2) && (Z.Get(i)>lim1-eps2) )
  {
   heaviside.Set(i,0.5);
  }
  // dentro da bolha
  if( (X.Get(i)<lim2-eps2) && (X.Get(i)>lim1+eps2) && 
      (Y.Get(i)<lim2-eps2) && (Y.Get(i)>lim1+eps2) && 
      (Z.Get(i)<lim2-eps2) && (Z.Get(i)>lim1+eps2) )
  {
   heaviside.Set(i,1.0);
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

 heaviside.Dim(numVerts);
 for( int i=0;i<numVerts;i++ )
 {
  // na interface
  if( (X.Get(i)<xlimSup+eps2) && (X.Get(i)>xlimInf-eps2) && 
      (Y.Get(i)<ylimSup+eps2) && (Y.Get(i)>ylimInf-eps2) && 
      (Z.Get(i)<zlimSup+eps2) && (Z.Get(i)>zlimInf-eps2) )
  {
   heaviside.Set(i,0.5);
  }
  // dentro da bolha
  if( (X.Get(i)<xlimSup-eps2) && (X.Get(i)>xlimInf+eps2) && 
      (Y.Get(i)<ylimSup-eps2) && (Y.Get(i)>ylimInf+eps2) && 
      (Z.Get(i)<zlimSup-eps2) && (Z.Get(i)>zlimInf+eps2) )
  {
   heaviside.Set(i,1.0);
  }
 }
}

void Model3D::setNuZDiskBC()
{
 real omega,aux,radius;
 rMax = Y.Max();

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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

// FSBC - Free Surface Boundary Condition
void Model3D::setDiskFSBC()
{
 real aux;
 heaviside.Dim(numVerts);
 rMax = Y.Max();

#if NUMGLEU == 5
 real numBCPoints = numVerts;
#else
 real numBCPoints = numNodes;
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

   real aux = 1.0;
   cc.Set(i,aux);
  }
 }
}

void Model3D::setAdimenDisk()
{
 real aux;
 real Red = 100;
 real factorz = 1.0/(Z.Max()-Z.Min());
 real xrMax = X.Max();
 real yrMax = Y.Max();

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
   aux = Z.Get(i)*( 1.0 + (X.Get(i)/Red)*factorz*0.3 );
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
void Model3D::setSphereToEllipsoid(real _factor)
{
 real aux;
 real dz = _factor;
 real dx = sqrt(1.0/dz);

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
void Model3D::setBiggerSphere(real _factor)
{
 real aux;
 real ds = _factor;

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
  real centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
  real centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
  real centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;
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

  real centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
  real centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
  real centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;

  X.Set(v5,centroidX);
  Y.Set(v5,centroidY);
  Z.Set(v5,centroidZ);
 }
}

void Model3D::edgeMidPointPositionCorrection()
{
 real xMid,yMid,zMid;

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
  real volume = getVolume(i);
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

// criando os nos na metade de cada aresta. Para isso eh necessario
// fazer uma lista de arestas e numera-las de forma que uma aresta comum
// tenha apenas 1 numero e seja compartilhada em todos os elementos que
// tem a aresta. 
/*
 * mapEdge.Set(edge,0,numVerts+edge); // numero da aresta
 * mapEdge.Set(edge,1,xMid ); // coordenada X do centro da aresta
 * mapEdge.Set(edge,2,yMid ); // coordenada Y do centro da aresta
 * mapEdge.Set(edge,3,zMid ); // coordenada Y do centro da aresta
 * mapEdge.Set(edge,4,faces[i].p1 ); // 1o noh
 * mapEdge.Set(edge,5,faces[i].p2 ); // 2o noh
 *
 * */
void Model3D::setMapEdge()
{
 int numFace = 6; // teraedro tem 6 arestas
 clVector faceaux(2);
 IFACE2D *faces = NULL;
 int listSize = numFace*numElems;
 faces = new IFACE2D[listSize];

 checkTetrahedronOrientation();

 for( int mele=0;mele<numElems;mele++ )
 {
  int v1 = (int) IEN.Get(mele,0);
  int v2 = (int) IEN.Get(mele,1);
  int v3 = (int) IEN.Get(mele,2);
  int v4 = (int) IEN.Get(mele,3);

  // 1st. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v2);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+0].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+0].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+0].p3 = mele;

  // 2nd. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+1].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+1].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+1].p3 = mele;

  // 3rd. edge
  faceaux.Set(0,v1);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+2].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+2].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+2].p3 = mele;

  // 4th. edge
  faceaux.Set(0,v2);
  faceaux.Set(1,v3);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+3].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+3].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+3].p3 = mele;

  // 5th. edge
  faceaux.Set(0,v2);
  faceaux.Set(1,v4);
  faceaux.Sort(); // para ordenar os vertices de uma aresta
  faces[numFace*mele+4].p1 = (int) faceaux.Get(0);
  faces[numFace*mele+4].p2 = (int) faceaux.Get(1);
  faces[numFace*mele+4].p3 = mele;

  // 6th. edge
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

//--------------------------------------------------
//  for( int i=0;i<listSize;i++ )
//   cout << faces[i].p1 << " "
//        << faces[i].p2 << " "
//        << faces[i].p3 << endl;
//-------------------------------------------------- 

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
 neighbourEdge.resize (listSize);
 mapEdge.Dim(listSize,6);

 // numeracao de arestas a partir de numVerts e associacao das arestas
 // aos elementos para inclusao na IEN
 // OBS.: the for loop is from 0 to listSize-1. This can be done because
 // the last line of faces (faces[listSize-1]) is always a repetion of
 // the previous one (faces[listSize-2]).
 for( int i=0;i<listSize-1;i++ )
 {
  real x1=X.Get(faces[i].p1);
  real y1=Y.Get(faces[i].p1);
  real z1=Z.Get(faces[i].p1);
  real x2=X.Get(faces[i].p2);
  real y2=Y.Get(faces[i].p2);
  real z2=Z.Get(faces[i].p2);

  real xMid = (x1+x2)*0.5;
  real yMid = (y1+y2)*0.5;
  real zMid = (z1+z2)*0.5;
  real length = vectorLength(x1-x2,y1-y2,z1-z2);

  // checking the minimum edge length
  if( length < minEdge )
   minEdge = length;

  mapEdge.Set(edge,0,numVerts+edge); // numero da aresta
  mapEdge.Set(edge,1,xMid ); // coordenada X do centro da aresta
  mapEdge.Set(edge,2,yMid ); // coordenada Y do centro da aresta
  mapEdge.Set(edge,3,zMid ); // coordenada Y do centro da aresta
  mapEdge.Set(edge,4,faces[i].p1 ); // 1o noh
  mapEdge.Set(edge,5,faces[i].p2 ); // 2o noh

  neighbourEdge.at(edge).push_back(faces[i].p3);

  // faces[listSize]
  while( (i<listSize-1 ) &&
         (faces[i].p1 == faces[i+1].p1) &&
         (faces[i].p2 == faces[i+1].p2) )
		 
  {
   neighbourEdge.at(edge).push_back(faces[i+1].p3);
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
  real xc    = mapEdge.Get(i,1);
  real yc    = mapEdge.Get(i,2);
  real zc    = mapEdge.Get(i,3);
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
 neighbourElem.resize (numNodes);

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
 clVector surfaceAux = heaviside==0.5;
 surface = surfaceAux.Find();
 clVector nonSurfaceAux = heaviside!=0.5;
 nonSurface = nonSurfaceAux.Find();
 
 // {x,y,z}Surface representam os valores de {X,Y,Z} dos nos da interface
 xSurface.Dim( surface.Dim() );
 ySurface.Dim( surface.Dim() );
 zSurface.Dim( surface.Dim() );
 for( int i=0;i<surface.Dim();i++ )
 {
  int surfaceNode = surface.Get(i);
  xSurface.Set(i,X.Get( surfaceNode ));
  ySurface.Set(i,Y.Get( surfaceNode ));
  zSurface.Set(i,Z.Get( surfaceNode ));
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
 boundaryVert.resize (0);

 for(int i=0;i<freeFace.DimI();i++ )
 {
  int v1 = freeFace.Get(i,2);
  int v2 = freeFace.Get(i,3);
  int v3 = freeFace.Get(i,4);
  boundaryVert.push_back(v1);
  boundaryVert.push_back(v2);
  boundaryVert.push_back(v3);
 }
 boundaryVert.sort();
 boundaryVert.unique();

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

  real hsum = heaviside.Get(v1)+heaviside.Get(v2)+
              heaviside.Get(v3)+heaviside.Get(v4);

  if( hsum > 1.5 )
   inElem.push_back(i);
  else
   outElem.push_back(i);
 }
}

// cria matrizes de mapeamentos de vizinhos opostos ao vertice em
// questao. Este metodo funciona com Semi-Lagrangian pois permite a
// procura de elementos considerando o vertice em questao. Uma descricao
// completa das matrizes eh definida dentro do proprio metodo.
void Model3D::setOFace()
{
 clMatrix mapViz(numElems,numElems);
 clMatrix faceFace(1000,NUMGLE+1);
 //clMatrix freeFace(1000,NUMGLE+1);
 freeFace.Dim(1000,NUMGLE+1);
 clMatrix faceFaceAux;
 clMatrix freeFaceAux;
 clMatrix mapVizAux(numElems*NUMGLE,NUMGLE);
 clMatrix mapVizAux2(numElems*NUMGLE,NUMGLE);
 clMatrix comb(NUMGLE,NUMGLE-1);  // triangulo = 2
 clVector verts(NUMGLE);           // tetraedro = 3
 clVector face(NUMGLE-1);
 clVector index(numElems*NUMGLE);
 clVector index2(numElems*NUMGLE);
 clVector index3(numElems*NUMGLE);
 clVector idcol(numElems*NUMGLE);
 clVector idcol2(numElems*NUMGLE);
 clVector idcol3(numElems*NUMGLE);
 clVector vect1(NUMGLE+1);
 clVector vert(NUMGLE-1);

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

 setNeighbour();

 int k = 0;

 for( int i=0;i<numElems;i++ )
 {
  verts.CopyRow(i,IEN);
  for( int j=0;j<NUMGLE;j++ )
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
// 	cout << v1 << " " << heaviside.Get(v1) << " " 
// 	     << v2 << " " << heaviside.Get(v2) << " "
// 	     << v3 << " " << heaviside.Get(v3) << " "
// 	     << v4 << " " <<  heaviside.Get(v4) << " " << *mele5 << endl;
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
	faceFace.Dim(iFace+faceFaceAux.DimI(),NUMGLE+1);
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
	freeFace.Dim(iFree+freeFaceAux.DimI(),NUMGLE+1);
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
	freeFace.Dim(iFree+freeFaceAux.DimI(),NUMGLE+1);
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
 faceFace.Dim(iFace,NUMGLE+1);
 freeFace.Dim(iFree,NUMGLE+1);
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
 
 oFace.Dim(numElems,NUMGLE);
 for( int i=0;i<numElems;i++ )
  for( int j=0;j<NUMGLE;j++ )
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

  for( int kk=0;kk<NUMGLE;kk++)
   if( (tetra.Get(kk) != vert.Get(0)) && 
	   (tetra.Get(kk) != vert.Get(1)) && 
	   (tetra.Get(kk) != vert.Get(2)) )
   {
	oFace.Set(elem1,kk,elem2);
	break;
   }

  tetra.CopyRow(elem2,IEN);
  for ( int kk=0;kk<NUMGLE;kk++ )
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
  for( int j=0;j<NUMGLE-1;j++ )
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
	  if( vec[m] != ii ) kFace = m; 
	 
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
 setVertNeighbour(); // neighbourVert (3D mesh)
 setInOutVert(); // inVert e boundaryVert
 setInOutElem(); // inElem e outElem
 setMapEdge(); 

 // update surface, edge matrix, surface neigh elems and points
 restoreMappingArrays();


 //setSurfaceTri(); // triang superficie - interfaceMesh
 //setConvexTri(); // triang parte externa do dominio - convexMesh
 //buildSurfMesh();

 setInterfaceDistance();
 setNormalAndKappa();

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

real Model3D::getAreaVert(int _v1,int _v2,int _v3)
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

 return 0.5*vectorLength(crossX,crossY,crossZ);
}

real Model3D::getAreaElem(int _elem)
{
 int v1=(int) surfMesh.IEN.Get(_elem,0);
 int v2=(int) surfMesh.IEN.Get(_elem,1);
 int v3=(int) surfMesh.IEN.Get(_elem,2);

 // vectors
 real x1 = surfMesh.X.Get(v2) - surfMesh.X.Get(v1);
 real y1 = surfMesh.Y.Get(v2) - surfMesh.Y.Get(v1);
 real z1 = surfMesh.Z.Get(v2) - surfMesh.Z.Get(v1);

 real x2 = surfMesh.X.Get(v3) - surfMesh.X.Get(v1);
 real y2 = surfMesh.Y.Get(v3) - surfMesh.Y.Get(v1);
 real z2 = surfMesh.Z.Get(v3) - surfMesh.Z.Get(v1);

 real crossX = (y1*z2)-(z1*y2);
 real crossY = -( (x1*z2)-(z1*x2) );
 real crossZ = (x1*y2)-(y1*x2);

 return 0.5*vectorLength(crossX,crossY,crossZ);
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

void Model3D::moveXPoints(clVector &_vec,real _dt)
{
 //X = X + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  real aux = X.Get(i)+(_vec.Get(i)*_dt);
  X.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  real aux = surfMesh.X.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.X.Set(i,aux);
 }
}

void Model3D::moveYPoints(clVector &_vec,real _dt)
{
 //Y = Y + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  real aux = Y.Get(i)+(_vec.Get(i)*_dt);
  Y.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  real aux = surfMesh.Y.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.Y.Set(i,aux);
 }
}

void Model3D::moveZPoints(clVector &_vec,real _dt)
{
 //Z = Z + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  real aux = Z.Get(i)+(_vec.Get(i)*_dt);
  Z.Set(i,aux);
 }

 // movimentando os pontos da malha de superficie (interface e convex) 
 // com velocidade _vec e _dt
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  real aux = surfMesh.Z.Get(i)+(_vec.Get(i)*_dt);
  surfMesh.Z.Set(i,aux);
 }
}

void Model3D::movePoints2ndOrder(clVector &_uSol,
                                 clVector &_vSol,
                                 clVector &_wSol,
								 real _dt,real _time)
{
 real T = 3.0;
 real time = _time + _dt/2.0;
 real pi = 3.14159265358;

 //X = X + _vec*_dt;
 for( int i=0;i<numVerts;i++ )
 {
  real Xp = X.Get(i)+(_uSol.Get(i)*_dt/2.0);
  real Yp = Y.Get(i)+(_vSol.Get(i)*_dt/2.0);
  real Zp = Z.Get(i)+(_wSol.Get(i)*_dt/2.0);

  real up = (2.0)*sin(pi*Xp)*
	              sin(pi*Xp)*
			      sin(2*pi*Yp)*
			      sin(2*pi*Zp)*
			      cos(pi*time/T);
  // ---> passar up em stepALE <---
  real xn = X.Get(i)+(up*_dt);
  X.Set(i,xn);

  real vp = (-1.0)*sin(2*pi*Xp)*
	               sin(pi*Yp)*
				   sin(pi*Yp)*
				   sin(2*pi*Zp)*
				   cos(pi*time/T);
  real yn = Y.Get(i)+(vp*_dt);
  Y.Set(i,yn);

  real wp = (-1.0)*sin(2*pi*Xp)*
	               sin(2*pi*Yp)*
				   sin(pi*Zp)*
				   sin(pi*Zp)*
				   cos(pi*time/T);
  real zn = Z.Get(i)+(wp*_dt);
  Z.Set(i,zn);
 }

 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  real Xps = surfMesh.X.Get(i)+(_uSol.Get(i)*_dt/2.0);
  real Yps = surfMesh.Y.Get(i)+(_vSol.Get(i)*_dt/2.0);
  real Zps = surfMesh.Z.Get(i)+(_wSol.Get(i)*_dt/2.0);

  real ups = (2.0)*sin(pi*Xps)*
	               sin(pi*Xps)*
		 	       sin(2*pi*Yps)*
		 	       sin(2*pi*Zps)*
			      cos(pi*time/T);
  real xns = surfMesh.X.Get(i)+(ups*_dt);
  surfMesh.X.Set(i,xns);

  real vps = (-1.0)*sin(2*pi*Xps)*
	                sin(pi*Yps)*
					sin(pi*Yps)*
					sin(2*pi*Zps)*
					cos(pi*time/T);
  real yns = surfMesh.Y.Get(i)+(vps*_dt);
  surfMesh.Y.Set(i,yns);

  real wps = (-1.0)*sin(2*pi*Xps)*
	                sin(2*pi*Yps)*
					sin(pi*Zps)*
					sin(pi*Zps)*
					cos(pi*time/T);
  real zns = surfMesh.Z.Get(i)+(wps*_dt);
  surfMesh.Z.Set(i,zns);
 }
}


SurfaceMesh* Model3D::getSurfMesh(){ return &surfMesh; }
SurfaceMesh* Model3D::getInterfaceMesh(){ return &interfaceMesh; }
SurfaceMesh* Model3D::getConvexMesh(){ return &convexMesh; }
Mesh3D* Model3D::getMesh3d(){ return &mesh3d; }
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
//clMatrix Model3D::getMapViz(){ return mapViz; }
//clMatrix Model3D::getFaceFace(){ return faceFace; }
clMatrix* Model3D::getOFace(){ return &oFace; }
real Model3D::getXCenter(){ return xCenter; }
real Model3D::getYCenter(){ return yCenter; }
real Model3D::getZCenter(){ return zCenter; }
clVector* Model3D::getSurface(){ return &surface; }
vector< list<int> >* Model3D::getNeighbourElem(){return &neighbourElem;}
vector< list<int> >* Model3D::getNeighbourVert(){return &neighbourVert;}
vector< list<int> >* Model3D::getNeighbourFace(){return &neighbourFace;}
vector< list<int> >* Model3D::getNeighbourPoint(){return &neighbourPoint;}
vector< list<int> >* Model3D::getFaceIEN(){return &faceIEN;}
list<int>* Model3D::getInVert(){return &inVert;}
list<int>* Model3D::getBoundaryVert(){return &boundaryVert;}
list<int>* Model3D::getInElem(){return &inElem;}
list<int>* Model3D::getOutElem(){return &outElem;}
real Model3D::getMinEdge(){return minEdge;}
vector<real> Model3D::getTriEdge(){return triEdge;}
void Model3D::setTetVol(vector< real > _tetVol){tetVol= _tetVol;}
vector<real> Model3D::getTetVol(){return tetVol;}
clVector* Model3D::getEdgeSize(){ return &edgeSize; }
void Model3D::setEdgeSize(clVector _edgeSize){ edgeSize = _edgeSize; }
vector<real> Model3D::getAverageTriLength(){ return averageTriLength; }
vector<real> Model3D::getAverageTriArea(){ return averageTriArea; }
vector<real> Model3D::getAverageTetVolume(){ return averageTetVolume; }

vector<real> Model3D::getInitSurfaceArea(){return initSurfaceArea;}
vector<real> Model3D::getInitSurfaceVolume(){return initSurfaceVolume;}
vector<real> Model3D::getSurfaceArea(){return surfaceArea;}
vector<real> Model3D::getSurfaceVolume(){return surfaceVolume;}
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
vector<real> Model3D::getMinArea(){return minArea;}
vector<real> Model3D::getMaxArea(){return maxArea;}
vector<real> Model3D::getMinLength(){return minLength;}
vector<real> Model3D::getMaxLength(){return maxLength;}
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
vector<real> Model3D::getMinVolume(){return minVolume;}
vector<real> Model3D::getMaxVolume(){return maxVolume;}
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
  nonSurface = _mRight.nonSurface;
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
  faceFace = _mRight.faceFace;
  freeFace = _mRight.freeFace;
  mapViz = _mRight.mapViz;
  oFace = _mRight.oFace;
  surfMesh = _mRight.surfMesh;
  interfaceMesh = _mRight.interfaceMesh;
  convexMesh = _mRight.convexMesh;
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
 real P0x = surfMesh.X.Get(_node);
 real P0y = surfMesh.Y.Get(_node);
 real P0z = surfMesh.Z.Get(_node);

 //int c1 = 0;
 real fx = 0;
 real fy = 0;
 real fz = 0;
 real sumArea = 0;
 real sumXCrossUnit = 0;
 real sumYCrossUnit = 0;
 real sumZCrossUnit = 0;

 int listSize = _myList.size();
 list<int>::iterator mele=_myList.begin();
 for( int i=0;i<listSize-1;i++ )
 {
  int v1 = *mele;++mele;
  int v2 = *mele;

  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);
  real P2x = surfMesh.X.Get(v2);
  real P2y = surfMesh.Y.Get(v2);
  real P2z = surfMesh.Z.Get(v2);

  // distance do ponto 0 ate metade do segmento 01
  real a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance do ponto 0 ate metade do segmento 02
  real b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // distance da metade do segmento 01 ate metade do segmento 02
  real c = distance(P1x,P1y,P1z,P2x,P2y,P2z);

  // vetores 
  real x1 = P1x-P0x;
  real y1 = P1y-P0y;
  real z1 = P1z-P0z;

  real x2 = P2x-P0x;
  real y2 = P2y-P0y;
  real z2 = P2z-P0z;

  real xReta = P2x-P1x;
  real yReta = P2y-P1y;
  real zReta = P2z-P1z;

  // vetores unitarios deslocados para origem do sistema (0,0,0)
  real x1Unit = x1/a;
  real y1Unit = y1/a;
  real z1Unit = z1/a;

  real x2Unit = x2/b;
  real y2Unit = y2/b;
  real z2Unit = z2/b;

  // normal to each triangular face
  clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

  // somatorio NAO ponderado pela area dos vetores unitarios normais 
  // aos triangulos encontrados na estrela do vertice
  sumXCrossUnit += cross.Get(0);
  sumYCrossUnit += cross.Get(1);
  sumZCrossUnit += cross.Get(2);

  // soma dos vetores 1Unit + 2Unit = resultante
  real xPlaneRes = x1Unit+x2Unit;
  real yPlaneRes = y1Unit+y2Unit;
  real zPlaneRes = z1Unit+z2Unit;

  clVector proj = projection(xPlaneRes,yPlaneRes,zPlaneRes,
	xReta,yReta,zReta);
  real xPlaneTang = proj.Get(0);
  real yPlaneTang = proj.Get(1);
  real zPlaneTang = proj.Get(2);

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
  real xPlaneNormal = xPlaneRes - xPlaneTang;
  real yPlaneNormal = yPlaneRes - yPlaneTang;
  real zPlaneNormal = zPlaneRes - zPlaneTang;

  real len = vectorLength(xPlaneNormal,yPlaneNormal,zPlaneNormal);

  // Unitario do vetor resultante situado no plano do triangulo
  // combinacao linear dos vetores unitarios das arestas do triangulo
  real xPlaneNormalUnit = xPlaneNormal/len;
  real yPlaneNormalUnit = yPlaneNormal/len;
  real zPlaneNormalUnit = zPlaneNormal/len;

  // normal ao plano integrada na distancia (MOD) dos 2 vertices medianos
  // force = resultante das componentes * tamanho da aresta que sera
  // usada como referencia no calculo da area do triangulo 
  real base = c/2.0;

  fx += xPlaneNormalUnit*base;
  fy += yPlaneNormalUnit*base;
  fz += zPlaneNormalUnit*base;

  // 1/3 of area P0-Pm01-Pm02
  real area = getArea(P0x,P0y,P0z,P1x,P1y,P1z,P2x,P2y,P2z);

  sumArea += (1.0/3.0)*area;

  // norb's correction
  //real fact = (a*b*c*c)/ (4*area*area);
  //sumArea += (1.0/4.0)*area*fact;
 }
 mele=_myList.end();

 real len = vectorLength(sumXCrossUnit,sumYCrossUnit,sumZCrossUnit);
 real xNormalUnit = sumXCrossUnit/len;
 real yNormalUnit = sumYCrossUnit/len;
 real zNormalUnit = sumZCrossUnit/len;

 // intensidade da forca resultante
 real force = -sqrt( (fx*fx)+(fy*fy)+(fz*fz) );

 // teste para saber o sentido correto de aplicacao da
 // pressao no noh.
 if( dotProd(fx,fy,fz,xNormalUnit,yNormalUnit,zNormalUnit) < 0.0 )
  force = -force;

 real pressure = force/sumArea;

 clVector vec(4);
 vec.Set(0,pressure);    // curvature
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

  clVector vec = getNormalAndKappaByDesbrun(node,
	              getNeighbourSurfacePoint(node));

  real pressure    = vec.Get(0);
  real xNormalUnit = vec.Get(1);
  real yNormalUnit = vec.Get(2);
  real zNormalUnit = vec.Get(3);

  surfMesh.xNormal.Set(node,xNormalUnit);
  surfMesh.yNormal.Set(node,yNormalUnit);
  surfMesh.zNormal.Set(node,zNormalUnit);
  surfMesh.curvature.Set(node,pressure);

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
//    cout << "pressure: " << pressure << endl;
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
  real aux = distance( X.Get(i),Y.Get(i),Z.Get(i),
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
real Model3D::getSurfaceVolume(int _region)
{
 real sumVolume = 0;
 //real sumCentroidX = 0;
 //real sumCentroidY = 0;
 //real sumCentroidZ = 0;
 for( int mele=0;mele<surfMesh.numElems;mele++ )
 {
  if( _region == surfMesh.elemIdRegion.Get(mele) )
  {
   // P1
   int v1 = surfMesh.IEN.Get(mele,0);
   real p1x = surfMesh.X.Get(v1);
   real p1y = surfMesh.Y.Get(v1);
   real p1z = surfMesh.Z.Get(v1);

   // P2
   int v2 = surfMesh.IEN.Get(mele,1);
   real p2x = surfMesh.X.Get(v2);
   real p2y = surfMesh.Y.Get(v2);
   real p2z = surfMesh.Z.Get(v2);

   // P3
   int v3 = surfMesh.IEN.Get(mele,2);
   real p3x = surfMesh.X.Get(v3);
   real p3y = surfMesh.Y.Get(v3);
   real p3z = surfMesh.Z.Get(v3);

   // element centroid
   real xCentroid = (p1x+p2x+p3x)/3.0;
   real yCentroid = (p1y+p2y+p3y)/3.0;
   real zCentroid = (p1z+p2z+p3z)/3.0;

   // distance from point 1 to 2
   real a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

   // distance from point 2 to 3
   real b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

   // unit vectors
   real x1Unit = (p2x-p1x)/a;
   real y1Unit = (p2y-p1y)/a;
   real z1Unit = (p2z-p1z)/a;

   real x2Unit = (p3x-p2x)/b;
   real y2Unit = (p3y-p2y)/b;
   real z2Unit = (p3z-p2z)/b;

   // calculando o produto vetorial de cada elemento triangular da superficie
   clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

   // somatorio ponderado pela area dos vetores unitarios normais 
   // aos triangulos encontrados na estrela do vertice
   real xNormalElem = cross.Get(0);
   real yNormalElem = cross.Get(1);
   real zNormalElem = cross.Get(2);

   real len = vectorLength(xNormalElem,yNormalElem,zNormalElem);

   real xNormalElemUnit = xNormalElem/len;
   real yNormalElemUnit = yNormalElem/len;
   real zNormalElemUnit = zNormalElem/len;

   real area = getArea(p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z);

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
 real vol = (1.0/3.0)*sumVolume;
 //real xc = sumCentroidX/(2*vol);
 //real yc = sumCentroidY/(2*vol);
 //real zc = sumCentroidZ/(2*vol);

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
real Model3D::getSurfaceVolumeTET(int _region)
{
 real sumVolume=0;
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
real Model3D::getSurfaceArea(int _region)
{
 real sumArea = 0;
 for(int tri=0;tri<surfMesh.numElems;tri++ )
  if( surfMesh.elemIdRegion.Get(tri) == _region )
   sumArea += getAreaElem(tri);
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
  real xc = centroid.Get(0);
  real yc = centroid.Get(1);
  real zc = centroid.Get(2);

  for( int elem=0;elem<surfMesh.numElems;elem++ )
  {
   if( surfMesh.elemIdRegion.Get(elem) == nb )
   {
	int v1 = surfMesh.IEN.Get(elem,0);
	real p1x = surfMesh.X.Get(v1);
	real p1y = surfMesh.Y.Get(v1);
	real p1z = surfMesh.Z.Get(v1);

	int v2 = surfMesh.IEN.Get(elem,1);
	real p2x = surfMesh.X.Get(v2);
	real p2y = surfMesh.Y.Get(v2);
	real p2z = surfMesh.Z.Get(v2);

	int v3 = surfMesh.IEN.Get(elem,2);
	real p3x = surfMesh.X.Get(v3);
	real p3y = surfMesh.Y.Get(v3);
	real p3z = surfMesh.Z.Get(v3);

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

	real vx = p1x-xc;
	real vy = p1y-yc;
	real vz = p1z-zc;

	real v1x = p2x-p1x;
	real v1y = p2y-p1y;
	real v1z = p2z-p1z;

	real v2x = p3x-p2x;
	real v2y = p3y-p2y;
	real v2z = p3z-p2z;

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
 real sumX = 0;
 real sumY = 0;
 real sumZ = 0;
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
 real xc = sumX/count;
 real yc = sumY/count;
 real zc = sumZ/count;

 clVector centroid(3);
 centroid.Set(0,xc);
 centroid.Set(1,yc);
 centroid.Set(2,zc);

 return centroid;
}

//--------------------------------------------------
// clVector Model3D::getNormalAndKappaByGauss(int _node,list<int> _myList)
// {
//  real P0x = surfMesh.X.Get(_node);
//  real P0y = surfMesh.Y.Get(_node);
//  real P0z = surfMesh.Z.Get(_node);
// 
//  // Normal Vector
//  real xNormal = 0;
//  real yNormal = 0;
//  real zNormal = 0;
// 
//  // Fundamental Form Coefficients
//  real E = 0;
//  real F = 0;
//  real G = 0;
//  real L = 0;
//  real M = 0;
//  real N = 0;
// 
//  int listSize = _myList.size();
//  list<int>::iterator mele=_myList.begin();
//  for( int i=0;i<listSize-1;i++ )
//  {
//   int v1 = *mele;++mele;
//   int v2 = *mele;
// 
//   real P1x = surfMesh.X.Get(v1);
//   real P1y = surfMesh.Y.Get(v1);
//   real P1z = surfMesh.Z.Get(v1);
//   real P2x = surfMesh.X.Get(v2);
//   real P2y = surfMesh.Y.Get(v2);
//   real P2z = surfMesh.Z.Get(v2);
// 
//   // distance do ponto 0 ate metade do segmento 01
//   real a = distance(P0x,P0y,P0z,P1x,P1y,P1z);
// 
//   // distance do ponto 0 ate metade do segmento 02
//   real b = distance(P0x,P0y,P0z,P2x,P2y,P2z);
// 
//   // vetores 
//   real x1 = P1x-P0x;
//   real y1 = P1y-P0y;
//   real z1 = P1z-P0z;
// 
//   real x2 = P2x-P0x;
//   real y2 = P2y-P0y;
//   real z2 = P2z-P0z;
// 
//   // vetores unitarios deslocados para origem do sistema (0,0,0)
//   real x1Unit = x1/a;
//   real y1Unit = y1/a;
//   real z1Unit = z1/a;
// 
//   real x2Unit = x2/b;
//   real y2Unit = y2/b;
//   real z2Unit = z2/b;
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
//   real ux = x1Unit - xNormal;
//   real uy = y1Unit - yNormal;
//   real uz = z1Unit - zNormal;
// 
//   // vector Rv
//   real vx = x2Unit - xNormal;
//   real vy = y2Unit - yNormal;
//   real vz = z2Unit - zNormal;
// 
//   // First Fundamental Form
//   E = ux*ux + uy*uy + uz*uz;
//   F = ux*vx + uy*vy + uz*vz;
//   G = vx*vx + vy*vy + vz*vz;
// 
//   // vector Ruu
//   real uux = 0;
//   real uuy = 0;
//   real uuz = 0;
// 
//   // vector Ruv
//   real uvx = 0;
//   real uvy = 0;
//   real uvz = 0;
// 
//   // vector Rvv
//   real vvx = 0;
//   real vvy = 0;
//   real vvz = 0;
//   
//   // Second Fundamental Form
//   L = uux*xNormal + uuy*yNormal + uuz*zNormal;
//   M = uvx*xNormal + uvy*yNormal + uvz*zNormal;
//   N = vvx*xNormal + vvy*yNormal + vvz*zNormal;
// 
//  }
//  mele=_myList.end();
// 
//  real len = vectorLength(xNormal,yNormal,zNormal);
//  real xNormalUnit = xNormal/len;
//  real yNormalUnit = yNormal/len;
//  real zNormalUnit = zNormal/len;
// 
//  real meanCurvature = (E*N + G*L - 2*F*M)/2*(E*G-F*F);
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
 real P0x = surfMesh.X.Get(_node);
 real P0y = surfMesh.Y.Get(_node);
 real P0z = surfMesh.Z.Get(_node);

 //int c1 = 0;
 real fx = 0;
 real fy = 0;
 real fz = 0;
 real sumVoronoiArea = 0;
 real sumXCrossUnit = 0;
 real sumYCrossUnit = 0;
 real sumZCrossUnit = 0;

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

  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);
  real P2x = surfMesh.X.Get(v2);
  real P2y = surfMesh.Y.Get(v2);
  real P2z = surfMesh.Z.Get(v2);
  real P3x = surfMesh.X.Get(v3);
  real P3y = surfMesh.Y.Get(v3);
  real P3z = surfMesh.Z.Get(v3);

  // distance from point 0 to 1  
  real a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance from point 0 to 2  
  real b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // distance from point 0 to 3  
  real c = distance(P0x,P0y,P0z,P3x,P3y,P3z);

  // distance from point 1 to 2  
  real d = distance(P1x,P1y,P1z,P2x,P2y,P2z);

  // distance from point 2 to 3  
  real e = distance(P2x,P2y,P2z,P3x,P3y,P3z);

  // vetores 
  real x1 = P1x-P0x;
  real y1 = P1y-P0y;
  real z1 = P1z-P0z;

  real x2 = P2x-P0x;
  real y2 = P2y-P0y;
  real z2 = P2z-P0z;

  // vetores unitarios deslocados para origem do sistema (0,0,0)
  real x1Unit = x1/a;
  real y1Unit = y1/a;
  real z1Unit = z1/a;

  real x2Unit = x2/b;
  real y2Unit = y2/b;
  real z2Unit = z2/b;

  // normal to each triangular face
  clVector cross = crossProd(x1Unit,y1Unit,z1Unit,x2Unit,y2Unit,z2Unit);

  // somatorio NAO ponderado pela area dos vetores unitarios normais 
  // aos triangulos encontrados na estrela do vertice
  sumXCrossUnit += cross.Get(0);
  sumYCrossUnit += cross.Get(1);
  sumZCrossUnit += cross.Get(2);

  // angles
  real alpha = acos( (a*a+d*d-b*b)/(2*a*d) );
  real beta = acos( (c*c+e*e-b*b)/(2*c*e) );

  real a1 = acos( (b*b+a*a-d*d)/(2*a*b) );
  real a2 = alpha;
  real a3 = acos( (b*b+d*d-a*a)/(2*b*d) );

//--------------------------------------------------
//   cout << a1 << " " << a2 << " " << " " << a3 << " " << a1+a2+a3 <<endl;
//   cout << a4 << " " << alpha << " " << beta <<endl;
//   cout << "----"<< endl;
//-------------------------------------------------- 

  fx += ( 1.0/tan(alpha)+1.0/tan(beta) )*x2;
  fy += ( 1.0/tan(alpha)+1.0/tan(beta) )*y2;
  fz += ( 1.0/tan(alpha)+1.0/tan(beta) )*z2;

  // area P0-P1-P2
  real area1 = getArea(P0x,P0y,P0z,P1x,P1y,P1z,P2x,P2y,P2z);

  // voronoi area
  if( a1 < 3.14/2.0 && a2 < 3.14/2.0 && a3 && 3.14/2.0 )
   sumVoronoiArea += (1.0/8.0)*( b*b*1.0/(tan(alpha)+a*a*1.0/(tan(a3))) );
  else
   if( a1 > 3.14/2.0 )
	sumVoronoiArea += area1/2.0;
   else
	sumVoronoiArea += area1/4.0;
 }
 mele=_myList.end();

 fx = fx/(2*sumVoronoiArea);
 fy = fy/(2*sumVoronoiArea);
 fz = fz/(2*sumVoronoiArea);

 real len = vectorLength(sumXCrossUnit,sumYCrossUnit,sumZCrossUnit);
 real xNormalUnit = sumXCrossUnit/len;
 real yNormalUnit = sumYCrossUnit/len;
 real zNormalUnit = sumZCrossUnit/len;

 // intensidade da forca resultante
 real force = -sqrt( (fx*fx)+(fy*fy)+(fz*fz) );

 // teste para saber o sentido correto de aplicacao da
 // pressao no noh.
 if( dotProd(fx,fy,fz,xNormalUnit,yNormalUnit,zNormalUnit) < 0.0 )
  force = -force;

 real pressure = force/2.0;

 clVector vec(4);
 vec.Set(0,pressure);
 vec.Set(1,xNormalUnit);
 vec.Set(2,yNormalUnit);
 vec.Set(3,zNormalUnit);

 return vec;
} // fecha metodo getNormalAndKappaByDesbrun

clVector Model3D::getNormalElem(int _elem)
{
 int v1 = surfMesh.IEN.Get(_elem,0);
 real p1x = surfMesh.X.Get(v1);
 real p1y = surfMesh.Y.Get(v1);
 real p1z = surfMesh.Z.Get(v1);

 int v2 = surfMesh.IEN.Get(_elem,1);
 real p2x = surfMesh.X.Get(v2);
 real p2y = surfMesh.Y.Get(v2);
 real p2z = surfMesh.Z.Get(v2);

 int v3 = surfMesh.IEN.Get(_elem,2);
 real p3x = surfMesh.X.Get(v3);
 real p3y = surfMesh.Y.Get(v3);
 real p3z = surfMesh.Z.Get(v3);
 
 // distance from point 1 to 2
 real a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

 // distance from point 2 to 3
 real b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

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
 real v1x = (p2x-p1x)/a;
 real v1y = (p2y-p1y)/a;
 real v1z = (p2z-p1z)/a;

 real v2x = (p3x-p2x)/b;
 real v2y = (p3y-p2y)/b;
 real v2z = (p3z-p2z)/b;

 clVector normal = crossProd(v1x,v1y,v1z,v2x,v2y,v2z);
 real length = vectorLength(normal.Get(0),normal.Get(1),normal.Get(2));

 // unit normal to each triangular face
 return normal/length;
} // fecha metodo getNormalElem

clVector Model3D::getNormalElem(int _v1,int _v2,int _v3)
{
 real p1x = surfMesh.X.Get(_v1);
 real p1y = surfMesh.Y.Get(_v1);
 real p1z = surfMesh.Z.Get(_v1);

 real p2x = surfMesh.X.Get(_v2);
 real p2y = surfMesh.Y.Get(_v2);
 real p2z = surfMesh.Z.Get(_v2);

 real p3x = surfMesh.X.Get(_v3);
 real p3y = surfMesh.Y.Get(_v3);
 real p3z = surfMesh.Z.Get(_v3);
 
 // distance from point 1 to 2
 real a = distance(p1x,p1y,p1z,p2x,p2y,p2z);

 // distance from point 2 to 3
 real b = distance(p2x,p2y,p2z,p3x,p3y,p3z);

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
 real v1x = (p2x-p1x)/a;
 real v1y = (p2y-p1y)/a;
 real v1z = (p2z-p1z)/a;

 real v2x = (p3x-p2x)/b;
 real v2y = (p3y-p2y)/b;
 real v2z = (p3z-p2z)/b;

 clVector normal = crossProd(v1x,v1y,v1z,v2x,v2y,v2z);
 real length = vectorLength(normal.Get(0),normal.Get(1),normal.Get(2));

 // unit normal to each triangular face
 return normal/length;
} // fecha metodo getNormalElem

real Model3D::triangleQualityMeasure(int _v1,int _v2, int _v3)
{
 real P0x = surfMesh.X.Get(_v1);
 real P0y = surfMesh.Y.Get(_v1);
 real P0z = surfMesh.Z.Get(_v1);

 real P1x = surfMesh.X.Get(_v2);
 real P1y = surfMesh.Y.Get(_v2);
 real P1z = surfMesh.Z.Get(_v2);

 real P2x = surfMesh.X.Get(_v3);
 real P2y = surfMesh.Y.Get(_v3);
 real P2z = surfMesh.Z.Get(_v3);

 // Triangle quality measure;
 real length12 = distance(P0x,P0y,P0z,P1x,P1y,P1z);
 real length13 = distance(P0x,P0y,P0z,P2x,P2y,P2z);
 real length23 = distance(P1x,P1y,P1z,P2x,P2y,P2z);

 real semiPerimeter = 0.5*(length12+length13+length23);
 real area = getArea(P0x,P0y,P0z,P1x,P1y,P1z,P2x,P2y,P2z);
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

 // triangle measure quality q = 3.4641*inRadius/h;
 return 3.4641*inRadius/h;
}

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
  real v1 = mapEdgeTri.Get(edge,1);
  real v2 = mapEdgeTri.Get(edge,2);
  real v3elem1 = mapEdgeTri.Get(edge,3);
  real v3elem2 = mapEdgeTri.Get(edge,4);
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



  real theta = angle3D( normalElem1.Get(0),
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
void Model3D::removePointByNeighbourCheck()
{
 for( int i=0;i<surfMesh.numVerts;i++ )
 {
  int vertID = surfMesh.vertIdRegion.Get(i);
  int elemListSize = neighbourSurfaceElem.at( i ).size();

  if( elemListSize < 3 )
  {
   removeSurfacePoint(i);
   cout << "------------- " << color(none,red,black) 
	<< "removing fake triangle: " << resetColor() 
	<< i << endl;
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
   removeSurfacePoint(i);

   cout << "------------- " << color(none,red,black) 
	<< "removing low-quality point cluster: " << resetColor() 
	<< i << endl;

   saveVTKSurface("./vtk/","surface",opersurf[vertID]);
   rspn[vertID]++;
   opersurf[vertID]++;
  }

//--------------------------------------------------
//   if( elemListSize == 4 && surfMesh.curvature.Get(i) > 50 )
//   {
//    removeSurfacePoint(i);
// 
//    cout << "------------- " << color(none,magenta,black) 
// 	    << "removing low-quality triangle with curvature (" 
// 		<< resetColor() << surfMesh.curvature.Get(i)
// 		<< color(none,magenta,black) 
// 		<< ") at (" 
// 		<< resetColor()
// 		<< surfMesh.vertIdRegion.Get(i)
// 		<< color(none,magenta,black) 
// 		<< "): "
// 		 << resetColor() << i << endl;
// 
//    saveVTKSurface("./vtk/","surface",opersurf[vertID]);
//    rspn[vertID]++;
//    opersurf[vertID]++;
//   }
//-------------------------------------------------- 
 }
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
  real vertID = surfMesh.vertIdRegion.Get(surfaceNode); 
  real curv = fabs(surfMesh.curvature.Get(surfaceNode));

  real maxEdgeLength = 0;
  int listSize = neighbourPoint.at(surfaceNode).size();
  list<int> plist = neighbourPoint.at(surfaceNode);
  list<int>::iterator vert=plist.begin();
  for( int i=0;i<listSize-1;i++ )
  {
   real P0x = surfMesh.X.Get(surfaceNode);
   real P0y = surfMesh.Y.Get(surfaceNode);
   real P0z = surfMesh.Z.Get(surfaceNode);

   int v1 = *vert;++vert;
   real P1x = surfMesh.X.Get(v1);
   real P1y = surfMesh.Y.Get(v1);
   real P1z = surfMesh.Z.Get(v1);

   real edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

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
 real aux;
 real xSum = 0.0;
 real ySum = 0.0;
 real zSum = 0.0;
 real distSum = 0.0;

 int listSize = neighbourPoint.at(_node).size();
 list<int> plist = neighbourPoint.at(_node);
 list<int>::iterator vert=plist.begin();
 for( int i=0;i<listSize-1;i++ )
 {
  real P0x = surfMesh.X.Get(_node);
  real P0y = surfMesh.Y.Get(_node);
  real P0z = surfMesh.Z.Get(_node);

  int v1 = *vert;++vert;
  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);

  real edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);

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
void Model3D::setTriEdge(vector< real > _triEdge)
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
 idMinVolume.resize(numSurface);  // ID of min tet volume
 idMaxVolume.resize(numSurface);  // ID of max tet volume
 minVolume.resize(numSurface);    // min tet volume
 maxVolume.resize(numSurface);    // max tet volume     
 idMaxArea.resize(numSurface);    // ID of min tri area   
 idMinArea.resize(numSurface);    // ID of max tri area                            
 maxArea.resize(numSurface);      // min triangle area    
 minArea.resize(numSurface);      // max triangle area                                     
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
 fill(maxVolume.begin(),maxVolume.end(),1.0E-20);
 fill(minVolume.begin(),minVolume.end(),1.0E+20);
 fill(idMaxVolume.begin(),idMaxVolume.end(),0);
 fill(idMinVolume.begin(),idMinVolume.end(),0);
 fill(maxArea.begin(),maxArea.end(),1.0E-20);
 fill(minArea.begin(),minArea.end(),1.0E+20);
 fill(idMaxArea.begin(),idMaxArea.end(),0);
 fill(idMinArea.begin(),idMinArea.end(),0);
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
  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);
  real P2x = surfMesh.X.Get(v2);
  real P2y = surfMesh.Y.Get(v2);
  real P2z = surfMesh.Z.Get(v2);
  real P3x = surfMesh.X.Get(v3);
  real P3y = surfMesh.Y.Get(v3);
  real P3z = surfMesh.Z.Get(v3);

  // centroid
  real xMid = (P1x+P2x+P3x)/3.0;
  real yMid = (P1y+P2y+P3y)/3.0;
  real zMid = (P1z+P2z+P3z)/3.0;

  // mid 12
  real xMid12 = (P1x+P2x)/2.0;
  real yMid12 = (P1y+P2y)/2.0;
  real zMid12 = (P1z+P2z)/2.0;

  // mid 13
  real xMid13 = (P1x+P3x)/2.0;
  real yMid13 = (P1y+P3y)/2.0;
  real zMid13 = (P1z+P3z)/2.0;

  // mid 23
  real xMid23 = (P2x+P3x)/2.0;
  real yMid23 = (P2y+P3y)/2.0;
  real zMid23 = (P2z+P3z)/2.0;

  // list of neighbouring points
  list<int> plist = neighbourVert.at(v1);
  for(list<int>::iterator vert=plist.begin(); vert != plist.end();++vert )
  {
   // *vert cannot be a surface mesh vertex
   if( *vert > surfMesh.numVerts )
   {
	int vertID = surfMesh.vertIdRegion.Get(v1);
	int vertID3d = vertIdRegion.Get(*vert);

	real Pxvert = X.Get(*vert);
	real Pyvert = Y.Get(*vert);
	real Pzvert = Z.Get(*vert);

	real height1 = distance(P1x,P1y,P1z,Pxvert,Pyvert,Pzvert);
	real height2 = distance(P2x,P2y,P2z,Pxvert,Pyvert,Pzvert);
	real height3 = distance(P3x,P3y,P3z,Pxvert,Pyvert,Pzvert);
	// centroid
	real height4 = distance(xMid,yMid,zMid,Pxvert,Pyvert,Pzvert);
	// xMid12 
	real height5 = distance(xMid12,yMid12,zMid12,Pxvert,Pyvert,Pzvert);
	// xMid13 
	real height6 = distance(xMid13,yMid13,zMid13,Pxvert,Pyvert,Pzvert);
	// xMid23 
	real height7 = distance(xMid23,yMid23,zMid23,Pxvert,Pyvert,Pzvert);

	real minHeight = min(height1,height2);
	minHeight = min(minHeight,height3);
	minHeight = min(minHeight,height4);
	minHeight = min(minHeight,height5);
	minHeight = min(minHeight,height6);
	minHeight = min(minHeight,height7);

	bool pressureTest=(surfMesh.phyBounds.at(v1) == "\"wallOutflow\"" ) ||
	                  (surfMesh.phyBounds.at(v2) == "\"wallOutflow\"" ) ||
	                  (surfMesh.phyBounds.at(v3) == "\"wallOutflow\"" );

	if( minHeight < 0.3*triEdge[vertID] && 
	//if( minHeight < 0.4*edgeSize.Get(v1)  && 
	    pressureTest 
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

  real P0x = surfMesh.X.Get(v0);
  real P0y = surfMesh.Y.Get(v0);
  real P0z = surfMesh.Z.Get(v0);

  real P1x = surfMesh.X.Get(v1);
  real P1y = surfMesh.Y.Get(v1);
  real P1z = surfMesh.Z.Get(v1);

  real P2x = surfMesh.X.Get(v2);
  real P2y = surfMesh.Y.Get(v2);
  real P2z = surfMesh.Z.Get(v2);

  // distance do ponto 0 ate ponto 1
  real a = distance(P0x,P0y,P0z,P1x,P1y,P1z);

  // distance do ponto 0 ate ponto 2
  real b = distance(P0x,P0y,P0z,P2x,P2y,P2z);

  // vetors
  real x1Unit = (P1x-P0x)/a;
  real y1Unit = (P1y-P0y)/a;

  real x2Unit = (P2x-P0x)/b;
  real y2Unit = (P2y-P0y)/b;

  real fx = x1Unit+x2Unit;
  real fy = y1Unit+y2Unit;

  // 1/2 of length P0-P1 and P0-P2
  real sumLength = (a+b)/2.0;

  /* 2D rotation of z = 90 degrees
   * x' = x*cos(z) - y*sin(z)
   * y' = x*sin(z) + y*cos(z)
   * */

  real xNormalUnit = +y1Unit*1 - y2Unit*1;
  real yNormalUnit = -x1Unit*1 + x2Unit*1;

  real len = vectorLength(xNormalUnit,yNormalUnit);

  xNormalUnit = xNormalUnit/len;
  yNormalUnit = yNormalUnit/len;

  // intensidade da forca resultante
  real force = sqrt( (fx*fx)+(fy*fy) );

  // aplicando o teste para saber o sentido correto de aplicacao da
  // pressao no noh.
  if( (fx*xNormalUnit+fy*yNormalUnit) > 0.0 )
   force = -force;

  real pressure = force/sumLength;
//--------------------------------------------------
//   cout << v0 << " " 
//        << v1 << " "
//        << v2 << " "
// 	   << fx << " "
// 	   << fy << " "
// 	   << xNormalUnit << " "
// 	   << yNormalUnit << " "
// 	   << pressure << endl;
//-------------------------------------------------- 

  surfMesh.curvature.Set(v0,fabs(pressure));
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
  real rMax = 1.0;
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

   //real aux = 0.0;
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

   real aux = 0.0;
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
  real rMax = 0.5;

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

clVector Model3D::computeConvexRegionCentroid2D(real _zPlane)
{
 real sumX = 0;
 real sumY = 0;
 //real sumZ = 0;
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

 real xc = sumX/count;
 real yc = sumY/count;
 //real zc = sumZ/count;

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
  real xc = centroid.Get(0);
  real yc = centroid.Get(1);

  for( int elem=0;elem<surfMesh.numElems;elem++ )
  {
   int v1 = surfMesh.IEN.Get(elem,0);
   real p1x = surfMesh.X.Get(v1);
   real p1y = surfMesh.Y.Get(v1);

   int v2 = surfMesh.IEN.Get(elem,1);
   real p2x = surfMesh.X.Get(v2);
   real p2y = surfMesh.Y.Get(v2);


   /*               
	*           o ----- o 
	*         v1         v2
	* */

   // reference vector
   real vx = p1x-xc;
   real vy = p1y-yc;

   // line edge vector
   real x1 = p1x-p2x;
   real y1 = p1y-p2y;

   // Normal
   /* 2D rotation of z = 90 degrees
	* x' = x*cos(z) - y*sin(z)
	* y' = x*sin(z) + y*cos(z)
	* */
   real v1x = -y1*1;
   real v1y = +x1*1;

   if( dotProd(vx,vy,v1x,v1y) > 0.0 )
   {
	surfMesh.IEN.Set(elem,0,v2);
	surfMesh.IEN.Set(elem,1,v1);
   }
  }
 }
}

//--------------------------------------------------
// void Model3D::applyBubbleVolumeCorrection()
// {
//  // surfMesh.elemIdRegion == 0 --> wall
//  // surfMesh.elemIdRegion == 1 --> bubble 1
//  // surfMesh.elemIdRegion == 2 --> bubble 2 , etc
//  for( int nb=1;nb<=surfMesh.elemIdRegion.Max();nb++ )
//  {
//   real aux = 0;
// 
//   real da = (initSurfaceArea[nb] - surfaceArea[nb]);
//   real dv = (initSurfaceVolume[nb] - surfaceVolume[nb]);
// 
//   real erroa = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
//   real errov = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);
// 
//   real TOL = initSurfaceVolume[nb]*0.00000001;
// 
//   int count = 0;
//   //while( fabs(dv) > 1E-06 )
//   //while( fabs(da) > 1E-06 )
//   //while( fabs(dv) > TOL && count < 30 )
//   while( fabs(errov) > TOL && count < 30 )
//   {
//    for( int i=0;i<surface.Dim();i++ )
//    {
// 	int surfaceNode = surface.Get(i);
// 
// 	if( surfMesh.vertIdRegion.Get(surfaceNode) == nb )
// 	{
// 	 aux = surfMesh.X.Get(surfaceNode) + 
// 	       surfMesh.xNormal.Get(surfaceNode)*triEdge[nb]*errov;
// 	       //surfMesh.xNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
// 	 X.Set(surfaceNode,aux);
// 	 surfMesh.X.Set(surfaceNode,aux);
// 
// 	 aux = surfMesh.Y.Get(surfaceNode) + 
// 	       surfMesh.yNormal.Get(surfaceNode)*triEdge[nb]*errov;
// 	       //surfMesh.yNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
// 	 Y.Set(surfaceNode,aux);
// 	 surfMesh.Y.Set(surfaceNode,aux);
// 
// 	 aux = surfMesh.Z.Get(surfaceNode) + 
// 	       surfMesh.zNormal.Get(surfaceNode)*triEdge[nb]*errov;
// 	       //surfMesh.zNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
// 	 Z.Set(surfaceNode,aux);
// 	 surfMesh.Z.Set(surfaceNode,aux);
// 	}
//    }
//    surfaceVolume[nb] = getSurfaceVolume(nb);
//    surfaceArea[nb] = getSurfaceArea(nb);
//    da = (initSurfaceArea[nb] - surfaceArea[nb]);
//    dv = (initSurfaceVolume[nb] - surfaceVolume[nb]);
//    erroa = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
//    errov = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);
//    count++;
//    //cout << nb << " " << dv << " " << initSurfaceVolume[nb] << " " << surfaceVolume[nb] << endl;
//   }
//   cout << endl;
//   cout << setw(20) << color(none,red,black) 
//                    << "|--------- VOLUME CORRECTION ---------|" << endl;
//   cout << setw(33) << color(none,white,black) << "|initial: " 
//                    << initSurfaceVolume[nb] << endl;
//   cout << setw(33) << color(none,white,black) 
//                    << "|final: " << surfaceVolume[nb] << endl;
//   cout << setw(33) << color(none,white,black) 
//                    << "|dv: " << dv << endl;
//   cout << setw(26) << color(none,white,black) 
//                    << "volume |error: " << fabs(errov) << endl;
//   cout << setw(21) << color(none,red,black) 
//                    << "     ---------------------------- " << endl;
//   cout << setw(33) << color(none,white,black) << "|initial: " 
//                    << initSurfaceArea[nb] << endl;
//   cout << setw(33) << color(none,white,black) 
//                    << "|final: " << surfaceArea[nb] << endl;
//   cout << setw(33) << color(none,white,black) 
//                    << "|da: " << da << endl;
//   cout << setw(26) << color(none,white,black) 
//                    << "  area |error: " << fabs(erroa) << endl;
//   cout << setw(21) << color(none,red,black) 
//                    << "     ---------------------------- " << endl;
//   cout << setw(28) << color(none,white,black) 
//                    << "number of iterations: " << count << endl;
//   cout << setw(20) << color(none,red,black) 
//                    << "|-------------------------------------|" << endl;
//   cout << resetColor() << endl;
//  }
// }
//-------------------------------------------------- 

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
                    << setw(18) << "errov" 
                    << setw(18) << "erroa" 
                    << setw(18) << "dv" 
                    << setw(18) << "da" 
     			    << endl;
  }

  real aux = 0;

  real da = (initSurfaceArea[nb] - surfaceArea[nb]);
  real dv = (initSurfaceVolume[nb] - surfaceVolume[nb]);

  real erroa = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
  real errov = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);

  real TOL = initSurfaceVolume[nb]*0.00000001;

  int count = 0;
  //while( fabs(dv) > 1E-06 )
  //while( fabs(da) > 1E-06 )
  //while( fabs(dv) > TOL && count < 30 )
  while( fabs(errov) > TOL && count < 30 )
  {
  file << setprecision(20) << scientific; 
  file << setw(10) << count << " " 
       << setw(17) << surfaceVolume[nb] << " " 
       << setw(17) << surfaceArea[nb] << " " 
       << setw(17) << errov << " " 
       << setw(17) << erroa << " " 
       << setw(17) << dv << " " 
       << setw(17) << da << " " 
       << setw(5) << setprecision(0) << fixed  
       << endl;

   for( int i=0;i<surface.Dim();i++ )
   {
	int surfaceNode = surface.Get(i);

	if( surfMesh.vertIdRegion.Get(surfaceNode) == nb )
	{
	 aux = surfMesh.X.Get(surfaceNode) + 
	       surfMesh.xNormal.Get(surfaceNode)*triEdge[nb]*errov;
	       //surfMesh.xNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
	 X.Set(surfaceNode,aux);
	 surfMesh.X.Set(surfaceNode,aux);

	 aux = surfMesh.Y.Get(surfaceNode) + 
	       surfMesh.yNormal.Get(surfaceNode)*triEdge[nb]*errov;
	       //surfMesh.yNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
	 Y.Set(surfaceNode,aux);
	 surfMesh.Y.Set(surfaceNode,aux);

	 aux = surfMesh.Z.Get(surfaceNode) + 
	       surfMesh.zNormal.Get(surfaceNode)*triEdge[nb]*errov;
	       //surfMesh.zNormal.Get(surfaceNode)*1.1*(dv/fabs(da));
	 Z.Set(surfaceNode,aux);
	 surfMesh.Z.Set(surfaceNode,aux);
	}
   }
   surfaceVolume[nb] = getSurfaceVolume(nb);
   surfaceArea[nb] = getSurfaceArea(nb);
   da = (initSurfaceArea[nb] - surfaceArea[nb]);
   dv = (initSurfaceVolume[nb] - surfaceVolume[nb]);
   erroa = (1.0 - surfaceArea[nb]/initSurfaceArea[nb]);
   errov = (1.0 - surfaceVolume[nb]/initSurfaceVolume[nb]);

   count++;
   //cout << nb << " " << dv << " " << initSurfaceVolume[nb] << " " << surfaceVolume[nb] << endl;
   
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
                   << "|dv: " << dv << endl;
  cout << setw(26) << color(none,white,black) 
                   << "volume |error: " << fabs(errov) << endl;
  cout << setw(21) << color(none,red,black) 
                   << "     ---------------------------- " << endl;
  cout << setw(33) << color(none,white,black) << "|initial: " 
                   << initSurfaceArea[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|final: " << surfaceArea[nb] << endl;
  cout << setw(33) << color(none,white,black) 
                   << "|da: " << da << endl;
  cout << setw(26) << color(none,white,black) 
                   << "  area |error: " << fabs(erroa) << endl;
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
 real sumUArea = 0;
 real sumArea = 0;
 for( int i=0;i<surfMesh.numElems;i++ )
 {
  int v1 = surfMesh.IEN.Get(i,0);
  int v2 = surfMesh.IEN.Get(i,1);
  int v3 = surfMesh.IEN.Get(i,2);

  real P1x = surfMesh.X.Get(v1);
  real P2x = surfMesh.X.Get(v2);
  real P3x = surfMesh.X.Get(v3); 

  if( P1x == X.Max() &&
      P2x == X.Max() && 
      P3x == X.Max() )
  {
   real area = getAreaElem(i);
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

  real vol = fabs(getVolume(elem));

  int maxVert = max(v1,v2);
  maxVert = max(maxVert,v3);
  maxVert = max(maxVert,v4);

  real edgeMean = ( edgeSize.Get(v1) +
                    edgeSize.Get(v2) +
					edgeSize.Get(v3) +
					edgeSize.Get(v4) )/4.0;

  real tet = edgeMean*edgeMean*edgeMean*sqrt(2.0)/12.0;

  if( vol > 4.0*tet )
  {
   int vAdd = numVerts; // aditional vertice
   real XvAdd = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )/4.0;
   real YvAdd = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )/4.0;
   real ZvAdd = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )/4.0;

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
