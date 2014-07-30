// =================================================================== //
// this is file Helmholtz3D.cpp, created at 21-Sep-2011                  //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Helmholtz3D.h"

Helmholtz3D::Helmholtz3D(){}

Helmholtz3D::Helmholtz3D( Model3D &_m )  
{
 getModel3DAttrib(_m);

 setSolver( new PCGSolver() );

 allocateMemoryToAttrib();

 k = 0.1;
}

Helmholtz3D::Helmholtz3D( Model3D &_m,Helmholtz3D &_d )  
{
 getModel3DAttrib(_m);

 setSolver( new PCGSolver() );

 allocateMemoryToAttrib();

 k = _d.getk();
 cSolOld = *_m.getEdgeSize();
}

void Helmholtz3D::init()
{
 cSolOld.CopyFrom( 0,cc );
}

void Helmholtz3D::initImposedField()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
 }

 /* loop at surfMesh.numVerts->numVerts: 
  *
  *
  * */
 clVector* vertIdRegion = m->getVertIdRegion();
 double minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double radius = 0.15; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   double factor = triEdge[vertIdRegion->Get(i)]/minEdge;
   if( interfaceDistance->Get(i) < 1.0*radius )
   {
	double aux = triEdge[vertIdRegion->Get(i)]/factor;
	convC.Set(i,aux);
   }
   else
   {
	double aux = triEdge[vertIdRegion->Get(i)]/(factor*0.2);
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   double aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::initThreeBubbles()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
 }

 clVector* vertIdRegion = m->getVertIdRegion();
 double minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double radius = 0.5; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   double factor = triEdge[vertIdRegion->Get(i)]/minEdge;
   if( interfaceDistance->Get(i) < 3.5*radius )
   {
	double aux = triEdge[vertIdRegion->Get(i)]/factor;
	convC.Set(i,aux);
   }
   else
   {
	double aux = triEdge[vertIdRegion->Get(i)]/(factor*0.2);
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   double aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::initRisingBubble()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
 }

 clVector* vertIdRegion = m->getVertIdRegion();
 double minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double radius = 0.5; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   double factor = triEdge[vertIdRegion->Get(i)]/minEdge;
   if( interfaceDistance->Get(i) < 2.5*radius )
   {
	double aux = triEdge[vertIdRegion->Get(i)]/factor;
	convC.Set(i,aux);
   }
   else
   {
	double aux = triEdge[vertIdRegion->Get(i)]/(factor*0.2);
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   double aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::initSessile()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
  
  // Adding bubble as boundary condition
  cc.Set(i,sumEdgeLength/listSize);
  if( surfMesh->Marker.Get(i) >= 0.5 )
   idbcc.AddItem(i);
 }

 clVector* vertIdRegion = m->getVertIdRegion();
 double minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double radius = 0.5; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   double factor = triEdge[vertIdRegion->Get(i)]/minEdge;
   if( interfaceDistance->Get(i) < 1.2*radius )
   {
	double aux = triEdge[vertIdRegion->Get(i)]/factor;
	convC.Set(i,aux);
   }
   else
   {
	double aux = triEdge[vertIdRegion->Get(i)];
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   double aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::initMicro()
{
 // init closer
 setCloserWall();

 // set liquid film thickness
 double thickness = 0.2;

 double triEdgeMin = *(min_element(triEdge.begin(),triEdge.end()));
 convC.Dim(numVerts);
 convC.SetAll(triEdgeMin);

 //double xMid = (X->Min()+X->Max())*0.5;
 double yMid = (Y->Min()+Y->Max())*0.5;
 double zMid = (Z->Min()+Z->Max())*0.5;
 //double diameter = ( (X->Max()-X->Min())+(Y->Max()-Y->Min()) )*0.5;
 //double diameter = ( (X->Max()-X->Min())+(Z->Max()-Z->Min()) )*0.5;
 double diameter = ( (Y->Max()-Y->Min())+(Z->Max()-Z->Min()) )*0.5;
 //double tol = diameter/3.0;
 double tol = diameter/2.0-thickness;
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  // middle and out of liquid/gas film
  // triEdgeMin*5 is the refinement level at circular.geo, where far
  // from the bubble the edge length is 5x bigger.
  if( //(X->Get(i) > xMid-epslocal && X->Get(i) < yMid+epslocal) &&
	  (Y->Get(i) > yMid-tol && Y->Get(i) < yMid+tol ) &&
	  (Z->Get(i) > zMid-tol && Z->Get(i) < zMid+tol ) &&
      ( interfaceDistance->Get(i) > 0.2*diameter ) )
   convC.Set(i,triEdgeMin*5);
  else
  {
   // add boundary condition value to convC. Note that the boundary
   // conditions is bases on the initial mesh, thus if the mesh is
   // refined, convC will be refined.
   convC.Set(i,cc.Get(closerWall.Get(i)));
   //convC.Set(i,triEdgeMin*0.7);
  }
 }
}

void Helmholtz3D::initSquareChannel()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
 }

 //double xMid = (X->Min()+X->Max())*0.5;
 double yMid = (Y->Min()+Y->Max())*0.5;
 double zMid = (Z->Min()+Z->Max())*0.5;
 //double diameter = ( (X->Max()-X->Min())+(Y->Max()-Y->Min()) )*0.5;
 //double diameter = ( (X->Max()-X->Min())+(Z->Max()-Z->Min()) )*0.5;
 double diameter = ( (Y->Max()-Y->Min())+(Z->Max()-Z->Min()) )*0.5;
 double epslocal = 0.1*diameter;
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  if( //(X->Get(i) > xMid-epslocal && X->Get(i) < yMid+epslocal) &&
	(Y->Get(i) > yMid-epslocal && Y->Get(i) < yMid+epslocal) &&
	(Z->Get(i) > zMid-epslocal && Z->Get(i) < zMid+epslocal) )
  {
   double aux = triEdge[0]*10;
   convC.Set(i,aux);
  }
  else
  {
   double aux = triEdge[0];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::initBackwardStep()
{
 init();

 // init closer
 setCloserWall();

 // set wall distance for each 3D vertex
 for( int i=0;i<numVerts;i++ )
 {
  double vert = closerWall.Get(i);
  double xCloser = X->Get(vert);
  double yCloser = Y->Get(vert);
  double zCloser = Z->Get(vert);
  double aux = distance( X->Get(i),Y->Get(i),Z->Get(i),
	                   xCloser,yCloser,zCloser );
  wallDistance.Set(i,aux);
 }
 wallDistance = wallDistance.adimensionalize();

 // considering that the smaller distance is equivalent to the wall edge
 // length (triEdge[0])
 //double minWallDistance = triEdge[0];

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->phyBounds.at(i).compare(5,6,"Normal") != 0 )
  {
   double P0x = surfMesh->X.Get(i);
   double P0y = surfMesh->Y.Get(i);
   double P0z = surfMesh->Z.Get(i);

   double sumEdgeLength = 0;
   int listSize = neighbourPoint->at(i).size();
   list<int> plist = neighbourPoint->at(i);
   for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
   {
	double P1x = surfMesh->X.Get(*vert);
	double P1y = surfMesh->Y.Get(*vert);
	double P1z = surfMesh->Z.Get(*vert);

	double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
	sumEdgeLength += edgeLength;
   }
   convC.Set(i,sumEdgeLength/listSize);
  }
  else // apply same strategy for non-boundary vertics if it is
       // NormalU,V,W c.c.
  {
  double ratio = triEdge[0]*wallDistance.Get(i)*2.5;
  double aux = triEdge[0]+ratio;
   convC.Set(i,aux);
  }
 }

 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double ratio = triEdge[0]*wallDistance.Get(i)*2.5;
  double aux = triEdge[0]+ratio;
  convC.Set(i,aux);
 }
}

void Helmholtz3D::init2Bubbles()
{
 init();

 /* loop at surfMesh: for each vertex, an average value is set based on
  * the umbrella operator (neighbors) for distance. Thus, each vertex
  * will have an associated distanced based on such an average distance.
  * */
 convC.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  double P0x = surfMesh->X.Get(i);
  double P0y = surfMesh->Y.Get(i);
  double P0z = surfMesh->Z.Get(i);

  double sumEdgeLength = 0;
  int listSize = neighbourPoint->at(i).size();
  list<int> plist = neighbourPoint->at(i);
  for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
  {
   double P1x = surfMesh->X.Get(*vert);
   double P1y = surfMesh->Y.Get(*vert);
   double P1z = surfMesh->Z.Get(*vert);

   double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
   sumEdgeLength += edgeLength;
  }
  convC.Set(i,sumEdgeLength/listSize);
 }

 clVector* vertIdRegion = m->getVertIdRegion();
 double minEdge = *min_element(triEdge.begin(),triEdge.end());
 for( int i=surfMesh->numVerts;i<numVerts;i++ )
 {
  double radius = 0.5; // sphere radius
  // outside mesh
  if( heaviside->Get(i) < 0.5 ) 
  {
   double factor = minEdge/triEdge[vertIdRegion->Get(i)];
   if( interfaceDistance->Get(i) < 0.6*radius )
   {
	double aux = triEdge[vertIdRegion->Get(i)]*factor*0.4;
	convC.Set(i,aux);
   }
   else
   {
	double aux = triEdge[vertIdRegion->Get(i)];
	convC.Set(i,aux);
   }
  }
  else                         // inside mesh
  {
   double aux = triEdge[vertIdRegion->Get(i)];
   convC.Set(i,aux);
  }
 }
}

void Helmholtz3D::assemble()
{
 int i,j,ii,jj;
 int v[NUMGLEU];
 double aux;
 clMatrix KcMat( numVerts,numVerts );
 clMatrix McMat( numVerts,numVerts );

 FEMLinElement3D linElem(*X,*Y,*Z);

 for( int mele=0;mele<numElems;mele++ )
 {
  for( int n=0;n<NUMGLEU;n++ )
   v[n] = (int) IEN->Get(mele,n);

  double dif = 1.0;

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

void Helmholtz3D::setCRHS()
{
 vcc = convC;
}

void Helmholtz3D::matMountC()
{
 /* k=5;       /\     more diffusion
  * k=4;      /||\
  * k=3;       ||     moderate diffusion 
  * k=1;      \||/
  * k=0.7;     \/     less diffusion
  * */
 for( int i=0;i<numVerts;i++ )
 {
  McLumped.Set( i,-1.0);
 }
 matc = McLumped - (k * Kc);

}

void Helmholtz3D::setBC()
{
 cc.Dim(numVerts);
 for( int i=0;i<surfMesh->numVerts;i++ )
 {
  if( surfMesh->Marker.Get(i) < 0.5 &&
	  surfMesh->phyBounds.at(i).compare(5,6,"Normal") != 0)
  {
   idbcc.AddItem(i);

   double P0x = surfMesh->X.Get(i);
   double P0y = surfMesh->Y.Get(i);
   double P0z = surfMesh->Z.Get(i);

   double sumEdgeLength = 0;
   int listSize = neighbourPoint->at(i).size();
   list<int> plist = neighbourPoint->at(i);
   for(list<int>::iterator vert=plist.begin();vert!=plist.end();++vert )
   {
	double P1x = surfMesh->X.Get(*vert);
	double P1y = surfMesh->Y.Get(*vert);
	double P1z = surfMesh->Z.Get(*vert);

	double edgeLength = distance(P0x,P0y,P0z,P1x,P1y,P1z);
	sumEdgeLength += edgeLength;
   }
   cc.Set(i,sumEdgeLength/listSize);
  }
 }
}

void Helmholtz3D::setUnCoupledCBC()
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

void Helmholtz3D::unCoupledC()
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

void Helmholtz3D::getModel3DAttrib(Model3D &_m)
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
 neighbourPoint = m->getNeighbourPoint();

 triEdge = m->getTriEdge();
 boundaryVert = m->getBoundaryVert();
 edgeSize = m->getEdgeSize();
}

void Helmholtz3D::allocateMemoryToAttrib()
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

 // closer
 closerWall.Dim( numVerts );
 wallDistance.Dim( numVerts );

 // auxiliar vectors
 ipc.Dim( numVerts,1 );

 // oldSol vectors
 cSolOld.Dim( numVerts );
}

clVector* Helmholtz3D::getCSol(){return &cSol;}
void Helmholtz3D::setSolver(Solver *s){solver = s;}
void Helmholtz3D::setModel3DEdgeSize(){ m->setEdgeSize(cSol); }

void Helmholtz3D::saveVTK( const char* _dir,const char* _filename, int _iter )
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
 vtkFile << setprecision(0) << fixed; 
 for( int i=0;i<numElems;i++ )
   vtkFile << "4 " << IEN->Get(i,0) << " "  
            	   << IEN->Get(i,1) << " " 
				   << IEN->Get(i,2) << " " 
				   << IEN->Get(i,3) << endl;
 vtkFile << endl;

 vtkFile <<  "CELL_TYPES " << numElems << endl;
 for( int i=0;i<numElems;i++ )
   vtkFile << "10 ";

 vtkFile << endl;

 vtkFile << "POINT_DATA " << numVerts << endl;
 vtkFile << "SCALARS boundary double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << cc.Get(i) << endl;

 vtkFile << endl;
 vtkFile << "SCALARS distance double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << wallDistance.Get(i) << endl;

 vtkFile << endl;
 vtkFile << "SCALARS RHS double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << convC.Get(i) << endl;

 vtkFile << endl;
 vtkFile << "SCALARS solution double" << endl;
 vtkFile << "LOOKUP_TABLE default"  << endl;

 vtkFile << setprecision(10) << scientific;
 for( int i=0;i<numVerts;i++ )
  vtkFile << cSol.Get(i) << endl;

 vtkFile << endl;

 vtkFile.close();
}

void Helmholtz3D::saveChordalEdge( const char* _dir,const char* _filename, int _iter )
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
 double nPoints = 1000;
 clVector xVert(nPoints);
 clVector yVert(nPoints);
 clVector zVert(nPoints);

 xVert.SetAll( (X->Max()+X->Min())/2.0 );
 yVert.SetAll( (Y->Max()+Y->Min())/2.0 );
 for( int i=0;i<nPoints;i++ )
 {
  double dz = i * ( (Z->Max()-Z->Min()) )/(nPoints-1);
  double pos = Z->Min()+dz;
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

void Helmholtz3D::setk(double _k){k = _k;}
double Helmholtz3D::getk(){return k;}


void Helmholtz3D::setCloserWall()
{
 // procurando vertices da parede que nao sao de c.c. de simetria
 // (NormalU,V,W)
 clVector boundary(0);
 for( int i=0;i<surfMesh->numVerts;i++ )
  if( surfMesh->Marker(i) < 0.5 && 
	  surfMesh->phyBounds.at(i).compare(5,6,"Normal") != 0)
   boundary.AddItem(i);

 clVector xBoundary( boundary.Dim() );
 clVector yBoundary( boundary.Dim() );
 clVector zBoundary( boundary.Dim() );
 for( int i=0;i<boundary.Dim();i++ )
 {
  int surfaceNode = boundary.Get(i);
  xBoundary.Set(i,X->Get( surfaceNode ));
  yBoundary.Set(i,Y->Get( surfaceNode ));
  zBoundary.Set(i,Z->Get( surfaceNode ));
 }

 clVector xVert(numVerts);
 clVector yVert(numVerts);
 clVector zVert(numVerts);
 X->CopyTo(0,xVert);
 Y->CopyTo(0,yVert);
 Z->CopyTo(0,zVert);

 closerWall = dsearchn(xBoundary,yBoundary,zBoundary,xVert,yVert,zVert);
 for( int i=0;i<closerWall.Dim();i++ )
 {
  int aux = closerWall.Get(i);
  closerWall.Set(i,boundary.Get(aux)); // alterando os valores de closer(i)
 }
}

bool Helmholtz3D::isInsideLiquidFilm(int _node, double _thickness)
{
 double x1 = X->Get(_node);
 double y1 = Y->Get(_node);
 double z1 = Z->Get(_node);

 int closer = closerWall.Get(_node);
 double x2 = X->Get(closer);
 double y2 = Y->Get(closer);
 double z2 = Z->Get(closer);
 
 return distance(x1,y1,z1,x2,y2,z2) < _thickness;
}



