// =================================================================== //
// this is file Model3D.cpp, created at 23-Ago-2007                    //
// maintained by Gustavo Rabello dos Anjos                             //
// e-mail: gustavo.rabello@gmail.com                                   //
// =================================================================== //


#include "Model3D.h"
#include "compare.h"

using namespace std;

Model3D::Model3D()
{
 IEN;       // matrix de triangulacao
 X;         // coordenada X de todos os pontos da malha
 Y;         // coordenada Y de todos os pontos da malha
 Z;         // coordenada Z de todos os pontos da malha
 uc;        // vetor velocidade na direcao U com condicao de contorno
 vc;        // vetor velocidade na direcao V com condicao de contorno
 wc;        // vetor velocidade na direcao W com condicao de contorno
 pc;        // vetor pressao P com condicao de contorno
 idbcu;     // vetor de indices de condicao de contorno para U
 idbcv;     // vetor de indices de condicao de contorno para V
 idbcw;     // vetor de indices de condicao de contorno para W
 idbcp;     // vetor de indices de condicao de contorno para P
 idbcc;     // vetor de indices de condicao de contorno para C
 outflow;   // vetor de indices de condicao de outflow para pressao
 numVerts;  // numero total de vertices na malha
 numNodes;  // numero total de nos na malha (inclui centroide)
 numElems;  // numero total de elementos da malha
 numGLEU;   // numero de graus de liberdade para velocidade
 numGLEP;   // numero de graus de liberdade para pressao
 neighbourElem; // lista de vizinhos de cada no
 //faceFace;  // matrix com lista de elementos vizinhos
 freeFace;  // matriz com lista de faces de contorno
 //mapViz;    // matriz com identificacao de vizinhos por elementos
 oFace;     // matriz com identificacao de face oposta ao vertice
}

Model3D::~Model3D(){}

bool Model3D::readVTK( const char* filename )
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

   IEN.Dim(numElems,4+1); // "+1" para o centroide
   numGLEP = 4; // triangulo linear
   numGLEU = 5; // elemento MINI
   numGLEC = 4; // elemento linear

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

 numNodes = numVerts + numElems;
 uc.Dim(numNodes);
 vc.Dim(numNodes);
 wc.Dim(numNodes);
 pc.Dim(numVerts);
 cc.Dim(numVerts);
 outflow.Dim(numNodes,1); // usado no metodo Galerkin

 return true;
} // fim do metodo vtkRead

/**
 * @brief metodo para leitura de arquivo do tipo BC para condicoes de
 * contorno. O arquivo a ser lido deve conter todos os nos de condicao
 * de contorno contendo colunas de vertices, colunas de u,v e p com
 * numeracao 1 para condicao do tipo Dirichlet e 2 para Newmann e seu
 * respectivo valor
 *
 * @return verdadeiro ou falso dependendo do sucesso da leitura
 **/
bool Model3D::readBC( const char* filename )
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
 return true;
} // fim do metodo readBC

void Model3D::clearBC()
{
 uc.SetAll(0.0);
 vc.SetAll(0.0);
 wc.SetAll(0.0);
 pc.SetAll(0.0);
}

void Model3D::setStep(int nX,int nY,int nZ)
{
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
	aux = j*dy;
	Y.Set(count,aux);
	aux = k*dz;
	Z.Set(count,aux);
	count++;
   }
  }
 }
}

void Model3D::setStepBC()
{
 for( int i=0;i<numVerts;i++ )
 {
  if( (X.Get(i)==X.Min()) || (Y.Get(i)==Y.Min()) || (Y.Get(i)==Y.Max()) )
  {
	idbcc.AddItem(i);
   idbcu.AddItem(i);
   idbcv.AddItem(i);

   uc.Set(i,0.0);
   vc.Set(i,0.0);
   if( (X.Get(i)==X.Min()) && (Y.Get(i)>(Y.Max()/2.0)) && (Y.Get(i)<Y.Max()) )
   {

	uc.Set(i,1.0);
	cc.Set(i,1.0);
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

void Model3D::setDisk(int nLados1Poli,int nCircMax,int nZ)
{
 real dr = 1;
 real  r = dr;
 int j = 0;
 real z = 0;
 real dl = ( (2*PI)/nLados1Poli)*dr;
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
   if( theta >= 2*PI ) break;
  }
  r=r+dr;
 }

 rMax = r - dr;
 j=0;
 z=0;
 for( int jz=1;jz<=nZ;jz++ )
 {
  //if( jz<=nZ/2 ) dz=0.1;
  //else 
  if( exp(z/nZ) <= exp(2.0/3.0) ) dz=exp(z/nZ);
  else dz = exp(2.0/3.0);
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
   aux = 1.0;
   cc.Set(i,aux);
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

   //uc.Set(i,radius*2.6505291e-02); // Sc = 2000 Z = 4
   //vc.Set(i,radius*3.6178208e-02); // Sc = 2000 Z = 4
   //wc.Set(i,-8.2427145e-01); // Sc = 2000 Z = 4

   uc.Set(i,radius*1.3690760e-04); // Sc = 2000 Z = 10
   vc.Set(i,radius*1.7819422e-04); // Sc = 2000 Z = 10
   wc.Set(i,-8.8528405e-01); // Sc = 2000 Z = 10
   
   //wc.Set(i,-8.8559326E-01); // Sc = 2000
   //wc.Set(i,-9.1044679e-01); // Sc = 10
   //wc.Set(i,-9.2281563e-01); // Sc = 5
  }

  if( Z.Get(i) == Z.Min() )
  {
   idbcu.AddItem(i);
   idbcv.AddItem(i);
   idbcw.AddItem(i);
   idbcc.AddItem(i);

   omega=1.0;

   aux = (-1.0)*Y.Get(i)*omega;
   uc.Set(i,aux);
   aux = X.Get(i)*omega;
   vc.Set(i,aux);
   aux = 0.0;
   wc.Set(i,aux);

   aux = 1.0;
   cc.Set(i,aux);
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

void Model3D::setBubbleBubbleBC()
{
 //real eps2 = 1E-5;
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
	  (Y.Get(i)-yCenter)*(Y.Get(i)-yCenter) +
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
  if( X.Get(i) == X.Max() )
  {
   idbcp.AddItem(i);
   pc.Set(i,0.0);
  }
 }
} // fecham metodo setBubbleBubbleBC

void Model3D::setCubeCubeBC(real param)
{
 real eps2 = 1.0E-10;
 real aux;
 xCenter = 1.5;
 yCenter = 1.5;
 zCenter = 1.5;
 cout << "xCenter = " << xCenter << endl;
 cout << "yCenter = " << yCenter << endl;
 cout << "zCenter = " << yCenter << endl;
 cout << "bubbleRadius = " << bubbleRadius << endl;

 // fora da bolha
 cc.SetAll(0.0);

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

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
  }
  // na interface
  if( (X.Get(i)<1.9+eps2) && (X.Get(i)>1.1-eps2) && 
      (Y.Get(i)<1.9+eps2) && (Y.Get(i)>1.1-eps2) && 
      (Z.Get(i)<1.9+eps2) && (Z.Get(i)>1.1-eps2) )
  {
   cc.Set(i,0.5);
  }
  // dentro da bolha
  if( (X.Get(i)<1.9-eps2) && (X.Get(i)>1.1+eps2) && 
      (Y.Get(i)<1.9-eps2) && (Y.Get(i)>1.1+eps2) && 
      (Z.Get(i)<1.9-eps2) && (Z.Get(i)>1.1+eps2) )
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


 string file = "./estadoBase/F_";
 file += _filename;
 file += ".dat";
 const char* filenameF = file.c_str();
 ifstream fFile( filenameF,ios::in );

 file = "./estadoBase/G_"; 
 file += _filename;
 file += ".dat";
 const char* filenameG = file.c_str();
 ifstream gFile( filenameG,ios::in );

 file = "./estadoBase/H_";
 file += _filename;
 file += ".dat";
 const char* filenameH =  file.c_str();
 ifstream hFile( filenameH,ios::in );

 file = "./estadoBase/dHdZ_";
 file += _filename;
 file += ".dat";
 const char* filenamedHdZ = file.c_str();
 ifstream dHdZFile( filenamedHdZ,ios::in );

 if( !fFile || !gFile || !hFile || !dHdZFile )
 {
  cerr << "Esta faltando algum arquivo do estado base para NuCte!" << endl;
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

   //uc.Set(i,radius*4.8979239E-02); // Z=4
   //vc.Set(i,radius*7.9371092E-02); // Z=4
   //wc.Set(i,-9.1936908E-01); // Z=4

   //uc.Set(i,radius*6.8574756E-03); // Z=6
   //vc.Set(i,radius*1.0335366E-02); // Z=6
   //wc.Set(i,-0.10196142E+01); // Z=6

   uc.Set(i,radius*0.11735664E-03); // Z=10
   vc.Set(i,radius*0.17501519E-03); // Z=10
   wc.Set(i,-0.10193840E+01); // Z=10
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
   idbcc.AddItem(i);

   aux = 0.0;
   uc.Set(i,aux);
   vc.Set(i,aux);
   wc.Set(i,aux);
   cc.Set(i,1.0);
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
   idbcp.AddItem(i);
   idbcu.AddItem(i);
   idbcv.AddItem(i);

   // cc = 0.5 -> noh da superficie
   cc.Set(i,0.5); // para funcionamento do ALE

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
 real Red = 100;
 real factorz = 1.0/(Z.Max()-Z.Min());

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
 for( int i=0;i<numNodes;i++ )
 {
  aux = (uc.Get(i)/rMax)*Red;
  uc.Set(i,aux);
  aux = (vc.Get(i)/rMax)*Red;
  vc.Set(i,aux);
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

void Model3D::setCentroid()
{
 real volume;
 V.Dim(numElems);
 real centroidX,centroidY,centroidZ;
 int v1,v2,v3,v4,v5,v[5];

 // realocando vetores...
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

 for( int i=0;i<numElems;i++ )
 {
  v[0]=v1=(int)IEN.Get(i,0);
  v[1]=v2=(int)IEN.Get(i,1);
  v[2]=v3=(int)IEN.Get(i,2);
  v[3]=v4=(int)IEN.Get(i,3);


  // este procedimento foi validado! Correto!
  volume = (-1.0/6.0) * (+1*( (X.Get(v2)*Y.Get(v3)*Z.Get(v4)) 
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
  V.Set(i,volume);

  if( fabs(volume)<1.0E-10)
   cerr << "tetraedro singular, verificar a qualidade da malha!" << endl;

  if( volume<0.0 )
  {
   IEN.Set(i,1,v[2]);
   IEN.Set(i,2,v[1]);
  };

  v[4]=v5=numVerts+i;

  IEN.Set(i,4,v5);
  centroidX = ( X.Get(v1)+X.Get(v2)+X.Get(v3)+X.Get(v4) )*0.25;
  centroidY = ( Y.Get(v1)+Y.Get(v2)+Y.Get(v3)+Y.Get(v4) )*0.25;
  centroidZ = ( Z.Get(v1)+Z.Get(v2)+Z.Get(v3)+Z.Get(v4) )*0.25;
  X.Set(v5,centroidX);
  Y.Set(v5,centroidY);
  Z.Set(v5,centroidZ);
 
 }
}

void Model3D::setNeighbour()
{
 neighbourElem.resize (numVerts);
 for( int i=0;i<numElems;i++ )
  for( int j= 0;j<numGLEP;j++ )
   neighbourElem.at( (int)IEN.Get(i,j) ).push_back(i);
}

void Model3D::setVertNeighbour()
{
 // cria lista de vizinhos para toda a malha
 int v1,v2,v3,v4;
 listElem plist;
 list<int>::iterator mele;
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
// xSurfaceViz -> coordenada x dos vizinhos
// ySurfaceViz -> coordenada y dos vizinhos
void Model3D::setSurface()
{
 int surfaceNode;
 listElem plist;
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

 /*        - nome: xSurfaceViz,ySurfaceViz
           - definicao: vetor de listas de coordenadas de vizinhos 
		                de vertices na interface 
		   - obs: no caso 2D n=2. no caso 3D n pode variar de vertice 
		          para vertice.

     ---   +---+----+-----+----+
      |    | b | a1 | ... | an |   a1...an = identificacao da coordenada  X e Y
      |    +---+----+-----+----+             dos vizinhos (lista)
    c |    | b | a1 | ... | an |   b = identificacao das coordenadas do   
      |    +---+----+-----+----+       do vertice de trabalho
      |    | b | a1 | ... | an |   c = numero de vertices na interface
     ---   +---+----+-----+----+

 */

 // procurando vertices da bolha
 clVector surfaceAux = cc==0.5;
 surface = surfaceAux.Find();
 clVector nonSurfaceAux = cc!=0.5;
 nonSurface = nonSurfaceAux.Find();

 // dimensionando vetores
 surfaceViz.resize ( surface.Dim() );
 xSurfaceViz.resize ( surface.Dim() );
 ySurfaceViz.resize ( surface.Dim() );
 zSurfaceViz.resize ( surface.Dim() );

 for( int i=0;i<surface.Dim();i++ )
 {
  surfaceNode = surface.Get(i);
  plist = neighbourVert.at(surfaceNode); 

  // adicionando no primeiro elemento da lista o valor do vertice de
  // trabalho
  surfaceViz.at( i ).push_back(surfaceNode); 
  xSurfaceViz.at( i ).push_back( X.Get(surfaceNode) );
  ySurfaceViz.at( i ).push_back( Y.Get(surfaceNode) );
  zSurfaceViz.at( i ).push_back( Z.Get(surfaceNode) );
  // procura dos vertices adjacentes na interface (um de cada lado)
  for( vert=plist.begin();vert!=plist.end();++vert )
  {
   if( cc.Get(*vert) == 0.5 )   
   {
	surfaceViz.at( i ).push_back(*vert);
	xSurfaceViz.at( i ).push_back( X.Get(*vert) );
	ySurfaceViz.at( i ).push_back( Y.Get(*vert) );
	zSurfaceViz.at( i ).push_back( Z.Get(*vert) );
   }
  }
 }
} // fecha metodo setSurface

void Model3D::setSurfaceFace()
{
 int v1,v2,v3,v4;
 listElem plist;
 list<int>::iterator mele;
 int surfaceNode;
 
 // mapeamento de faces da interface com numeracao
 elemSurface.resize ( numVerts ); 
 // mapeamento de vertices das faces numeradas por elemSurface
 // super dimensionado!!!
 neighbourFaceVert.resize ( 5*numNodes );
 // contador de faces na superficie
 int count = 0; 

 for( int i=0;i<surface.Dim();i++ )
 {
  surfaceNode = surface.Get(i);
  plist = neighbourElem.at(surfaceNode);
  for( mele=plist.begin(); mele != plist.end();++mele )
  {
   v1 = (int) IEN.Get(*mele,0);
   v2 = (int) IEN.Get(*mele,1);
   v3 = (int) IEN.Get(*mele,2);
   v4 = (int) IEN.Get(*mele,3);
   // pegando o elemento com 3 vertices na superfice e 1 dentro da bolha
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	if( v1 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v1);
	if( v2 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v2);
	if( v3 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v3);

	//neighbourFaceVert.at( count ).sort();
    count = count + 1;
   }
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	if( v1 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v1);
	if( v2 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v2);
	if( v4 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v4);

	//neighbourFaceVert.at( count ).sort();
    count = count + 1;
   }
   if( cc.Get(v1)==0.5 && cc.Get(v2)==0 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	if( v1 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v1);
	if( v3 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v3);
	if( v4 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v4);

	//neighbourFaceVert.at( count ).sort();
    count = count + 1;
   }
   if( cc.Get(v1)==0 && cc.Get(v2)==0.5 && 
	   cc.Get(v3)==0.5 && cc.Get(v4)==0.5 )
   {
	elemSurface.at( surfaceNode ).push_back(count);
	if( v2 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v2);
	if( v3 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v3);
	if( v4 != surfaceNode )
	 neighbourFaceVert.at( count ).push_back(v4);

	//neighbourFaceVert.at( count ).sort();
    count = count + 1;
   }
  }
  //cout << "---------" << surfaceNode << "------------" << endl;
  //std::ostream_iterator< int > output( cout, " " );
  //std::copy( elemSurface.at(surfaceNode).begin(), elemSurface.at(surfaceNode).end(), output );
  //std::copy( neighbourFaceVert.at(surfaceNode).begin(), neighbourFaceVert.at(surfaceNode).end(), output );
  //cout << endl;
 }
 neighbourFaceVert.resize (count); // trim do vector para numero real de itens
}

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
 neighbourElem.resize (numVerts);
 for( int i=0;i<numElems;i++ )
 {
  for( int j=0;j<numGLEP;j++ )
   neighbourElem.at( (int)IEN.Get(i,j) ).push_back(i);

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
 
 //        - nome: faceFree
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
 int aresta;
 int jAresta;
 int vertOp1,vertOp2;
 int kFace;
 int nele;
 int count;
 listElem plist;
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
 delete faces;

 setVertNeighbour();
 setSurface();
 setSurfaceFace();
}

bool Model3D::testFace(int v1, int v2, int v3, int v4)
{
 real V,Ax1,Ax2,Ay1,Ay2,Az1,Az2;
 real prodEsc;

  V = (-1.0/6.0) * (+1*( (X.Get(v2)*Y.Get(v3)*Z.Get(v4)) 
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

// atraves da matriz freeFace, este metodo cria vetor de vertices
// ordenados e nao repeditos para imposicao das condicoes de contorno
// eh necessario habilitar a alocacao global da matriz freeFace
/*
clVector Model3D::sortFreeVector()
{
 clVector aux1,aux2;
 clVector vertFree,elemFree;

 for( int i=0;i<freeFace.DimI();i++ )
 {
  //aux1.AddItem( freeFace.Get(i,1) );
  aux2.AddItem( (int)freeFace.Get(i,2) );
  aux2.AddItem( (int)freeFace.Get(i,3) );
  aux2.AddItem( (int)freeFace.Get(i,4) );
 }
 //elemFree = aux1.Unique();
 vertFree = aux2.Unique();
 return vertFree;
}
*/

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

void Model3D::getNonZeros()
{
}

clVector Model3D::getX(){ return X; }
clVector* Model3D::getPointerX(){ return &X; }
real Model3D::getMaxX(){ return X.Max(); }
real Model3D::getMinX(){ return X.Min(); }
void Model3D::setX(clVector _X){ X = _X; }
clVector Model3D::getY(){ return Y; }
clVector* Model3D::getPointerY(){ return &Y; }
real Model3D::getMinY(){ return Y.Min(); }
real Model3D::getMaxY(){ return Y.Max(); }
void Model3D::setY(clVector _Y){ Y = _Y; }
clVector Model3D::getZ(){ return Z; }
real Model3D::getMaxZ(){ return Z.Max(); }
real Model3D::getMinZ(){ return Z.Min(); }
clVector* Model3D::getPointerZ(){ return &Z; }
void Model3D::setZ(clVector _Z){ Z = _Z; }
clVector Model3D::getUC(){ return uc; }
clVector* Model3D::getPointerUC(){ return &uc; }
clVector Model3D::getVC(){ return vc; }
clVector* Model3D::getPointerVC(){ return &vc; }
clVector Model3D::getWC(){ return wc; }
clVector* Model3D::getPointerWC(){ return &wc; }
clVector Model3D::getPC(){ return pc; }
clVector* Model3D::getPointerPC(){ return &pc; }
clVector Model3D::getCC(){ return cc; }
clVector* Model3D::getPointerCC(){ return &cc; }
clVector Model3D::getOutflow(){ return outflow; }
clVector* Model3D::getPointerOutflow(){ return &outflow; }
clVector Model3D::getIdbcu(){ return idbcu; }
clVector* Model3D::getPointerIdbcu(){ return &idbcu; }
clVector Model3D::getIdbcv(){ return idbcv; }
clVector* Model3D::getPointerIdbcv(){ return &idbcv; }
clVector Model3D::getIdbcw(){ return idbcw; }
clVector* Model3D::getPointerIdbcw(){ return &idbcw; }
clVector Model3D::getIdbcp(){ return idbcp; }
clVector* Model3D::getPointerIdbcp(){ return &idbcp; }
clVector Model3D::getIdbcc(){ return idbcc; }
clVector* Model3D::getPointerIdbcc(){ return &idbcc; }
clMatrix Model3D::getIEN(){ return IEN; }
clMatrix* Model3D::getPointerIEN(){ return &IEN; }
int Model3D::getNumVerts(){ return numVerts; }
int Model3D::getNumNodes(){ return numNodes; }
int Model3D::getNumElems(){ return numElems; }
int Model3D::getNumGLEU(){ return numGLEU; }
int Model3D::getNumGLEP(){ return numGLEP; }
int Model3D::getNumGLEC(){ return numGLEC; }
//clMatrix Model3D::getMapViz(){ return mapViz; }
//clMatrix Model3D::getFaceFace(){ return faceFace; }
clMatrix Model3D::getFreeFace(){ return freeFace; }
clMatrix Model3D::getOFace(){ return oFace; }
clMatrix* Model3D::getPointerOFace(){ return &oFace; }
real Model3D::getXCenter(){ return xCenter; }
real Model3D::getYCenter(){ return yCenter; }
real Model3D::getZCenter(){ return zCenter; }
real Model3D::getBubbleRadius(){ return bubbleRadius; }
clVector Model3D::getSurface(){ return surface; }
clVector* Model3D::getPointerSurface(){ return &surface; }

