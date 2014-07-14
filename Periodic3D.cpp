/**
 * \file:    Periodic3D.cpp
 * 
 * \author:  Gustavo Charles P. de Oliveira (tavolesliv@gmail.com)
 * \author:  Norberto Mangiavacchi (norberto@uerj.br)
 *
 * \date     April 9th, 2013, 17:43 AM
 *
 * \brief    Class liable to set PBC (Periodic Boundary Conditions)
 *           over boundaries of domains of simulation accompanying
 *           \link <Periodic3Dh> \endlink.
 *
 * \version  1.0
 *
 * \detailed Inclusion of periodicity will require a mesh periodic
 * 			 geometrically.
 * 			 
 */


#include "Periodic3D.h"

using namespace std;


Periodic3D::Periodic3D()
{
    nyPointsL = 0;
    nyPointsR = 0;
    VecXMin.Dim(0);
    VecXMax.Dim(0);

	MasterIndices.resize(0);
	SlaveIndices.resize(0);

}


Periodic3D::Periodic3D(Model3D &_M3D)
{
    M3DPtr = &_M3D;
    NumVerts = M3DPtr->getNumVerts();
    NumNodes = M3DPtr->getNumNodes();
    XPtr = M3DPtr->getX();
    YPtr = M3DPtr->getY();
    ZPtr = M3DPtr->getZ();
	MasterIndicesPtr = M3DPtr->getPbcIndicesLeft();
	SlaveIndicesPtr = M3DPtr->getPbcIndicesRight();
}


Periodic3D::~Periodic3D() {}

/** \brief It reads VTK mesh file and saves two vectors
 *  containing the pairing indices. Also, it reorders the pairs,
 *  otherwise.
 *
 *  \return{ void }
 *  
 */
void Periodic3D::MountPeriodicVectorsNew()
{
   nyPointsL = MasterIndicesPtr->size();
   nyPointsR = SlaveIndicesPtr->size();
   SetIndicesVector(MasterIndicesPtr,SlaveIndicesPtr);
 
   vector<int> xAux(nyPointsL);
   int ibR = 0;
   double YRight = 0.0;
   double ZRight = 0.0;
   /** Test of dimension:
	*  
	* \remark{It is already done in setGenericBC(). However, it only
	* checks the spatial correspondence based on simple extrusion.}
	*
	**/
   if ( nyPointsL != nyPointsR )
   { 
	   cout << "Error: vectors for PBC don't match dimensions!" << endl;
	   cout << nyPointsL << "!=" << nyPointsR << endl;
	   cerr << "PBC implementation is invalid! Stopping..." << endl;
	   exit(1);
   }
   else /* If test of dimension is OK, mounts. */
   {
	   cout << "Mounting vectors of periodic nodes...\n" << endl;
  
	   for ( int i = 0; i < nyPointsL; i++ )
	   {
		   int ibL = MasterIndices.at(i);

		   double YLeft = YPtr->Get(ibL);
		   double ZLeft = ZPtr->Get(ibL);
		   
		   for ( int j = 0; j < nyPointsR; j++ )
		   {
			   ibR = SlaveIndices.at(j);
			   YRight = YPtr->Get(ibR);
			   ZRight = ZPtr->Get(ibR);
			   
			   double deltaY = fabs( YLeft - YRight);
			   double deltaZ = fabs( ZLeft - ZRight);

			   // if finds pair, jumps out
			   if ( ( deltaY < 1E-10 )
			     && ( deltaZ < 1E-10 ) )
				 break;

			   
		   }

		   xAux.at(i) = ibR; // reorders after 'break'
	   }
   }

   SlaveIndices = xAux;
   
   /*Printing pairs */
   cout << "\t >>>> Periodic Pairing Mounted <<<<"  << endl;
   for (int i = 0; i < nyPointsL; ++i)
   {
	  int ibL = MasterIndices.at(i);
	  int ibR = SlaveIndices.at(i);
	  
	  cout << "(" << i << ")\t Index ibL: " << ibL << "\t pairs with \t Index ibR: " << ibR << endl;
   }

} /* End of function */


void Periodic3D::SetIndicesVector(vector<int>* _master, vector<int>* _slave)
{
	MasterIndices.resize(nyPointsL);
	SlaveIndices.resize(nyPointsL);

	for ( int i = 0; i != nyPointsL; ++i )
	{
		MasterIndices.at(i) = _master->at(i);
		SlaveIndices.at(i) = _slave->at(i);
	}

}


/** \brief Sets the Periodic Boundary Condition for scalar through
 * copy process. 
 * 
 *  \param[in] & _Scalar
 *  \param[in] & _VecXMin;
 *  \param[in] & _VecXMax;
 *  \param[in] int L: number of boundary points;
 *  \param[in] string direction: direction of copy
 *
 *  \note See method \link "::"<SetVelocityPBC> \endlink. 
 *
 */
void Periodic3D::SetPureScalarPBC(clVector &_Scalar, clVector &_VecXMin, clVector &_VecXMax, int L, string direction)
{
    if ( direction == "RL" ) 
	{
	 	cout << "Copying scalar field from RIGHT to LEFT..." << endl;
    	
		int right;
		double sRight;
		
		for (right = 0; right < L; right++)
		{
		  int index2 = _VecXMax.Get(right);
		  int left = right;
		  int index = _VecXMin.Get(left);
		  sRight = _Scalar.Get(index2);
		  _Scalar.Set(index,sRight);
		}
	}
	else
	{	
	 	cout << "Copying scalar field from LEFT to RIGHT..." << endl;
	    
		int left;
    	double sRight;
    
    	for (left = 0; left < L; left++) 
    	{
          int index = _VecXMin.Get(left); 
          int right = left;
          int index2 = _VecXMax.Get(right); 
          sRight = _Scalar.Get(index2); 
          _Scalar.Set(index,sRight); 
    	}
	}

} /* End of function */


void Periodic3D::SetPureScalarPBCNew(clVector &_Scalar, vector<int>* _master, vector<int>* _slave, int L, string direction)
{
    if ( direction == "RL" ) 
	{
	 	cout << "Copying scalar field from RIGHT to LEFT..." << endl;
    	
		int right;
		double sRight;
		
		for (right = 0; right < L; right++)
		{
		  int index2 = _slave->at(right);
		  int left = right;
		  int index = _master->at(left);
		  sRight = _Scalar.Get(index2);
		  _Scalar.Set(index,sRight);
		}
	}
	else
	{	
	 	cout << "Copying scalar field from LEFT to RIGHT..." << endl;
	    
		int left;
    	double sRight;
    
    	for (left = 0; left < L; left++) 
    	{
          int index = _master->at(left); 
          int right = left;
          int index2 = _slave->at(right); 
          sRight = _Scalar.Get(index2); 
          _Scalar.Set(index,sRight); 
    	}
	}

} /* End of function */


/** \brief Sets the PBC for velocity through copy process.
 *  \param[in] & _uVelocity;
 *  \param[in] & _vVelocity;
 *  \param[in] & _wVelocity;
 *  \param[in] & _VecXMin;
 *  \param[in] & _VecXMax;
 *  \param[in] int L: number of boundary points;
 *  \param[in] string direction: if equal to "RL", attributes  V(X_R):= V(X_L);
 *                       otherwise, attributes V(X_L):=V(X_R);
 *                       where V = (u,v,w);
 *                       X_L = (x_L,y_L,z_L) in Gamma_{left} and
 *                       X_R = (x_R,y_R,z_R) in Gamma_{right} 
 *
 */
void Periodic3D::SetVelocityPBC(clVector &_uVelocity, clVector &_vVelocity,
                                clVector &_wVelocity, clVector &_VecXMin,
                                clVector &_VecXMax, int L, string direction)
{
	if ( direction  == "RL" ) /* copy from right to left */
	{
		cout << "Copying velocity field from RIGHT to LEFT..." << endl;
        
		int right;
        double uVelRight, vVelRight, wVelRight;
        
        for (right = 0; right < L; right++) 
        {
            int index2 = _VecXMax.Get(right); 
            int left = right;
            int index = _VecXMin.Get(left); 
            uVelRight = _uVelocity.Get(index2); 
            vVelRight = _vVelocity.Get(index2); 
            wVelRight = _wVelocity.Get(index2); 
            _uVelocity.Set(index,uVelRight); 
            _vVelocity.Set(index,vVelRight); 
            _wVelocity.Set(index,wVelRight); 

        }
	}
	else /* copy from left to right */
    {
        cout << "Copying velocity field from LEFT to RIGHT..." << endl;
        
        int left;
        double uVelLeft, vVelLeft, wVelLeft;
        
        for (left = 0; left < L; left++) 
        {
            int index = _VecXMin.Get(left); 
            int right = left;
            int index2 = _VecXMax.Get(right); 
            uVelLeft = _uVelocity.Get(index); 
            vVelLeft = _vVelocity.Get(index); 
            wVelLeft = _wVelocity.Get(index); 
            _uVelocity.Set(index2,uVelLeft); 
            _vVelocity.Set(index2,vVelLeft); 
            _wVelocity.Set(index2,wVelLeft); 
        }
	}
} /* End of function */


/* Idem, but with container <vector> */
void Periodic3D::SetVelocityPBCNew(clVector &_uVelocity, clVector &_vVelocity,
                                clVector &_wVelocity, vector<int>* _master,
                                vector<int>* _slave, int L, string direction)
{
	if ( direction  == "RL" ) /* copy from right to left */
	{
		cout << "Copying velocity field from RIGHT to LEFT..." << endl;
        
		int right;
        double uVelRight, vVelRight, wVelRight;
        
        for (right = 0; right < L; right++) 
        {
            int index2 = _slave->at(right); 
            int left = right;
            int index = _master->at(left); 
            uVelRight = _uVelocity.Get(index2); 
            vVelRight = _vVelocity.Get(index2); 
            wVelRight = _wVelocity.Get(index2); 
            _uVelocity.Set(index,uVelRight); 
            _vVelocity.Set(index,vVelRight); 
            _wVelocity.Set(index,wVelRight); 

        }
	}
	else /* copy from left to right */
    {
        cout << "Copying velocity field from LEFT to RIGHT..." << endl;
        
        int left;
        double uVelLeft, vVelLeft, wVelLeft;
        
        for (left = 0; left < L; left++) 
        {
            int index = _master->at(left); 
            int right = left;
            int index2 = _slave->at(right); 
            uVelLeft = _uVelocity.Get(index); 
            vVelLeft = _vVelocity.Get(index); 
            wVelLeft = _wVelocity.Get(index); 
            _uVelocity.Set(index2,uVelLeft); 
            _vVelocity.Set(index2,vVelLeft); 
            _wVelocity.Set(index2,wVelLeft); 
        }
	}
} /* End of function */


/** \brief Sets the PBC for pressure through copy process.
 *
 * \param[in] & _Pressure
 * \param[in] & _VecXMin;
 * \param[in] & _VecXMax;
 * \param[in] int L: number of boundary points;
 * \param[in] string direction: direction of copy
 *
 * \note Method must be discussed.
 *
 * \note See method \link "::"<SetVelocityPBC> \endlink. 
 *
 */
void Periodic3D::SetPurePressurePBC(clVector &_Pressure, clVector &_VecXMin, clVector &_VecXMax, int L, string direction)
{

 	if ( direction == "RL" )
	{
	 	cout << "Copying pressure field from RIGHT to LEFT..." << endl;
	 	
		int right;
		double pRight;

		for (right = 0; right < L; right++)
		{
		   int index2 = _VecXMax.Get(right);
		   int left = right;
		   int index = _VecXMin.Get(left);
		   pRight = _Pressure.Get(index2);
		   _Pressure.Set(index,pRight);
		}
	}

	else
	{
	 	cout << "Copying pressure field from LEFT to RIGHT..." << endl;

	 	int left;
    	double pRight;
    
    	for (left = 0; left < L; left++) 
    	{
          int index = _VecXMin.Get(left); 
          int right = left;
          int index2 = _VecXMax.Get(right);
          pRight = _Pressure.Get(index2); 
          _Pressure.Set(index,pRight); // fills entries of p(X_L) indices with p(X_R)
       }
	}

} /* End of function */


/* Idem, but with container <vector> */
void Periodic3D::SetPurePressurePBCNew(clVector &_Pressure, vector<int>* _master, vector<int>* _slave, 
                                                                                 int L, string direction)
{

 	if ( direction == "RL" )
	{
	 	cout << "Copying pressure field from RIGHT to LEFT..." << endl;
	 	
		int right;
		double pRight;

		for (right = 0; right < L; right++)
		{
		   int index2 = _slave->at(right);
		   int left = right;
		   int index = _master->at(left);
		   pRight = _Pressure.Get(index2);
		   _Pressure.Set(index,pRight);
		}
	}

	else
	{
	 	cout << "Copying pressure field from LEFT to RIGHT..." << endl;

	 	int left;
    	double pRight;
    
    	for (left = 0; left < L; left++) 
    	{
          int index = _master->at(left); 
          int right = left;
          int index2 = _slave->at(right);
          pRight = _Pressure.Get(index2); 
          _Pressure.Set(index,pRight); // fills entries of p(X_L) indices with p(X_R)
       }
	}

} /* End of function */


/** \brief Sets PBC for pressure as p(X_L) = p(X_R) + c, where c is a scalar
 * 		   (jump), X_L = (x_L,y_L) in Gamma_{left} and 
 * 		   X_R = (x_R,y_R) in Gamma_{right} 
 *
 *  \note Method must be discussed.
 *
 */
void Periodic3D::SetJumpPressurePBC(clVector &_Pressure, clVector &_VecXMin, clVector &_VecXMax, int L, double jump)
{
    int left;
    double pRight;
    
    for (left = 0; left < L; left++) 
    {
        int index = _VecXMin.Get(left); 
        int right = left;
        int index2 = _VecXMax.Get(right); 
        pRight = _Pressure.Get(index2); 
        
        double aux = pRight + jump;
        _Pressure.Set(index,aux); // set pressure jump between X_L,X_R
     
    }
    
} /* End of function */


/** Get blocks */
int Periodic3D::GetNyPointsL() { return nyPointsL; };
int Periodic3D::GetNyPointsR() { return nyPointsR; };

clVector* Periodic3D::GetVecXMin() { return &VecXMin; };
clVector* Periodic3D::GetVecXMax() { return &VecXMax; };

vector<int>* Periodic3D::GetMasterIndices() { return &MasterIndices; };
vector<int>* Periodic3D::GetSlaveIndices() { return &SlaveIndices; };

