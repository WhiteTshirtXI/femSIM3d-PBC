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
    nyPointsM = 0;
    YLeftBoundaryVector.Dim(0);
    YRightBoundaryVector.Dim(0);
    ZLeftBoundaryVector.Dim(0);
    ZRightBoundaryVector.Dim(0);
    VecXMin.Dim(0);
    VecXMax.Dim(0);
    VecXMid.Dim(0);

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
 *  \attention Mesh extrusion must be over axis X. 
 *  
 */
void Periodic3D::MountPeriodicVectors(Model3D &_M3D)
{
    const double XMin = XPtr->Min();
    const double XMax = XPtr->Max();
    
    VecXMin = XPtr->FindValue(XMin); // returns vector with indices xMin
    VecXMax = XPtr->FindValue(XMax); // returns vector with indices xMax
    VecXMid = XPtr->FindComplementaryValues(XMin,XMax);
    nyPointsL = VecXMin.Dim(); 
    nyPointsR = VecXMax.Dim(); 
    nyPointsM = VecXMid.Dim();
   

    /* Test of dimension */
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
        		
		// Eliminating corner points of periodicity
		for (int i = 0; i < 4; ++i)
		{
		  VecXMin.Delete(0);
		  VecXMax.Delete(0);
		  nyPointsL--;
		}
		
		YLeftBoundaryVector.Dim(nyPointsL);
        YRightBoundaryVector.Dim(nyPointsL);
        ZLeftBoundaryVector.Dim(nyPointsL);
        ZRightBoundaryVector.Dim(nyPointsL);

        clVector XAux;
        XAux = VecXMax;
              
        for ( int i = 0; i < nyPointsL; i++ )
        {
            int ibL = VecXMin.Get(i);

			double YLeft = YPtr->Get(ibL);
            YLeftBoundaryVector.Set(i,YLeft);
            double ZLeft = ZPtr->Get(ibL);
            ZLeftBoundaryVector.Set(i,ZLeft);
            
			double deltaYOld = 1000.0;
			double deltaZOld = 1000.0;

            for ( int j = 0; j < nyPointsL; j++ )
            {
                int ibR = XAux.Get(j);
                double YRight = YPtr->Get(ibR);
                double ZRight = ZPtr->Get(ibR);
                
				double deltaY = fabs( YLeft - YRight );
				double deltaZ = fabs( ZLeft - ZRight );

                if ( ( deltaY < deltaYOld ) &&
				     ( deltaZ < deltaZOld ) ) // reorders pairing
                {
					VecXMax.Set(i,ibR);
					deltaYOld = deltaY;
					deltaZOld = deltaZ;
                    YRightBoundaryVector.Set(i,YRight);
                    ZRightBoundaryVector.Set(i,ZRight);
                }
            
            }
            
        }
    
    }
    
    /* Printing pairs */
	cout << "\t >>>> Periodic Pairing Mounted <<<<" << endl;
	for (int i = 0; i < nyPointsL; ++i)
	{
		int ibL = VecXMin.Get(i);
		int ibR = VecXMax.Get(i);
		
		cout << "(" << i << ")\t Index ibL: " << ibL << "\t pairs with \t Index ibR: " << ibR << endl;

	}	

    //VerifyParallelismY(YLeftBoundaryVector,YRightBoundaryVector);
    //VerifyParallelismZ(ZLeftBoundaryVector,ZRightBoundaryVector);
    
} /* End of function */


/** \brief It reads VTK mesh file and saves two vectors
 *  containing the pairing indices. Also, it reorders the pairs,
 *  otherwise.
 *
 *  \attention Mesh extrusion must be over axis X. 
 *  
 */
void Periodic3D::MountPeriodicVectorsNew(Model3D &_M3D)
{
    nyPointsL = MasterIndicesPtr->size();
    nyPointsR = SlaveIndicesPtr->size();
	SetIndicesVector(MasterIndicesPtr,SlaveIndicesPtr);

    /* Test of dimension */
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
        
		vector<int> aux (nyPointsL);
		aux = SlaveIndices;

        for ( int i = 0; i < nyPointsL; i++ )
        {
            int ibL = MasterIndices.at(i);

			double YLeft = YPtr->Get(ibL);
            double ZLeft = ZPtr->Get(ibL);
            
			double deltaYOld = 1000.0;
			double deltaZOld = 1000.0;

            for ( int j = 0; j < nyPointsL; j++ )
            {
                int ibR = aux.at(j);
                double YRight = YPtr->Get(ibR);
                double ZRight = ZPtr->Get(ibR);
                
				double deltaY = fabs( YLeft - YRight);
				double deltaZ = fabs( ZLeft - ZRight);

				/* Might work with overwriting from setGenericBC() */
                if ( ( deltaY < deltaYOld ) &&
				     ( deltaZ < deltaZOld ) ) // reorders pairing
                {
				 	//SlaveIndices.at(i) = ibR;
					deltaYOld = deltaY;
					deltaZOld = deltaZ;
                } 
            }
        }
    }
    
    /* Printing pairs */
	cout << "\t >>>> Periodic Pairing Mounted <<<<" << endl;
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

	for ( size_t i = 0; i != nyPointsL; ++i )
	{
		MasterIndices.at(i) = _master->at(i);
		SlaveIndices.at(i) = _slave->at(i);
	}

}


/** \brief Given two vectors u,v in \Rn, it verifies if their
 *  y-coordinates are equal, for the purposes of this class. 
 *  On the other hand, it is equivalent to verify if u,v are equal. 
 *
 */
bool Periodic3D::VerifyParallelismY(clVector _u, clVector _v)
{
   int j;

    string Warning = "Warning! Left and right y-components are not equal. \n"
    "PBC should not be implemented.";
    
    for ( int i = 0; i < _u.Dim(); i++ )
    {
        double yU = _u.Get(i);
        double yV = _v.Get(i);
        
        if ( yU != yV )
        {
            cout << Warning << endl;
            cout << "Entry: " << i << endl;
            cout << "yU = " << setprecision(16) << yU << "; yV = " << setprecision(16) << yV << endl;
            
			j = i;
			ForcedParallelismY(j,_u,_v);
        	
			return false;
		}

    }
    
	 cout << "Test was successfully completed for y! PBC can be implemented." << endl;
    
	 return true;
    
} /* End of function */


/** \brief Forces Y-parallelism not reached by \link "::"<VerifyParallelismY>
 * \endlink, if any.
 *
 *  \attention Remark: It was observed that meshes generated by GMSH
 *  				   software have limited precision. Hence, the test 
 *  				   of paired coordinates might fail when the numbers
 *  				   differs from such digits. Hence, the enforcement
 *  				   was suggested.
 *
 */
void Periodic3D::ForcedParallelismY(int j, clVector _u, clVector _v)
{
 	cout << "Forcing paralellism on y..." << endl;

	double yU = _u.Get(j);
 	_v.Set(j,yU); 
	double yV2 = _v.Get(j);
            
	cout << "yU = " << setprecision(16) << yU << "; yV = " << setprecision(16) << yV2 << endl;

	 cout << "Test was successfully completed for y! PBC can be implemented." << endl;

} /* End of function */


/** \brief Given two vectors u,v in \Rn, it verifies if their
 *  y-coordinates are equal, for the purposes of this class. 
 *  On the other hand, it is equivalent to verify if u,v are equal. 
 *
 */
bool Periodic3D::VerifyParallelismZ(clVector _u, clVector _v)
{
	int j;
    
    string Warning = "Warning! Left and right z-components are not equal. \n"
    "PBC should not be implemented.";
    
    for ( int i = 0; i < _u.Dim(); i++ )
    {
        double zU = _u.Get(i);
        double zV = _v.Get(i);
        
        if ( zU != zV )
        {
            cout << Warning << endl;
            cout << "Entry: " << i << endl;
            cout << "zU = " << setprecision(16) << zU << "; zV = " << setprecision(16) << zV << endl;
			
			j = i;
			ForcedParallelismZ(j,_u,_v);
            
			return false;
        }
        
    }
    
    cout << "Test was successfully completed for z! PBC can be implemented." << endl;
    
    return true;
    
} /* End of function */


/** \brief Forces Z-parallelism not reached by \link <VerifyParallelismZ>
 * \endlink.
 *
 *  \attention Remark: It was observed that meshes generated by GMSH
 *  				   software have limited precision. Hence, the test 
 *  				   of paired coordinates might fail when the numbers
 *  				   differs from such digits. Hence, the enforcement
 *  				   was suggested.
 *
 */
void Periodic3D::ForcedParallelismZ(int j, clVector _u, clVector _v)
{
 	cout << "Forcing paralellism on z..." << endl;

	double zU = _u.Get(j);
 	_v.Set(j,zU); 
	double zV2 = _v.Get(j);
            
	cout << "zU = " << setprecision(16) << zU << "; zV = " << setprecision(16) << zV2 << endl;

	 cout << "Test was successfully completed for z! PBC can be implemented." << endl;

} /* End of function */



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


void Periodic3D::SetPureScalarPBCVector(clVector &_Scalar, vector<int> &_master, vector<int> &_slave, int L, string direction)
{
    if ( direction == "RL" ) 
	{
	 	cout << "Copying scalar field from RIGHT to LEFT..." << endl;
    	
		int right;
		double sRight;
		
		for (right = 0; right < L; right++)
		{
		  int index2 = _slave.at(right);
		  int left = right;
		  int index = _master.at(left);
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
          int index = _master.at(left); 
          int right = left;
          int index2 = _slave.at(right); 
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
void Periodic3D::SetVelocityPBCVector(clVector &_uVelocity, clVector &_vVelocity,
                                clVector &_wVelocity, vector<int> &_master,
                                vector<int> &_slave, int L, string direction)
{
	if ( direction  == "RL" ) /* copy from right to left */
	{
		cout << "Copying velocity field from RIGHT to LEFT..." << endl;
        
		int right;
        double uVelRight, vVelRight, wVelRight;
        
        for (right = 0; right < L; right++) 
        {
            int index2 = _slave.at(right); 
            int left = right;
            int index = _master.at(left); 
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
            int index = _master.at(left); 
            int right = left;
            int index2 = _slave.at(right); 
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
void Periodic3D::SetPurePressurePBCVector(clVector &_Pressure, vector<int> &_master, vector<int> &_slave, 
                                                                                 int L, string direction)
{

 	if ( direction == "RL" )
	{
	 	cout << "Copying pressure field from RIGHT to LEFT..." << endl;
	 	
		int right;
		double pRight;

		for (right = 0; right < L; right++)
		{
		   int index2 = _slave.at(right);
		   int left = right;
		   int index = _master.at(left);
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
          int index = _master.at(left); 
          int right = left;
          int index2 = _slave.at(right);
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


/**  \brief It extracts mesh vertices (without centroid) over
 *   \Omega - \Gammas_{left,right} (but includes top, bottom) before
 *   setting centroids. 
 *   
 *   \note Approach to be used when considering future implementation
 *         of PBC with full reassembling of the global matrices. 
 *
 */
void Periodic3D::ExtractMiddleVerts(const char *MeshFileName)
{
 	char auxstr[255];
	int i;
	double coords[3];
	clVector XAux;
    
	ifstream vtkFile( MeshFileName, ios::in);
    
	if ( !vtkFile )
	{
		cerr << "Missing 'file'.vtk. Nothing to read..." << endl;
		exit(1);
	}
    
	while ( ( !vtkFile.eof() && strcmp( auxstr,"POINTS" ) != 0 ) )
	{
	 	vtkFile >> auxstr;
	}
    
	if ( !vtkFile.eof() )
	{
		vtkFile >> NumVerts;
		vtkFile >> auxstr;
        
		XAux.Dim(NumVerts);
		
		for ( i = 0; i < NumVerts; i++ )
		{
            vtkFile >> coords[0];
            XAux.Set(i,coords[0]);
		}
    
    }
    
    for ( i = 0; i < NumVerts; i++ )
    {
	  	if ( XAux.Get(i) > XAux.Min() && XAux.Get(i) < XAux.Max() )
		{
		 	VecXMidVerts.Append(i);
		}
    }
} /* End of function */


/** Get blocks */
int Periodic3D::GetNyPointsL() { return nyPointsL; };
int Periodic3D::GetNyPointsR() { return nyPointsR; };
int Periodic3D::GetNyPointsM() { return nyPointsM; };
int Periodic3D::GetNumVertsMid() { return NumVertsMid; };

clVector* Periodic3D::GetVecXMin() { return &VecXMin; };
clVector* Periodic3D::GetVecXMax() { return &VecXMax; };
clVector* Periodic3D::GetVecXMid() { return &VecXMid; };
clVector* Periodic3D::GetVecXMidVerts() { return &VecXMidVerts; };

vector<int> Periodic3D::GetMasterIndices() { return MasterIndices; };
vector<int> Periodic3D::GetSlaveIndices() { return SlaveIndices; };

