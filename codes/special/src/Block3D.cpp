/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2025 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlkMesh.h"
#include "Block3D.h"
#include "MLine.h"
#include "MDomain.h"

#include "Prj.h"
#include "Dimension.h"
#include "BlockFaceSolver.h"
#include "Transfinite.h"
#include "StrGrid.h"
#include "NodeMesh.h"
#include <fstream>
#include <iomanip>


BeginNameSpace( ONEFLOW )

Block3D::Block3D()
{
    int nMDomain = 6;
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = new MDomain();
        mDomain->pos = iMDomain;
        mDomain->coorMap = & this->coorMap;
        mDomainList.push_back( mDomain );
    }

    this->AddLocalPt( 1, 4, 8, 5 );
    this->AddLocalPt( 2, 3, 7, 6 );
    this->AddLocalPt( 1, 2, 6, 5 );
    this->AddLocalPt( 4, 3, 7, 8 );
    this->AddLocalPt( 1, 2, 3, 4 );
    this->AddLocalPt( 5, 6, 7, 8 );
}

Block3D::~Block3D()
{
    DeletePointer( mDomainList );
    DeletePointer( facelist );
}

void Block3D::Alloc()
{
    AllocateVector( x3d, ni, nj, nk );
    AllocateVector( y3d, ni, nj, nk );
    AllocateVector( z3d, ni, nj, nk );
}

int Block3D::GetNSubDomain()
{
    int nSubDomain = 0;
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        nSubDomain += mDomain->GetNsubDomain();
    }
    return nSubDomain;
}

void Block3D::ConstructTopo()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->ConstructMultiDomainTopo();
    }
    std::map< int, IntSet > p2dMap;
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        int nSize = mDomain->candidate_bcpoints.size();
        for ( int i = 0; i < nSize; ++ i )
        {
            int pt = mDomain->candidate_bcpoints[ i ];
            std::map< int, IntSet >::iterator iter;
            iter = p2dMap.find( pt );
            if ( iter == p2dMap.end() )
            {
                IntSet iset;
                iset.insert( iMDomain );
                p2dMap.insert( std::pair< int, IntSet >( pt, iset ) );
            }
            else
            {
                iter->second.insert( iMDomain );
            }
        }
        int kkk = 1;
    }

    IntField ctrl_points;
    std::map< int, IntSet >::iterator iter;
    for ( iter = p2dMap.begin(); iter != p2dMap.end(); ++ iter )
    {
        if ( iter->second.size() == 3 )
        {
            ctrl_points.push_back( iter->first );
        }
    }

    //for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    //{
    //    MDomain * mDomain = mDomainList[ iMDomain ];
    //    mDomain->CalcDomainCtrlPoints( ctrl_points );
    //}

    int p1, p2, p3, p4, p5, p6, p7, p8;

    GetCornerPoint( p1, 0, 2, 4 );
    GetCornerPoint( p2, 1, 2, 4 );
    GetCornerPoint( p3, 1, 3, 4 );
    GetCornerPoint( p4, 0, 3, 4 );
    GetCornerPoint( p5, 0, 2, 5 );
    GetCornerPoint( p6, 1, 2, 5 );
    GetCornerPoint( p7, 1, 3, 5 );
    GetCornerPoint( p8, 0, 3, 5 );

    this->controlpoints.push_back( p1 );
    this->controlpoints.push_back( p2 );
    this->controlpoints.push_back( p3 );
    this->controlpoints.push_back( p4 );
    this->controlpoints.push_back( p5 );
    this->controlpoints.push_back( p6 );
    this->controlpoints.push_back( p7 );
    this->controlpoints.push_back( p8 );

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDomainCtrlPoints( this->controlpoints, this->localpt[iMDomain ] );
    }
    this->CalcBlkDim();
}

void Block3D::GetCornerPoint( int & pt, int id1, int id2, int id3 )
{
    MDomain * d1 = mDomainList[ id1 ];
    MDomain * d2 = mDomainList[ id2 ];
    MDomain * d3 = mDomainList[ id3 ];

    int nSize = d1->candidate_ctrlpoints.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        int ip = d1->candidate_ctrlpoints[ i ];
        bool flag1 = InArray( ip, d2->candidate_ctrlpoints );
        if ( ! flag1 ) continue;
        bool flag2 = InArray( ip, d3->candidate_ctrlpoints );
        if ( ! flag2 ) continue;
        pt = ip;
        break;
    }
}

void Block3D::SetInterfaceBc()
{
    int nFaces = this->facelist.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        Face2D * face2d = this->facelist[ iFace ];
        int domain_id = face2d->face_id;
        if ( face2d->bcType == -1 )
        {
            face2d->t = new Face2D();

            BlkF2C & face_struct = blkFaceSolver.face2Block[ domain_id ];
            int n_neibor = face_struct.cellList.size();
            int blk1 = face_struct.cellList[ 0 ] - 1;
            int blk2 = face_struct.cellList[ 1 ] - 1;

            int tblk = -1;
            if ( this->blk_id == blk1 )
            {
                tblk = blk2;
            }
            else
            {
                tblk = blk1;
            }
            Face2D * facet = blkFaceSolver.GetBlkFace( tblk, domain_id );
            face2d->t->bcType = tblk + 1;
            face2d->t->st = facet->st;
            face2d->t->ed = facet->ed;
        }
    }
}

void Block3D::CalcBlkDim()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDim2D();
    }

    MDomain * d0 = mDomainList[ 0 ];
    MDomain * d2 = mDomainList[ 2 ];
    //d0: nj,nk
    //d2: ni,nk

    this->ni = d2->ni;
    this->nj = d0->ni;
    this->nk = d0->nj;

    IntField iList, jList, kList;
    Add( iList, jList, kList, 1, 1, 1 );
    Add( iList, jList, kList, ni, 1, 1 );
    Add( iList, jList, kList, ni, nj, 1 );
    Add( iList, jList, kList, 1, nj, 1 );
    Add( iList, jList, kList, 1, 1, nk );
    Add( iList, jList, kList, ni, 1, nk );
    Add( iList, jList, kList, ni, nj, nk );
    Add( iList, jList, kList, 1, nj, nk );

    for ( int iPoint = 0; iPoint < iList.size(); ++ iPoint )
    {
        int pt = this->controlpoints[ iPoint ];
        int i = iList[ iPoint ];
        int j = jList[ iPoint ];
        int k = kList[ iPoint ];
        CalcCoor c;
        c.SetCoor( i, j, k );
        this->coorMap.insert( std::pair<int, CalcCoor>( pt, c ) );
    }

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcCoor();
    }

    CreateFaceList();

    int kkk = 1;
}

void Block3D::CreateFaceList()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CreateInpFaceList( facelist );
    }

    int kkk = 1;

}


void Block3D::CreateBlockMesh()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->SetBlkBcMesh( this );
    }

    this->GenerateBlockMesh();
}

void Block3D::GenerateBlockMesh()
{
    std::fstream file;
    Prj::OpenPrjFile( file, "grid/blkfaceplot.dat", std::ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ 0 ] << " " << y3d[ i ][ j ][ 0 ] << " " << z3d[ i ][ j ][ 0 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ nk-1 ] << " " << y3d[ i ][ j ][ nk-1 ] << " " << z3d[ i ][ j ][ nk-1 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ 0 ][ k ] << " " << y3d[ i ][ 0 ][ k ] << " " << z3d[ i ][ 0 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ nj - 1 ][ k ] << " " << y3d[ i ][ nj - 1 ][ k ] << " " << z3d[ i ][ nj - 1 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ 0 ][ j ][ k ] << " " << y3d[ 0 ][ j ][ k ] << " " << z3d[ 0 ][ j ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ ni - 1 ][ j ][ k ] << " " << y3d[ ni-1 ][ j ][ k ] << " " << z3d[ ni-1 ][ j ][ k ] << "\n";
        }
    }
    Prj::CloseFile( file );
    TransfiniteInterpolation( x3d, ni, nj, nk );
    TransfiniteInterpolation( y3d, ni, nj, nk );
    TransfiniteInterpolation( z3d, ni, nj, nk );
    //for ( int k = 1; k < nk - 1; ++ k )
    //{
    //    for ( int j = 1; j < nj - 1; ++ j )
    //    {
    //        for ( int i = 1; i < ni - 1; ++ i )
    //        {
    //            Real d = ( x3d[ i ][ j ][ nk - 1 ] - x3d[ i ][ j ][ 0 ] ) / ( nk - 1 );
    //            x3d[ i ][ j ][ k ] = x3d[ i ][ j ][ 0 ] + k * d;
    //            y3d[ i ][ j ][ k ] = y3d[ i ][ j ][ 0 ] + k * d;
    //            z3d[ i ][ j ][ k ] = z3d[ i ][ j ][ 0 ] + k * d;
    //        }
    //    }
    //}


    Prj::OpenPrjFile( file, "grid/blkplot.dat", std::ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << ", K = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            for ( int i = 0; i < ni; ++ i )
            {
                file << x3d[ i ][ j ][ k ] << " " << y3d[ i ][ j ][ k ] << " " << z3d[ i ][ j ][ k ] << "\n";
            }
        }
    }
    Prj::CloseFile( file );

    Prj::OpenPrjFile( file, "grid/blkfaceplot111.dat", std::ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ 0 ] << " " << y3d[ i ][ j ][ 0 ] << " " << z3d[ i ][ j ][ 0 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ nk-1 ] << " " << y3d[ i ][ j ][ nk-1 ] << " " << z3d[ i ][ j ][ nk-1 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ 0 ][ k ] << " " << y3d[ i ][ 0 ][ k ] << " " << z3d[ i ][ 0 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ nj - 1 ][ k ] << " " << y3d[ i ][ nj - 1 ][ k ] << " " << z3d[ i ][ nj - 1 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ 0 ][ j ][ k ] << " " << y3d[ 0 ][ j ][ k ] << " " << z3d[ 0 ][ j ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ ni - 1 ][ j ][ k ] << " " << y3d[ ni-1 ][ j ][ k ] << " " << z3d[ ni-1 ][ j ][ k ] << "\n";
        }
    }
    Prj::CloseFile( file );
    int kkk = 1;
}

void Block3D::FillStrGrid( Grid * gridIn, int iZone )
{
    StrGrid * grid = StrGridCast( gridIn );
    int ni = this->ni;
    int nj = this->nj;
    int nk = this->nk;

    grid->id = iZone;
    grid->ni = ni;
    grid->nj = nj;
    grid->nk = nk;
    grid->SetBasicDimension();
    grid->nodeMesh->CreateNodes( grid->nNodes );

    grid->SetLayout();

    SetGridXYZ3D( grid, x3d, y3d, z3d );
}

EndNameSpace
