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


#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )
class NodeMesh
{
public:
    NodeMesh();
    NodeMesh( const NodeMesh & rhs );
    ~NodeMesh();
public:
    RealField xN, yN, zN;
    RealField pmin, pmax;
    bool boxFlag;
public:
    HXSize_t GetNumberOfNodes() { return xN.size(); }
    void CreateNodes( int numberOfNodes );
    void CalcMinMaxBox();
    void AddPoint( Real xp, Real yp, Real zp );
};


EndNameSpace
