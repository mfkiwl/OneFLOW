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
#include "Constant.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <map>


BeginNameSpace( ONEFLOW )

#define IMPLEMENT_GRID_CLONE( TYPE ) \
Grid * Clone() const { return new TYPE(); }
//Grid * Clone() const { return new TYPE( * this ); }

#define REGISTER_GRID( TYPE ) \
    Grid * TYPE ## _myClass = \
        Grid::Register( #TYPE, new TYPE() );

class DataBook;
class NodeMesh;
class InterFace;
class SlipFace;
class DataBase;
class IFaceLink;

class Grid
{
public:
    Grid();
    virtual ~Grid();
public:
    virtual Grid * Clone() const = 0;
public:
    static Grid * SafeClone( const std::string & type );
    static Grid * Register( const std::string & type, Grid * clone );
    static std::map < std::string, Grid * > * classMap;
public:
    std::string name;
    int dimension;
    int type, level;
    int id, localId;
    int nNodes;
    int nFaces, nCells;
    int nBFaces;
    int nIFaces;
    int volBcType;
    NodeMesh * nodeMesh;
    InterFace * interFace;
    SlipFace * slipFace;
    DataBase * dataBase;
public:
    DataBase * GetDataBase() { return dataBase; };
public:
    void BasicInit();
    void Free();
    virtual void Init();
public:
    bool IsOneD();
    bool IsTwoD();
    bool IsThreeD();
public:
    virtual void ReadGrid ( std::fstream & file ) {};
    virtual void WriteGrid( std::fstream & file ) {};
    virtual void Decode( DataBook * databook ){};
    virtual void Encode( DataBook * databook ){};
    virtual void ReadGrid( DataBook * databook ){};
    virtual void WriteGrid( DataBook * databook ){};
    virtual void ModifyBcType( int bcType1, int bcType2 ) {};
    virtual void GenerateLgMapping( IFaceLink * iFaceLink ){};
    virtual void ReGenerateLgMapping( IFaceLink * iFaceLink ){};
    virtual void UpdateOtherTopologyTerm( IFaceLink * iFaceLink ){};
public:
    virtual void GetMinMaxDistance( Real & dismin, Real & dismax ) {};
    virtual void CalcMetrics() {};
};

EndNameSpace
