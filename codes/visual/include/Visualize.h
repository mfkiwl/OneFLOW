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
#include "HXArray.h"
#include <sstream>


BeginNameSpace( ONEFLOW )

class Visualize
{
public:
    Visualize();
    virtual ~Visualize();
public:
    virtual void Visual(){};
};

class Plot
{
public:
    Plot();
    ~Plot();
public:
    static std::ostringstream * oss;
    static int nWords;
    static void DumpField( RealField & field );
    static void DumpField( IntField & l2g, RealField & x );
    static void DumpFaceNodeNumber( LinkField & f2n );
    static void DumpFaceNodeLink( LinkField & f2n );
    static void DumpFaceElementLink( IntField & elementId, int nElem );
};

int GetTotalNumFaceNodes( LinkField & f2n );

EndNameSpace
