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
#include "HXClone.h"
#include "Task.h"
#include "Plate.h"
#include <sstream>
#include <fstream>


BeginNameSpace( ONEFLOW )

DEFINE_DATA_CLASS( CreateTurbPlateTask );
void SetTurbPlateTask();

class UnsGrid;

class TurbVelCut : public CuttingClass
{
public:
    TurbVelCut();
    ~TurbVelCut();
public:
    void Dump();
    void Dump( LamData * lamData, std::fstream & file, int axis );
    void DumpDetail( LamData * lamData, std::fstream & file, int axis );
public:
    void DumpNormal();
    void DumpDetail();
};

class TurbFriCut : public CuttingClass
{
public:
    TurbFriCut();
    ~TurbFriCut();
public:
    void Dump();
    void Dump( LamData * lamData, std::fstream & file, int axis );
};

class TurbFlatPlateTask : public Task
{
public:
    TurbFlatPlateTask();
    ~TurbFlatPlateTask() override;
public:
    CuttingClass * velCut;
    CuttingClass * friCut;
public:
    void Run() override;
    void OutProfile( CuttingClass * cut );
};

EndNameSpace
