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

#include "UnsteadyTaskReg.h"
#include "InterField.h"
#include "ActionState.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Restart.h"
#include "Unsteady.h"
#include "UnsteadyImp.h"
#include "Update.h"
#include "FieldWrap.h"
#include "FieldAlloc.h"
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtil.h"
#include "FieldRecord.h"
#include "UVisualize.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"
#include <map>
#include <iostream>


BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterUnsteadyTask )

void RegisterUnsteadyTask()
{
    REGISTER_DATA_CLASS( UpdateUnsteadyFlow );
    REGISTER_DATA_CLASS( CalcUnsteadyCriterion );
}


void UpdateUnsteadyFlow( StringField & data )
{
    int sTid = SolverState::tid;

    Unsteady * unsteady = CreateUnsteady( sTid );
    unsteady->UpdateUnsteady( sTid );
    delete unsteady;
}

void CalcUnsteadyCriterion( StringField & data )
{
    int sTid = SolverState::tid;
    Unsteady * unsteady = CreateUnsteady( sTid );
    unsteady->CalcUnsteadyCriterion();
    delete unsteady;
}

EndNameSpace
