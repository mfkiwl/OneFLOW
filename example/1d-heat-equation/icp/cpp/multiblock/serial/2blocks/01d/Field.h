#pragma once
#include <vector>
#include "Vec1d.h"

class Zone;
class Grid;
class Region;

class Field
{
public:
    Field() {}
    virtual ~Field() {};
public:
    virtual void Init( Grid * grid ) {}
    virtual void UpdateOldField() {}
    virtual void FTCS( Zone * zone ) {}
    virtual void CN( Zone * zone ) {}
    virtual void ICP( Zone * zone ) {}
    virtual void PhysicalBoundary( Zone * zone ) {}
    virtual void InterfaceBoundary( Zone * zone ) {}
    virtual void PostProcess( Grid * grid ) {}
    virtual void Rhs( Vec1d & u, Vec1d & r ) {};
public:
    void RungeKutta( Zone * zone, int nStage, int istage );
    void RungeKutta1( Zone * zone, int istage );
    void RungeKutta3( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
public:
    Vec1d u_e;
    Vec1d u, un;
    Vec1d res;
    int ni;
    double dt;
};

