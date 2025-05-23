#pragma once
#include <vector>

class Zone;
class Grid;
class Region;

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un;
    std::vector<double> r;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid );
public:
    void FTCS( Zone * zone );
    void CN( Zone * zone );
    void ICP( Zone * zone );
    void RungeKutta( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
public:
    void Rhs( std::vector<double> & u, std::vector<double> & r );

    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );

    void PhysicalBoundary( Zone * zone );

    void PostProcess( Grid * grid );
    void UpdateOldField();
};
