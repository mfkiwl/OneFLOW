#include "Field.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void FieldSub::Init( Grid * grid )
{
    this->InitHeatEquation( grid );
}

void FieldSub::InitHeatEquation( Grid * grid )
{
    this->ni = grid->x.size();
    std::cout << "ni = " << ni << "\n";

    Vec1d & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = dx / 10.0;
    this->total_time = 1.0;
    this->nt = std::round( total_time / dt );

    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->total_time  = " << this->total_time << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    Global::nt = nt;

    this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
    this->beta = this->alpha * dt / ( dx * dx );

    int nghost = 2;
    int ni_total = ni + nghost;

    int N = this->ni;

    int ighost = 1;
    int iist = 0 - ighost;
    int iied = N - 1 + ighost;
    u_e.Allocate( iist, iied );
    u.Allocate( iist, iied );
    un.Allocate( iist, iied );
    res.Allocate( 0, N, 0 ); //N+1

    for ( int i = 0; i < ni; ++ i )
    {
        double xm = x[ i ];
        u_e[ i ] = - std::exp( -total_time ) * std::sin( std::numbers::pi * xm ); //theory solution
        u[ i ] = - std::sin( std::numbers::pi * xm ); //initial condition @ t=0
    }
    int kkk = 1;
}

void FieldSub::FTCS( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = u[ i ] + dt * res[ i ];
    }
}

void FieldSub::CN( Zone * zone )
{
    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    std::vector<double> a( ni );//0:ni-1
    std::vector<double> b( ni );//0:ni-1
    std::vector<double> c( ni );//0:ni-1
    std::vector<double> d( ni );//0:ni-1

    for ( int i = 0; i < ni; ++ i )
    {
        a[ i ] = - rr;
        b[ i ] = 1.0 + 2.0 * rr;
        c[ i ] = - rr;
    }

    for ( int i = 0; i < ni; ++ i )
    {
        d[ i ] = rr * u[ i - 1 ] + ( 1.0 - 2.0 * rr ) * u[ i ] + rr * u[ i + 1 ];
    }

    a[ 0 ] = 0;
    c[ ni - 1 ] = 0;

    d[ 0 ] -= ( - rr ) * u[ 0 ];
    d[ ni - 1 ] -= ( - rr ) * u[ ni + 1 ];

    std::vector<double> values( d.size() );

    thomas_algorithm( a, b, c, d, values );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = values[ i ];
    }
}

void FieldSub::ICP( Zone * zone )
{
    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    std::vector<double> a( ni );//0:ni-1
    std::vector<double> b( ni );//0:ni-1
    std::vector<double> c( ni );//0:ni-1
    std::vector<double> d( ni );//0:ni-1

    for ( int i = 0; i < ni; ++ i )
    {
        a[ i ] = 1.0 / 12.0 - rr;
        b[ i ] = 10.0 / 12.0 + 2.0 * rr;
        c[ i ] = 1.0 / 12.0 - rr;
    }

    a[ 0 ] = 0;
    b[ 0 ] = 1;
    c[ 0 ] = 0;

    a[ ni - 1 ] = 0;
    b[ ni - 1 ] = 1;
    c[ ni - 1 ] = 0;

    for ( int i = 0; i < ni; ++ i )
    {
        double aa = 1.0 / 12.0 + rr;
        double bb = 10.0 / 12.0 - 2.0 * rr;
        double cc = 1.0 / 12.0 + rr;
        d[ i ] = aa * u[ i - 1 ] + bb * u[ i ] + cc * u[ i + 1 ];
    }

    d[ 0 ] = 0;
    d[ ni - 1 ] = 0;

    std::vector<double> values( d.size() );

    thomas_algorithm( a, b, c, d, values );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = values[ i ];
    }
}

void FieldSub::RungeKutta( Zone * zone, int istage )
{
    if ( istage == 0 )
    {
        this->RungeKutta3Stage0( zone );
        return;
    }

    if ( istage == 1 )
    {
        this->RungeKutta3Stage1( zone );
        return;
    }

    if ( istage == 2 )
    {
        this->RungeKutta3Stage2( zone );
        return;
    }
}

void FieldSub::UpdateOldField()
{
    this->un = this->u;
}

void FieldSub::RungeKutta3Stage0( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = u[ i ] + dt * res[ i ];
    }
}

void FieldSub::RungeKutta3Stage1( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u[ i ] + 0.25 * dt * res[ i ];
    }
}

void FieldSub::RungeKutta3Stage2( Zone * zone )
{
    this->Rhs( this->u, this->res );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * res[ i ];
    }
}

void FieldSub::Rhs( Vec1d & u, Vec1d & r )
{
    double coef = this->alpha / ( dx * dx );

    for ( int i = 0; i < ni; ++ i )
    {
        r[ i ] = coef * ( u[ i + 1 ] - 2.0 * u[ i ] + u[ i - 1 ] );
    }
}

void FieldSub::PhysicalBoundary( Zone * zone )
{
    int nbccos = zone->bccos.size();
    for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
    {
        ZoneBc * zonebc = zone->bccos[ ibcco ];
        Region region;
        region.SetRegion( zonebc->pnts );
        Boundary( region, zonebc->bcType );
    }
}

void FieldSub::Boundary( Region &region, int bcType )
{
    if ( bcType == BCInflow )
    {
        this->InflowBc( region );
    }
    else if ( bcType == BCExtrapolate || bcType == BCOutflow )
    {
        this->OutflowBc( region );
    }
}

void FieldSub::InflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        if ( i == 1 )
        {
            idir = -1;
        }
        int ii = i - 1; //index from 0
        int ighost = ii + idir;
        int iinner = ii - idir;
        //this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
        this->u[ ighost ] = - this->u[ iinner ];
        int kkk = 1;
    }
}

void FieldSub::OutflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        if ( i == 1 )
        {
            idir = -1;
        }
        int ii = i - 1; //index from 0
        int ighost = ii + idir;
        int iinner = ii - idir;
        //this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
        this->u[ ighost ] = - this->u[ iinner ];
    }
}
