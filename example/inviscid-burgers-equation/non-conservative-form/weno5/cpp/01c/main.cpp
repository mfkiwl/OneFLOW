#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args));
}

class vec1d
{
public:
    std::vector<double> data;
    int ist = 0;
public:
    void Allocate( int ist, int ied, double value = 0 )
    {
        int nelement = ied - ist + 1;
        this->data.resize( nelement, value );
        this->ist = ist;
    }
    std::size_t size()
    {
        return this->data.size();
    }

    double operator [] ( int i ) const
    {
        return data[ i - ist ];
    }

    double & operator [] ( int i )
    {
        return data[ i - ist ];
    }
    
    vec1d & operator = ( const vec1d & rhs )
    {
        if ( this == & rhs ) return * this;
        this->data = rhs.data;

        return * this;
    }
};

void DumpCsvFile( const std::string & filename, vec1d & x, std::vector<vec1d> & u );
double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int nx, vec1d & u, vec1d & f );
void wenoR( int nx, vec1d & u, vec1d & f );
void rhs( int nx, double dx, vec1d & u, vec1d & r );
void numerical( int nx, int ns, int nt, double dx, double dt );

double wcL( double v1, double v2, double v3, double v4, double v5 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 1.0e-1 / ( SQR( eps + s1 ) );
    double c2 = 6.0e-1 / ( SQR( eps + s2 ) );
    double c3 = 3.0e-1 / ( SQR( eps + s3 ) );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    // candiate stencils
    double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
    double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
    double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;

    // reconstructed value at interface
    double f = ( w1 * q1 + w2 * q2 + w3 * q3 );

    return f;
}

double wcR( double v1, double v2, double v3, double v4, double v5 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 3.0e-1 / SQR( eps + s1 );
    double c2 = 6.0e-1 / SQR( eps + s2 );
    double c3 = 1.0e-1 / SQR( eps + s3 );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    // candiate stencils;
    double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0;
    double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
    double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;

    // reconstructed value at interface
    double f = ( w1 * q1 + w2 * q2 + w3 * q3 );

    return f;
}

//-----------------------------------------------------------------------------
// WENO reconstruction for upwind direction (positive; left to right)
// u(i): solution values at finite difference grid nodes i = 1,...,N+1
// f(j): reconstructed values at nodes j = i+1/2; j = 1,...,N
//-----------------------------------------------------------------------------
void wenoL( int nx, vec1d & u, vec1d & f )
{
    for ( int i = 0; i <= nx - 1; ++ i )
    {
        double v1 = u[ i - 2 ];
        double v2 = u[ i - 1 ];
        double v3 = u[ i ];
        double v4 = u[ i + 1 ];
        double v5 = u[ i + 2 ];
        f[ i ] = wcL( v1, v2, v3, v4, v5 );
    }
}

//-----------------------------------------------------------------------------
// CRWENO reconstruction for downwind direction (negative; right to left)
// u(i): solution values at finite difference grid nodes i = 1,...,N+1
// f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
//-----------------------------------------------------------------------------
void wenoR( int nx, vec1d & u, vec1d & f )
{
    for ( int i = 1; i <= nx; ++ i )
    {
        double v1 = u[ i - 2 ];
        double v2 = u[ i - 1 ];
        double v3 = u[ i ];
        double v4 = u[ i + 1 ];
        double v5 = u[ i + 2 ];
        f[ i ] = wcR( v1, v2, v3, v4, v5 );
    }
}

void boundary( int nx, vec1d & u )
{
    //left bc
    int i = 0;
    u[ i ] = 0.0;
    u[ i - 2 ] = 3.0 * u[ i ] - 2.0 * u[ i + 1 ];
    u[ i - 1 ] = 2.0 * u[ i ] - u[ i + 1 ];
    //right bc
    i = nx;
    u[ i ] = 0.0;
    u[ i + 1 ] = 2.0 * u[ i ] - u[ i - 1 ];
    u[ i + 2 ] = 3.0 * u[ i ] - 2.0 * u[ i - 1 ];
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs( int nx, double dx, vec1d & u, vec1d & r )
{
    vec1d uL;
    uL.Allocate( 0, nx - 1, 0 );

    vec1d uR;
    uR.Allocate( 0, nx, 0 );

    wenoL( nx, u, uL );
    wenoR( nx, u, uR );

    for ( int i = 1; i < nx; ++ i )
    {
        if ( u[ i ] >= 0.0 )
        {
            r[ i ] = - u[ i ] * ( uL[ i ] - uL[ i - 1 ] ) / dx;
        }
        else
        {
            r[ i ] = - u[ i ] * ( uR[ i + 1 ] - uR[ i ] ) / dx;
        }
    }
}

//-----------------------------------------------------------------------------
// Compute numerical solution
//   - Time integration using Runge-Kutta third order
//   - 5th-order Compact WENO scheme for spatial terms
//-----------------------------------------------------------------------------
void numerical( int nx, int ns, int nt, double dx, double dt )
{
    vec1d x;
    x.Allocate( 0, nx, 0 ); // npoints = N = nx + 1

    int ighost = 2;
    //int ighost = 0;
    int ist = 0 - ighost;
    int ied = nx + ighost;

    std::vector<vec1d> u( ns + 1 );

    for ( int i = 0; i < u.size(); ++ i )
    {
        u[ i ].Allocate( ist, ied, 0 );
    }

    vec1d un;
    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step

    vec1d ut;
    ut.Allocate( ist, ied, 0 ); // temporary array during RK3 integration

    vec1d r;
    r.Allocate( 0, nx - 1, 0 );

    int k = 0; // record index
    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );

    for ( int i = 0; i < nx + 1; ++ i )
    {
        x[ i ] = dx * ( i );
        un[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
    }

    boundary( nx, un );
    u[ k ] = un; //store solution at t = 0;

    for ( int it = 1; it <= nt; ++ it )
    {
        rhs( nx, dx, un, r );
        for ( int i = 1; i < nx; ++ i )
        {
            ut[ i ] = un[ i ] + dt * r[ i ];
        }
        boundary( nx, ut );
        rhs( nx, dx, ut, r );

        for ( int i = 1; i < nx; ++ i )
        {
            ut[ i ] = 0.75 * un[ i ] + 0.25 * ut[ i ] + 0.25 * dt * r[ i ];
        }
        boundary( nx, ut );
        rhs( nx, dx, ut, r );

        for ( int i = 1; i < nx; ++ i )
        {
            un[ i ] = ( 1.0 / 3.0 ) * un[ i ] + ( 2.0 / 3.0 ) * ut[ i ] + ( 2.0 / 3.0 ) * dt * r[ i ];
        }
        boundary( nx, un );

        if ( it % freq == 0 )
        {
            k = k + 1;
            u[ k ] = un; 
            std::print( "k={}, ns={}\n", k, ns );
        }
    }

    DumpCsvFile( "field_final.csv", x, u );
}

void DumpCsvFile( const std::string &filename, vec1d &x, std::vector<vec1d> &u )
{
    std::fstream file;
    file.open( filename.c_str(), std::fstream::out );
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::string str = {};
        str += std::format( "{:.16f} ", x[ i ] );
        int nj = u.size();
        for ( int j = 0; j < nj; ++ j )
        {
            double value = u[ j ][ i ];
            str += std::format( "{:.16f}", value );
            if ( j != nj - 1 )
            {
                str += " ";
            }
        }
        std::format_to(std::ostream_iterator<char>(file), "{}\n", str );
    }
    file.close();
}


int main( int argc, char ** argv )
{
    int nx = 200;
    int ns = 10;
    double dt = 0.0001;
    double tm = 0.25;

    double dx = 1.0 / nx;
    int nt = std::round( tm / dt );
    double ds = tm / ns;

    std::print( "nx={}\n", nx );
    std::print( "ns={}\n", ns );
    std::print( "dt={}\n", dt );
    std::print( "tm={}\n", tm );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );

    numerical( nx, ns, nt, dx, dt );

    return 0;
}
