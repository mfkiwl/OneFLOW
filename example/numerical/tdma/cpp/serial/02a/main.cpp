import std;

void Print( std::vector<double> & x, const std::string & name = "vector" )
{
    std::print( "{} = ", name );
    for ( auto v: x )
    {
        std::print( "{} ", v );
    }
    std::println();
}

void ThomasAlgorithm( int N, double * b, double * a, double * c, double * x, double * q )
{
    double *l,*u,*d,*y;
    l = new double[ N ];
    u = new double[ N ];
    d = new double[ N ];
    y = new double[ N ];
    /* LU Decomposition */
    d[ 0 ] = a[ 0 ];
    u[ 0 ] = c[ 0 ];
    for ( int i = 1; i < N - 1; i++ )
    {
        l[ i ] = b[ i ] / d[ i - 1 ];
        d[ i ] = a[ i ] - l[ i ] * u[ i - 1 ];
        u[ i ] = c[ i ];
    }
    l[ N - 1 ] = b[ N - 1 ] / d[ N - 2 ];
    d[ N - 1 ] = a[ N - 1 ] - l[ N - 1 ] * u[ N - 2 ];
    /* Forward Substitution [L][y] = [q] */
    y[ 0 ] = q[ 0 ];
    for ( int i = 1; i < N; i++ )
    {
        y[ i ] = q[ i ] - l[ i ] * y[ i - 1 ];
    }
    /* Backward Substitution [U][x] = [y] */
    x[ N - 1 ] = y[ N - 1 ] / d[ N - 1 ];
    for ( int i = N - 2; i >= 0; i-- )
    {
        x[ i ] = ( y[ i ] - u[ i ] * x[ i + 1 ] ) / d[ i ];
    }
    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    return;
}

void ThomasAlgorithmVersion1( int N, double * b, double * a, double * c, double * x, double * q )
{
    double *l,*u,*d,*y;
    l = new double[ N ];
    u = new double[ N ];
    d = new double[ N ];
    y = new double[ N ];
    /* LU Decomposition */
    d[ 0 ] = a[ 0 ];
    u[ 0 ] = c[ 0 ];

    for ( int i = 1; i < N - 1; i++ )
    {
        l[ i - 1 ] = b[ i - 1  ] / d[ i - 1  ];
        d[ i ] = a[ i ] - l[ i - 1 ] * u[ i - 1 ];
        u[ i ] = c[ i ];
    }
    l[ N - 2 ] = b[ N - 2 ] / d[ N - 2 ];
    d[ N - 1 ] = a[ N - 1 ] - l[ N - 2 ] * u[ N - 2 ];
    /* Forward Substitution [L][y] = [q] */
    y[ 0 ] = q[ 0 ];
    for ( int i = 1; i < N; i++ )
    {
        y[ i ] = q[ i ] - l[ i - 1 ] * y[ i - 1 ];
    }
    /* Backward Substitution [U][x] = [y] */
    x[ N - 1 ] = y[ N - 1 ] / d[ N - 1 ];
    for ( int i = N - 2; i >= 0; i-- )
    {
        x[ i ] = ( y[ i ] - u[ i ] * x[ i + 1 ] ) / d[ i ];
    }
    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    return;
}

void ThomasAlgorithmBAK( int N, double * b, double * a, double * c, double * x, double * q )
{
    double *l,*u,*d,*y;
    l = new double[ N ];
    u = new double[ N ];
    d = new double[ N ];
    y = new double[ N ];
    /* LU Decomposition */
    d[ 0 ] = a[ 0 ];
    u[ 0 ] = c[ 0 ];
    for ( int i = 0; i < N - 2; i++ )
    {
        l[ i ] = b[ i ] / d[ i ];
        d[ i + 1 ] = a[ i + 1 ] - l[ i ] * u[ i ];
        u[ i + 1 ] = c[ i + 1 ];
    }
    l[ N - 2 ] = b[ N - 2 ] / d[ N - 2 ];
    d[ N - 1 ] = a[ N - 1 ] - l[ N - 2 ] * u[ N - 2 ];
    /* Forward Substitution [L][y] = [q] */
    y[ 0 ] = q[ 0 ];
    for ( int i = 1; i < N; i++ )
    {
        y[ i ] = q[ i ] - l[ i - 1 ] * y[ i - 1 ];
    }
    /* Backward Substitution [U][x] = [y] */
    x[ N - 1 ] = y[ N - 1 ] / d[ N - 1 ];
    for ( int i = N - 2; i >= 0; i-- )
    {
        x[ i ] = ( y[ i ] - u[ i ] * x[ i + 1 ] ) / d[ i ];
    }
    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    return;
}

void thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    std::vector<double> & d )
{
    size_t N = d.size();

    std::vector<double> c_star( N, 0.0 );
    std::vector<double> d_star( N, 0.0 );

    c_star[ 0 ] = c[ 0 ] / b[ 0 ];
    d_star[ 0 ] = d[ 0 ] / b[ 0 ];

    for ( int i = 1; i < N; ++ i )
    {
        double r = 1.0 / ( b[ i ] - a[ i ] * c_star[ i - 1 ] );
        d_star[ i ] = r * ( d[ i ] - a[ i ] * d_star[ i - 1 ] );
        c_star[ i ] = r * c[ i ];
    }

    d[ N - 1 ] = d_star[ N - 1 ];

    for ( int i = N - 2; i >= 0; -- i )
    {
        d[ i ] = d_star[ i ] - c_star[ i ] * d[ i + 1 ];
    }
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    std::vector<double> a{ 0, -1, -1, -1, -1 };
    std::vector<double> b{ 2,  2,  2,  2,  2 };
    std::vector<double> c{ -1, -1, -1, -1, 0 };
    std::vector<double> d{ 1.0, 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> q = d;

    std::vector<double> x( q.size() );

    Print( d, "d" );                                                                                                                                                                      
    thomas_algorithm( a, b, c, d );
    Print( d, "results" );

    int N = a.size();
    ThomasAlgorithm( N, a.data(),b.data(), c.data(), x.data(), q.data() );
    Print( x, "x" );
    Print( q, "q" );

    ThomasAlgorithmVersion1( N, a.data(), b.data(), c.data(), x.data(), q.data() );
    Print( x, "xVersion1" );

    return 0;
}
