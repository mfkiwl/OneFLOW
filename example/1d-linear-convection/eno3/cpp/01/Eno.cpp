#include "Eno.h"

//std::vector<std::vector<double>> enocoef = {
//    { 11.0/6, -7.0/6,  1.0/3 }£¬
//    { 1.0/3,  5.0/6,  -1.0/6 },
//    { -1.0/6, 5.0/6,  1.0/3 },
//    { 1.0/3,  -7.0/6,  11.0/6 }
//};

std::vector<std::vector<double>> enocoef = {
    { 3.0/2, -1.0/2 },
    { 1.0/2,  1.0/2 },
    { -1.0/2, 3.0/2 }
};

void enoL( int N, Vec1d & u, Vec1d & f )
{
    int ighost = 3;
    int iorder = 2;
    int ist = 0 - ighost;
    int ied = N - 1 + ighost;

    VecWrap dd;
    dd.Allocate( iorder, ist, ied, 0 );

    std::vector<int> ir( N + 1 );

    for ( int j = -ighost; j < N + ighost; ++j )
    {
        dd[ 0 ][ j ] = u[ j ];
    }

    for ( int i = 1; i < iorder; ++ i )
    {
        for ( int j = -ighost; j < N + ighost - 1; ++j )
        {
            dd[ i ][ j ] = dd[ i - 1 ][ j + 1 ] - dd[ i - 1 ][ j ];
        }
    }

    for ( int j = 0; j <= N; ++ j )
    {
        ir[ j ] = j;
        for ( int i = 1; i < iorder; ++ i )
        {
            if ( std::abs( dd[ i ][ ir[ j ] - 1 ] ) <= std::abs( dd[ i ][ ir[ j ] ] ) )
            {
                ir[ j ] = ir[ j ] - 1;
            }
        }
    }

    // reconstruction u(j+1_2)
    for ( int i = 0; i <= N; ++ i )
    {
        int kk = ir[ i ];
        int l = i - kk;
        int L = l + 1;
        f[ i ] = 0;
        for ( int m = 0; m < iorder; ++ m )
        {
            f[ i ] += u[ kk + m ] * enocoef[ L ][ m ];
        }
    }
}

void enoR( int N, Vec1d & u, Vec1d & f )
{
    int ighost = 3;
    int iorder = 2;
    int ist = 0 - ighost;
    int ied = N - 1 + ighost;

    VecWrap dd;
    dd.Allocate( iorder, ist, ied, 0 );

    std::vector<int> ir( N + 1 );

    for ( int j = -ighost; j < N + ighost; ++j )
    {
        dd[ 0 ][ j ] = u[ j ];
    }

    for ( int i = 1; i < iorder; ++ i )
    {
        for ( int j = -ighost; j < N + ighost - 1; ++j )
        {
            dd[ i ][ j ] = dd[ i - 1 ][ j + 1 ] - dd[ i - 1 ][ j ];
        }
    }

    for ( int j = 0; j <= N; ++ j )
    {
        ir[ j ] = j + 1;
        for ( int i = 1; i < iorder; ++ i )
        {
            if ( std::abs( dd[ i ][ ir[ j ] - 1 ] ) <= std::abs( dd[ i ][ ir[ j ] ] ) )
            {
                ir[ j ] = ir[ j ] - 1;
            }
        }
    }

    // reconstruction u(j+1_2)
    for ( int i = 0; i <= N; ++ i )
    {
        int kk = ir[ i ];
        int l = i - kk;
        int L = l + 1;
        f[ i ] = 0;
        for ( int m = 0; m < iorder; ++ m )
        {
            f[ i ] += u[ kk + m ] * enocoef[ L ][ m ];
        }
    }
}

void enoL( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        enoL( ni, u.vec( m ), f.vec( m ) );
    }
}

void enoR( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        enoR( ni, u.vec( m ), f.vec( m ) );
    }
}

