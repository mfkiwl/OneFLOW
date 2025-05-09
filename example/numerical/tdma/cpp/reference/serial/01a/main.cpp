import std;

const int size = 15;
                                                                                                                                                                                   
int main( int argc, char ** argv )
{
    //Memory allocation and generation of matrix
    double * x = new double[ size ];
    for ( int i = 0; i < size; ++ i )
    {
        x[ i ] = 0.0;
    }
    double * F = new double[ size ];
    double ** A = new double * [ size ];
    for ( int i = 0; i < size; ++ i )
    {
        A[ i ] = new double[ size ];
        for ( int j = 0; j < size; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
        F[ i ] = (double)i;
    }
    A[ 0 ][ 0 ] = -2.0;
    A[ 0 ][ 1 ] = 1.0;
    A[ size - 1 ][ size - 2 ] = 1.0;
    A[ size - 1 ][ size - 1 ] = -2.0;
    for ( int i = 1; i < size - 1; ++ i )
    {
        A[ i ][ i ] = -2.0;
        A[ i ][ i - 1 ] = 1.0;
        A[ i ][ i + 1 ] = 1.0;
    }

    //Cyclic reduction
    for ( int i = 0; i < std::log2( size + 1 ) - 1; ++ i )
    {
        std::print( "i={}\n", i );
        for ( int j = std::pow( 2, i + 1 ) - 1; j < size; j = j + std::pow( 2, i + 1 ) )
        {
            int offset = std::pow( 2, i );
            int index1 = j - offset;
            int index2 = j + offset;
            double alpha = A[ j ][ index1 ] / A[ index1 ][ index1 ];
            double gamma = A[ j ][ index2 ] / A[ index2 ][ index2 ];
            std::print( "index1,j,index2 = ({},{},{})\n", index1, j, index2 );
            for ( int k = 0; k < size; ++ k )
            {
                A[ j ][ k ] -= ( alpha * A[ index1 ][ k ] + gamma * A[ index2 ][ k ] );
                std::print( "({},{})", j, k );
            }
            std::println();
            F[ j ] -= ( alpha * F[ index1 ] + gamma * F[ index2 ] );
        }
    }

    delete [] x;
    delete [] F;
    for ( int i = 0; i < size; ++ i )
    {
        delete [] A[ i ];
    }
    delete [] A;

    return 0;
}
