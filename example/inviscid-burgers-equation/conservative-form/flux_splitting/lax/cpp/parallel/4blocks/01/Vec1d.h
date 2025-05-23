#pragma once
#include <vector>

class Vec1d
{
public:
    std::vector<double> data;
    int ist = 0;
    int ied;
public:
    void Allocate( int ist, int ied, double value = 0 )
    {
        int nelement = ied - ist + 1;
        this->data.resize( nelement, value );
        this->ist = ist;
        this->ied = ied;
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

    Vec1d & operator = ( const Vec1d & rhs )
    {
        if ( this == & rhs ) return * this;
        this->data = rhs.data;
        this->ist = rhs.ist;
        this->ied = rhs.ied;
        return * this;
    }

    Vec1d & operator = ( const double value )
    {
        for ( int i = 0; i < data.size(); ++ i )
        {
            data[ i ] = value;
        }
        return * this;
    }
};
