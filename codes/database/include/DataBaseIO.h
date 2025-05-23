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
#include "HXDefine.h"
#include "HXArray.h"
#include "DataBook.h"
#include <string>
#include <fstream>

BeginNameSpace( ONEFLOW )

class DataBook;

template < typename TIO, typename T >
void HXRead( TIO * tio, T & value );

template < typename T >
void HXRead( std::fstream * file, T & value );

template < typename TIO, typename T >
void HXRead( TIO * tio, T * field, int nElement );

template < typename T >
void HXRead( std::fstream * file, T * field, int nElement );

template < typename TIO, typename T >
void HXReadVector( TIO * tio, std::vector< T > & field );

template < typename TIO, typename T >
void HXRead( TIO * tio, std::vector< T > & field );

template < typename TIO, typename T >
void HXRead( TIO * tio, HXVector< T > & field );

template < typename T >
void HXReadVector( std::fstream * file, std::vector< T > & field );

template < typename T >
void HXRead( DataBook * dataBook, std::vector< std::vector< T > > & field2D );

template < typename T >
void HXRead( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

template < typename vect2D >
void HXReadVector2D( DataBook * dataBook, vect2D & field2D );

template < typename T >
void HXRead( std::fstream * file, std::vector< T > & field );

template < typename T >
void HXRead( std::fstream * file, HXVector< T > & field );

template < typename TIO, typename T >
void HXWrite( TIO * tio, T & value );

template < typename T >
void HXWrite( std::fstream * file, T & value );

template < typename TIO, typename T >
void HXWrite( TIO * tio, T * field, int nElement );

template < typename T >
void HXWrite( std::fstream * file, T * field, int nElement );

template < typename TIO, typename T >
void HXWriteVector( TIO * tio, std::vector< T > & field );

template < typename TIO, typename T >
void HXWrite( TIO * tio, std::vector< T > & field );

template < typename TIO, typename T >
void HXWrite( TIO * tio, HXVector< T > & field );

template < typename T >
void HXWriteVector( std::fstream * file, std::vector< T > & field );

template < typename T >
void HXWrite( std::fstream * file, std::vector< T > & field );

template < typename T >
void HXWrite( std::fstream * file, HXVector< T > & field );

template < typename vect2D >
void HXWriteVector2D( DataBook * dataBook, vect2D & field2D );

template < typename T >
void HXWrite( DataBook * dataBook, std::vector< std::vector< T > > & field2D );

template < typename T >
void HXWrite( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

template < typename T >
void HXAppend( DataBook * dataBook, T & value );

template < typename T >
void HXAppend( DataBook * dataBook, T * field, int nElement );

template < typename T >
void HXAppendVector( DataBook * dataBook, std::vector< T > & field );

template < typename T >
void HXAppend( DataBook * dataBook, std::vector< T > & field );

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< T > & field );

template < typename vect2D >
void HXAppendVector2D( DataBook * dataBook, vect2D & field2D );

template < typename T >
void HXAppend( DataBook * dataBook, std::vector< std::vector< T > > & field2D );

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

void HXRead( DataBook * dataBook, std::string & cs );
void HXWrite( DataBook * dataBook, std::string & cs );

void HXRead( DataBook * dataBook, MRField * field );
void HXWrite( DataBook * dataBook, MRField * field );

template < typename TIO, typename T >
void HXRead( TIO * tio, T & value )
{
    tio->Read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename T >
void HXRead( std::fstream * file, T & value )
{
    file->read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename TIO, typename T >
void HXRead( TIO * tio, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Read( field, nElement * sizeof( T ) );
}

template < typename T >
void HXRead( std::fstream * file, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->read( reinterpret_cast< char * >( field ), nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXReadVector( TIO * tio, std::vector< T > & field )
{
    size_t nElement = field.size();
    if ( nElement <= 0 ) return;
    tio->Read( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXRead( TIO * tio, std::vector< T > & field )
{
    HXReadVector( tio, field );
}

template < typename TIO, typename T >
void HXRead( TIO * tio, HXVector< T > & field )
{
    HXReadVector( tio, field );
}

template < typename T >
void HXReadVector( std::fstream * file, std::vector< T > & field )
{
    int nElement = field.size();
    HXRead( file, & field[ 0 ], nElement );
}

template < typename T >
void HXRead( std::fstream * file, std::vector< T > & field )
{
    HXReadVector( file, field );
}

template < typename T >
void HXRead( std::fstream * file, HXVector< T > & field )
{
    HXReadVector( file, field );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, T & value )
{
    tio->Write( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename T >
void HXWrite( std::fstream * file, T & value )
{
    file->write( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Write( field, nElement * sizeof( T ) );
}

template < typename T >
void HXWrite( std::fstream * file, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->write( reinterpret_cast< char * >( field ), nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXWriteVector( TIO * tio, std::vector< T > & field )
{
    int nElement = field.size();
    if ( nElement <= 0 ) return;
    tio->Write( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, std::vector< T > & field )
{
    HXWriteVector( tio, field );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, HXVector< T > & field )
{
    HXWriteVector( tio, field );
}

template < typename T >
void HXWriteVector( std::fstream * file, std::vector< T > & field )
{
    HXSize_t nElement = static_cast<int> (field.size());
    HXWrite( file, & field[ 0 ], nElement );
}

template < typename T >
void HXWrite( std::fstream * file, std::vector< T > & field )
{
    HXWriteVector( file, field );
}

template < typename T >
void HXWrite( std::fstream * file, HXVector< T > & field )
{
    HXWriteVector( file, field );
}

template < typename vect2D >
void HXReadVector2D( DataBook * dataBook, vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        auto & field = field2D[ iElem ];

        int nSubElem = 0;
        HXRead( dataBook, nSubElem );

        field.resize( nSubElem );
        HXRead( dataBook, field );
    }
}

template < typename T >
void HXRead( DataBook * dataBook, std::vector< std::vector< T > > & field2D )
{
    HXReadVector2D( dataBook, field2D );
}

template < typename T >
void HXRead( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    HXReadVector2D( dataBook, field2D );
}

template < typename vect2D >
void HXWriteVector2D( DataBook * dataBook, vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        auto & field = field2D[ iElem ];

        int nSubElem = field.size();
        HXWrite( dataBook, nSubElem );

        HXWrite( dataBook, field );
    }
}

template < typename T >
void HXWrite( DataBook * dataBook, std::vector< std::vector< T > > & field2D )
{
    HXWriteVector2D( dataBook, field2D );
}

template < typename T >
void HXWrite( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    HXWriteVector2D( dataBook, field2D );
}

template < typename T >
void HXAppend( DataBook * dataBook, T & value )
{
    dataBook->Append( & value, sizeof( T ) );
}

template < typename T >
void HXAppend( DataBook * dataBook, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    dataBook->Append( field, nElement * sizeof( T ) );
}

template < typename T >
void HXAppendVector( DataBook * dataBook, std::vector< T > & field )
{
    HXSize_t nElement = field.size();
    if ( nElement <= 0 ) return;
    dataBook->Append( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename T >
void HXAppend( DataBook * dataBook, std::vector< T > & field )
{
    HXAppendVector( dataBook, field );
}

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< T > & field )
{
    HXAppendVector( dataBook, field );
}

template < typename vect2D >
void HXAppendVector2D( DataBook * dataBook, vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        auto & field = field2D[ iElem ];

        HXSize_t nSubElem = field.size();
        HXAppend( dataBook, nSubElem );

        HXAppend( dataBook, field );
    }
}

template < typename T >
void HXAppend( DataBook * dataBook, std::vector< std::vector< T > > & field2D )
{
    HXAppendVector2D( dataBook, field2D );
}

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    HXAppendVector2D( dataBook, field2D );
}

EndNameSpace
