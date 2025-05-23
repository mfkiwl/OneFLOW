#include "cgnslib.h"
#include <iostream>
#include <vector>
#include <numbers>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <string> 

class Grid;
void ReadCgnsGrid( const std::string & filename );

class Grid
{
public:
    int zoneIndex;
    std::vector<double> x;
};

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid )
    {
        this->ni = grid->x.size();
        std::cout << "ni = " << ni << "\n";

        std::vector<double> & x = grid->x;
        this->dx = x[ 1 ] - x[ 0 ];
        this->dt = dx / 10.0;
        this->t = 1.0;
        int nt = static_cast<int>( t / dt ) + 1;
        std::cout << "nt = " << nt << "\n";

        this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
        this->beta = this->alpha * dt / ( dx * dx );
        std::cout << "alpha = " << std::setprecision(15) << this->alpha << "\n";
        std::cout << "beta = " << std::setprecision(15) << this->beta << "\n";

        u_e.resize( ni );
        u.resize( ni );
        un.resize( ni );

        for ( int i = 0; i < ni; ++ i )
        {
            u_e[ i ] = - std::exp( -t ) * std::sin( std::numbers::pi * x[i] ); //theory solution
            un[ i ] = - std::sin( std::numbers::pi * x[ i ] ); //initial condition @ t=0
        }
        un[ 0 ] = 0.0;
        un[ ni - 1 ] = 0.0;
    }

    void Solve( Grid * grid )
    {
        for ( int it = 0; it < nt; ++ it )
        {
            for ( int i = 1; i < ni - 1; ++ i )
            {
                u[ i ] = un[ i ] + beta * ( un[ i + 1 ] - 2.0 * un[ i ] + un[ i - 1 ] );
            }
            //boundary
            u[ 0 ] = 0.0; // boundary condition at x = -1
            u[ ni - 1 ] = 0.0; // boundary condition at x = 1
            this->update( un, u );
        }
    }

    void PostProcess( Grid * grid )
    {
        std::vector<double> & x = grid->x;
        //compute L2 norm of the error
        std::vector<double> u_error( ni );
        for ( int i = 0; i < ni; ++ i )
        {
            u_error[ i ] = un[ i ] - u_e[ i ];
        }

        this->DumpErrorDetails( u_error );

        std::string csvname = "field_final" + std::to_string( grid->zoneIndex ) + ".csv";
        this->DumpCsvFile( csvname, x, u_e, un, u_error );
    }

    void update( std::vector<double> &un, std::vector<double> &u )
    {
        for ( int i = 0; i < u.size(); ++ i )
        {
            un[ i ] = u[ i ];
        }
    }

    void DumpErrorDetails( std::vector<double> &u_error )
    {
        int ni = u_error.size();
        double rms_error = compute_l2norm( ni, u_error );
        double max_error = compute_max_error( ni, u_error );
        std::cout << "max_error = " << std::setprecision(15) << max_error << "\n";
        //create output file for L2-norm
        std::fstream file;
        file.open("output.txt", std::fstream::out);
        std::format_to(std::ostream_iterator<char>(file), "Error details: \n");
        std::format_to(std::ostream_iterator<char>(file), "L-2 Norm = {0}\n", rms_error);
        std::format_to(std::ostream_iterator<char>(file), "Maximum Norm = {0}\n", max_error);
        file.close();
    }

    void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &uerror )
    {
        std::fstream file;
        file.open(filename.c_str(), std::fstream::out);
        std::format_to(std::ostream_iterator<char>(file), "x ue un uerror\n");
        for ( int i = 0; i < x.size(); ++ i )
        {
            std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
        }
        file.close();
    }

    double compute_l2norm( int ni, std::vector<double> & r )
    {
        double rms = 0.0;
        for ( int i = 1; i < ni - 1; ++ i )
        {
            rms += r[ i ] * r[ i ];
        }
        rms = std::sqrt( rms / ( ni - 2 ) );
        return rms;
    }

    double compute_max_error( int ni, std::vector<double> & u_error )
    {
        double val_max = -1;
        int ipos = -1;
        for ( int i = 1; i < ni - 1; ++ i )
        {
            //val_max = std::max( val_max, std::abs( u_error[ i ] ) );
            if ( val_max < std::abs( u_error[ i ] ) )
            {
                ipos = i;
                val_max = std::abs( u_error[ i ] );
            }
        }
        std::cout << " ipos = " << ipos << "\n";
        return val_max;
    }
};

class Global
{
public:
    static std::vector<Grid *> grids;
    static std::vector<Field *> fields;
};

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;

class Solver
{
public:
    std::vector<double> u_e;
    std::vector<double> u,un;
    std::vector<double> error;
public:
    void Run()
    {
        this->ReadGrid();
        this->InitFields();
        this->SolveMultiZones();
        this->PostProcess();
    }

    void ReadGrid()
    {
        ReadCgnsGrid( "../heat1d2blocks.cgns" );
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            grid->zoneIndex = iZone;
        }
    }

    void InitFields()
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = new Field();
            Global::fields.push_back( field );
            field->Init( grid );
        }
    }

    void SolveMultiZones()
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = Global::fields[ iZone ];
            field->Solve( grid );
        }
    }

    void PostProcess()
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = Global::fields[ iZone ];
            field->PostProcess( grid );
        }
    }

    void PrintField( std::vector<double> &f )
    {
        int icount = 0;
        for ( int i = 0; i < f.size(); ++ i )
        {
            std::cout << std::setprecision(15) << f[ i ] << " ";
            icount ++;
            if ( icount % 5 == 0 )
            {
                std::cout << "\n";
            }
        }
        std::cout << "\n";
        std::cout << "\n";
    }
};

void ReadCgnsGrid( const std::string & filename )
{
    int fileId = -1;
    cg_open( filename.c_str(), CG_MODE_READ, &fileId);
    std::cout << "fileId = " << fileId << "\n";
    int nbases = -1;
    cg_nbases( fileId, &nbases );
    std::cout << "nbases = " << nbases << "\n";
    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";
        int nzones = -1;
        cg_nzones(fileId, baseId, &nzones);
        std::cout << "nzones = " << nzones << "\n";
        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            cg_index_dim( fileId, baseId, zoneId, &index_dim );
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize(index_dim*3);

            char zonename[33];
            cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "i = " << i << " isize = " << isize[i] << "\n";
            }

            std::vector<cgsize_t> irmin(index_dim);
            std::vector<cgsize_t> irmax(index_dim);
            int nNodes = 1;
            for ( int m = 0; m < index_dim; ++ m )
            {
                /* lower range index */
                irmin[ m ] = 1;
                /* upper range index of vertices */
                irmax[ m ] = isize[ m ];
                nNodes *= irmax[ m ];
            }
            std::cout << "nNodes = " << nNodes << "\n";

            ZoneType_t zoneType;
            cg_zone_type( fileId, baseId, zoneId, &zoneType );
            std::cout << "zoneType = " << zoneType << " ZoneTypeName = " << ZoneTypeName[zoneType] << "\n";
            int ncoords = -1;
            cg_ncoords( fileId, baseId, zoneId, &ncoords );
            std::cout << "ncoords = " << ncoords << "\n";

            Grid * grid = new Grid();
            Global::grids.push_back( grid );

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = icoord + 1;
                DataType_t dataType;
                char coordname[33];
                cg_coord_info( fileId, baseId, zoneId, coorId, &dataType, coordname );
                std::cout << "coordname = " << coordname << "\n";
                std::cout << "dataType = " << dataType << " DataTypeName = " <<  DataTypeName[dataType] << "\n";
                //std::vector<double> coord( nNodes );
                std::vector<char> coord( nNodes * sizeof(double) );
                grid->x.resize( nNodes );

                cg_coord_read( fileId, baseId, zoneId, coordname, dataType, irmin.data(), irmax.data(), coord.data() );
                double * xd = reinterpret_cast<double*>(const_cast<char *>(coord.data()));
                for ( int i = 0; i < nNodes; ++ i )
                {
                    //std::cout << coord[i] << " ";
                    std::cout << xd[i] << " ";
                    grid->x[ i ] = xd[ i ];
                    if( (i+1)%5 == 0 ) std::cout << "\n";
                }
                std::cout << "\n";
            }

            int nbocos = -1;

            cg_nbocos( fileId, baseId, zoneId, &nbocos );
            std::cout << "nbocos = " << nbocos << "\n";
            for ( int iboco = 0; iboco < nbocos; ++ iboco )
            {
                int bccoId = iboco + 1;
                GridLocation_t location;
                cg_boco_gridlocation_read( fileId, baseId, zoneId, bccoId, &location );
                std::cout << "iboco = " << iboco <<  " location = " << location << " GridLocationName = " << GridLocationName[location] << "\n";

                char boconame[ 33 ];
                BCType_t bocotype;
                PointSetType_t ptset_type;
                cgsize_t npnts = 0;
                std::vector<int> normalIndex(index_dim,-1);
                cgsize_t normalListSize = 0;
                DataType_t normalDataType;
                int ndataset = -1;

                cg_boco_info( fileId, baseId, zoneId, bccoId, boconame, &bocotype, &ptset_type, 
                    &npnts, normalIndex.data(), &normalListSize, &normalDataType, &ndataset);
                std::cout << "boconame = " << boconame <<  " bocotype = " << bocotype << " BCTypeName = " << BCTypeName[bocotype] << "\n";
                std::cout << "ptset_type = " << ptset_type <<  " PointSetTypeName = " << PointSetTypeName[ptset_type] << "\n";
                std::cout << "npnts = " << npnts << "\n";
                std::cout << "normalIndex = ";
                for ( int i = 0; i < index_dim; ++ i )
                {
                    std::cout << normalIndex[i] << " ";
                }
                std::cout << "\n";
                std::cout << "normalListSize = " << normalListSize << "\n";
                std::cout << "normalDataType = " << normalDataType << " DataTypeName = " << DataTypeName[normalDataType] << "\n";
                std::cout << "ndataset = " << ndataset << "\n";

                std::vector<char> normalList(nNodes*iphysdim*sizeof(double));

                std::vector<cgsize_t> pnts(npnts*index_dim);
                cg_boco_read( fileId, baseId, zoneId, bccoId, pnts.data(), normalList.data() );
                std::cout << "pnts = ";
                for ( int i = 0; i < pnts.size(); ++ i )
                {
                    std::cout << pnts[ i ] << " ";
                }
                std::cout << "\n";

                double * normal_d = reinterpret_cast<double*>(const_cast<char *>(normalList.data()));

                //std::cout << "normalList = ";
                //for ( int i = 0; i < nNodes*iphysdim; ++ i )
                //{
                //    std::cout << normal_d[ i ] << " ";
                //}
                //std::cout << "\n";
            }
        }
    }
    cg_close(fileId);
}

int main(int argc, char **argv)
{
    Solver solver;
    solver.Run();
    //solver.solve();
    return 0;
}
