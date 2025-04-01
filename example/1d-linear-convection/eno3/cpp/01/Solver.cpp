#include "Solver.h"
#include "Field.h"
#include "HeatField.h"
#include "ConvectionField.h"
#include "BurgersField.h"
#include "EulerField.h"
#include "CgnsUtil.h"
#include "Parallel.h"
#include "Post.h"
#include "Weno.h"
#include "ZoneState.h"
#include "global.h"
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <print> 
#include <set>
#include <unordered_map>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

Solver::Solver()
{
}

Solver::~Solver()
{
    Parallel::Finalize();
}

void Solver::Init()
{
    Parallel::Init();

    if ( Parallel::IsServer() )
    {
        std::ifstream f( "../linearconvection.json" );
        json data = json::parse( f );
        std::cout << "data=" << data.dump( 4 ) << std::endl;
        Global::istart = data[ "istart" ];
        if ( Global::istart == 1 )
        {
            Read_iter();
        }

        std::string equation = data[ "equation" ];
        if ( equation == "heat" )
        {
            Global::governing_equation = GoverningEquation::Heat;
            Global::nequ = 1;
        }
        else if ( equation == "linearconvection" )
        {
            Global::governing_equation = GoverningEquation::LinearConvection;
            Global::nequ = 1;
        }
        else if ( equation == "burgers" )
        {
            Global::governing_equation = GoverningEquation::Burgers;
            Global::nequ = 1;
        }
        else
        {
            Global::governing_equation = GoverningEquation::Euler;
            Global::nequ = 3;
        }
        Global::iconservation = data[ "iconservation" ];
        Global::iviscous = data[ "iviscous" ];
        Global::nsave = data[ "nsave" ];
        Global::idump_initial_field = data[ "idump_initial_field" ];
        Global::ifinite_volume = data[ "ifinite_volume" ];
        Global::total_time = data[ "total_time" ];
        Global::dt = data[ "dt" ];
        std::cout << "Global::total_time  = " << Global::total_time << "\n";
        std::cout << "Global::dt  = " << Global::dt << "\n";

        json &s = data[ "scheme" ];

        std::cout << "s=" << s.dump( 4 ) << std::endl;
        Global::scheme.read( s );

        if ( Global::scheme.reconstruction == to_int( BasicScheme::ENO ) ||
            Global::scheme.reconstruction == to_int( BasicScheme::WENO ) ||
            Global::scheme.reconstruction == to_int( BasicScheme::CRWENO )
            )
        {
            Global::nghost = 3;
        }
        else
        {
            Global::nghost = 1;
        }

        this->gridfile = data[ "grid" ];
    }
    HXBcastData( &Global::istart, 1, Parallel::serverid );
    HXBcastData( &Global::governing_equation, 1, Parallel::serverid );
    HXBcastData( &Global::iconservation, 1, Parallel::serverid );
    HXBcastData( &Global::iviscous, 1, Parallel::serverid );
    HXBcastData( &Global::nsave, 1, Parallel::serverid );
    HXBcastData( &Global::idump_initial_field, 1, Parallel::serverid );
    HXBcastData( &Global::ifinite_volume, 1, Parallel::serverid );
    HXBcastData( &Global::total_time, 1, Parallel::serverid );
    HXBcastData( &Global::dt, 1, Parallel::serverid );
    HXBcastData( &Global::scheme, 1, Parallel::serverid );
    HXBcastData( &Global::nghost, 1, Parallel::serverid );
    HXBcastData( &Global::nequ, 1, Parallel::serverid );
    HXBcastString( this->gridfile, Parallel::serverid );

    PrintPidHeader();
    std::cout << "Global::istart = " << static_cast<int>( Global::istart ) << "\n";
    PrintPidHeader();
    std::cout << "Global::governing_equation = " << static_cast<int>( Global::governing_equation ) << "\n";
    PrintPidHeader();
    std::cout << "Global::iconservation = " << Global::iconservation << "\n";
    PrintPidHeader();
    std::cout << "Global::iviscous = " << Global::iviscous << "\n";
    PrintPidHeader();
    std::cout << "Global::nsave = " << Global::nsave << "\n";
    PrintPidHeader();
    std::cout << "Global::idump_initial_field = " << Global::idump_initial_field << "\n";
    PrintPidHeader();
    std::cout << "Global::ifinite_volume = " << Global::ifinite_volume << "\n";
    PrintPidHeader();
    std::cout << "Global::total_time = " << Global::total_time << "\n";
    PrintPidHeader();
    std::cout << "Global::dt = " << Global::dt << "\n";
    PrintPidHeader();
    std::cout << "Global::scheme.inviscid = " << Global::scheme.inviscid << "\n";
    PrintPidHeader();
    std::cout << "Global::scheme.viscous = " << Global::scheme.viscous << "\n";
    PrintPidHeader();
    std::cout << "Global::scheme.time_scheme = " << Global::scheme.time_scheme << "\n";
    PrintPidHeader();
    std::cout << "Global::scheme.reconstruction = " << Global::scheme.reconstruction << "\n";
    PrintPidHeader();
    std::cout << "Global::nghost = " << Global::nghost << "\n";
    PrintPidHeader();
    std::cout << "this->gridfile = " << this->gridfile << "\n";
 }

void Solver::Read_iter()
{
    std::ifstream f( "iter.json" );
    json data = json::parse( f );

    Global::iter_start = data[ "iter" ];
}

void Solver::Dump_iter()
{
    //std::print( "Global::iter={}", Global::iter + 1 );
    json j;
    // add a number that is stored as double (note the implicit conversion of j to an object)
    j["iter"] = Global::iter + 1;

    // ��JSON����д���ļ�  
    std::ofstream file("iter.json"); 
    if (file.is_open()) {  
        file << j.dump( 4 ); // 4��ʾ��ʽ������  
        file.close();  
        std::cout << "JSON������д��iter.json�ļ�" << std::endl;  
    } else {  
        std::cerr << "�޷����ļ�" << std::endl;  
    }  
}

void Solver::Run()
{
    this->Init();
    this->ReadGrid();
    this->InitTopo();
    this->InitFields();
    this->SolveFields();
    this->PostProcess();
}

void Solver::ReadGrid()
{
    ReadCgnsGridBaseZone( this->gridfile );
    ReadCgnsGrid( this->gridfile );
}

void Solver::CreateField()
{
    Field * field = nullptr;
    if ( Global::governing_equation == GoverningEquation::Heat )
    {
        field = new HeatField();
    }
    else if ( Global::governing_equation == GoverningEquation::LinearConvection )
    {
        field = new ConvectionField();
    }
    else if ( Global::governing_equation == GoverningEquation::Burgers )
    {
        field = new BurgersField();
    }
    else
    {
        field = new EulerField();
    }
    Global::fields.push_back( field );
}

void Solver::InitFields()
{
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Solver::InitFields() ZoneState::nZones = " << ZoneState::nZones << "\n";
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        CreateField();
    }

    InitFieldCommon();

    if ( Global::istart == 0 )
    {
        InitFieldAsRestart();
    }
    else
    {
        ReadFlowField();
    }
    
    this->Boundary();
    this->UpdateOldField();
}

void Solver::InitFieldAsRestart()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = Global::fields[ iZone ];
        field->InitFieldAsRestart( grid );
    }
}

void Solver::ReadFlowField()
{
    std::string filename = "field_final.csv";
    std::fstream file;
    if ( Parallel::pid == Parallel::serverid )
    {
        file.open( filename.c_str(), std::fstream::in );
    }
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int send_pid = ZoneState::pids[ iZone ];
        int recv_pid = Parallel::serverid;

        if ( Parallel::pid == send_pid )
        {
            int local_zoneid = ZoneState::g2lzoneids[ iZone ];

            Grid * grid = Global::grids[ local_zoneid ];
            Field * field = Global::fields[ local_zoneid ];
            field->ReadFlowField( file, grid );
        }

        //HXSendRecvString( Global::file_string, send_pid, Parallel::serverid );

        //if ( Parallel::pid == Parallel::serverid )
        //{
        //    total_string += Global::file_string;
        //}
    }
    if ( Parallel::pid == Parallel::serverid )
    {
        file.close();
    }
}

void Solver::InitFieldCommon()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        grid->zoneIndex = iZone;
        Field * field = Global::fields[ iZone ];
        field->InitFieldCommon( grid );
    }
}

void Solver::DumpInitialFields()
{
    this->DumpField();
}

void Solver::TimeIntegral()
{
    BasicScheme time_scheme = to_BasicScheme( Global::scheme.time_scheme );
    switch ( time_scheme ) {
    case BasicScheme::RK1:
        this->RungeKutta( 1 );
        break;
    case BasicScheme::RK2:
        this->RungeKutta( 2 );
        break;
    case BasicScheme::RK3:
        this->RungeKutta( 3 );
        break;
    default:
        this->CrankNicolsonSeries();
    }
}

void Solver::SolveFields()
{
    if ( Global::idump_initial_field == 1 )
    {
        this->DumpInitialFields();
    }

    for ( int it = Global::iter_start; it < Global::nt; ++ it )
    {
        Global::iter = it;
        this->TimeIntegral();

        if ( ( Global::iter + 1 ) % Global::nsave == 0 )
        {
            std::print( "it = {} nt = {}\n", Global::iter + 1, Global::nt );
            this->DumpField();
        }
    }
}

void Solver::FTCS()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->FTCS( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::CrankNicolsonSeries()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->CrankNicolsonSeries( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::ICP()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->ICP( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::RungeKutta( int nStage )
{
    for ( int istage = 0; istage < nStage; ++ istage )
    {
        this->RungeKutta( nStage, istage );
    }
    this->UpdateOldField();
}

void Solver::RungeKutta( int nStage, int istage )
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->RungeKutta( zone, nStage, istage );
    }
    this->Boundary();
}

void Solver::Boundary()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->PhysicalBoundary( zone );
    }
    ExchangeInterfaceField();
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->InterfaceBoundary( zone );
    }
}

void Solver::UpdateOldField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->UpdateOldField();
    }
}

void Solver::UploadInterfaceField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        Field * field = Global::fields[ iZone ];

        int nsend_zones = interface->send_to_zones.size();
        for ( int iSend = 0; iSend < nsend_zones; ++ iSend )
        {
            int zone_to_send = interface->send_to_zones[ iSend ];
            std::vector<int> & donorfaces_for_send = interface->donorfaces_for_send[ iSend ];
            std::vector<int> & donorijk_for_send = interface->donorijk_for_send[ iSend ];
            std::vector<double> & donordata_for_send = interface->donordata_for_send[ iSend ];

            int nInterFaces = donorfaces_for_send.size();
            int index_dim = 1;
            int ngsize = Global::nghost + 1 - Global::ifinite_volume;
            for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
            {
                int ijkpos = index_dim * iFace * ngsize;
                int data_pos = iFace * ngsize * Global::nequ;
  
                for ( int ig = Global::ifinite_volume; ig <= Global::nghost; ++ ig )
                {
                    int iig = ig - Global::ifinite_volume;
                    int id_cell = donorijk_for_send[ ijkpos + iig ] - 1;
                    for ( int iequ = 0; iequ < Global::nequ; ++ iequ )
                    {
                        Vec1d & u = field->u.vec( iequ );
                        donordata_for_send[ data_pos ++ ] = u[ id_cell ];
                    }
                }
            }
        }
    }
}

void Solver::UpdateInterfaceField()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int ndonor_zones = Global::interfaceTopo.linkmap[ iZone ].size();
        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donor_zone = Global::interfaceTopo.linkmap[ iZone ][ iNei ];
            int send_pid = ZoneState::pids[ iZone ];
            int recv_pid = ZoneState::pids[ donor_zone ];
            int nsend = -1;
            std::vector<double> donordata;
            if ( Parallel::pid != send_pid && Parallel::pid != recv_pid ) continue;
            if ( Parallel::pid == send_pid )
            {
                int local_zoneid = ZoneState::g2lzoneids[ iZone ];
                Interface * interface = Global::interfaces[ local_zoneid ];
                donordata = interface->donordata_for_send[ iNei ];
                nsend = donordata.size();
            }
            HXSendRecvData( &nsend, 1, send_pid, recv_pid );

            if ( Parallel::pid == recv_pid )
            {
                donordata.resize( nsend );
            }
            HXSendRecvData( donordata.data(), donordata.size(), send_pid, recv_pid );

            if ( Parallel::pid == recv_pid )
            {
                int local_donor_zoneid = ZoneState::g2lzoneids[ donor_zone ];
                Interface * donor_interface = Global::interfaces[ local_donor_zoneid ];
                int nSize = donor_interface->neighbor_donor_zones.size();
                int ipos = -1;
                for ( int i = 0; i < nSize; ++ i )
                {
                    int nei_zone = donor_interface->neighbor_donor_zones[ i ];
                    if ( nei_zone == iZone )
                    {
                        ipos = i;
                        break;
                    }
                }

                std::vector<int> & neighbor_donorfaces = donor_interface->neighbor_donorfaces[ ipos ];
                std::vector<int> & sub_local_faceids = donor_interface->sub_local_faceids[ ipos ];
                for ( int i = 0; i < neighbor_donorfaces.size(); ++ i )
                {
                    int local_faceid = sub_local_faceids[ i ];

                    int index_dim = 1;
                    int ngsize =  Global::nghost + 1 - Global::ifinite_volume;
                    int ijkpos = index_dim * local_faceid * ngsize;
                    int data_pos = local_faceid * ngsize * Global::nequ;
                    int donor_data_pos = i * ngsize * Global::nequ;

                    for ( int ig = Global::ifinite_volume; ig <= Global::nghost; ++ ig )
                    {
                        for ( int iequ = 0; iequ < Global::nequ; ++ iequ )
                        {
                            double donor_value = donordata[ donor_data_pos ++ ];
                            donor_interface->data_recv[ data_pos ++ ] = donor_value;
                        }
                    }
                }
            }
        }
    }
}

void Solver::DownloadInterfaceField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];
        Field * field = Global::fields[ iZone ];

        int nInterFaces = interface->zoneList.size();

        int index_dim = 1;
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int ngsize =  Global::nghost + 1 - Global::ifinite_volume;
            int ijkpos = index_dim * iFace * ngsize;
            int data_pos = iFace * ngsize * Global::nequ;
            for ( int ig = Global::ifinite_volume; ig <= Global::nghost; ++ ig )
            {
                int ig_cell = interface->ijk_ghosts[ ijkpos ++ ] - 1;
                if ( ig == 0 )
                {
                    for ( int iequ = 0; iequ < Global::nequ; ++ iequ )
                    {
                        Vec1d & u = field->u.vec( iequ );
                        double donor_value = interface->data_recv[ data_pos ++ ];
                        double valueb = 0.5 * ( u[ ig_cell ] + donor_value );
                        u[ ig_cell ] = valueb;
                    }
                }
                else
                {
                    for ( int iequ = 0; iequ < Global::nequ; ++ iequ )
                    {
                        Vec1d & u = field->u.vec( iequ );
                        double donor_value = interface->data_recv[ data_pos ++ ];
                        u[ ig_cell ] = donor_value;
                    }
                }
            }
        }
    }
}

void Solver::ExchangeInterfaceField()
{
    this->UploadInterfaceField();
    this->UpdateInterfaceField();
    this->DownloadInterfaceField();
}

void Solver::DumpField()
{
    std::string filename = std::format( "field_final{}.csv", Global::iter+1 );
    this->DumpField( filename );
}

void Solver::DumpField( const std::string & filename )
{
    std::string total_string = {};
    std::fstream file;
    if ( Parallel::pid == Parallel::serverid )
    {
        file.open( filename.c_str(), std::fstream::out );
    }
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int send_pid = ZoneState::pids[ iZone ];
        int recv_pid = Parallel::serverid;

        Global::file_string = {};

        if ( Parallel::pid == send_pid )
        {
            int local_zoneid = ZoneState::g2lzoneids[ iZone ];

            Grid * grid = Global::grids[ local_zoneid ];
            Field * field = Global::fields[ local_zoneid ];
            field->DumpField( grid );
        }

        HXSendRecvString( Global::file_string, send_pid, Parallel::serverid );

        if ( Parallel::pid == Parallel::serverid )
        {
            total_string += Global::file_string;
        }
    }
    if ( Parallel::pid == Parallel::serverid )
    {
        std::format_to( std::ostream_iterator<char>( file ), "{}", total_string );
        file.close();
    }
}

void Solver::PostProcess()
{
    std::string filename = "field_final.csv";
    this->DumpField( filename );
    this->Dump_iter();
}

void Solver::PrintField( std::vector<double> &f )
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

void Solver::InitTopo()
{
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Solver::InitTopo() " << "\n";

    Global::donor_zone_sets.resize( LocalZone::nZones );
    Global::donor_zones.resize( LocalZone::nZones );
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        int global_zoneid = LocalZone::global_zoneids[ iZone ];
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "iZone = " << iZone << " global_zoneid = " << global_zoneid << "\n";

        Interface * interface = new Interface();
        interface->zoneid = global_zoneid;
        Global::interfaces.push_back( interface );
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Zone * zone = Global::zones[ iZone ];
        Interface * interface = Global::interfaces[ iZone ];

        int nbc1to1s = zone->bc1to1s.size();
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout <<  "nbc1to1s = " << nbc1to1s << "\n";

        for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
        {
            ZoneBc1To1 * bc1to1 = zone->bc1to1s[ ibc1to1 ];
            int zoneid = bc1to1->zoneid;
            int donor_zoneid = bc1to1->donor_zoneid;
            std::cout << "Parallel::pid = " << Parallel::pid << " ";
            std::cout << "zoneid = " << zoneid << " donor_zoneid = " << donor_zoneid << "\n";

            Region region;
            region.SetRegion( bc1to1->pnts );

            Region donor_region;
            donor_region.SetRegion( bc1to1->donor_pnts );

            Transform transform;
            transform.begin1 = region.start;
            transform.begin2 = donor_region.start;
            transform.transform = bc1to1->transform;
            transform.Init();

            interface->CalcInterface( &transform, region.start, region.end, donor_zoneid );
        }
        int nInterfaces = interface->zoneList.size();
        int ngsize = Global::nghost + 1 - Global::ifinite_volume;
        int nData = nInterfaces * ngsize * Global::nequ;
        interface->data_recv.resize( nData );
        interface->data_send.resize( nData );
    }

    for ( int iProc = 0; iProc < Parallel::nProc; ++ iProc )
    {
        int nSize = -1;
        if ( iProc == Parallel::pid )
        {
            nSize = Global::facePairList.size();
        }
        HXBcastData( &nSize, 1, iProc );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "nSize = " << nSize << "\n";
        std::vector<FacePair> tmp;
        if ( iProc == Parallel::pid )
        {
            tmp = Global::facePairList;
        }
        else
        {
            tmp.resize( nSize );
        }

        HXBcastData( tmp.data(), tmp.size(), iProc );
        Global::AddFacePairList( Global::mpi_facePairList, tmp );
    }

    for ( int i = 0; i < Global::mpi_facePairList.size(); ++ i )
    {
        FacePair &facePair = Global::mpi_facePairList[ i ];
        Global::InsertFacePairMap( facePair );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        facePair.Print();
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Zone * zone = Global::zones[ iZone ];

        Interface * interface = Global::interfaces[ iZone ];
        int nInterfaces = interface->local_faceids.size();
        for ( int iInterface = 0; iInterface < nInterfaces; ++ iInterface )
        {
            int local_faceid = interface->local_faceids[ iInterface ];
            int proc_global_faceid = interface->proc_global_faceids[ iInterface ];
            FacePair & facePair = Global::facePairList[ proc_global_faceid ];
            int global_faceid = Global::InsertFacePairMap( facePair );
            interface->global_faceids.push_back( global_faceid );
            interface->global_local_face_map.insert( std::make_pair( global_faceid, local_faceid ) );
        }
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];
        int nInterFaces = interface->zoneList.size();
        std::set<int> &donor_zoneSet = Global::donor_zone_sets[ iZone ];
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            donor_zoneSet.insert( donor_zoneid );
        }

        std::vector<int> &donor_zones = Global::donor_zones[ iZone ];
        for ( std::set<int>::iterator iter = donor_zoneSet.begin(); iter != donor_zoneSet.end(); ++ iter )
        {
            donor_zones.push_back( *iter );
        }

        interface->neighbor_donor_zones = donor_zones;

        std::unordered_map<int, int> donor_zonelocal;

        for ( int idonor = 0; idonor < donor_zones.size(); ++ idonor )
        {
            int donor_zone = donor_zones[ idonor ];
            donor_zonelocal.insert( std::make_pair( donor_zone, idonor ) );
        }
        int ndonors = donor_zones.size();
        std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;
        neighbor_donorfaces.resize( ndonors );

        std::vector<std::vector<int>> & sub_local_faceids = interface->sub_local_faceids;
        sub_local_faceids.resize( ndonors );

        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            int ineighbor = donor_zonelocal[ donor_zoneid ];
            std::vector<int> &donorfaces = neighbor_donorfaces[ ineighbor ];
            int global_faceid = interface->global_faceids[ iFace ];
            donorfaces.push_back( global_faceid );
            int local_faceid = interface->local_faceids[ iFace ];

            std::vector<int> & sub_local_faces = sub_local_faceids[ ineighbor ];
            sub_local_faces.push_back( local_faceid );
        }
    }

    Global::interfaceTopo.InitNeighborInfo();
    Global::interfaceTopo.SwapNeighborInfo();
}