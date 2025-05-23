#pragma once
#include <vector>
#include "global.h"
#include "CgnsGrid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

enum class Scheme
{
    FTCS = 0,
    CN, //Crank�CNicolson
    RungeKutta
};

class Solver
{
public:
    Solver();
    ~Solver();
public:
    Scheme scheme;
public:
    void Run();
    void ReadGrid();
    void InitFields();
    void InitTopo();
    void SolveFields();
public:
    void Boundary();
    void UpdateOldField();
    void UpdateRungeKuttaOldField( int istage );
    void UploadInterfaceField();
    void UpdateInterfaceField();
    void DownloadInterfaceField();
    void ExchangeInterfaceField();
    void PrintField( std::vector<double> & f );
public:
    void PostProcess();
public:
    void FTCS();
    void CN();
    void RungeKutta();
    void RungeKutta( int istage );
};

class Post
{
public:
    std::vector<int> zoneids;
    std::vector<double> x;
    std::vector<double> u_e;
    std::vector<double> un;
    std::vector<double> u;
public:
    void Process();
    void ReorderZones();
    void GatherField();
public:
    void AddVector( std::vector<double> & a, std::vector<double> & b );
    void DumpField();
    void DumpErrorDetails( std::vector<double> & u_error );
    void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );
};