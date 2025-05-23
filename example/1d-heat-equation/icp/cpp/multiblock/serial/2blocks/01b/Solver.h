#pragma once
#include <vector>
#include "Global.h"
#include "Grid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

class Solver
{
public:
    Solver();
    ~Solver();
public:
    //Scheme scheme;
    int nghost;
public:
    void Init();
    void Run();
    void ReadGrid();
    void InitFields();
    void InitTopo();
    void SolveFields();
    void TimeIntegral();
    void DumpInitialFields();
public:
    void Boundary();
    void UpdateOldField();
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
    void ICP();
public:
    void RungeKutta( int nStage );
    void RungeKutta( int nStage, int istage );
};
