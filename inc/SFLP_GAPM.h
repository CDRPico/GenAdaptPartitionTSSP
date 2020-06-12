// Created by CDRPico
// 11/06/2020 13:47

#pragma once

#ifndef SFLP_GAPM_H
#define SFLP_GAPM_H

using namespace std;

#include"InstanceSFLP.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<cfloat>

ILOSTLBEGIN

//Useful Definitions
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;

struct Master_CPX {
    //Variables
    IloNumVarArray x;
    IloNumVarArray3 y;
    //Objective
    IloExpr objective;
    //Constraints
    IloRangeArray dem_satisfaction;
    IloRangeArray linking_constraints;
    IloRangeArray Feasibility; 
};

struct SubProblem_CPX {
    //Variables
    IloNumVarArray2 y;
    //Objective
    IloExpr objective;
    //Constraints
    IloRangeArray dem_constraints;
    IloRangeArray cap_constraints;
};

class SFLP_GAPM : public InstanceSFLP 
{
public:
    vector<vector<size_t>> partition;
    vector<double> obj_opt_sps;

    //COnstructor takes the name of the instance and the name of the stoch instance
    SFLP_GAPM(string &inst_name, string &stoch_inst_name);

    //Compute expectation of stochastic parameters
    vector<double> expected_demand(vector<size_t> &element);

    //Create master problem given a certain partition
    void MasterProblemCreation ();

    // Solve master problem given a certain partition
    void MasterProblemSolution (double &LB, const double &TL);

    //Create the subproblem
    void SPProblemCreation ();

    //Modify the subproblem to solve a different scenario
    void SPProblemModification(vector<size_t> &element, const bool& = false);

    //Solve the subproblem
    void SPProbleSolution(vector<double> &stoch, vector<double> &lambda, double &obj, vector<size_t> &element);

    // Given the current x_bar and solution of subprobs, the UB is computed
    double compute_UB();

private:
    //Master env
    IloEnv SFLP = IloEnv();

    //Master model
    IloModel master = IloModel(SFLP);

    //Master entities
    Master_CPX master_entities;

    //SubProblem env
    IloEnv SFLP_sp = IloEnv();

    //Subproblem model
    IloModel subprob = IloModel(SFLP_sp);

    //Subproblem entities
    SubProblem_CPX sp_entities;
};

#endif