// Created by CDRPico
// 20/06/2020 00:48

#pragma once

#ifndef BendersAPM_H
#define BendersAPM_H

#include"SFLP_GAPM.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"

ILOSTLBEGIN

struct entitiesBendSFLP {
    //Variables
	IloNumVarArray x;
	IloNumVarArray theta;
	IloNumVarArray2 y;
	//Objective
	IloExpr objective;
	//Constraints
	IloRangeArray Feasibility;
};

class BendersSFLP : public SFLP_GAPM
{
public:
    // Stochatic parameters
	vector<vector<double>> stoch;

	// dual multipliers
	vector<solution_sps> sp_info;

	//Stochastic parameters aggregated for a certain partition
	vector<vector<double>> stoch_agg;

	//dual multipliers agg
	vector<solution_sps> sp_info_agg;

    //Instantiate using constructor of base class
    BendersSFLP (string &inst_name, string &stoch_inst_name) : SFLP_GAPM(inst_name, stoch_inst_name) {}
    BendersSFLP() = default;

    //Create Master Problem
    void CreateMaster(const char &algo);

	//Recovering solution information
	void RecoverSolution(IloCplex &master_cpx);

	//Printing the solution
	void PrintSolution(solFeat &org);

    //compute left hand size of a cut (single or aggregated) constant
    vector<double> lhs;

    //Entities of the problem
    entitiesBendSFLP ent_sflp;

    //Master Environment
    IloEnv Mast_Bend = IloEnv();

	//Master model
	IloModel Mast_mod = IloModel(Mast_Bend);

	//Solution info
	double obj_fin;
	double LB;
	double GAP;
	size_t status;
	size_t exploredNodes;
	double cpx_runtime;
};

#endif