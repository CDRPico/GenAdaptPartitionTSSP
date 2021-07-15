// Created by CDRPico
// 14/10/2020 17:35

#pragma once

#ifndef BendersAPM_CP_H
#define BendersAPM_CP_H

#include"CapPlan_GAPM.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include"UsfFunctions.h"

ILOSTLBEGIN

class BendersCP : public CP_GAPM 
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

	//compute left hand size of a cut (single or aggregated) constant
	vector<double> lhs;

	//Instantiate using constructor of base class
	BendersCP(string &file_ce_core, string &file_ce_stoch, string &file_ce_time, string &testname) : 
		CP_GAPM(file_ce_core, file_ce_stoch, file_ce_time, testname) {};

	//Run Benders
	double runBenders(const char &algo, vector<vector<size_t>> &part, vector<double> &xb);

	// Mixed Benders APM + Single cut
	double mixedBenders(const char &algo, vector<vector<size_t>> &part, vector<double> &xb);
	
	//Create Master Problem
	void CreateMaster(const char &algo, vector<vector<size_t>> &part, vector<double> &xb);

	//Update Master problem
	void updateMaster(solFeat &org, size_t &nc, bool &violated);
	void disaggPartition(solFeat &org, bool &violated);

	//Aggregate cuts right after disaggregation
	void cutsindiasgg(solFeat &org);

	//Solve Master
	void solveMaster(solFeat &org, size_t &iterations);

	//Entities of the problem
	entitiesBendCP ent_sflp;

	//Master Environment
	IloEnv Mast_Bend = IloEnv();

	//Master model
	IloModel Mast_mod = IloModel(Mast_Bend);

	//Recovering solution information
	void RecoverSolution(IloCplex &master_cpx, solFeat &org);

	//Printing the solution
	void PrintSolution(solFeat &org);

	//Solution info
	vector<double> avg_scenarios;
	double obj_fin;
	double LB;
	double GAP;
	size_t status;
	size_t exploredNodes;
	double cpx_runtime;
	vector<double> current_theta;
};

#endif