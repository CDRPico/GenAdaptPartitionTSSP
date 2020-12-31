// Created by CDRPico
// 09/07/2020 05:42

#pragma once

#ifndef SFCMFP_GAPM_H
#define SFCMFP_GAPM_H


#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<gurobi_c++.h>
#include<cfloat>
#include"UsfFunctions.h"
#include"InstancesSFCMFP.h"

using namespace std;

ILOSTLBEGIN

//Definitions for Gurobi
struct SubProblemMFP_GRB {
	//Variables
	GRBVar **y = NULL;

	GRBConstr *flow_constraints = NULL;

	GRBConstr *cap_constraints = NULL;

	SubProblemMFP_GRB() = default;
};

struct Master_CPX_SFCMFP {
	//Variables
	IloNumVarArray x;
	IloNumVarArray3 y;
	//Objective
	IloExpr objective;
	//Constraints
	IloRangeArray flow_satisfaction;
	IloRangeArray cap_constraints;
};


class SFCMFP_GAPM : public InstanceSFCMFP
{
public:
	vector<vector<size_t>> partition;
	vector<double> obj_opt_sps;
	//Max demand, worst case scenario
	double max_dem;

	//COnstructor takes the name of the instance and the name of the stoch instance
	SFCMFP_GAPM(string &inst_name, string &stoch_inst_name);

	//Compute expectation of stochastic parameters
	vector<double> expected_demand(const vector<size_t> &scenarios);

	//Create master problem given a certain partition
	IloModel MasterProblemCreation(const char &algo, bool = false);

	// Solve master problem given a certain partition
	void MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL);

	//Create the subproblem CPLEX
	/*void SPProblemCreation_GRB();

	//Modify the subproblem to solve a different scenario CPLEX
	void SPProblemModification_GRB(vector<size_t> &element, bool = false);

	//Solve the subproblem
	void SPProbleSolution_GRB(vector<double> &stoch, solution_sps *sp_info, bool = false);

	// Given the current x_bar and solution of subprobs, the UB is computed
	//double compute_UB(vector<solution_sps> &sp_info, double solution_sps::*obj);

	// Labeling second stage variables
	void Labeling_y();

	//Labeling demand constraints
	void Label_flow_constr();
	//Labeling capacity constraints
	void Label_capacity_constr();*/

protected:
	//Master entities
	Master_CPX_SFCMFP master_entities;

	//Gurobi env
	GRBEnv SFLP_sp_grb = GRBEnv();

	//Gurobi Model
	GRBModel subprob_grb = GRBModel(SFLP_sp_grb);

	//Subproblem entities Gurobi
	SubProblemMFP_GRB sp_entities_grb;

	//Names second stage variables
	vector<vector<string>> Label_y;

	//Demand constraint label
	vector<string> Label_demconst;
	//Capacity constraint label
	vector<string> Label_capconst;
};

#endif