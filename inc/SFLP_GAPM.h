// Created by CDRPico
// 11/06/2020 13:47

#pragma once

#ifndef SFLP_GAPM_H
#define SFLP_GAPM_H

using namespace std;

#include"InstanceSFLP.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<gurobi_c++.h>
#include<cfloat>
#include"UsfFunctions.h"

ILOSTLBEGIN

//Definitions for Gurobi
struct SubProblem_GRB {
	//Variables
	GRBVar **y = NULL;

	GRBConstr *dem_constraints = NULL;

	GRBConstr *cap_constraints = NULL;

	SubProblem_GRB() = default;
};

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
	//Max demand, worst case scenario
	double max_dem;
	vector<double> cp;

	//COnstructor takes the name of the instance and the name of the stoch instance
	SFLP_GAPM(string &inst_name, string &stoch_inst_name);
	SFLP_GAPM(const size_t &nFac, const size_t &nCli, const vector<double> &fc, const vector<vector<double>> &dc,
		const vector<double> &fcap, const vector<vector<double>> &sto_p, const size_t &nSc, const vector<double> &pr,
		const vector<double> &xb, const vector<vector<size_t>> &par);

	//Compute expectation of stochastic parameters
	vector<double> expected_demand(vector<size_t> &element);

	//Create master problem given a certain partition
	IloModel MasterProblemCreation(const char &algo, bool = false);

	// Solve master problem given a certain partition
	void MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL);

	//Create the subproblem CPLEX
	void SPProblemCreation_CPX();
	//Create the subproblem CPLEX
	void SPProblemCreation_GRB();

	//Modify the subproblem to solve a different scenario CPLEX
	void SPProblemModification_CPX(vector<size_t> &element, bool = false);
	//Modify the subproblem to solve a different scenario CPLEX
	void SPProblemModification_GRB(vector<size_t> &element, bool = false);

	//Solve the subproblem
	void SPProbleSolution_CPX(vector<double> &stoch, solution_sps *sp_info, bool = false);
	//Solve the subproblem
	void SPProbleSolution_GRB(vector<double> &stoch, solution_sps *sp_info, bool = false);

	//function declarations only to allow outerbenders running correctly templated
	void MasterProblemModification(IloCplex *cplex_master, IloModel &master, const char &algo, bool = false) {};

	//Check active cuts
	void check_active_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool = false) {};

	//New iteration cuts
	void new_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool &violated, bool = false) {};

	// Given the current x_bar and solution of subprobs, the UB is computed
	//double compute_UB(vector<solution_sps> &sp_info, double solution_sps::*obj);

	// Labeling second stage variables
	void Labeling_y();

	//Labeling demand constraints
	void Label_demand_constr();
	//Labeling capacity constraints
	void Label_capacity_constr();

	vector<double> part_prob;

	//Type of GAPM (scenarios or full problem)
	char gapm_type;

protected:
	//Master entities
	Master_CPX master_entities;

	//SubProblem env
	IloEnv SFLP_sp_cpx = IloEnv();

	//Subproblem model
	IloModel subprob_cpx = IloModel(SFLP_sp_cpx);

	//Subproblem entities CPLEX
	SubProblem_CPX sp_entities_cpx;

	//Gurobi env
	GRBEnv SFLP_sp_grb = GRBEnv();

	//Gurobi Model
	GRBModel subprob_grb = GRBModel(SFLP_sp_grb);

	//Subproblem entities Gurobi
	SubProblem_GRB sp_entities_grb;

	//Names second stage variables
	vector<vector<string>> Label_y;

	//Demand constraint label
	vector<string> Label_demconst;
	//Capacity constraint label
	vector<string> Label_capconst;
};

#endif