// Created by CDRPico
// 07/10/2020 04:35

#pragma once

#ifndef CP_GAPM_H
#define CP_GAPM_H

using namespace std;

#include<list>
#include"InstanceSFLP.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<gurobi_c++.h>
#include<cfloat>
#include"UsfFunctions.h"
#include"SMPS_ElecPlan.h"


ILOSTLBEGIN


class CP_GAPM : public Inst_ElecPlan {
public:
	struct Master_CPX {
		//Variables
		IloNumVarArray x;
		IloNumVarArray2 y;
		//Objective
		IloExpr objective;
		//Constraints
		IloRangeArray ax_b;
		IloRangeArray linking_constraints;
		IloRangeArray Feasibility;
	};

	//Definitions for Gurobi
	struct SubProblem_GRB {
		//Variables
		GRBVar *y = NULL;

		GRBConstr *link_constr = NULL;

		SubProblem_GRB() = default;
	};

	vector<vector<size_t>> partition;
	vector<double> obj_opt_sps;
	//Max demand, worst case scenario
	double max_dem;
	vector<double> cp;
	size_t nFacilities;

	//COnstructor takes the name of the instance and the name of the stoch instance
	CP_GAPM(string &inst_name, string &stoch_name, string &time_name, string &stoch_rnd);

	//Compute expectation of stochastic parameters
	vector<double> expected_demand(vector<size_t> &element);

	//Create master problem given a certain partition
	IloModel MasterProblemCreation(const char &algo, bool = false);

	// Solve master problem given a certain partition
	void MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL);

	//Create the subproblem CPLEX
	void SPProblemCreation_CPX();
	//Create the subproblem Gurobi
	void SPProblemCreation_GRB();

	//Modify the subproblem to solve a different scenario CPLEX
	void SPProblemModification_CPX(vector<size_t> &element, bool = false);
	void SPProblemModification_GRB(vector<size_t> &element, bool = false);

	//Solve the subproblem
	void SPProbleSolution_CPX(vector<double> &stoch, solution_sps *sp_info, bool = false);
	void SPProbleSolution_GRB(vector<double> &stoch, solution_sps *sp_info, bool = false);

	//function declarations only to allow outerbenders running correctly templated
	void MasterProblemModification(IloCplex *cplex_master, IloModel &master, const char &algo, bool = false) {};

	//Check active cuts
	void check_active_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool = false) {};

	//New iteration cuts
	void new_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool &violated, bool = false) {};

	//Identify the conex component of a solution
	//Store origins and destinations which appear in a certain solution
	vector<size_t> ori_sol;
	vector<size_t> dest_sol;
	//Origin and destination nodes for the trees into the conex component
	vector<vector<size_t>> ori_comcon_sol;
	vector<vector<size_t>> dest_comcon_sol;
	//Store demand and capacity assigned on a certain solution for a tree in the conex component
	vector<double> flows_comcon_sol;
	vector<double> flows_orig_comcon_sol;
	vector<vector<double>> caps_sol;

	//Order of demands in constraints
	vector<size_t> dem_constr;
	void link_dem_con();

	//Generate the conex component of a given solution
	//This function depends on s, because the conex component might be different foir each scenario
	void LookCompCon(const size_t &s);

	//Flows in the conex component
	//Function to compute the total flow in a tree of the comp conex
	void FlowsComCon(const size_t &s);
	//Capacity of plants in the same tree of the conex component
	void CapacityComCon(const size_t &s);

	//Check the cases in the conex component
	//demand > capacity
	//demand < capacity
	//demand = capacity
	void CheckCasesComcon(const size_t &s);

	//to compute partition prob based on the number of scenraios in an element
	vector<double> part_prob;

	//Type of GAPM (scenarios or full problem)
	char gapm_type;


	vector<vector<vector<vector<size_t>>>> subpartitions;
	vector<vector<size_t>> subpart_clients;
	vector<vector<vector<double>>> prob_sp_scen;
	vector<vector<double>> prob_sp_acum;
	vector<vector<vector<double>>> expected_subpart;

	void GenSubpartitions();

	void GenScen(vector<size_t> client, vector<size_t> &sc, const size_t &size_v, const size_t &k, double &new_prob);

	vector<vector<double>> master_scens_sp;
	vector<double> prob_scens;
	vector<size_t> clients_attended;
	void GenExpectedScen();
	void join_scen(vector<double> new_scen, size_t k, double prob);

	vector<vector<double>> sp_scenarios_full;
	void FinalScenariosSubparts();

	//Identify total demand-offer in a solution tree
	double difDemOf(size_t &nTree, vector<size_t> &offnodes, vector<size_t> demnodes);

	//Full problem for benchmark
	void FullCP();


protected:
	//Master entities
	Master_CPX master_entities;

	//Gurobi env
	GRBEnv SFLP_sp_grb = GRBEnv();

	//Gurobi Model
	GRBModel subprob_grb = GRBModel(SFLP_sp_grb);

	//Subproblem entities Gurobi
	SubProblem_GRB sp_entities_grb;

};


#endif