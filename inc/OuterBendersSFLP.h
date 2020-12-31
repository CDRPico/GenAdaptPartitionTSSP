// Created by CDRPico
// 23/06/2020 01:14

#pragma once

#ifndef OuterBendersGAPM
#define OuterBendersGAPM

#include"SFLP_GAPM.h"
#include"UsfFunctions.h"

struct cut_pool {
	//CPLEX cuts preserved from one iteration to the next one
	//Scenarios for a certain cut
	vector<size_t> elem_cut;
	//cut name
	string cut_name;
	//is active
	bool active;
	//x coefficients
	vector<double> x_coefs;
	//theta coefficients
	vector<double> theta_coefs;
	IloRange cut_pool_constr;
	//constructor
	cut_pool() {}
	cut_pool(const IloEnv *env, const IloExpr *cut_pool_lhs, double rhs, vector<size_t> elem_cut, const string &cut_name, bool active, vector<double> x_coefs, vector<double> theta_coefs) :
		elem_cut(elem_cut), cut_name(cut_name), active(active), x_coefs(x_coefs), theta_coefs(theta_coefs), cut_pool_constr(*env, *cut_pool_lhs, rhs) {
		//IloRange cut_pool_constr = IloRange(*env, *cut_pool_lhs, rhs);
	}
};

class OuterBendersSFLP : public SFLP_GAPM
{
public:
	//entities of the master problem
	entitiesBendSFLP ent_sflp;

	//struct to store pool of cuts
	vector<cut_pool> cp;

	//to store current solution
	//x_bar inherited from base class
	vector<double> current_theta;

	//store old number of cuts (useful to remove old inactive cuts)
	size_t num_old_cuts;

	//Instantiate using constructor of base class
	OuterBendersSFLP(string &inst_name, string &stoch_inst_name);
	~OuterBendersSFLP() = default;

	//Create master problem given a certain partition
	IloModel MasterProblemCreation(const char &algo, bool = false);

	//Modify master problem using pool cuts
	void MasterProblemModification(IloCplex *cplex_master, IloModel &master, const char &algo, bool = false);

	void MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL);

	//Check active cuts
	void check_active_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool = false);

	//New iteration cuts
	void new_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool &violated, bool = false);

	//Remove inactive constraints
	void remove_inactive();


	//Master Environment
	IloEnv Mast_Bend = IloEnv();

	//Master model
	IloModel Mast_mod = IloModel(Mast_Bend);

	//Master cplex
	IloCplex Mast_cplex = IloCplex(Mast_Bend);

	vector<double> avg_scenarios;
};

#endif