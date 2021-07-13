// Created by CDRPico
// 07/10/2020 04:35

#pragma once

#ifndef CP_GAPM_H
#define CP_GAPM_H


#include<list>
#include"InstanceSFLP.h"
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<gurobi_c++.h>
#include<cfloat>
#include"UsfFunctions.h"
#include"SMPS_ElecPlan.h"
#include<Eigen/Dense>
#include<Eigen/Sparse>

//This structure will keep track of a matrix which charaterizes a subpartition
struct eigen_characterisation {
	//A sparse matrix which multiply the scenarios and then check which of them belong to this partition
	Eigen::SparseMatrix<double> subpart_matrix;
	//A vector which contains the values to check the inequalities
	Eigen::MatrixXd rhs_checkpart;
};

//A class to manipulate the entire partition (it includes many matrices which are multiplied by the scenarios)
class eigen_partition {
public:
	//Many matrices
	vector<eigen_characterisation> pr;

	//A constructor for each eigen_characterisation
	//it takes the estimated number of nonzeros into the matrix and the size of the rhs vector
	void instantiate_matrix(const size_t &nnz_reserve, const size_t &nineq);

	//method to add a new subpartition to control_subregions_new
	void add_subregion_store(Eigen::MatrixXd &break_values, Eigen::SparseMatrix<double> &subpart_matrix);

	//copy of the scenarios
	Eigen::MatrixXd cp_scen;
	vector<double> scenarios;
	
	//Copy scenarios from vector of vectors to dense matrix
	void copy_scenarios(vector<vector<double>> &full_scenarios);

	//A method to check which scenarios belong to each subpartition
	//Here the matrix multiplication is carried out
	vector<vector<size_t>> scenario_classif();
};

using namespace std;

ILOSTLBEGIN

struct part_characterisation {
	//store clients of a given tree
	vector<vector<size_t>> clients;
	//store signs of inequalities
	//L stands for smaller or equal
	//G stands for greater
	vector<vector<string>> ineq_sign;
	//store demand breakpoints
	vector<vector<double>> demand_breakpoints;
};

struct part_matrix_char {
	//for this implementation I will keep clients that belong to the current partition in a vector
	//A matrix that will contain the coefficients of the inequalities 
	vector<vector<double>> part_ineq_coef;
	//a vector of string containing the inequalities sign
	vector<string> ineq_sign;
	//the breakpoints of the demand
	vector<double> demand_breakpoints;
};

class whole_partition {
public:

	vector<part_characterisation> pr;

	//method to add a new subpartition to control_subregions_new
	void add_subregion_store(vector<vector<string>> &ineqs, vector<vector<double>> &break_values, vector<vector<size_t>> &clients);
};


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

	//Create a vector to map the position where each plant is allocated in terminals vector
	vector<size_t> ordered_terminals;
	//Create a vector to map the position where each client is allocated in dem_constr vector
	vector<size_t> ordered_dem_constr;
	void index_plants_clients();

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

	//store the total demand into a tree of the conex component
	vector<double> total_demand_tree;
	void TotalDemandTree(const size_t &s);

	//Check the cases in the conex component
	//demand > capacity
	//demand < capacity
	//demand = capacity
	//This method returns vector of vectors where each tree must be partitioned
	void CheckCasesComcon(const size_t &s);
	//Demand is integer but split points might be float, so we need to transform the result to integer splitting points


	//to compute partition prob based on the number of scenraios in an element
	vector<double> part_prob;

	//Type of GAPM (scenarios or full problem)
	char gapm_type;

	//To control the partition and breakpoints
	//vector of inequelity signs is not necessary here
	//part_characterisation general_tree_part;
	//Used to return the features a new partition for each subpartition in the current iteration
	part_characterisation comcon_part_tree;
	//vector of part_characterisation
	//each position will contain conditions that one subregion must satisfy
	whole_partition control_subregions_old;
	//at the beginning of each iteration new becomes old
	//from old we generate new subregions for each given old region
	whole_partition control_subregions_new;

	//-------------- THIS IS ALTERNATIVE TO MANAGE THE PARTITION USING EIGEN LIBRARY ------------------------------
	eigen_partition matrix_partition_old;
	eigen_partition matrix_partition_new;

	//using the vector of vector from CheckCasesComcon Method
	//we can update subregions
	void update_subregions(const size_t &iteration, const size_t &s);

	void update_subregions_eigen(const size_t &iteration, const size_t &s);

	//at each iteration the control_subregions_new will be initialized
	//either because it is the first iteration
	//or because s = 0 and first subregions must be added
	void initialize_subregions_it(whole_partition &control_tobe_updated);
	//the same method but using eigen matrices
	void initialize_subregions_eigen(eigen_partition &control_tobe_updated);

	//recursive method to generate all combinations of subregions
	void recursive_gen_subregions(vector<vector<string>> &ineqs, vector<vector<double>> &break_values, size_t &counting, whole_partition &control_tobe_updated);
	void recursive_gen_subregions_eigen(Eigen::MatrixXd &break_values, Eigen::SparseMatrix<double> &subpart_matrix, size_t &counting, eigen_partition &control_tobe_updated);

	//A method to update subregion by pushing back the new constraints of the current iterations
	void update_subregions_it(const size_t &s, whole_partition &new_control, whole_partition &control_tobe_updated);

	void update_subregions_it_eigen(const size_t &s, eigen_partition &new_control, eigen_partition &control_tobe_updated);

	//classify scenarios after creating subregions
	void classify_scenarios();
	bool classify_one_scenario(vector<double> &demand_tree, size_t &i);

	//check if a subregion already exists in the array of subregions
	bool sr_exists(part_characterisation &new_sr);

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

	//vector of vectors to store where trees of conex compnent will be partitioned
	// 1 position for greater and smaller or equal
	// 2 positions for smaller, between [], greater
	vector<vector<double>> where_part_tree;


	//Full problem for benchmark
	void FullCP();

	//An inner method to run the especial version of the GAPM
	void Inner_GAPM(const char &algo);

protected:
	//Master entities
	Master_CPX master_entities;

	//Gurobi env
	GRBEnv SFLP_sp_grb = GRBEnv();

	//Gurobi Model
	GRBModel subprob_grb = GRBModel(SFLP_sp_grb);

	//Subproblem entities Gurobi
	SubProblem_GRB sp_entities_grb;

private:

	//I'm gonna create an inner GAPM, especialized for subdividing regions
	//as the problems suggests
	double inner_lb_gapm;
	double inner_ub_gapm;
	double inner_gap;
	MyClock inner_gampm_runningtime;
	size_t inner_gapm_iterations;


};


#endif