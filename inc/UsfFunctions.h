// Created by CDRPico
// 11/06/2020 14:05



#ifndef USF_H
#define USF_H

#include<vector>
#include<iostream>
#include<random>
#include<sstream>
#include<chrono>
#include<unordered_set>
#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<utility> 
#include<Eigen/Dense>
#include<Eigen/Sparse>

using namespace std;

#define tolabscuts		1e-10
#define tolrelcuts		1e-4
#define timelimit		21600
#define tol_diff_disag	1e-4
#define timeslot		60
#define GAP_threshold           1e-3
#define expectation_threshold   1e-3
#define tol_difference 1e-6

//Useful Definitions CPLEX
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;


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

struct entitiesBendCP {
	//Variables
	IloNumVarArray x;
	IloNumVarArray theta;
	IloNumVarArray y;
	//Objective
	IloExpr objective;
	//Constraints
	IloRangeArray Feasibility;
};

class MyClock
{
public:
	chrono::steady_clock sc;
	chrono::steady_clock::time_point start;
	chrono::steady_clock::time_point end;
	MyClock();
};

struct solFeat {
	char algo;
	double feasCuts = 0;
	double optCuts = 0;
	size_t cnode = 0;
	size_t depth = 0;
	double user_feasCuts = 0;
	double user_optCuts = 0;
	vector<double> x_prev;
	bool part_modified = false;
	bool x_modified = false;

	solFeat() = default;
};

//Subproblem info solution
struct solution_sps {
	//Vector of scenarios summarized
	vector<size_t> scen;
	//Bool indicating the status of the solution
	// 1 will be optimal
	// 2 infeasible
	// 9 timelimit reached
	size_t F;
	//Objective value
	double obj;
	//dual multipliers of key constraints
	vector<double> lambda;
};

class disag_procedure 
{
public:
	disag_procedure() = default;
	~disag_procedure();
	//functions designed to refine the current partition
	void disaggregation(vector<solution_sps> &sp_info, vector<vector<size_t>> &partition, const double &nScenarios);
	//function which returns the new partition
	vector<vector<size_t>> refine(vector<solution_sps> &sp_info, vector<vector<size_t>> &partition, const double &nScenarios);
	//function to refine one element into the partition
	vector<vector<size_t>> refine_element(vector<solution_sps> &sp_info,vector<size_t> &element, const double &nScenarios);
	//function to campare scenarios pairwise (their duals)
	bool compare_duals(solution_sps &sp_info_s1, solution_sps &sp_info_s2, const size_t &s1, const size_t &s2);
	//Compute the partition element probs
	vector<double> compute_part_prob(vector<vector<size_t>> &partition, const double &nScenarios);
};

//Read instance data
template <typename T>
void read_instance_data(T &ProblemInstance, string &inst_name, string &stoch_inst_name);

vector<size_t> create_partition(const size_t &scenarios);

//Computing taxicab norm
double taxicab(vector<double> &vect);

//Computing norm-2
double norm2(vector<double> &vect);

//Compute worst case demand SFLP
double max_demand(vector<vector<double>> &stoch_dem);

bool smallerWithTolerance(double smaller, double larger);

bool compareSolutions(vector<double> x_prev, vector<double> cur_x);


template<class T>
void AddVarsMaster(T &BendersProb, const char &algo);

template<class T>
void AddVarsMaster(T &BendersProb, const char &algo, const char &linear);

template<class T>
void ValidInequalities(T &BendersProb, const char &algo);

template<class T>
void FeasibilityConstraint(T &BendersProb);

template <class T>
void DeleteAll(vector<T> &data, const vector<size_t>& deleteIndices);

class whole_partition;
struct part_characterisation;
struct eigen_characterisation;
class eigen_partition;
void DeleteAllWP(whole_partition & data, const vector<size_t>& deleteIndices);

void DeleteAllWP(eigen_partition & data, const vector<size_t>& deleteIndices);

bool is_integer(float k);

void remove_duplicates(std::vector<size_t> &v);

//Function to compare to vector of strings and check if boot are equal or not
template<class T>
bool compareVectors(vector<T> &base_case, vector<T> &compare);

template<class T>
pair<bool, size_t > findInVector(vector<T>  & vecOfElements, const T  & element);

pair <bool, size_t> findInVector(vector<double> & vecOfElements, double & element);

void removeColumn(Eigen::MatrixXd& matrix, size_t colToRemove);

template<class T>
void subtractScalar(vector<T> &vec, T &scalar);

#endif
