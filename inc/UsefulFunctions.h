// Created by CDRPico
// 11/06/2020 14:05

using namespace std;
#include<vector>
#include<iostream>
#include<random>
#include<sstream>
#include"SFLP_GAPM.h"

#define tolabscuts		1e-10
#define tolrelcuts		1e-4
#define timelimit		21600
#define tol_diff_disag	1e-4

#ifndef USF_H
#define USF_H

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

class SFLP_GAPM;

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

#endif