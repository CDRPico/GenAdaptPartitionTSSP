// Created by CDRPico
// 11/06/2020 12:47

#pragma once

#ifndef GAPM_H
#define GAPM_H

#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<vector>
#include"../inc/UsfFunctions.h"

using namespace std;

class GAPM
{
private:
	/* data */
public:
	// Data to store solution information
	size_t num_stoch_param;
	size_t termination;
	double execution_time;
	size_t iterations;
	double LB;
	double UB;
	double GAP;
	vector<vector<size_t>> partition;

	// Stochatic parameters
	vector<vector<double>> stoch;

	// dual multipliers
	vector<solution_sps> sp_info;

	//Stochastic parameters aggregated for a certain partition
	vector<vector<double>> stoch_agg;

	//dual multipliers agg
	vector<solution_sps> sp_info_agg;

	//Function to run the entire algorithm
	template<typename T>
	void gapm_algorithm(T &ProblemInstance, const char &algo, const size_t &nsp);

	//Function to do an iteration of the algorithm
	template<class T>
	size_t body_gapm(T &ProblemInstance, const char &algo);

	template<class T>
	size_t body_gapm(T &ProblemInstance, const char &algo, vector<double> &prev_x);

	//outer benders gapm
	template <class T>
	size_t body_gapm(T &ProblemInstance, IloCplex *cplex_master, IloModel &master, const char &algo, vector<double> &prev_x, double &prev_lb);

	//Compute the current solution GAP
	void compute_gap();

	bool stopping_criteria(const size_t &nScenarios, vector<double> &prob, const size_t &stoch_param_num);

};

#endif // GAPM_H
