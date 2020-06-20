// Created by CDRPico
// 11/06/2020 12:47

#pragma once

#ifndef GAPM_H
#define GAPM_H

#include"ilcplex/ilocplex.h"
#include"ilconcert/iloiterator.h"
#include<vector>
#include"../inc/UsefulFunctions.h"

using namespace std;

#define timelimit               21600
#define tol_diff_disag          1e-4
#define GAP_threshold           1e-3
#define expectation_threshold   1e-3

class GAPM
{
private:
	/* data */
public:
	// Data to store solution information
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
	void gapm_algorithm(T &ProblemInstance);

	//Function to do an iteration of the algorithm
	template<typename T>
	size_t body_gapm(T &ProblemInstance);

	//functions designed to refine the current partition
	void disaggregation(const double &nScenarios);
	//function which returns the new partition
	vector<vector<size_t>> refine(const double &nScenarios);
	//function to refine one element into the partition
	vector<vector<size_t>> refine_element(vector<size_t> &element, const double &nScenarios);
	//function to campare scenarios pairwise (their duals)
	bool compare_duals(const size_t &s1, const size_t &s2);

	//Compute the current solution GAP
	void compute_gap();

	bool stopping_criteria(const size_t &nScenarios, vector<double> &prob);

};

#endif // GAPM_H
