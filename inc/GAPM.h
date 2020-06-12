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

#define timelimit 21600
#define tol_diff_disag 1e-4
#define GAP_threshold 1e-4

class GAPM
{
private:
    /* data */
public:
    // Data to store solution information
    double LB;
    double UB;
    vector<vector<size_t>> partition;

    // Stochatic parameters
    vector<vector<double>> stoch;

    // dual multipliers
    vector<vector<double>> lambda;

    //Function to do an iteration of the algorithm
    template<typename T>
    void body_gapm(T &ProblemInstance);

    //functions designed to refine the current partition
    void disaggregation(const double &nScenarios, vector<vector<double>> &lambda);
    //function which returns the new partition
    vector<vector<size_t>> refine(const double &nScenarios, vector<vector<double>> &lambda);
    //function to refine one element into the partition
    vector<vector<size_t>> refine_element(vector<size_t> &element, const double &nScenarios, vector<vector<double>> &lambda);
    //function to campare scenarios pairwise (their duals)
    bool compare_duals (const vector<vector<double>> &lambda, const size_t &s1, const size_t &s2);

    size_t stopping_criteria(const vector<double> &RHS, vector<vector<double>> &lambda);

};

#endif // GAPM_H