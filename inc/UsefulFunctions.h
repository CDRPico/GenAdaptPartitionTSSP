// Created by CDRPico
// 11/06/2020 14:05

using namespace std;
#include<vector>
#include<iostream>
#include<random>
#include<sstream>
#include"SFLP_GAPM.h"

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
//class SFLP_GAPM;

//Read instance data
template <typename T>
void read_instance_data(T &ProblemInstance, string &inst_name, string &stoch_inst_name);

vector<size_t> create_partition(const size_t &scenarios);

//Computing taxicab norm
double taxicab(vector<double> &vect);

//Computing norm-2
double norm2(vector<double> &vect);