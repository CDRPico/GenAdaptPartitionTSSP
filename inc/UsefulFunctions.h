// Created by CDRPico
// 11/06/2020 14:05

using namespace std;
#include<vector>
#include<iostream>
#include<random>
#include<sstream>

class SFLP_GAPM;

//Read instance data
template <typename T>
string read_instance_data(T &ProblemInstance, string &inst_name, string &stoch_inst_name);

vector<size_t> create_partition(const size_t &scenarios);

//Computing taxicab norm
double taxicab(vector<double> &vect);

//Computing norm-2
double norm2(vector<double> &vect);