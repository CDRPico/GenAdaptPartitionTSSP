// Created by CDRPico
// 11/06/2020 14:06

#include"../inc/UsefulFunctions.h"

template <typename T>
void read_instance_data(T &ProblemInstance, string &inst_name, string &stoch_inst_name){
    const char* fln = inst_name.c_str();
    ifstream file(fln);

    //read the general data of the instance
    ProblemInstance.read_instance(file);

    const char* fln_stoch = stoch_inst_name.c_str();
    ifstream file_stoch(fln_stoch);

    //read the scenarios info
    ProblemInstance.read_stoch_data(file_stoch);
}
template void read_instance_data(SFLP_GAPM &ProblemInstance, string &inst_name, string &stoch_inst_name);

//suitable function to create the initial partition
vector<size_t> create_partition(const size_t &Scenarios) {
	vector<size_t> partition;
	for (size_t i = 0; i < Scenarios; i++) {
		partition.push_back(i + 1);
	}
	return partition;
}


//taxicab norm
double taxicab(vector<double> &vect) {
	double taxicab = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		taxicab += fabs(vect[i]);
	}
	return taxicab;
}

//norm-2
double norm2(vector<double> &vect) {
	double norm = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		norm += pow(vect[i],2);
	}
	return norm;
}