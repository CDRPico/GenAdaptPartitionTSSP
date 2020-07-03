// Created by CDRPico
// 11/06/2020 13:05

#include"../inc/GAPM.h"
#include"../inc/SFLP_GAPM.h"
#include"../inc/OuterBendersSFLP.h"
#include <chrono>

template<typename T>
void GAPM::gapm_algorithm(T &ProblemInstance, const char &algo) {
	//initialize size of stochastic parameters and duals
	partition = ProblemInstance.partition;
	//Creating the clock to control the execution time
	chrono::steady_clock sc;
	auto start = sc.now();
	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>>(end - start);
	execution_time = time_span.count();
	iterations = 0;
	GAP = 1;
	UB = DBL_MAX;
	LB = 0.0;
	termination = 9;

	//create master for outer benders
	IloModel master = ProblemInstance.MasterProblemCreation(removeobjs, algo);
	IloCplex cplex_master(master.getEnv());

	//subproblem creation
	ProblemInstance.SPProblemCreation_GRB();

	//main loop of the algorithm
	while (GAP > GAP_threshold && execution_time < timelimit) {
		if (algo == 'n') {
			termination = body_gapm(ProblemInstance, algo);
		} else {
			termination = body_gapm_outben(ProblemInstance, &cplex_master, master, algo);
		}
		
		//update execution time
		end = sc.now();
		auto time_span = static_cast<chrono::duration<double>>(end - start);
		execution_time = time_span.count();
		if (termination == 2) {
			goto outwhile;
		}
		//restart vectors
		sp_info.clear();
		stoch.clear();
		sp_info_agg.clear();
		stoch_agg.clear();
		//sp_info.lambda.resize(ProblemInstance.nScenarios);
		stoch.resize(ProblemInstance.nScenarios);
	}
outwhile:
	/* DISPLAY */
	cout << "CRP: "
		<< termination << " "
		<< execution_time << " "
		<< UB << " "
		<< LB << " "
		<< (UB - LB) / (1e-10 + LB) << " "
		<< iterations << " " //iterations rather than explored nodes
		<< 0 << " " //No pending nodes to be explored
		<< partition.size()
		<< endl;
}
template void GAPM::gapm_algorithm(SFLP_GAPM &ProblemInstance, const char &algo);
template void GAPM::gapm_algorithm(OuterBendersSFLP &ProblemInstance, const char &algo);

//This method makes an iteration of the algorithm for a given Master and subproblem
template <typename T>
size_t GAPM::body_gapm(T &ProblemInstance, const char &algo) {
	//Update number of iterations 
	iterations++;
	sp_info.resize(ProblemInstance.nScenarios);
	stoch.resize(ProblemInstance.nScenarios);
	//Here we need to create the master problem
	// Notice, this need s to be created after each iteration because we change both
	// number of constraints and variables, cplex does not have some of these features
	cout << "current iteration is " << iterations << endl << endl;
	bool removeobjs = false;
	if (iterations > 1) removeobjs = true;

	IloModel master = ProblemInstance.MasterProblemCreation(removeobjs, algo);
	IloCplex cplex_master(master.getEnv());
	//IloCplex* cplex_pointer = &cple

	//Solve the master model
	ProblemInstance.MasterProblemSolution(&cplex_master, master, LB, timelimit - execution_time);

	//Solve every subproblem
	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		vector<size_t> el(1, s);
		sp_info[s].scen = el;
		if (s == 0) {
			ProblemInstance.SPProblemModification_GRB(sp_info[s].scen, true);
		}
		else {
			ProblemInstance.SPProblemModification_GRB(sp_info[s].scen);
		}
		double obj = 0.0;
		ProblemInstance.SPProbleSolution_GRB(stoch[s], &sp_info[s], true);
	}

	//Compute the upper bound
	UB = 0;
	for (size_t i = 0; i < ProblemInstance.nFacilities; i++) {
		UB += ProblemInstance.x_bar[i] * ProblemInstance.fixed_costs[i];
	}

	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		UB += ProblemInstance.probability[s] * sp_info[s].obj;
	}
	//UB = ProblemInstance.compute_UB(sp_info);

	//Update solution Gap
	compute_gap();

	double continue_algorithm = 9;
	cout << "THE CURRENT GAP IS " << GAP << endl;
	if (GAP > GAP_threshold) {
		//run disaggregation procedure locally
		disag_procedure to_dis;
		to_dis.disaggregation(sp_info, partition, ProblemInstance.nScenarios);
		//reset stoch params agg and dual mult agg
		stoch_agg.clear();
		stoch_agg.resize(partition.size());
		sp_info_agg.clear();
		sp_info_agg.resize(partition.size());
		ProblemInstance.partition = partition;

		for (size_t s = 0; s < partition.size(); s++) {
			sp_info_agg[s].scen = partition[s];
			if (s == 0) {
				ProblemInstance.SPProblemModification_GRB(partition[s], true);
			}
			else {
				ProblemInstance.SPProblemModification_GRB(partition[s]);
			}
			double obj = 0.0;
			ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s]);
		}

		//run checking stopping criteria
		bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability);
		if (optimality == true) continue_algorithm = 2;
	}
	else {
		continue_algorithm = 1;
	}
	return continue_algorithm;
}
template size_t GAPM::body_gapm(SFLP_GAPM &ProblemInstance, const char &algo);

//This method makes an iteration of the algorithm for a given Master and subproblem
//Outer Benders
template <typename T>
size_t GAPM::body_gapm_outben(T &ProblemInstance,IloCplex *cplex_master, IloModel &master, const char &algo) {
	//Update number of iterations 
	iterations++;
	sp_info.resize(ProblemInstance.nScenarios);
	stoch.resize(ProblemInstance.nScenarios);
	//Here we need to create the master problem
	// Notice, this need s to be created after each iteration because we change both
	// number of constraints and variables, cplex does not have some of these features
	cout << "current iteration is " << iterations << endl << endl;
	bool removeobjs = false;
	if (iterations > 1) removeobjs = true;

	//Modify and Solve the master model
	ProblemInstance.MasterProblemModification(&cplex_master, master, algo);
	ProblemInstance.MasterProblemSolution(&cplex_master, master, LB, timelimit - execution_time);

	bool violated = false;
	//Solve every subproblem
	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		vector<size_t> el(1, s);
		sp_info[s].scen = el;
		if (s == 0) {
			ProblemInstance.SPProblemModification_GRB(sp_info[s].scen, true);
		}
		else {
			ProblemInstance.SPProblemModification_GRB(sp_info[s].scen);
		}
		double obj = 0.0;
		ProblemInstance.SPProbleSolution_GRB(stoch[s], &sp_info[s], true);
	}

	//Compute the upper bound
	UB = 0;
	for (size_t i = 0; i < ProblemInstance.nFacilities; i++) {
		UB += ProblemInstance.x_bar[i] * ProblemInstance.fixed_costs[i];
	}

	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		UB += ProblemInstance.probability[s] * sp_info[s].obj;
	}
	//UB = ProblemInstance.compute_UB(sp_info);

	//Update solution Gap
	compute_gap();

	double continue_algorithm = 9;
	cout << "THE CURRENT GAP IS " << GAP << endl;
	if (GAP > GAP_threshold) {
		//reset stoch params agg and dual mult agg
		stoch_agg.clear();
		stoch_agg.resize(partition.size());
		sp_info_agg.clear();
		sp_info_agg.resize(partition.size());
		ProblemInstance.partition = partition;

		for (size_t s = 0; s < partition.size(); s++) {
			sp_info_agg[s].scen = partition[s];
			if (s == 0) {
				ProblemInstance.SPProblemModification_GRB(partition[s], true);
			}
			else {
				ProblemInstance.SPProblemModification_GRB(partition[s]);
			}
			double obj = 0.0;
			ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s]);
		}

		ProblemInstance.new_cuts(stoch_agg, &sp_info_agg, violated);
		if (violated == false) {
			//run disaggregation procedure locally
			disag_procedure to_dis;
			to_dis.disaggregation(sp_info, partition, ProblemInstance.nScenarios);
		}

		//check current active cuts
		ProblemInstance.check_active_cuts(stoch, &sp_info);
		

		//run checking stopping criteria
		bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability);
		if (optimality == true) continue_algorithm = 2;
	}
	else {
		continue_algorithm = 1;
	}
	return continue_algorithm;
}
template size_t GAPM::body_gapm(OuterBendersSFLP &ProblemInstance, const char &algo);

void GAPM::compute_gap() {
	if (UB >= 1e300) {
		GAP = 1;
	}
	else {
		GAP = (UB - LB) / (1e-10 + LB);
	}
}

bool GAPM::stopping_criteria(const size_t &nScenarios, vector<double> &prob) {
	bool check_op = true;
	//aggregated expectation and atomized expectations for each element into the partition are computed
	for (size_t p = 0; p < partition.size(); p++) {
		double atom_expectation = 0.0;
		double element_prob = partition[p].size() / nScenarios;
		double agg_expectation = 0.0;
		for (size_t j = 0; j < sp_info_agg[p].lambda.size(); j++) {
			agg_expectation += (sp_info_agg[p].lambda[j] * stoch_agg[p][j]);
			for (size_t s = 0; s < partition[p].size(); s++) {
				atom_expectation += prob[partition[p][s]] / element_prob * (sp_info[partition[p][s]].lambda[j] * stoch[partition[p][s]][j]);
			}
		}

		//check equality of expectations, this is done for each element, when one element does not have the same value, procces is aborted
		if (fabs(agg_expectation - atom_expectation) > expectation_threshold) {
			check_op = false;
			break;
		}
	}
	return check_op;
}