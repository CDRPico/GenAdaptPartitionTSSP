// Created by CDRPico
// 11/06/2020 13:05

#include"../inc/GAPM.h"
#include"../inc/SFLP_GAPM.h"
#include"../inc/OuterBendersSFLP.h"
#include"../inc/BendersAPM_SFLP.h"
#include"../inc/CapPlan_GAPM.h"
#include <chrono>

template<typename T>
void GAPM::gapm_algorithm(T &ProblemInstance, const char &algo, const size_t &nsp) {
	//initialize size of stochastic parameters and duals
	partition = ProblemInstance.partition;
	num_stoch_param = nsp;
	//First partition prob
	ProblemInstance.part_prob.resize(1);
	ProblemInstance.part_prob[0] = 1.0;
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
	IloModel master = ProblemInstance.MasterProblemCreation(false, algo);
	IloCplex cplex_master(master.getEnv());
	cplex_master.extract(master);
	cplex_master.exportModel("masterproblem_testCP.lp");

	//subproblem creation
	ProblemInstance.SPProblemCreation_GRB();
	/*stoch_agg.clear();
	stoch_agg.resize(partition.size());
	sp_info_agg.clear();
	sp_info_agg.resize(partition.size());
	ProblemInstance.SPProblemModification_GRB(partition[0], true);
	ProblemInstance.SPProbleSolution_GRB(stoch_agg[0], &sp_info_agg[0]);*/
	vector<double> prev_x(ProblemInstance.x_bar.size());
	double prev_lb = 0.0;

	//main loop of the algorithm
	while (GAP > GAP_threshold && execution_time < timelimit) {
		if (algo == 'n') {
			termination = body_gapm(ProblemInstance, algo);
		}
		/*else if (algo == 'o') {
			termination = body_gapm(ProblemInstance, &cplex_master, master, algo, prev_x, prev_lb);
			prev_x = ProblemInstance.x_bar;
			if (LB > prev_lb) prev_lb = LB;
		}
		else {
			termination = body_gapm(ProblemInstance, algo, prev_x);
			prev_x = ProblemInstance.x_bar;
		}*/

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
		<< partition.size() << " ";

	if (algo == 'o') {
		cout << ProblemInstance.cp.size() << endl;
	}
	else {
		cout << endl;
	}
}
template void GAPM::gapm_algorithm(SFLP_GAPM &ProblemInstance, const char &algo, const size_t &nsp);
template void GAPM::gapm_algorithm(OuterBendersSFLP &ProblemInstance, const char &algo, const size_t &nsp);
template void GAPM::gapm_algorithm(CP_GAPM &ProblemInstance, const char &algo, const size_t &nsp);

//This method makes an iteration of the algorithm for a given Master and subproblem
template <class T>
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
	cplex_master.extract(master);
	cplex_master.exportModel("masterproblem_testCP.lp");
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

	ofstream file;
	file.open("multiplicadores.csv");
	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		string scen;
		for (size_t i = 0; i < 70; i++) {
			if (i < 50) {
				scen = scen + " " + to_string(sp_info[s].lambda[i]) + " " + to_string(ProblemInstance.stoch_param[i][s]);
			}
			else {
				scen = scen + " " + to_string(sp_info[s].lambda[i]);
			}
		}
		scen = scen + "\n";
		file << scen;
	}
	file.close();

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
		if (ProblemInstance.gapm_type == 's') {
			ProblemInstance.part_prob = to_dis.compute_part_prob(partition, ProblemInstance.nScenarios);
		}
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
		bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability, num_stoch_param);
		if (optimality == true) continue_algorithm = 2;
	}
	else {
		continue_algorithm = 1;
	}
	return continue_algorithm;
}
template size_t GAPM::body_gapm(SFLP_GAPM &ProblemInstance, const char &algo);
template size_t GAPM::body_gapm(CP_GAPM &ProblemInstance, const char &algo);

//This method makes an iteration of the algorithm for a given Benders Master 
/*template <class T>
size_t GAPM::body_gapm(T &ProblemInstance, const char &algo, vector<double> &prev_x) {
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

	//Initialize benders
	//BendersSFLP MasterBenders();
	BendersSFLP MasterBenders(ProblemInstance.nFacilities,
						 ProblemInstance.nClients,
						 ProblemInstance.fixed_costs,
						 ProblemInstance.dist_costs,
						 ProblemInstance.facil_capacities,
						 ProblemInstance.stoch_param,
						 ProblemInstance.nScenarios,
						 ProblemInstance.probability,
						 ProblemInstance.x_bar,
						 partition);
	//Create and run procedure
	vector<vector<size_t>> part;
	part.push_back(vector<size_t>(create_partition(ProblemInstance.nScenarios)));
	double dummy;
	dummy = MasterBenders.CreateMaster(algo, part, ProblemInstance.x_bar);

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

	//check current active cuts
	ProblemInstance.check_active_cuts(stoch, &sp_info);

	//Compute the upper bound
	UB = 0;
	for (size_t i = 0; i < ProblemInstance.nFacilities; i++) {
		UB += ProblemInstance.x_bar[i] * ProblemInstance.fixed_costs[i];
	}

	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		UB += ProblemInstance.probability[s] * sp_info[s].obj;
	}
	//UB = ProblemInstance.compute_UB(sp_info);
	cout << "The current Upper Bound is " << UB << endl << endl;

	//Update solution Gap
	compute_gap();

	double continue_algorithm = 9;
	cout << "THE CURRENT GAP IS " << GAP << endl;
	if (GAP > GAP_threshold) {
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
			ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s], true);
		}

		bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability, ProblemInstance.nClients);
		if (optimality == true) return 2;

		//reset stoch params agg and dual mult agg
		ProblemInstance.new_cuts(stoch_agg, &sp_info_agg, violated);
		if (violated == false) {
			bool signal_remove = false;
			for (size_t i = 0; i < ProblemInstance.x_bar.size(); i++) {
				if (ProblemInstance.x_bar[i] != prev_x[i]) {
					signal_remove = true;
					break;
				}
			}
			//run disaggregation procedure locally
			disag_procedure to_dis;
			to_dis.disaggregation(sp_info, partition, ProblemInstance.nScenarios);
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
				ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s], true);
			}
			ProblemInstance.partition.clear();
			ProblemInstance.partition = partition;
			ProblemInstance.new_cuts(stoch_agg, &sp_info_agg, violated);
		}


		//run checking stopping criteria
		optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability, ProblemInstance.nClients);
		if (optimality == true) continue_algorithm = 2;
	}
	else {
		continue_algorithm = 1;
	}
	return continue_algorithm;
}
template size_t GAPM::body_gapm(OuterBendersSFLP &ProblemInstance, const char &algo, vector<double> &prev_x);*/

//This method makes an iteration of the algorithm for a given Master and subproblem
//Outer Benders
template <class T>
size_t GAPM::body_gapm(T &ProblemInstance, IloCplex *cplex_master, IloModel &master, const char &algo, vector<double> &prev_x, double &prev_lb) {
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
	ProblemInstance.MasterProblemModification(cplex_master, master, algo);
	ProblemInstance.MasterProblemSolution(cplex_master, master, LB, timelimit - execution_time);

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

	//check current active cuts
	ProblemInstance.check_active_cuts(stoch, &sp_info);

	//Compute the upper bound
	UB = 0;
	for (size_t i = 0; i < ProblemInstance.nFacilities; i++) {
		UB += ProblemInstance.x_bar[i] * ProblemInstance.fixed_costs[i];
	}

	for (size_t s = 0; s < ProblemInstance.nScenarios; s++) {
		UB += ProblemInstance.probability[s] * sp_info[s].obj;
	}
	//UB = ProblemInstance.compute_UB(sp_info);
	cout << "The current Upper Bound is " << UB << endl << endl;

	//Update solution Gap
	compute_gap();

	double continue_algorithm = 9;
	cout << "THE CURRENT GAP IS " << GAP << endl;
	if (GAP > GAP_threshold) {
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
			ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s], true);
		}

		bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability, ProblemInstance.nClients);
		if (optimality == true) return 2;

		//reset stoch params agg and dual mult agg
		ProblemInstance.new_cuts(stoch_agg, &sp_info_agg, violated);
		if (violated == false) {
			bool signal_remove = false;
			for (size_t i = 0; i < ProblemInstance.x_bar.size(); i++) {
				if (ProblemInstance.x_bar[i] != prev_x[i]) {
					signal_remove = true;
					break;
				}
			}
			if (signal_remove == true && LB > prev_lb) {
				ProblemInstance.MasterProblemModification(cplex_master, master, algo, true);
			}
			//run disaggregation procedure locally
			disag_procedure to_dis;
			to_dis.disaggregation(sp_info, partition, ProblemInstance.nScenarios);
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
				ProblemInstance.SPProbleSolution_GRB(stoch_agg[s], &sp_info_agg[s], true);
			}
			ProblemInstance.partition.clear();
			ProblemInstance.partition = partition;
			ProblemInstance.new_cuts(stoch_agg, &sp_info_agg, violated);
		}


		//run checking stopping criteria
		optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability, ProblemInstance.nClients);
		if (optimality == true) continue_algorithm = 2;
	}
	else {
		continue_algorithm = 1;
	}
	return continue_algorithm;
}
template size_t GAPM::body_gapm(OuterBendersSFLP &ProblemInstance, IloCplex *cplex_master, IloModel &master, const char &algo, vector<double> &prev_x, double &prev_lb);

void GAPM::compute_gap() {
	if (UB >= 1e300) {
		GAP = 1;
	}
	else {
		GAP = fabs(UB - LB) / (1e-10 + fabs(LB));
	}
}

bool GAPM::stopping_criteria(const size_t &nScenarios, vector<double> &prob, const size_t &stoch_param_num) {
	bool check_op = true;
	//aggregated expectation and atomized expectations for each element into the partition are computed
	for (size_t p = 0; p < partition.size(); p++) {
		double atom_expectation = 0.0;
		double element_prob = double(partition[p].size()) / double(nScenarios);
		double agg_expectation = 0.0;
		for (size_t j = 0; j < stoch_param_num; j++) {
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