// Created by CDRPico
// 14/10/2020 17:35

#include"../inc/BendersAPM_CP.h"
#include<chrono>

void BendersCP::updateMaster(solFeat &org, size_t &nc, bool &violated)
{
	//Callback execution time
	chrono::steady_clock sc_cut;
	auto start_cut = sc_cut.now();
	auto end_cut = sc_cut.now();

	//Clear elements to store information of dual multipliers and subproblem solutions
	sp_info.clear();
	stoch.clear();
	sp_info.resize(nScenarios);
	stoch.resize(nScenarios);
	size_t scenarios = 0;
	//The number of scenarios will depend on the Benders implementation
	if (org.algo == 's' || org.algo == 'r') {
		scenarios = 1;
	}
	else {
		scenarios = nScenarios;
	}

	//Important elements for dual values
	stoch_agg.clear();
	stoch_agg.resize(partition.size());
	sp_info_agg.clear();
	sp_info_agg.resize(partition.size());
	lhs.clear();
	lhs.resize(partition.size());

	//Elements for single cut implementation
	double single_lhs = 0.0;
	IloExpr new_single_cut(Mast_Bend);
	double const_part_single = 0.0;

	//Checking if any cut should added
	for (size_t p = 0; p < partition.size(); p++) {
		sp_info_agg[p].scen = partition[p];
		//Solve the aggregated subproblem
		SPProblemModification_GRB(partition[p]);
		SPProbleSolution_GRB(stoch_agg[p], &sp_info_agg[p], true);

		//Compute the current rhs to check if the cut is violated
		double element_prob = 0.0;
		double theta_ac = 0.0;
		if (org.algo == 'a' || org.algo == 'm') {
			for (size_t s = 0; s < partition[p].size(); s++) {
				element_prob += probability[partition[p][s]];
				theta_ac += current_theta[partition[p][s]] * probability[partition[p][s]];
			}
		}
		else if (org.algo == 's') {
			element_prob = probability[partition[p][0]];
		}
		else if (org.algo == 'r'){
			for (size_t s = 0; s < partition[p].size(); s++) {
				element_prob += probability[partition[p][s]];
			}
		}

		//constant of current first stage multiplied by duals of capacity constarints
		double const_fs = 0.0;
		for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
			const_fs += element_prob * x_bar[i] * sp_info_agg[p].lambda[stoch_constr.size() + i];
		}
		double rhs = theta_ac;
		//compute the constant part of the cut
		lhs[p] = 0.0;
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			lhs[p] += element_prob * stoch_agg[p][j] * sp_info_agg[p].lambda[j];
		}
		lhs[p] += const_fs;

		//Adding the cut
		if (org.algo == 'a' || org.algo == 'm') {
			if (smallerWithTolerance(rhs, lhs[p])) {
				violated = true;
				IloExpr new_cut(Mast_Bend);
				for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
					new_cut += element_prob * ent_sflp.x[i] * sp_info_agg[p].lambda[stoch_constr.size() + i];
				}
				for (size_t s = 0; s < partition[p].size(); s++) {
					new_cut -= ent_sflp.theta[partition[p][s]] * probability[partition[p][s]];
				}
				double derecho = -1.0 * (lhs[p] - const_fs);
				//cout << new_cut << " <= " << derecho << endl;
				Mast_mod.add(new_cut <= derecho);
				new_cut.end();
				nc += 1;
				org.user_optCuts += 1;
			}
		}
		else {
			// building single cut
			single_lhs += lhs[p];
			const_part_single += lhs[p] - const_fs;
			for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
				new_single_cut += element_prob * ent_sflp.x[i] * sp_info_agg[p].lambda[stoch_constr.size() + i];
			}
		}
	}
	if (smallerWithTolerance(current_theta[0], single_lhs)) {
		if (org.algo == 's' || org.algo == 'r') {
			violated = true;
			new_single_cut -= ent_sflp.theta[0];
			Mast_mod.add(new_single_cut <= -const_part_single);
			nc += 1;
			org.user_optCuts += 1;
		}		
	}

	
	//Callback execution time
	end_cut = sc_cut.now();
	auto time_span_cut = static_cast<chrono::duration<double>>(end_cut - start_cut);
	double runtime_cut = time_span_cut.count();
	cout << "it takes " << runtime_cut << "to execute the lazy callback " << endl;
}

void BendersCP::disaggPartition(solFeat &org, bool &violated) {
	//Here we redefine org.par_modified to know in the future if the partition was modiied or not
	org.part_modified = false;
	if (violated == false && (org.algo == 'a' || org.algo == 'r') && partition.size() < nScenarios) {// && wherefrom == CPX_CALLBACK_MIP_INCUMBENT_NODESOLN) {
		//chrono::steady_clock sc;
		//auto start = sc.now();
		//auto end = sc.now();
		//partition is modified
		size_t cur_part_size = partition.size();
		for (size_t s = 0; s < nScenarios; s++) {
			vector<size_t> el(1, s);
			sp_info[s].scen = el;
			if (s == 0) {
				SPProblemModification_GRB(sp_info[s].scen, true);
			}
			else {
				SPProblemModification_GRB(sp_info[s].scen);
			}
			double obj = 0.0;
			SPProbleSolution_GRB(stoch[s], &sp_info[s], true);
		}
		//end = sc.now();
		//auto time_span = static_cast<chrono::duration<double>>(end - start);
		//double runtime = time_span.count();
		//cout << "it takes " << runtime << "to solve " << nScenarios << " subproblems" << endl;;
		disag_procedure to_dis;
		to_dis.disaggregation(sp_info, partition, nScenarios);
		if (cur_part_size != partition.size()) {
			org.part_modified = true;
		}
		//Add cuts for the new partition
		cutsindiasgg(org);
	}
}

void BendersCP::cutsindiasgg(solFeat &org) {
	double sum_theta = 0.0;
	for (size_t s = 0; s < nScenarios; s++) {
		sum_theta += current_theta[s];
	}

	for (size_t p = 0; p < partition.size(); p++) {
		double element_prob = 0.0;
		vector<double> elem_demand = expected_demand(partition[p]);
		for (size_t s = 0; s < partition[p].size(); s++) {
			element_prob += probability[partition[p][s]];
		}
		double theta_ac = element_prob * sum_theta;//element_prob * 

		//constant of current first stage multiplied by duals of capacity constarints
		double const_fs = 0.0;
		for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
			const_fs += element_prob * x_bar[i] * sp_info[partition[p][0]].lambda[stoch_constr.size() + i];
		}
		double rhs = theta_ac;
		//compute the constant part of the cut
		double lhs_p = 0.0;
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			lhs_p += element_prob * elem_demand[j] * sp_info[partition[p][0]].lambda[j];
		}
		lhs_p += const_fs;

		if (smallerWithTolerance(rhs, lhs_p)) {
			IloExpr new_cut(Mast_Bend);
			for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
				new_cut += element_prob * ent_sflp.x[i] * sp_info[partition[p][0]].lambda[stoch_constr.size() + i];
			}
			for (size_t s = 0; s < partition[p].size(); s++) {
				new_cut -= ent_sflp.theta[partition[p][s]] * probability[partition[p][s]];
			}
			double derecho = -1.0 * (lhs_p - const_fs);
			//cout << new_cut << " <= " << derecho << endl;
			Mast_mod.add(new_cut <= derecho);
			new_cut.end();
			org.user_optCuts += 1;
		}
	}

}


void BendersCP::CreateMaster(const char &algo, vector<vector<size_t>> &part, vector<double> &xb) {
	partition = part;
	//Add variables
	AddVarsMaster(*this, algo, 'l');

	//first stage constraints
	//Adding the first stage constraints
	IloRangeArray ax_b(Mast_Bend);
	for (size_t i = 1; i < first_st_const.size(); i++) {
		IloExpr constr(Mast_Bend);
		for (size_t j = 0; j < first_st_var.size(); j++) {
			constr += first_st_coefs[j][i] * ent_sflp.x[j];
		}
		if (first_st_constsense[i] == "E") {
			ax_b.add(constr == first_st_rhs[i]);
		}
		else if (first_st_constsense[i] == "L") {
			ax_b.add(constr <= first_st_rhs[i]);
		}
		else if (first_st_constsense[i] == "G") {
			ax_b.add(constr >= first_st_rhs[i]);
		}
		constr.end();
	}
	Mast_mod.add(ax_b);

	//Add valid inequalities
	ent_sflp.y = IloNumVarArray(Mast_Bend, second_st_var.size());
	for (size_t i = 0; i < second_st_var.size(); i++) {
		ent_sflp.y[i] = IloNumVar(Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
		Mast_mod.add(ent_sflp.y[i]);
	}

	size_t counter = 0;
	//Compute the average demand
	IloRangeArray sec_const(Mast_Bend);
	vector<vector<size_t>> part_start;
	part_start.push_back(vector<size_t>(create_partition(nScenarios)));
	vector<double> elem_demand = expected_demand(part_start[0]);
	part_start.clear();
	for (size_t i = first_st_const.size(); i < second_st_coefs[0].size(); i++) {
		IloExpr constr(Mast_Bend);
		for (size_t j = 0; j < first_st_var.size(); j++) {
			constr += first_st_coefs[j][i] * ent_sflp.x[j];
		}
		for (size_t j = 0; j < second_st_var.size(); j++) {
			constr += second_st_coefs[j][i] * ent_sflp.y[j];
		}
		if (counter < second_st_const.size() - stoch_constr.size()) {
			if (second_st_constsense[i - first_st_const.size()] == "E") {
				sec_const.add(constr == second_st_rhs[i - first_st_const.size()]);
			}
			else if (second_st_constsense[i - first_st_const.size()] == "L") {
				sec_const.add(constr <= second_st_rhs[i - first_st_const.size()]);
			}
			else if (second_st_constsense[i - first_st_const.size()] == "G") {
				sec_const.add(constr >= second_st_rhs[i - first_st_const.size()]);
			}
		}
		else {
			if (second_st_constsense[i - first_st_const.size()] == "E") {
				sec_const.add(constr == elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
			}
			else if (second_st_constsense[i - first_st_const.size()] == "L") {
				sec_const.add(constr <= elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
			}
			else if (second_st_constsense[i - first_st_const.size()] == "G") {
				sec_const.add(constr >= elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
			}
		}
		counter += 1;
		constr.end();
	}
	Mast_mod.add(sec_const);

	//Adding the inequality ensuring the sum of thetas is greater than the cost computed
	IloRangeArray VI(Mast_Bend);
	IloExpr sec_st_cost(Mast_Bend);
	for (size_t i = 0; i < second_st_var.size(); i++) {
		sec_st_cost -= second_st_coefs[i][0] * ent_sflp.y[i];
	}
	if (algo == 's' || algo == 'r') {
		sec_st_cost += ent_sflp.theta[0];
	}
	else {
		for (size_t s = 0; s < nScenarios; s++) {
			sec_st_cost += probability[s] * ent_sflp.theta[s];
		}
	}
	VI.add(sec_st_cost >= 0);
	sec_st_cost.end();
	Mast_mod.add(VI);

}

void BendersCP::solveMaster(solFeat &org) {
	/* Create cplex enviroment to solve the problem */
	
	//Setting parameters to solve master problem
	IloCplex cplex(Mast_Bend);
	cplex.setParam(IloCplex::Param::Preprocessing::Reduce, CPX_PREREDUCE_NOPRIMALORDUAL);
	cplex.setParam(IloCplex::Param::Preprocessing::Linear, IloFalse);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
	cplex.setParam(IloCplex::Param::TimeLimit, timelimit);

	cplex.extract(Mast_mod);
	cplex.exportModel("masterBenders_CP.lp");

	//cplex.use(BendersCallback_CP(this->Mast_Bend, *this, org));

	//Solve the model
	cplex.solve();

	//Recovering the solution
	RecoverSolution(cplex, org);
}

void BendersCP::RecoverSolution(IloCplex &master_cpx, solFeat &org) {
	try {
		x_bar.resize(first_st_var.size());
		for (size_t i = 0; i < first_st_var.size(); i++) {
			x_bar[i] = master_cpx.getValue(ent_sflp.x[i]);
			if (x_bar[i] > 0.0)
				cout << "Plant " << first_st_var[i] << " has capacity " << x_bar[i] << " assigned." << endl;
		}

		size_t scenarios;
		if (org.algo == 's' || org.algo == 'r') {
			scenarios = 1;
		}
		else {
			scenarios = nScenarios;
		}

		current_theta.resize(scenarios);
		for (size_t s = 0; s < scenarios; s++) {
			current_theta[s] = master_cpx.getValue(ent_sflp.theta[s]);
		}

		//statistics
		obj_fin = master_cpx.getObjValue();
		LB = master_cpx.getBestObjValue();
		GAP = master_cpx.getMIPRelativeGap();
		status = master_cpx.getStatus();
		exploredNodes = master_cpx.getNnodes();
		cpx_runtime = master_cpx.getTime();
	}
	catch (IloException& ex) {
		cerr << " ERROR: " << ex << endl;
	}
}

void BendersCP::PrintSolution(solFeat &org) {
	cout << "Final partition size is: " << partition.size() << endl;
	cout << "cap_plan:"
		<< status << " "
		<< cpx_runtime << " "
		<< obj_fin << " "
		<< LB << " "
		<< GAP << " "
		<< exploredNodes << " "
		<< 0 << " "//0 nodes pending to be explored
		<< partition.size() << " "
		<< org.optCuts << " "
		<< org.user_optCuts << " "
		<< org.feasCuts
		<< endl;
}

double BendersCP::runBenders(const char &algo, vector<vector<size_t>> &part, vector<double> &xb) {
	//Benders executon time execution time
	chrono::steady_clock sc_cut;
	auto start_cut = sc_cut.now();
	auto end_cut = sc_cut.now();
	
	
	CreateMaster(algo, part, xb);

	//Subproblem Creation
	SPProblemCreation_GRB();

	//Gathering solution info
	solFeat org;
	org.algo = algo;
	org.feasCuts = 1;
	org.optCuts = 0;
	org.user_optCuts = 0;
	org.x_prev.resize(nFacilities, 0.0);

	//solve first master
	solveMaster(org);

	if (algo == 'm') {
		sp_info.clear();
		stoch.clear();
		sp_info.resize(nScenarios);
		stoch.resize(nScenarios);
		for (size_t s = 0; s < nScenarios; s++) {
			vector<size_t> el(1, s);
			sp_info[s].scen = el;
			if (s == 0) {
				SPProblemModification_GRB(sp_info[s].scen, true);
			}
			else {
				SPProblemModification_GRB(sp_info[s].scen);
			}
			double obj = 0.0;
			SPProbleSolution_GRB(stoch[s], &sp_info[s], true);
		}

		cutsindiasgg(org);
		solveMaster(org);
	}/**/
	/*else if (algo == 'a') {
		sp_info.clear();
		stoch.clear();
		sp_info.resize(nScenarios);
		stoch.resize(nScenarios);
		for (size_t s = 0; s < nScenarios; s++) {
			vector<size_t> el(1, s);
			sp_info[s].scen = el;
			if (s == 0) {
				SPProblemModification_GRB(sp_info[s].scen, true);
			}
			else {
				SPProblemModification_GRB(sp_info[s].scen);
			}
			double obj = 0.0;
			SPProbleSolution_GRB(stoch[s], &sp_info[s], true);
		}

		double sum_theta = 0.0;
		for (size_t s = 0; s < nScenarios; s++) {
			sum_theta += current_theta[s];
		}


		vector<vector<size_t>> partition1;
		for (size_t s = 0; s < nScenarios; s++) {
			partition1.push_back({ s });
		}

		for (size_t p = 0; p < partition1.size(); p++) {
			double element_prob = 0.0;
			vector<double> elem_demand = expected_demand(partition1[p]);
			for (size_t s = 0; s < partition1[p].size(); s++) {
				element_prob += probability[partition1[p][s]];
			}
			double theta_ac = element_prob * sum_theta;//element_prob * 

			//constant of current first stage multiplied by duals of capacity constarints
			double const_fs = 0.0;
			for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
				const_fs += element_prob * x_bar[i] * sp_info[partition1[p][0]].lambda[stoch_constr.size() + i];
			}
			double rhs = theta_ac;
			//compute the constant part of the cut
			double lhs_p = 0.0;
			for (size_t j = 0; j < stoch_constr.size(); j++) {
				lhs_p += element_prob * elem_demand[j] * sp_info[partition1[p][0]].lambda[j];
			}
			lhs_p += const_fs;

			if (smallerWithTolerance(rhs, lhs_p)) {
				IloExpr new_cut(Mast_Bend);
				for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
					new_cut += element_prob * ent_sflp.x[i] * sp_info[partition1[p][0]].lambda[stoch_constr.size() + i];
				}
				for (size_t s = 0; s < partition1[p].size(); s++) {
					new_cut -= ent_sflp.theta[partition1[p][s]] * probability[partition1[p][s]];
				}
				double derecho = -1.0 * (lhs_p - const_fs);
				//cout << new_cut << " <= " << derecho << endl;
				Mast_mod.add(new_cut <= derecho);
				new_cut.end();
				org.user_optCuts += 1;
			}
		}

		solveMaster(org);
	}*/

	double prev_opt = -DBL_MAX;

	GAP = fabs(prev_opt - obj_fin) / (1e-10 + fabs(prev_opt));

	size_t nc = 10;
	while (nc != 0 || org.part_modified == true) {
		//Add cuts
		prev_opt = obj_fin;
		nc = 0;
		bool violated = false;
		org.part_modified = false;
		if (GAP > 1e-8) {
			updateMaster(org, nc, violated);
			//update partition
			disaggPartition(org, violated);
			//solve the master updated
			solveMaster(org);
			org.feasCuts += 1;
		}
		GAP = fabs(prev_opt - obj_fin) / (1e-10 + fabs(prev_opt));
	}

	for (size_t i = 0; i < first_st_var.size(); i++) {
		if (x_bar[i] > 0.0)
			cout << "Plant " << first_st_var[i] << " has capacity " << x_bar[i] << " assigned." << endl;
	}

	//Printing the solution
	PrintSolution(org);
	part = partition;


	end_cut = sc_cut.now();
	auto time_span_cut = static_cast<chrono::duration<double>>(end_cut - start_cut);
	double runtime_cut = time_span_cut.count();
	cout << "it takes " << runtime_cut << "to execute the benders algorithm " << endl;

	return LB;
}


double BendersCP::mixedBenders(const char &algo, vector<vector<size_t>> &part, vector<double> &xb) {
	//Benders executon time execution time
	chrono::steady_clock sc_cut;
	auto start_cut = sc_cut.now();
	auto end_cut = sc_cut.now();


	CreateMaster(algo, part, xb);

	//Subproblem Creation
	SPProblemCreation_GRB();

	//Gathering solution info
	solFeat org;
	org.algo = algo;
	org.feasCuts = 0;
	org.optCuts = 0;
	org.user_optCuts = 0;
	org.x_prev.resize(nFacilities, 0.0);

	//solve first master
	solveMaster(org);

	size_t nc = 10;
	while (nc != 0 || org.part_modified == true) {
		//Add cuts
		nc = 0;
		bool violated = false;
		updateMaster(org, nc, violated);

		org.x_prev.resize(nFacilities, 0.0);
		for (size_t i = 0; i < nFacilities; i++) {
			org.x_prev[i] = x_bar[i];
		}
		//solve the master updated
		double prevobj = obj_fin;
		solveMaster(org);
		

		if (compareSolutions(org.x_prev, x_bar)) {
			//update partition
			disaggPartition(org, violated);
			violated = true;
			//run single cut
			const char aux_algo = 's';
			vector<vector<size_t>> cur_part;
			for (size_t p = 0; p < partition.size(); p++) {
				cur_part.push_back({ partition[p] });
			}
			partition.clear();
			for (size_t s = 0; s < nScenarios; s++) {
				partition.push_back(vector<size_t>({ s }));
			}

			solFeat aux_org;
			aux_org.algo = aux_algo;
			aux_org.feasCuts = 0;
			aux_org.optCuts = 0;
			aux_org.user_optCuts = 0;
			aux_org.x_prev.resize(nFacilities, 0.0);
			for (size_t i = 0; i < nFacilities; i++) {
				aux_org.x_prev[i] = x_bar[i];
			}
			

			while (compareSolutions(aux_org.x_prev, x_bar) && fabs(prevobj - obj_fin)>tolabscuts){
				size_t nc1 = 0;
				bool violated1 = false;
				updateMaster(aux_org, nc1, violated1);
				//solve the master updated
				double prevobj = obj_fin;
				solveMaster(org);
			}

			partition.clear();
			for (size_t p = 0; p < cur_part.size(); p++) {
				partition.push_back({ cur_part[p] });
			}
			cur_part.clear();
		}

		
	}

	for (size_t i = 0; i < first_st_var.size(); i++) {
		if (x_bar[i] > 0.0)
			cout << "Plant " << first_st_var[i] << " has capacity " << x_bar[i] << " assigned." << endl;
	}

	//Printing the solution
	PrintSolution(org);
	part = partition;


	end_cut = sc_cut.now();
	auto time_span_cut = static_cast<chrono::duration<double>>(end_cut - start_cut);
	double runtime_cut = time_span_cut.count();
	cout << "it takes " << runtime_cut << "to execute the benders algorithm " << endl;

	return LB;
}