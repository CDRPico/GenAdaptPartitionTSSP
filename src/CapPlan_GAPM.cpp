// Created by CDRPico
// 07/10/2020 04:35

#include"../inc/CapPlan_GAPM.h"
#include<algorithm>

CP_GAPM::CP_GAPM(string &inst_name, string &stoch_name, string &time_name, string &stoch_rnd) {
	//read files
	const char* fln = time_name.c_str();
	ifstream filest(fln);
	read_time(filest);

	const char* fln1 = inst_name.c_str();
	ifstream filesc(fln1);
	read_core(filesc);

	const char* fln2 = stoch_name.c_str();
	ifstream filess(fln2);
	read_stoch(filess);
	cummalitve_prob();

	//List number of client attended on each constraint
	link_dem_con();
	
	/*size_t nnn = 100000;
	generate_instances(nnn);
	string testname = "capplan_large_100000.txt";
	write_st(testname);*/

	const char* flr = stoch_rnd.c_str();
	ifstream stochfile(flr);
	read_st(stochfile);

	//compute expected demands
	compute_expec_dem();

	//initialize partition
	partition.push_back(vector<size_t>(create_partition(nScenarios)));
	part_prob = probability;
	//initialize values of obj sps
	obj_opt_sps.resize(nScenarios, 0.0);
	x_bar.resize(first_st_var.size(), 0.0);
	nFacilities = first_st_var.size();
	GenOrgDest();
	fixed_costs.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		fixed_costs[i] = first_st_coefs[i][0];
	}
	//instatiate subproblem
	//SPProblemCreation();
	//Worst case scenario demand
	//max_dem = max_demand(stoch_param);

	if (gapm_type != 's') {
		part_prob.resize(1);
		part_prob[0] = 1.0;
		sp_scenarios_full.resize(1);
		sp_scenarios_full[0].resize(stoch_constr.size());
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			sp_scenarios_full[0][j] = expected_dem[j];
		}
	}
	
}

vector<double> CP_GAPM::expected_demand(vector<size_t> &element) {
	vector<double> agg_demand_elem(stoch_constr.size(), 0.0);
	double partProb = 0.0;
	for (size_t s = 0; s < element.size(); s++) {
		partProb += probability[s];
	}
	//Computing average demands
	for (size_t j = 0; j < stoch_constr.size(); j++) {
		for (size_t s = 0; s < element.size(); s++) {
			agg_demand_elem[j] += stoch_param[j][element[s]] * probability[s] / partProb;
		}
	}
	return agg_demand_elem;
}

IloModel CP_GAPM::MasterProblemCreation(const char &algo, bool removeobjs)
{
	size_t total_elements;
	if (gapm_type == 's') {
		total_elements = partition.size();
	}
	else {
		total_elements = sp_scenarios_full.size();
	}
	
	//Master env
	IloEnv CP = IloEnv();

	//Master model
	IloModel master = IloModel(CP);
	try {
		// Creating first stage variables
		// Facilities to be open
		//Populate OF
		master_entities.objective = IloExpr(CP);
		master_entities.x = IloNumVarArray(CP, first_st_var.size());
		for (size_t i = 0; i < first_st_var.size(); i++) {
			master_entities.x[i] = IloNumVar(CP, 0.0, DBL_MAX, ILOFLOAT);
			master.add(master_entities.x[i]);
			master_entities.objective += first_st_coefs[i][0] * master_entities.x[i];
		}
		


		// Second stage decision variables (distribution of clients)
		master_entities.y = IloNumVarArray2(CP, second_st_var.size());
		for (size_t i = 0; i < second_st_var.size(); i++) {
			//Completing second stage variables
			master_entities.y[i] = IloNumVarArray(CP, total_elements);
			for (size_t s = 0; s < total_elements; s++) {
				master_entities.y[i][s] = IloNumVar(CP, 0.0, DBL_MAX, ILOFLOAT);
				master.add(master_entities.y[i][s]);
				double elementProb = part_prob[s];
				master_entities.objective += (elementProb * second_st_coefs[i][0] * master_entities.y[i][s]);
			}
		}
		
		//Add objective function
		master.add(IloMinimize(CP, master_entities.objective));

		//Create the arrays corresponding to the constraints that will be added
		master_entities.ax_b = IloRangeArray(CP);
		master_entities.linking_constraints = IloRangeArray(CP);


		//Adding the first stage constraints
		for (size_t i = 1; i < first_st_const.size(); i++) {
			IloExpr constr(CP);
			for (size_t j = 0; j < first_st_var.size(); j++) {
				constr += first_st_coefs[j][i] * master_entities.x[j];
			}
			if (first_st_constsense[i] == "E") {
				master_entities.ax_b.add(constr == first_st_rhs[i]);
			}
			else if (first_st_constsense[i] == "L") {
				master_entities.ax_b.add(constr <= first_st_rhs[i]);
			}
			else if (first_st_constsense[i] == "G") {
				master_entities.ax_b.add(constr >= first_st_rhs[i]);
			}
			constr.end();
		}
		master.add(master_entities.ax_b);

		//Linking constraints
		for (size_t s = 0; s < total_elements; s++) {
			size_t counter = 0;
			vector<double> elem_demand;
			if (gapm_type == 's') {
				elem_demand = expected_demand(partition[s]);
			}
			else {
				elem_demand = sp_scenarios_full[s];
			}
			
			for (size_t i = first_st_const.size(); i < second_st_coefs[0].size(); i++) {
				IloExpr constr(CP);
				for (size_t j = 0; j < first_st_var.size(); j++) {
					constr += first_st_coefs[j][i] * master_entities.x[j];
				}
				for (size_t j = 0; j < second_st_var.size(); j++) {
					constr += second_st_coefs[j][i] * master_entities.y[j][s];
				}
				if (counter < second_st_const.size() - stoch_constr.size()) {
					if (second_st_constsense[i - first_st_const.size()] == "E") {
						master_entities.linking_constraints.add(constr == second_st_rhs[i - first_st_const.size()]);
					}
					else if (second_st_constsense[i - first_st_const.size()] == "L") {
						master_entities.linking_constraints.add(constr <= second_st_rhs[i - first_st_const.size()]);
					}
					else if (second_st_constsense[i - first_st_const.size()] == "G") {
						master_entities.linking_constraints.add(constr >= second_st_rhs[i - first_st_const.size()]);
					}
				}
				else {
					if (second_st_constsense[i - first_st_const.size()] == "E") {
						master_entities.linking_constraints.add(constr == elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
					}
					else if (second_st_constsense[i - first_st_const.size()] == "L") {
						master_entities.linking_constraints.add(constr <= elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
					}
					else if (second_st_constsense[i - first_st_const.size()] == "G") {
						master_entities.linking_constraints.add(constr >= elem_demand[counter - (second_st_const.size() - stoch_constr.size())]);
					}
				}
				counter += 1;
				constr.end();
			}
		}
		master.add(master_entities.linking_constraints);
	}
	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception" << endl;
	}
	return master;

}


void CP_GAPM::MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL)
{
	//Setting up cplex
	cplex_master->setParam(IloCplex::Param::Threads, 1);
	cplex_master->setParam(IloCplex::Param::TimeLimit, TL);
	cplex_master->setParam(IloCplex::Param::Benders::Strategy, 3);

	//Solving
	cplex_master->extract(master);
	cplex_master->exportModel("master_capacityplanning.lp");
	cplex_master->solve();

	// Solution recovery
	x_bar.resize(first_st_var.size());
	for (size_t i = 0; i < first_st_var.size(); i++) {
		x_bar[i] = cplex_master->getValue(master_entities.x[i]);
		/*if (x_bar[i] > 0) {
			cout << "Plant " << first_st_var[i] << " has capacity " << x_bar[i] << " assigned." << endl;
		}*/
	}

	size_t total_elements;
	if (gapm_type == 's') {
		total_elements = partition.size();
	}
	else {
		total_elements = sp_scenarios_full.size();
	}
	
	y_bar.resize(second_st_var.size());
	for (size_t i = 0; i < second_st_var.size(); i++) {
		y_bar[i].resize(total_elements);
		for (size_t s = 0; s < total_elements; s++) {
			y_bar[i][s] = cplex_master->getValue(master_entities.y[i][s]);
		}
	}
	IloNumArray lowerlim(cplex_master->getEnv(), master_entities.linking_constraints.getSize());
	IloNumArray upperlim(cplex_master->getEnv(), master_entities.linking_constraints.getSize());
	cplex_master->getRHSSA(lowerlim, upperlim, master_entities.linking_constraints);
	/*for (size_t s = 0; s < master_entities.linking_constraints.getSize(); s++) {
		cout << " Sensitivity analysis constraint " << s+1 << "[" << lowerlim[s] << " - " << upperlim[s] << "] " << cplex_master->getDual(master_entities.linking_constraints[s]) << endl;
	}*/
	LB = cplex_master->getObjValue();
	LookCompCon(0);
	CapacityComCon(0);
	FlowsComCon(0);
	CheckCasesComcon(0);
	//TODO: Compare flows and decide according cases
	//Check latex file from Eduardo to split according the case
	//GenSubpartitions();
	//GenExpectedScen();
	//FinalScenariosSubparts();
	double tot = 0.0;
	for (size_t s = 0; s < part_prob.size(); s++) {
		tot += part_prob[s];
	}
	LB = cplex_master->getObjValue();

}

void CP_GAPM::SPProblemCreation_GRB() {
	try {
		//Create the variables
		sp_entities_grb.y = new GRBVar[second_st_var.size()];
		for (size_t j = 0; j < second_st_var.size(); j++) {
			sp_entities_grb.y[j] = subprob_grb.addVar(0.0, GRB_INFINITY, second_st_coefs[j][0],
				GRB_CONTINUOUS, second_st_var[j]);
		}

		size_t counter = 0;
		for (size_t i = first_st_const.size(); i < second_st_coefs[0].size(); i++) {
			GRBLinExpr constr;
			for (size_t j = 0; j < second_st_var.size(); j++) {
				constr += second_st_coefs[j][i] * sp_entities_grb.y[j];
			}
			double rhs_byx = 0.0;
			for (size_t j = 0; j < first_st_var.size(); j++) {
				rhs_byx += first_st_coefs[j][i] * x_bar[j];
			}
			rhs_byx = rhs_byx * -1.0;
			if (counter < second_st_const.size() - stoch_constr.size()) {
				if (second_st_constsense[i - first_st_const.size()] == "E") {
					subprob_grb.addConstr(constr == second_st_rhs[i - first_st_const.size()] + rhs_byx, second_st_const[counter]);
				}
				else if (second_st_constsense[i - first_st_const.size()] == "L") {
					subprob_grb.addConstr(constr <= second_st_rhs[i - first_st_const.size()] + rhs_byx, second_st_const[counter]);
				}
				else if (second_st_constsense[i - first_st_const.size()] == "G") {
					subprob_grb.addConstr(constr >= second_st_rhs[i - first_st_const.size()] + rhs_byx, second_st_const[counter]);
				}
			}
			else {
				if (second_st_constsense[i - first_st_const.size()] == "E") {
					subprob_grb.addConstr(constr == stoch_param[counter - (second_st_const.size() - stoch_constr.size())][0] + rhs_byx, second_st_const[counter]);
				}
				else if (second_st_constsense[i - first_st_const.size()] == "L") {
					subprob_grb.addConstr(constr <= stoch_param[counter - (second_st_const.size() - stoch_constr.size())][0] + rhs_byx, second_st_const[counter]);
				}
				else if (second_st_constsense[i - first_st_const.size()] == "G") {
					subprob_grb.addConstr(constr >= stoch_param[counter - (second_st_const.size() - stoch_constr.size())][0] + rhs_byx, second_st_const[counter]);
				}
			}
			counter += 1;
			constr.clear();
		}
		subprob_grb.update();
		subprob_grb.write("subprob_capplan.lp");
	}
	catch (GRBException e) {
		cerr << "Error number: " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error GRB" << endl;
	}
}


void CP_GAPM::SPProblemModification_GRB(vector<size_t> &element, bool mod_x) {
	// s is used to modify RHS
	//modifying demands
	//TODO: compute expected demand
	vector<double> elem_demand;
	if (gapm_type == 's') {
		elem_demand = expected_demand(element);
	}
	else {
		elem_demand = sp_scenarios_full[0];
	}
	for (size_t j = 0; j < stoch_constr.size(); j++) {
		GRBConstr dem_con = subprob_grb.getConstrByName(stoch_constr[j]);
		//Modify the value of RHS
		dem_con.set(GRB_DoubleAttr_RHS, elem_demand[j]);
	}
	size_t counter = 0;
	for (size_t i = first_st_const.size(); i < second_st_coefs[0].size(); i++) {
		double rhs_byx = 0.0;
		for (size_t j = 0; j < first_st_var.size(); j++) {
			rhs_byx += first_st_coefs[j][i] * x_bar[j];
		}
		rhs_byx = rhs_byx * -1.0;
		GRBConstr dem_con = subprob_grb.getConstrByName(second_st_const[i - first_st_const.size()]);
		if (counter < second_st_const.size() - stoch_constr.size()) {
			dem_con.set(GRB_DoubleAttr_RHS, second_st_rhs[i - first_st_const.size()] + rhs_byx);
		}
		else {
			dem_con.set(GRB_DoubleAttr_RHS, elem_demand[counter - (second_st_const.size() - stoch_constr.size())] + rhs_byx);
		}
		counter += 1;
	}
	//Updating the model
	subprob_grb.update();
	subprob_grb.write("subprob_capplan.lp");
}

//Solving the subproblem via Gurobi
void CP_GAPM::SPProbleSolution_GRB(vector<double> &stoch, solution_sps *sp_info, bool benders) {
	//Setting up the configuration of gurobi
	subprob_grb.set(GRB_IntParam_Method, 0);
	subprob_grb.set(GRB_IntParam_Presolve, 0);
	subprob_grb.set(GRB_IntParam_OutputFlag, 0);
	subprob_grb.set(GRB_IntParam_InfUnbdInfo, 1);
	subprob_grb.optimize();

	//stoch
	stoch.clear();
	sp_info->lambda.clear();

	//Status of the solution
	sp_info->F = subprob_grb.get(GRB_IntAttr_Status);
	//If optimal, get duals
	if (sp_info->F == 2) {
		//cout << "Dual multipliers of demand constr ";
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			vector<double> elem_demand = expected_demand(sp_info->scen);
			stoch.push_back(elem_demand[j]);
			sp_info->lambda.push_back(subprob_grb.getConstrByName(stoch_constr[j].c_str()).get(GRB_DoubleAttr_Pi));
			//cout << "min " << subprob_grb.getConstrByName(stoch_constr[j].c_str()).get(GRB_DoubleAttr_SARHSLow) << " max " << subprob_grb.getConstrByName(stoch_constr[j].c_str()).get(GRB_DoubleAttr_SARHSUp) << endl;
			//cout << sp_info->lambda.back() << " ";
		}
		//cout << endl;
		if (benders) {
			//cout << "Dual multipliers of cap constr: ";
			for (size_t i = 0; i < second_st_const.size() - stoch_constr.size(); i++) {
				sp_info->lambda.push_back(subprob_grb.getConstrByName(second_st_const[i].c_str()).get(GRB_DoubleAttr_Pi));
				//cout << sp_info->lambda.back() << " ";
			}
			//cout << endl;
		}
		sp_info->obj = subprob_grb.get(GRB_DoubleAttr_ObjVal);
	}
	//If infeasible, get extreme ray
	else if (sp_info->F == 3) {
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			vector<double> elem_demand = expected_demand(sp_info->scen);
			stoch.push_back(elem_demand[j]);
			sp_info->lambda.push_back(-1.0 * subprob_grb.getConstrByName(stoch_constr[j].c_str()).get(GRB_DoubleAttr_FarkasDual));
		}
		sp_info->obj = DBL_MAX;
	}

	/*y_bar.resize(second_st_var.size());
	for (size_t i = 0; i < second_st_var.size(); i++) {
		y_bar[i] = sp_entities_grb.y[i].get(GRB_DoubleAttr_X);
		if (y_bar[i] > 0.0) {
			cout << "Var " << second_st_var[i] << " has value " << y_bar[i] << endl;
		}
	}*/
}

//Function to build the conex component of a given solution
void CP_GAPM::LookCompCon(const size_t &s) {
	//obtain all origins and destinations served from plants for each element into the partition
	dest_sol.clear();
	ori_sol.clear();
	for (size_t i = 0; i < second_st_var.size(); i++) {
		if (y_bar[i][s] > 0.0) {
			ori_sol.push_back(origins[i]);
			dest_sol.push_back(destinations[i]);
		}
	}
	//Generate the conex component
	ori_comcon_sol.clear();
	dest_comcon_sol.clear();
	while (ori_sol.size() > 0) {
		vector<size_t> to_check_or;
		vector<size_t> to_check_de;
		to_check_or.push_back(ori_sol[0]);
		to_check_de.push_back(dest_sol[0]);
		ori_comcon_sol.push_back({ ori_sol[0] });
		dest_comcon_sol.push_back({ dest_sol[0] });
		//delete elements from original vectors
		vector<size_t>::iterator it1;
		vector<size_t>::iterator it2;
		it1 = ori_sol.begin();
		it2 = dest_sol.begin();
		ori_sol.erase(it1);
		dest_sol.erase(it2);
		while (to_check_or.size() > 0 || to_check_de.size() > 0) {
			vector<size_t> matchesDelete;
			if (to_check_de.size() > 0) {
				for (size_t s = 0; s < dest_sol.size(); s++) {
					if (dest_sol[s] == to_check_de[0]) {
						to_check_or.push_back(ori_sol[s]);
						matchesDelete.push_back(s);
						//add arc to the current tree
						ori_comcon_sol.back().push_back(ori_sol[s]);
						dest_comcon_sol.back().push_back(dest_sol[s]);
					}
				}
				DeleteAll(ori_sol, matchesDelete);
				DeleteAll(dest_sol, matchesDelete);
				matchesDelete.clear();
				vector<size_t>::iterator it;
				it = to_check_de.begin();
				to_check_de.erase(it);
			}

			if (to_check_or.size() > 0) {
				for (size_t s = 0; s < ori_sol.size(); s++) {
					if (ori_sol[s] == to_check_or[0]) {
						to_check_de.push_back(ori_sol[s]);
						matchesDelete.push_back(s);
						//add arc to the current tree
						ori_comcon_sol.back().push_back(ori_sol[s]);
						dest_comcon_sol.back().push_back(dest_sol[s]);
					}
				}
				DeleteAll(ori_sol, matchesDelete);
				DeleteAll(dest_sol, matchesDelete);
				matchesDelete.clear();
				vector<size_t>::iterator it;
				it = to_check_or.begin();
				to_check_or.erase(it);
			}
		}
	}
}

void CP_GAPM::FlowsComCon(const size_t &s) {
	flows_comcon_sol.clear();
	flows_comcon_sol.resize(ori_comcon_sol.size());
	for (size_t k = 0; k < ori_comcon_sol.size(); k++) {
		flows_comcon_sol[k] = 0.0;
		for (size_t f = 0; f < ori_comcon_sol[k].size(); f++) {
			size_t pos = 0;
			for (size_t a = 0; a < second_st_var.size(); a++) {
				if (origins[a] == ori_comcon_sol[k][f] && destinations[a] == dest_comcon_sol[k][f])
					break;
				pos++;
			}
			flows_comcon_sol[k] += y_bar[pos][s];
		}
	}
}

void CP_GAPM::CapacityComCon(const size_t &s) {
	flows_orig_comcon_sol.clear();
	flows_orig_comcon_sol.resize(ori_comcon_sol.size(), 0.0);
	for (size_t i = 0; i < ori_comcon_sol.size(); i++) {
		vector<size_t> or_nodes;
		or_nodes = ori_comcon_sol[i];
		remove_duplicates(or_nodes);
		for (size_t j = 0; j < or_nodes.size(); j++) {
			vector<size_t>::iterator it = find(terminals.begin(), terminals.end(), or_nodes[j]);
			size_t pos = distance(terminals.begin(), it);
			flows_orig_comcon_sol[i] += x_bar[pos];
		}
	}
}



void CP_GAPM::GenScen(vector<size_t> client, vector<size_t> &sc, const size_t &size_v, const size_t &k, double &new_prob) {
	if (sc.size() == size_v) {
		//Classify the scenario according to the total demand
		double sum_up = 0.0;
		for (size_t a = 0; a < sc.size(); a++) {
			sum_up += sc[a];
		}
		size_t positions;
		if (is_integer(flows_comcon_sol[k])) {
			positions = 3;
		}
		else {
			positions = 2;
		}
		if (positions == 3 && sum_up == flows_comcon_sol[k]) {
			subpartitions[k][1].push_back(sc);
			prob_sp_scen[k][1].push_back(new_prob);
		}
		else if (sum_up < flows_comcon_sol[k]) {
			subpartitions[k][0].push_back(sc);
			prob_sp_scen[k][0].push_back(new_prob);
		}
		else if (positions == 3 && sum_up > flows_comcon_sol[k]) {
			subpartitions[k][2].push_back(sc);
			prob_sp_scen[k][2].push_back(new_prob);
		}
		else if (positions == 2 && sum_up > flows_comcon_sol[k]) {
			subpartitions[k][1].push_back(sc);
			prob_sp_scen[k][1].push_back(new_prob);
		}
	}
	else {
		vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), client[sc.size()]);
		size_t pos = distance(dem_constr.begin(), it);
		for (size_t s = 0; s < indep_rhs_dist[pos].size(); s++) {
			sc.push_back(indep_rhs_dist[pos][s]);
			double pr = new_prob;
			new_prob = new_prob * marginal_prob[pos][s];
			GenScen(client, sc, client.size(), k, new_prob);
			sc.pop_back();
			if (marginal_prob[pos][s] < 1e-6) {
				new_prob = pr;
			}
			else {
				new_prob = new_prob / marginal_prob[pos][s];
			}
		}
	}
}

void CP_GAPM::link_dem_con() {
	dem_constr.resize(stoch_constr.size());
	for (size_t t = 0; t < stoch_constr.size(); t++) {
		string cons = stoch_constr[t].substr(4, stoch_constr[t].size());
		sscanf(cons.c_str(), "%zu", &dem_constr[t]);
	}
}


void CP_GAPM::GenExpectedScen() {
	size_t k = 0;
	part_prob.clear();
	for (size_t s1 = 0; s1 < expected_subpart[k].size(); s1++) {
		vector<double> new_scen(expected_subpart[k][s1]);
		double prob = prob_sp_acum[k][s1];
		k++;
		join_scen(new_scen, k, prob);
		k--;
	}
	for (size_t kk = 0; kk < subpart_clients.size(); kk++) {
		clients_attended.insert(clients_attended.end(), subpart_clients[kk].begin(), subpart_clients[kk].end());
	}
}


void CP_GAPM::join_scen(vector<double> new_scen, size_t k, double prob) {
	if (k == expected_subpart.size()) {
		master_scens_sp.push_back(new_scen);
		part_prob.push_back(prob);
	}
	else {
		for (size_t s1 = 0; s1 < expected_subpart[k].size(); s1++) {
			new_scen.insert(new_scen.end(), expected_subpart[k][s1].begin(), expected_subpart[k][s1].end());
			double prev = prob;
			prob = prob * prob_sp_acum[k][s1];
			k++;
			join_scen(new_scen, k, prob);
			k--;
			if (prob_sp_acum[k][s1] < 1e-8) {
				prob = prev;
			}
			else {
				prob = prob / prob_sp_acum[k][s1];
			}
			new_scen.erase(new_scen.end() - expected_subpart[k][s1].size(), new_scen.end());
		}
	}
}


void CP_GAPM::FinalScenariosSubparts() {
	sp_scenarios_full.resize(part_prob.size());
	for (size_t s = 0; s < part_prob.size(); s++) {
		sp_scenarios_full[s].resize(stoch_constr.size(), 0.0);
		for (size_t j = 0; j < stoch_constr.size(); j++) {
			size_t cl = dem_constr[j];
			vector<size_t>::iterator it = find(clients_attended.begin(), clients_attended.end(), cl);
			if (it == clients_attended.end()) {
				sp_scenarios_full[s][j] = expected_dem[j];
			}
			else {
				size_t index = distance(clients_attended.begin(), it);
				sp_scenarios_full[s][j] = master_scens_sp[s][index];
			}
		}
	}
}


void CP_GAPM::FullCP() {
	//Full problem executon time
	chrono::steady_clock sc_cut;
	auto start_cut = sc_cut.now();
	auto end_cut = sc_cut.now();

	partition.clear();
	for (size_t s = 0; s < nScenarios; s++) {
		partition.push_back({ s });
	}
	part_prob = probability;


	IloModel master = MasterProblemCreation(false, 's');
	IloCplex cplex_master(master.getEnv());
	cplex_master.extract(master);

	double LB;
	MasterProblemSolution(&cplex_master, master, LB, timelimit);

	end_cut = sc_cut.now();
	auto time_span_cut = static_cast<chrono::duration<double>>(end_cut - start_cut);
	double runtime_cut = time_span_cut.count();
	cout << "it takes " << runtime_cut << "to execute the benders algorithm " << endl;

	cout << "CRP: "
		<< cplex_master.getStatus() << " "
		<< runtime_cut << " "
		<< "UB_NA" << " "
		<< LB << " "
		<< 0 << " "
		<< 0 << " " //iterations rather than explored nodes
		<< 0 << " " //No pending nodes to be explored
		<< partition.size() << " ";
}


//Check the cases in the conex component
//demand > capacity
//demand < capacity
//demand = capacity
//This method returns vector of vectors where each tree must be partitioned 
void CP_GAPM::CheckCasesComcon(const size_t &s) {
	//ori_comcon_sol to get the current number of trees in the component
	size_t ntrees = ori_comcon_sol.size();
	//identify case for each tree into the conex component
	for (size_t i = 0; i < ntrees; i++) {
		//case 1 demand > capacity
		if (flows_orig_comcon_sol[i] < flows_comcon_sol[i]) {
			size_t viol_cli;
			//first we need to find the client who is not being fully served
			//for each client check the total demand served
			double total_demand_tree = 0.0;
			for (size_t j = 0; dest_comcon_sol[i].size(); j++) {
				//dem sat is the total dem sent to client
				double dem_sat = 0.0;
				for (size_t a = 0; a < second_st_var.size(); a++) {
					if (destinations[a] == dest_comcon_sol[i][j])
						dem_sat += y_bar[a][s];
				}
				//Now look for the total demand of client in this scenario
				size_t cli = dest_comcon_sol[i][j];
				vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), cli);
				size_t pos = distance(dem_constr.begin(), it);
				double cli_dem = sp_scenarios_full[s][pos];
				if (dem_sat < cli_dem){
					//obtain the index of the client who is not being served completely
					viol_cli = pos;
				}
				total_demand_tree += sp_scenarios_full[s][pos];
			}
			//get total demand and substract the demand of unsatisfied client
			double total_minus_uns = total_demand_tree - sp_scenarios_full[s][viol_cli];
			vector<double> split_at{ total_minus_uns, total_demand_tree };
			// First, total tree demand minus demand of client unsatisfied
			// secondly, the total flow in the tree minus the demand of unsatisfied client
			where_part_tree.push_back(split_at);
		//case demand < capacity
		} else if (flows_orig_comcon_sol[i] > flows_comcon_sol[i]) {
			size_t viol_pl;
			double total_tree_offer = 0.0;
			//We need to find the plant having surplus
			//for each plant check the total product delivered
			//compare it against the x_bar
			for (size_t j = 0; j < ori_comcon_sol.size(); j++) {
				//total delivered by the plant
				double pr_del = 0.0;
				for (size_t a = 0; a < second_st_var.size(); a++) {
					if (origins[a] == ori_comcon_sol[i][j])
						pr_del += y_bar[a][s];
				}
				//Now look for the total capacity of the plant, using x_bar
				size_t pl = ori_comcon_sol[i][j];
				vector<size_t>::iterator it = find(terminals.begin(), terminals.end(), pl);
				size_t pos = distance(terminals.begin(), it);
				double pl_of = x_bar[pos];
				//
				if (pl_of > pr_del) {
					viol_pl = pos;
					break;
				}
				total_tree_offer += x_bar[pos];
			}
			double total_minus_surplus = total_tree_offer - x_bar[viol_pl];
			vector<double> split_at{ total_minus_surplus, total_tree_offer };
			where_part_tree.push_back(split_at);
		} else {
			//TODO: run the separation already made in GenSubpartitions method
			//This is the total demand served in this tree flows_comcon_sol[k])
			//Using this we can separte
			vector<double> split_at{ flows_comcon_sol[i] };
			where_part_tree.push_back(split_at);
		}
	}
}

void CP_GAPM::GenSubpartitions() {
	subpartitions.resize(flows_comcon_sol.size());
	subpart_clients.resize(flows_comcon_sol.size());
	prob_sp_scen.resize(flows_comcon_sol.size());
	prob_sp_acum.resize(flows_comcon_sol.size());
	expected_subpart.resize(flows_comcon_sol.size());

	size_t positions = 0;
	for (size_t k = 0; k < ori_comcon_sol.size(); k++) {
		sort(dest_comcon_sol[k].begin(), dest_comcon_sol[k].end());
		subpart_clients[k] = dest_comcon_sol[k];
		vector<size_t>::iterator ct = unique(subpart_clients[k].begin(), subpart_clients[k].end());
		subpart_clients[k].resize(distance(subpart_clients[k].begin(), ct));

		if (is_integer(flows_comcon_sol[k])) {
			positions = 3;
		}
		else {
			positions = 2;
		}
		subpartitions[k].resize(positions);
		prob_sp_scen[k].resize(positions);
		//expected_subpart[k].resize(positions);

		vector<size_t> new_scen;
		double new_prob;
		vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), subpart_clients[k][0]);
		size_t pos = distance(dem_constr.begin(), it);
		for (size_t s = 0; s < indep_rhs_dist[pos].size(); s++) {
			new_scen.push_back(indep_rhs_dist[pos][s]);
			new_prob = marginal_prob[pos][s];
			GenScen(subpart_clients[k], new_scen, subpart_clients[k].size(), k, new_prob);
			new_scen.pop_back();
		}

		for (size_t s1 = 0; s1 < prob_sp_scen[k].size(); s1++) {
			double prob_part = 0.0;
			vector<double> expected_sp_part(subpartitions[k][s1][0].size(), 0.0);
			for (size_t s2 = 0; s2 < prob_sp_scen[k][s1].size(); s2++) {
				prob_part += prob_sp_scen[k][s1][s2];
				for (size_t s3 = 0; s3 < subpartitions[k][s1][s2].size(); s3++) {
					expected_sp_part[s3] += prob_sp_scen[k][s1][s2] * subpartitions[k][s1][s2][s3];
					if (s2 == prob_sp_scen[k][s1].size() - 1) {
						expected_sp_part[s3] = expected_sp_part[s3] / prob_part;
					}
				}
			}
			expected_subpart[k].push_back(expected_sp_part);
			prob_sp_acum[k].push_back(prob_part);
		}
	}

}