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
	index_plants_clients();
	//instatiate subproblem
	//SPProblemCreation();
	//Worst case scenario demand
	//max_dem = max_demand(stoch_param);

	//gapm_type = 's';
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
			std::string varName = "x("+std::to_string(i)+")";
			master_entities.x[i] = IloNumVar(CP, 0.0, DBL_MAX, ILOFLOAT, varName.c_str());
			master.add(master_entities.x[i]);
			master_entities.objective += first_st_coefs[i][0] * master_entities.x[i];
		}
		


		// Second stage decision variables (distribution of clients)
		master_entities.y = IloNumVarArray2(CP, second_st_var.size());
		for (size_t i = 0; i < second_st_var.size(); i++) {
			//Completing second stage variables
			master_entities.y[i] = IloNumVarArray(CP, total_elements);
			for (size_t s = 0; s < total_elements; s++) {
				std::string varName = "y("+std::to_string(i)+","+std::to_string(s)+")";
				master_entities.y[i][s] = IloNumVar(CP, 0.0, DBL_MAX, ILOFLOAT, varName.c_str());
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
			/*if (y_bar[i][s] > 0.0) {
				cout << y_bar[i][s] << endl;
			}*/
		}
	}
	IloNumArray lowerlim(cplex_master->getEnv(), master_entities.linking_constraints.getSize());
	IloNumArray upperlim(cplex_master->getEnv(), master_entities.linking_constraints.getSize());
	cplex_master->getRHSSA(lowerlim, upperlim, master_entities.linking_constraints);
	/*for (size_t s = 0; s < master_entities.linking_constraints.getSize(); s++) {
		cout << " Sensitivity analysis constraint " << s+1 << "[" << lowerlim[s] << " - " << upperlim[s] << "] " << cplex_master->getDual(master_entities.linking_constraints[s]) << endl;
	}*/
	

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

void CP_GAPM::TotalDemandTree(const size_t &s) {
	total_demand_tree.clear();
	total_demand_tree.resize(ori_comcon_sol.size(), 0.0);
	for (size_t i = 0; i < ori_comcon_sol.size(); i++) {
		for (size_t j = 0; j < dest_comcon_sol[i].size(); j++) {
			size_t cli = dest_comcon_sol[i][j];
			vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), cli);
			size_t pos = distance(dem_constr.begin(), it);
			if (gapm_type == 's') {
				vector<double> elem_demand = expected_demand(partition[s]);
				total_demand_tree[i] = elem_demand[pos];
			}
			else {
				total_demand_tree[i] += sp_scenarios_full[s][pos];
			}
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
	comcon_part_tree.clients.clear();
	comcon_part_tree.ineq_sign.clear();
	comcon_part_tree.demand_breakpoints.clear();
	where_part_tree.clear();
	vector<double> elem_demand = expected_demand(partition[s]);
	//identify case for each tree into the conex component
	for (size_t i = 0; i < ntrees; i++) {
		//case 1 demand > capacity
		//flujo total no puede superar la capacidad
		//TODO:: cambiar flows_comcon_sol por demanda total en el arbol
		vector<size_t> unique_dest_tree = dest_comcon_sol[i];
		sort(unique_dest_tree.begin(), unique_dest_tree.end());
		vector<size_t>::iterator it;
		it = unique(unique_dest_tree.begin(), unique_dest_tree.end());
		unique_dest_tree.resize(distance(unique_dest_tree.begin(), it));
		//unique element of origins
		vector<size_t> unique_orig_tree = ori_comcon_sol[i];
		sort(unique_orig_tree.begin(), unique_orig_tree.end());
		it = unique(unique_orig_tree.begin(), unique_orig_tree.end());
		unique_orig_tree.resize(distance(unique_orig_tree.begin(), it));
		if (total_demand_tree[i] - flows_orig_comcon_sol[i] > tol_difference && unique_dest_tree.size() > 1) {
			size_t viol_cli = -1;
			//first we need to find the client who is not being fully served
			//for each client check the total demand served
			for (size_t j = 0; j< unique_dest_tree.size(); j++) {
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
				double cli_dem = 0.0;
				if (gapm_type == 's') {
					cli_dem = elem_demand[pos];
				}
				else {
					cli_dem += sp_scenarios_full[s][pos];
				}
				if (dem_sat < cli_dem){
					//obtain the index of the client who is not being served completely
					viol_cli = pos;
					break;
				}
			}
			//get total demand and substract the demand of unsatisfied client
			assert(viol_cli > -1);
			double total_minus_uns = 0.0; 
			if (gapm_type == 's') {
				total_minus_uns = total_demand_tree[i] - elem_demand[viol_cli];
			}
			else {
				total_minus_uns = total_demand_tree[i] - sp_scenarios_full[s][viol_cli];
			}
			vector<double> split_at{ total_minus_uns, total_demand_tree[i] };
			// First, total tree demand minus demand of client unsatisfied
			// secondly, the total flow in the tree minus the demand of unsatisfied client
			where_part_tree.push_back(split_at);
		//case demand < capacity
		} else if (flows_orig_comcon_sol[i] - flows_comcon_sol[i] > tol_difference && unique_orig_tree.size()>1) {
			size_t viol_pl;
			double total_tree_offer = 0.0;
			//We need to find the plant having surplus
			//for each plant check the total product delivered
			//compare it against the x_bar
			for (size_t j = 0; j < unique_orig_tree.size(); j++) {
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
	comcon_part_tree.clients = dest_comcon_sol;
	comcon_part_tree.demand_breakpoints = where_part_tree;
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

//We take the old region and for each 'scenario' (subpartition)
//a new set of subregions is generated
void CP_GAPM::update_subregions(const size_t &iteration, const size_t &s) {
	//if the old subregions is empty we are at the first iteration
	//In that case, the coming comcon_part_tree is used to build the current subregions
	if (iteration <= 1) { //modify this condition to first iteration
		//call method to create first subregions based on comcon_part_tree
		//reset control_subregions_new to empty
		control_subregions_new.pr.clear();
		initialize_subregions_it(control_subregions_new);
	}
	else {
		//set the new as the old and then modify the new one when s=0
		if (s == 0) {
			control_subregions_old.pr.clear();
			control_subregions_old.pr = control_subregions_new.pr;
			control_subregions_new.pr.clear();
		}
		whole_partition new_control;
		update_subregions_it(s, new_control, control_subregions_old);
	}
}

void CP_GAPM::update_subregions_eigen(const size_t &iteration, const size_t &s) {
	//if the old subregions is empty we are at the first iteration
	//In that case, the coming comcon_part_tree is used to build the current subregions
	if (iteration <= 1) {
		//call method to create first subregions based on comcon_part_tree
		//reset control_subregions_new to empty
		matrix_partition_new.pr.clear();
		initialize_subregions_eigen(matrix_partition_new);
	}
	else {
		//set the new as the old and then modify the new one when s=0
		if (s == 0) {
			matrix_partition_old.pr.clear();
			matrix_partition_old.pr = matrix_partition_new.pr;
			matrix_partition_new.pr.clear();
		}
		eigen_partition new_control;
		update_subregions_it_eigen(s, new_control, matrix_partition_old);
	}
}

void CP_GAPM::update_subregions_it(const size_t &s, whole_partition &new_control, whole_partition &control_tobe_updated) {
	//First we run the initialization with s and new_control
	initialize_subregions_it(new_control);
	//now with s and new_control we put alltogether in the new control
	for (size_t p = 0; p < new_control.pr.size(); p++) {
		//create a new entity to store info of the new subregion
		part_characterisation new_sr;
		new_sr.clients = control_tobe_updated.pr[s].clients;
		new_sr.ineq_sign = control_tobe_updated.pr[s].ineq_sign;
		new_sr.demand_breakpoints = control_tobe_updated.pr[s].demand_breakpoints;

		//control_subregions_new.pr.resize(control_subregions_new.pr.size() + 1);
		//plug the older information
		//control_subregions_new.pr[control_subregions_new.pr.size() - 1].clients = control_tobe_updated.pr[s].clients;
		//control_subregions_new.pr[control_subregions_new.pr.size() - 1].ineq_sign = control_tobe_updated.pr[s].ineq_sign;
		//control_subregions_new.pr[control_subregions_new.pr.size() - 1].demand_breakpoints = control_tobe_updated.pr[s].demand_breakpoints;
		//push back all the new information
		for (size_t k = 0; k < new_control.pr[p].clients.size(); k++) {
			new_sr.clients.push_back(new_control.pr[p].clients[k]);
			new_sr.ineq_sign.push_back(new_control.pr[p].ineq_sign[k]);
			new_sr.demand_breakpoints.push_back(new_control.pr[p].demand_breakpoints[k]);
			//control_subregions_new.pr[control_subregions_new.pr.size() - 1].clients.push_back(new_control.pr[p].clients[k]);
			//control_subregions_new.pr[control_subregions_new.pr.size() - 1].ineq_sign.push_back(new_control.pr[p].ineq_sign[k]);
			//control_subregions_new.pr[control_subregions_new.pr.size() - 1].demand_breakpoints.push_back(new_control.pr[p].demand_breakpoints[k]);
		}
		//check if new_sr already exists in control_subregions_new.pr
		if (s > 0) {
			bool addnewsr = sr_exists(new_sr);
			if (!addnewsr) {
				control_subregions_new.add_subregion_store(new_sr.ineq_sign, new_sr.demand_breakpoints, new_sr.clients);
			}
		}
		else {
			control_subregions_new.add_subregion_store(new_sr.ineq_sign, new_sr.demand_breakpoints, new_sr.clients);
		}
		
	}
	
}


void CP_GAPM::update_subregions_it_eigen(const size_t &s, eigen_partition &new_control, eigen_partition &control_tobe_updated) {
	//First we run the initialization with s and new_control
	initialize_subregions_eigen(new_control);
	//now with s and new_control we put alltogether in the new control
	for (size_t p = 0; p < new_control.pr.size(); p++) {
		//Append new control with new elements got from new subregions generated
		eigen_characterisation new_sr;
		new_sr.subpart_matrix = control_tobe_updated.pr[s].subpart_matrix;
		new_sr.rhs_checkpart = control_tobe_updated.pr[s].rhs_checkpart;

		size_t cur_nr = new_sr.subpart_matrix.rows();
		new_sr.subpart_matrix.conservativeResize(cur_nr + new_control.pr[p].subpart_matrix.rows(), new_control.pr[p].subpart_matrix.cols());
		size_t cur_nr_rhs = new_sr.rhs_checkpart.rows();
		new_sr.rhs_checkpart.conservativeResize(cur_nr_rhs + new_control.pr[p].rhs_checkpart.rows(), 1);

		//cout << new_control.pr[p].subpart_matrix << endl;
		
		/*for (Eigen::Index c = 0; c < new_sr.subpart_matrix.cols(); ++c) {
			new_sr.subpart_matrix.startVec(c);
			for (Eigen::SparseMatrix<double>::InnerIterator itC(new_control.pr[p].subpart_matrix, c); itC; ++itC)
				new_sr.subpart_matrix.insertBack(itC.row() + cur_nr, c) = itC.value();
		}*/		
		std::vector<Eigen::Triplet<double> > tripletList;
		for (int k = 0; k < new_sr.subpart_matrix.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(new_sr.subpart_matrix, k); it; ++it)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
			}
		}
		for (int k = 0; k < new_control.pr[p].subpart_matrix.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(new_control.pr[p].subpart_matrix, k); it; ++it)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(it.row() + cur_nr, it.col(), it.value()));
			}
		}
		new_sr.subpart_matrix.resize(cur_nr + new_control.pr[p].subpart_matrix.rows(), new_control.pr[p].subpart_matrix.cols());
		new_sr.subpart_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

		//cout << new_sr.subpart_matrix << endl;
		for (size_t k = 0; k < new_control.pr[p].rhs_checkpart.rows(); k++) {
			new_sr.rhs_checkpart(cur_nr_rhs + k, 0) = new_control.pr[p].rhs_checkpart(k, 0);
		}

		//cout << new_sr.rhs_checkpart << endl;
		matrix_partition_new.add_subregion_store(new_sr.rhs_checkpart, new_sr.subpart_matrix);
	}
}

bool CP_GAPM::sr_exists(part_characterisation &new_sr) {
	for (size_t p = 0; p < control_subregions_new.pr.size(); p++) {
		bool vector_found = false;
		//variable to count number of element in the subregion that are found in the current array of subregions
		size_t cont = 0;
		size_t k_new = 0;
		vector<vector<size_t>> copy_clients = control_subregions_new.pr[p].clients;
		while (k_new < new_sr.clients.size()) {
		find_new_element:
			if (k_new >= new_sr.clients.size()) {
				break;
			}
			vector<size_t> tree_cl = new_sr.clients[k_new];
			for (size_t k = 0; k < copy_clients.size(); k++) {
				bool v_equals = false;
				if (copy_clients[k].size() == tree_cl.size())
					v_equals = compareVectors(copy_clients[k], tree_cl);
				if (v_equals) {
					cont += 1;
					k_new += 1;
					vector_found = true;
					copy_clients.erase(copy_clients.begin() + k);
					goto find_new_element;
				}
			}
			k_new += 1;
		}

		//If it gets here, check inequality signs
		if (cont == new_sr.clients.size()) {
			cont = 0;
			k_new = 0;
			vector<vector<string>> copy_ineqsign = control_subregions_new.pr[p].ineq_sign;
			while (k_new < new_sr.ineq_sign.size()) {
			find_new_el:
				if (k_new >= new_sr.ineq_sign.size()) {
					break;
				}
				vector<string> tree_insig = new_sr.ineq_sign[k_new];
				for (size_t k = 0; k < copy_ineqsign.size(); k++) {
					bool v_equals = false;
					if (copy_ineqsign.size() == tree_insig.size())
						v_equals = compareVectors(copy_ineqsign[k], tree_insig);
					if (v_equals) {
						cont += 1;
						k_new += 1;
						vector_found = true;
						copy_ineqsign.erase(copy_ineqsign.begin() + k);
						goto find_new_el;
					}
				}
				k_new += 1;
			}
		}

		//If it gets here, check demand_breakpoints
		if (cont == new_sr.ineq_sign.size()) {
			cont = 0;
			k_new = 0;
			vector<vector<double>> copy_brkpnt = control_subregions_new.pr[p].demand_breakpoints;
			while (k_new < new_sr.demand_breakpoints.size()) {
			find_new_d:
				if (k_new >= new_sr.demand_breakpoints.size()) {
					if (cont == new_sr.demand_breakpoints.size()) {
						return true;
					}
					break;
				}
				vector<double> tree_brk = new_sr.demand_breakpoints[k_new];
				for (size_t k = 0; k < copy_brkpnt.size(); k++) {
					bool v_equals = false;
					if (copy_brkpnt.size() == tree_brk.size())
						v_equals = compareVectors(copy_brkpnt[k], tree_brk);
					if (v_equals) {
						cont += 1;
						k_new += 1;
						vector_found = true;
						copy_brkpnt.erase(copy_brkpnt.begin() + k);
						goto find_new_d;
					}
				}
				k_new += 1;
			}
		}
	}
	return false;
}

void CP_GAPM::initialize_subregions_it(whole_partition &control_tobe_updated) {
	//using comcon_part_tree we create an initial set of subregions
	for (size_t k = 0; k < comcon_part_tree.demand_breakpoints[0].size() + 1; k++) {
		//create a vector to store the inequality sign of the current subpartition that will be created
		vector<vector<string>> ineqs;
		vector<vector<double>> break_values;
		if (comcon_part_tree.demand_breakpoints[0].size() >= 1 && k == 0) {
			ineqs.push_back({ "leq" });
			break_values.push_back({ comcon_part_tree.demand_breakpoints[0][0] });
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 1 && k == 1) {
			ineqs.push_back({ "g" });
			break_values.push_back({ comcon_part_tree.demand_breakpoints[0][0] });
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 2 && k == 1) {
			ineqs.push_back({ "g", "leq" });
			break_values.push_back({ comcon_part_tree.demand_breakpoints[0][0], comcon_part_tree.demand_breakpoints[0][1] });
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 2 && k == 2) {
			ineqs.push_back({ "g" });
			break_values.push_back({ comcon_part_tree.demand_breakpoints[0][1] });
		}
		size_t counting = 1;
		recursive_gen_subregions(ineqs, break_values, counting, control_tobe_updated);
	}
}

void CP_GAPM::initialize_subregions_eigen(eigen_partition &control_tobe_updated) {
	//count the number of client present in the conex component
	size_t tc = 0;
	for (size_t k = 0; k < comcon_part_tree.clients.size(); k++) {
		tc += comcon_part_tree.clients[k].size();
		if (comcon_part_tree.demand_breakpoints[k].size() > 1) {
			tc += comcon_part_tree.clients[k].size();
		}
	}
	for (size_t k = 0; k < comcon_part_tree.demand_breakpoints[0].size() + 1; k++) {
		Eigen::SparseMatrix<double> subpart_matrix(1, dem_constr.size());
		subpart_matrix.reserve(tc);
		Eigen::MatrixXd rhs_checkpart;
		double inse;
		double new_rows_matrix = 0;
		size_t counting = 0;
		if (comcon_part_tree.demand_breakpoints[0].size() >= 1 && k == 0) {
			//Positive coefficients in matrix
			inse = 1.0;
			new_rows_matrix = 1;
			//Resize matrix (not necessary here because it was defined with 1 row)
			//Resize rhs matrix
			size_t nr_rhs = rhs_checkpart.rows();
			rhs_checkpart.conservativeResize(nr_rhs + 1, 1);
			//Populate rhs
			rhs_checkpart(nr_rhs, 0) = comcon_part_tree.demand_breakpoints[0][0];
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 1 && k == 1) {
			//Negative coefficients in matrix
			inse = -1.0;
			new_rows_matrix = 1;
			//Resize matrix (not necessary here because it was defined with 1 row)
			//Resize rhs matrix
			size_t nr_rhs = rhs_checkpart.rows();
			rhs_checkpart.conservativeResize(nr_rhs + 1, 1);
			//Populate rhs
			rhs_checkpart(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[0][0];
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 2 && k == 1) {
			//Positive and Negative coefficients in matrix
			inse = -1.0;
			new_rows_matrix = 2;
			//Resize matrix (Only 1 new row because we already have 1)
			size_t nr = subpart_matrix.rows();
			subpart_matrix.conservativeResize(nr + 1, dem_constr.size());
			//Resize rhs matrix
			size_t nr_rhs = rhs_checkpart.rows();
			rhs_checkpart.conservativeResize(nr_rhs + 2, 1);
			//Populate rhs
			rhs_checkpart(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[0][0];
			rhs_checkpart(nr_rhs + 1, 0) = comcon_part_tree.demand_breakpoints[0][1];
		}
		else if (comcon_part_tree.demand_breakpoints[0].size() == 2 && k == 2) {
			//Negative coefficients in matrix
			inse = -1.0;
			new_rows_matrix = 1;
			//Resize matrix (not necessary here because it was defined with 1 row)
			//Resize rhs matrix
			size_t nr_rhs = rhs_checkpart.rows();
			rhs_checkpart.conservativeResize(nr_rhs + 1, 1);
			//Populate rhs
			rhs_checkpart(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[0][1];
		}
		//Populate the matrix
		size_t nr = subpart_matrix.rows();
		if (new_rows_matrix == 1) {
			for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
				subpart_matrix.insert(nr - 1, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = inse;
			}
		} else if (new_rows_matrix == 2) {
			for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
				subpart_matrix.insert(nr - 2, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = inse;
			}
			for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
				subpart_matrix.insert(nr - 1, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = -1.0*inse;
			}
		}
		counting++;
		recursive_gen_subregions_eigen(rhs_checkpart, subpart_matrix, counting, control_tobe_updated);
	}
	
}


//recursively creating all possible combinations of subregions
void CP_GAPM::recursive_gen_subregions(vector<vector<string>> &ineqs, vector<vector<double>> &break_values, size_t &counting, whole_partition &control_tobe_updated) {
	//add elements until the last position of comcon_part_tree is reached
	if (counting < comcon_part_tree.clients.size()) {
		for (size_t k = 0; k < comcon_part_tree.demand_breakpoints[counting].size() + 1; k++) {
			if (comcon_part_tree.demand_breakpoints[counting].size() >= 1 && k == 0) {
				ineqs.push_back({ "leq" });
				break_values.push_back({ comcon_part_tree.demand_breakpoints[counting][0] });
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 1 && k == 1) {
				ineqs.push_back({ "g" });
				break_values.push_back({ comcon_part_tree.demand_breakpoints[counting][0] });
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 2 && k == 1) {
				ineqs.push_back({ "g", "leq" });
				break_values.push_back({ comcon_part_tree.demand_breakpoints[counting][0], comcon_part_tree.demand_breakpoints[counting][1] });
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 2 && k == 2) {
				ineqs.push_back({ "g" });
				break_values.push_back({ comcon_part_tree.demand_breakpoints[counting][1] });
			}
			counting += 1;
			recursive_gen_subregions(ineqs, break_values, counting, control_tobe_updated);
			ineqs.pop_back();
			break_values.pop_back();
			counting -= 1;
		}
	}
	else {
		//a method to push back a new subpartition
		control_tobe_updated.add_subregion_store(ineqs, break_values, comcon_part_tree.clients);
	}
}

void CP_GAPM::recursive_gen_subregions_eigen(Eigen::MatrixXd &break_values, Eigen::SparseMatrix<double> &subpart_matrix, size_t &counting, eigen_partition &control_tobe_updated) {
	if (counting < comcon_part_tree.clients.size()) {
		for (size_t k = 0; k < comcon_part_tree.demand_breakpoints[counting].size() + 1; k++) {
			double inse;
			double new_rows_matrix = 0;
			if (comcon_part_tree.demand_breakpoints[counting].size() >= 1 && k == 0) {
				//Positive coefficients in matrix
				inse = 1.0;
				new_rows_matrix = 1;
				//Resize matrix (1 new position)
				size_t nr = subpart_matrix.rows();
				subpart_matrix.conservativeResize(nr + 1, dem_constr.size());
				//Resize rhs matrix
				size_t nr_rhs = break_values.rows();
				break_values.conservativeResize(nr_rhs + 1, 1);
				//Populate rhs
				break_values(nr_rhs, 0) = comcon_part_tree.demand_breakpoints[counting][0];
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 1 && k == 1) {
				//Negative coefficients in matrix
				inse = -1.0;
				new_rows_matrix = 1;
				//Resize matrix (1 new position)
				size_t nr = subpart_matrix.rows();
				subpart_matrix.conservativeResize(nr + 1, dem_constr.size());
				//Resize rhs matrix
				size_t nr_rhs = break_values.rows();
				break_values.conservativeResize(nr_rhs + 1, 1);
				//Populate rhs
				break_values(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[counting][0];
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 2 && k == 1) {
				//Positive and Negative coefficients in matrix
				inse = -1.0;
				new_rows_matrix = 2;
				//Resize matrix (2 new positions)
				size_t nr = subpart_matrix.rows();
				subpart_matrix.conservativeResize(nr + 2, dem_constr.size());
				//Resize rhs matrix
				size_t nr_rhs = break_values.rows();
				break_values.conservativeResize(nr_rhs + 2, 1);
				//Populate rhs
				break_values(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[counting][0];
				break_values(nr_rhs + 1, 0) = comcon_part_tree.demand_breakpoints[counting][1];
			}
			else if (comcon_part_tree.demand_breakpoints[counting].size() == 2 && k == 2) {
				//Negative coefficients in matrix
				inse = -1.0;
				new_rows_matrix = 1;
				//Resize matrix (1 new position)
				size_t nr = subpart_matrix.rows();
				subpart_matrix.conservativeResize(nr + 1, dem_constr.size());
				//Resize rhs matrix
				size_t nr_rhs = break_values.rows();
				break_values.conservativeResize(nr_rhs + 1, 1);
				//Populate rhs
				break_values(nr_rhs, 0) = -1.0 * comcon_part_tree.demand_breakpoints[counting][1];
			}
			//cout << break_values << endl;
			//Populate the matrix
			size_t nr = subpart_matrix.rows();
			if (new_rows_matrix == 1) {
				for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
					subpart_matrix.insert(nr - 1, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = inse;
				}
			}
			else if (new_rows_matrix == 2) {
				for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
					subpart_matrix.insert(nr - 2, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = inse;
				}
				for (size_t i = 0; i < comcon_part_tree.clients[counting].size(); i++) {
					subpart_matrix.insert(nr - 1, dem_constr[ordered_dem_constr[comcon_part_tree.clients[counting][i]]]) = -1.0*inse;
				}
			}
			counting++;
			recursive_gen_subregions_eigen(break_values, subpart_matrix, counting, control_tobe_updated);
			subpart_matrix.conservativeResize(subpart_matrix.rows() - new_rows_matrix, dem_constr.size());
			break_values.conservativeResize(break_values.rows() - new_rows_matrix, 1);
			counting -= 1;
		}
	}
	else {
		//a method to push back a new subpartition
		//cout << subpart_matrix << endl;
		control_tobe_updated.add_subregion_store(break_values, subpart_matrix);
	}
}

void CP_GAPM::classify_scenarios() {
	partition.clear();
	partition.resize(control_subregions_new.pr.size());
	for (size_t s = 0; s < nScenarios; s++) {
		//check where should be this scenario added
		//compute scenario demand of each tree into conex component
		for (size_t i = 0; i < control_subregions_new.pr.size(); i++) {
			vector<double> demand_tree (control_subregions_new.pr[i].clients.size(), 0.0);
			for (size_t j = 0; j < control_subregions_new.pr[i].clients.size(); j++) {
				for (size_t k = 0; k < control_subregions_new.pr[i].clients[j].size(); k++) {
					size_t cl = control_subregions_new.pr[i].clients[j][k];
					//Now look for the index of the client to get its demands
					vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), cl);
					size_t pos = distance(dem_constr.begin(), it);
					demand_tree[j] += stoch_param[pos][s];
				}
			}
			//we now have the total demand for each tree for a given scenario s
			//Now we need to check where this scenario should be added
			//we get a bool value from classify_one_scenario to know id the scenario belongs to the current i subregion
			bool belongto = classify_one_scenario(demand_tree, i);
			//add the scenario to the element into the parition
			if (belongto == true) {
				partition[i].push_back(s);
			}
		}
	}
}

bool CP_GAPM::classify_one_scenario(vector<double> &demand_tree, size_t &i) {
	//to totalize the number of conditions satisfied
	size_t ncond_satisfied = 0;
	for (size_t j = 0; j < control_subregions_new.pr[i].ineq_sign.size(); j++) {
		if (control_subregions_new.pr[i].ineq_sign[j].size() == 1) {
			//there is only one condition to be checked as inequality
			if (control_subregions_new.pr[i].ineq_sign[j][0] == "leq") {
				if (! (demand_tree[j] <= control_subregions_new.pr[i].demand_breakpoints[j][0]))
					return false;
			} else if (control_subregions_new.pr[i].ineq_sign[j][0] == "g") {
				if (!(demand_tree[j] > control_subregions_new.pr[i].demand_breakpoints[j][0]) )
					return false;
			}
		}
		else if (control_subregions_new.pr[i].ineq_sign[j].size() == 2) {
			//there are 2 conditions to be checked as inequalities
			if (control_subregions_new.pr[i].ineq_sign[j][0] == "g" && control_subregions_new.pr[i].ineq_sign[j][1] == "leq") {
				if (!((demand_tree[j] > control_subregions_new.pr[i].demand_breakpoints[j][0]) && (demand_tree[j] <= control_subregions_new.pr[i].demand_breakpoints[j][1])))
					return false;
			}
		}
	}
	return true;
}

void CP_GAPM::Inner_GAPM(const char &algo){
	//Initialize control variables
	auto time_span = static_cast<chrono::duration<double>>(inner_gampm_runningtime.end - inner_gampm_runningtime.start);
	double run_time = time_span.count();
	inner_gapm_iterations = 0;
	inner_gap = 1;
	inner_lb_gapm = -DBL_MAX;
	inner_ub_gapm = DBL_MAX;

	//create the master problem
	/*IloModel master = MasterProblemCreation(algo);
	IloCplex cplex_master(master.getEnv());
	cplex_master.extract(master);
	master.end();
	cplex_master.end();*/

	//subproblem creation
	SPProblemCreation_GRB();

	//previous x solution
	//MAYBE useful to evaluate if my solution changes
	vector<double> prev_x(x_bar.size());

	bool cont = true;
	while (cont) {
		//update number of iterations
		inner_gapm_iterations++;
		//create master again
		IloModel master = MasterProblemCreation(algo);
		IloCplex cplex_master(master.getEnv());
		cplex_master.extract(master);
		MasterProblemSolution(&cplex_master, master, inner_lb_gapm, timelimit - run_time);

		vector<vector<size_t>> cp_part = partition;
		//runn all functions to divide the outcome space
		for (size_t s = 0; s < cp_part.size(); s++) {
			LookCompCon(s);
			CapacityComCon(s);
			FlowsComCon(s);
			TotalDemandTree(s);
			CheckCasesComcon(s);
			//update_subregions(inner_gapm_iterations, s);
			update_subregions_eigen(inner_gapm_iterations, s);

			//GenSubpartitions();
			//GenExpectedScen();
			//FinalScenariosSubparts();

			
			
		}

		matrix_partition_new.copy_scenarios(stoch_param);
		partition = matrix_partition_new.scenario_classif();

		//classify_scenarios();

		//Delete empty subregions
		vector<size_t> empty_index;
		for (size_t s = 0; s < partition.size(); s++) {
			if (partition[s].size() == 0) {
				empty_index.push_back(s);
			}
		}
		DeleteAll(partition, empty_index);
		//DeleteAllWP(control_subregions_new, empty_index);
		DeleteAllWP(matrix_partition_new, empty_index);
		inner_lb_gapm = cplex_master.getObjValue();
		part_prob.clear();
		part_prob.resize(partition.size());
		for (size_t s = 0; s < partition.size(); s++) {
			part_prob[s] = partition[s].size() / double(nScenarios);
		}
		
		size_t tot_s = 0;
		vector<size_t> prpr;
		for (size_t s = 0; s < partition.size(); s++) {
			tot_s += partition[s].size();
			for (size_t p = 0; p < partition[s].size(); p++) {
				prpr.push_back(partition[s][p]);
			}
		}
		vector<size_t> unique_par = prpr;
		sort(unique_par.begin(), unique_par.end());
		vector<size_t>::iterator it;
		it = unique(unique_par.begin(), unique_par.end());
		unique_par.resize(distance(unique_par.begin(), it));

		cout << "current lower bound is: " << inner_lb_gapm << endl;
	}
}

void CP_GAPM::index_plants_clients() {
	//resizing vectors
	ordered_terminals.resize(terminals.size());
	ordered_dem_constr.resize(dem_constr.size());

	//find and allocate the index positions
	for (size_t i = 0; i < ordered_terminals.size(); i++) {
		vector<size_t>::iterator it = find(terminals.begin(), terminals.end(), (i + 1));
		size_t pos = distance(terminals.begin(), it);
		ordered_terminals[i] = pos;
	}

	for (size_t i = 0; i < ordered_dem_constr.size(); i++) {
		vector<size_t>::iterator it = find(dem_constr.begin(), dem_constr.end(), (i + 1));
		size_t pos = distance(dem_constr.begin(), it);
		ordered_dem_constr[i] = pos;
	}
}


void whole_partition::add_subregion_store(vector<vector<string>> &ineqs, vector<vector<double>> &break_values, vector<vector<size_t>> &clients) {
	pr.resize(pr.size() + 1);
	//add elements to the last position of the subregions control
	pr[pr.size() - 1].clients = clients;
	pr[pr.size() - 1].ineq_sign = ineqs;
	pr[pr.size() - 1].demand_breakpoints = break_values;
}


void eigen_partition::instantiate_matrix(const size_t &nnz_reserve, const size_t &nineq) {
	pr.resize(pr.size() + 1);

	pr.back().subpart_matrix.reserve(nnz_reserve);
	pr.back().rhs_checkpart.resize(nineq,1);
}

void eigen_partition::add_subregion_store(Eigen::MatrixXd &break_values, Eigen::SparseMatrix<double> &subpart_matrix) {
	pr.resize(pr.size() + 1);

	pr[pr.size() - 1].subpart_matrix = subpart_matrix;
	pr[pr.size() - 1].rhs_checkpart = break_values;
}

void eigen_partition::copy_scenarios(vector<vector<double>> &full_scenarios) {
	size_t nr = full_scenarios.size();
	size_t nc = full_scenarios[0].size();

	cp_scen.resize(nr, nc);

	for (size_t i = 0; i < nr; i++) {
		for (size_t j = 0; j < nc; j++) {
			cp_scen(i, j) = full_scenarios[i][j];
		}
	}

	for (size_t j = 0; j < nc; j++) {
		scenarios.push_back(j);
	}
}

vector<vector<size_t>> eigen_partition::scenario_classif() {
	vector<vector<size_t>> partition;

	for (size_t p = 0; p < pr.size(); p++) {
		//Multiply the subpart matrix by the remaining scenarios and check which of them belong to this subpartition
		//Remove those which belong to this sp from the current setlist of scenarios
		partition.push_back({});

		//vector which will contain the elements that should be removed from current set of scenarios
		vector<size_t> to_rem;

		Eigen::MatrixXd result;

		result = pr[p].subpart_matrix * cp_scen;

		//cout << result << endl;
		
		//Now check each column from result ato check which scenarios should belong to this subregion p
		for (size_t j = 0; j < result.cols(); j++) {
			bool belong = true;
			for (size_t i = 0; i < result.rows(); i++) {
				if (result(i, j) > pr[p].rhs_checkpart(i, 0)) {
					belong = false;
					break;
				}
			}
			if (belong == true) {
				partition[p].push_back(scenarios[j]);
				to_rem.push_back(j);
			}
		}


		while (to_rem.size() > 0) {
			removeColumn(cp_scen, to_rem[0]);
			scenarios.erase(scenarios.begin() + to_rem[0]);
			to_rem.erase(to_rem.begin());
			size_t su = 1;
			subtractScalar(to_rem, su);
		}
	}

	return partition;
}