// Created by CDRPico
// 23/06/2020 01:14

#include"../inc/OuterBendersSFLP.h"
#include <cctype>
#include <algorithm>


ILOLAZYCONSTRAINTCALLBACK4(OLazyOptCuts,
	OuterBendersSFLP&, Instance,
	solFeat&, org,
	bool&, first_out,
	const char&, algo)
{
	//bool reaggregated = false;
	//size_t cont = 0;
//comeback:
	//Callback execution time
	size_t scenarios = 0;
	if (org.algo == 's') {
		scenarios = 1;
	}
	else {
		scenarios = Instance.nScenarios;
	}

	//get current solution
	for (size_t i = 0; i < Instance.nFacilities; i++) {
		Instance.x_bar[i] = getValue(Instance.ent_sflp.x[i]);
	}
	vector<double> current_theta(scenarios);
	for (size_t s = 0; s < scenarios; s++) {
		current_theta[s] = getValue(Instance.ent_sflp.theta[s]);
	}

	//Important elements for dual values
	vector<vector<double>> stoch_agg;
	stoch_agg.resize(Instance.partition.size());
	vector<solution_sps> sp_info_agg;
	sp_info_agg.resize(Instance.partition.size());
	vector<double> lhs;
	lhs.resize(Instance.partition.size());

	//Checking if any cut should added
	bool violated = false;
	double single_lhs = 0.0;
	IloExpr new_single_cut(getEnv());
	double const_part_single = 0.0;

	for (size_t p = 0; p < Instance.partition.size(); p++) {
		sp_info_agg[p].scen = Instance.partition[p];
		if (p == 0) {
			Instance.SPProblemModification_GRB(Instance.partition[p], true);
		}
		else {
			Instance.SPProblemModification_GRB(Instance.partition[p]);
		}
		Instance.SPProbleSolution_GRB(stoch_agg[p], &sp_info_agg[p], true);
		//Compute the current rhs to check if the cut is violated
		double element_prob = 0.0;
		double theta_ac = 0.0;
		if (org.algo != 's') {
			for (size_t s = 0; s < Instance.partition[p].size(); s++) {
				element_prob += Instance.probability[Instance.partition[p][s]];
				theta_ac += current_theta[Instance.partition[p][s]] * Instance.probability[Instance.partition[p][s]];
			}
		}
		else {
			element_prob = Instance.probability[Instance.partition[p][0]];
		}

		//constant of current first stage multiplied by duals of capacity constarints
		double const_fs = 0.0;
		for (size_t i = 0; i < Instance.nFacilities; i++) {
			const_fs += element_prob * Instance.facil_capacities[i] * Instance.x_bar[i] * sp_info_agg[p].lambda[Instance.nClients + i];
		}
		double rhs = theta_ac;
		//compute the constant part of the cut
		lhs[p] = 0.0;
		for (size_t j = 0; j < Instance.nClients; j++) {
			lhs[p] += element_prob * stoch_agg[p][j] * sp_info_agg[p].lambda[j];
		}
		lhs[p] += const_fs;
		if (org.algo != 's') {
			if (smallerWithTolerance(rhs, lhs[p])) {
				violated = true;
				IloExpr new_cut(getEnv());
				vector<double> x_coefs;
				vector<double> theta_coefs;
				for (size_t i = 0; i < Instance.nFacilities; i++) {
					x_coefs.push_back(element_prob * Instance.facil_capacities[i] * sp_info_agg[p].lambda[Instance.nClients + i]);
					new_cut += x_coefs.back() * Instance.ent_sflp.x[i];
				}
				for (size_t s = 0; s < Instance.partition[p].size(); s++) {
					theta_coefs.push_back(Instance.probability[Instance.partition[p][s]]);
					new_cut -= Instance.ent_sflp.theta[Instance.partition[p][s]] * theta_coefs.back();
				}
				double derecho = -1.0 * (lhs[p] - const_fs);
				//cout << new_cut << " <= " << derecho << endl;
				string cut_name = "test_" + to_string(p + 1);
				Instance.cp.emplace_back(&Instance.Mast_Bend, &new_cut, derecho, Instance.partition[p], cut_name, true, x_coefs, theta_coefs);
				add(new_cut <= derecho).end();
				org.optCuts += 1;
				new_cut.end();
			}
		}
		else {
			// do single cut
			single_lhs += lhs[p];
			const_part_single += lhs[p] - const_fs;
			for (size_t i = 0; i < Instance.nFacilities; i++) {
				new_single_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * sp_info_agg[p].lambda[Instance.nClients + i];
			}
		}
	}
	if (smallerWithTolerance(current_theta[0], single_lhs) && org.algo == 's') {
		violated = true;
		new_single_cut -= Instance.ent_sflp.theta[0];
		//cout << new_single_cut << " <= " << const_part_single << endl;
		add(new_single_cut <= -const_part_single).end();
		org.optCuts += 1;
	}
	
}

ILOUSERCUTCALLBACK2(OUserOptCuts,
	OuterBendersSFLP&, Instance,
	solFeat&, org)
{
	if (org.cnode == 0) {
		size_t scenarios = 0;
		if (org.algo == 's') {
			scenarios = 1;
		}
		else {
			scenarios = Instance.nScenarios;
		}

		//get current solution
		for (size_t i = 0; i < Instance.nFacilities; i++) {
			Instance.x_bar[i] = getValue(Instance.ent_sflp.x[i]);
		}
		vector<double> current_theta(scenarios);
		for (size_t s = 0; s < scenarios; s++) {
			current_theta[s] = getValue(Instance.ent_sflp.theta[s]);
		}

		//Important elements for dual values
		vector<vector<double>> stoch_agg;
		stoch_agg.resize(Instance.partition.size());
		vector<solution_sps> sp_info_agg;
		sp_info_agg.resize(Instance.partition.size());
		vector<double> lhs;
		lhs.resize(Instance.partition.size());

		//Checking if any cut should added
		bool violated = false;
		double single_lhs = 0.0;
		IloExpr new_single_cut(getEnv());
		double const_part_single = 0.0;

		for (size_t p = 0; p < Instance.partition.size(); p++) {
			sp_info_agg[p].scen = Instance.partition[p];
			if (p == 0) {
				Instance.SPProblemModification_GRB(Instance.partition[p], true);
			}
			else {
				Instance.SPProblemModification_GRB(Instance.partition[p]);
			}
			Instance.SPProbleSolution_GRB(stoch_agg[p], &sp_info_agg[p], true);
			//Compute the current rhs to check if the cut is violated
			double element_prob = 0.0;
			double theta_ac = 0.0;
			if (org.algo != 's') {
				for (size_t s = 0; s < Instance.partition[p].size(); s++) {
					element_prob += Instance.probability[Instance.partition[p][s]];
					theta_ac += current_theta[Instance.partition[p][s]] * Instance.probability[Instance.partition[p][s]];
				}
			}
			else {
				element_prob = Instance.probability[Instance.partition[p][0]];
			}

			//constant of current first stage multiplied by duals of capacity constarints
			double const_fs = 0.0;
			for (size_t i = 0; i < Instance.nFacilities; i++) {
				const_fs += element_prob * Instance.facil_capacities[i] * Instance.x_bar[i] * sp_info_agg[p].lambda[Instance.nClients + i];
			}
			double rhs = theta_ac;
			//compute the constant part of the cut
			lhs[p] = 0.0;
			for (size_t j = 0; j < Instance.nClients; j++) {
				lhs[p] += element_prob * stoch_agg[p][j] * sp_info_agg[p].lambda[j];
			}
			lhs[p] += const_fs;
			if (org.algo != 's') {
				if (smallerWithTolerance(rhs, lhs[p])) {
					violated = true;
					IloExpr new_cut(getEnv());
					vector<double> x_coefs;
					vector<double> theta_coefs;
					for (size_t i = 0; i < Instance.nFacilities; i++) {
						x_coefs.push_back(element_prob * Instance.facil_capacities[i] * sp_info_agg[p].lambda[Instance.nClients + i]);
						new_cut += x_coefs.back() * Instance.ent_sflp.x[i];
					}
					for (size_t s = 0; s < Instance.partition[p].size(); s++) {
						theta_coefs.push_back(Instance.probability[Instance.partition[p][s]]);
						new_cut -= Instance.ent_sflp.theta[Instance.partition[p][s]] * theta_coefs.back();
					}
					double derecho = -1.0 * (lhs[p] - const_fs);
					//cout << new_cut << " <= " << derecho << endl;
					string cut_name = "test_" + to_string(p + 1);
					Instance.cp.emplace_back(&Instance.Mast_Bend, &new_cut, derecho, Instance.partition[p], cut_name, true, x_coefs, theta_coefs);
					add(new_cut <= derecho).end();
					org.optCuts += 1;
					new_cut.end();
				}
			}
			else {
				// do single cut
				single_lhs += lhs[p];
				const_part_single += lhs[p] - const_fs;
				for (size_t i = 0; i < Instance.nFacilities; i++) {
					new_single_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * sp_info_agg[p].lambda[Instance.nClients + i];
				}
			}
		}
		if (smallerWithTolerance(current_theta[0], single_lhs) && org.algo == 's') {
			violated = true;
			new_single_cut -= Instance.ent_sflp.theta[0];
			//cout << new_single_cut << " <= " << const_part_single << endl;
			add(new_single_cut <= -const_part_single).end();
			org.optCuts += 1;
		}
	}
}


ILONODECALLBACK2(Onodecallback_test, size_t&, cnode, size_t&, depth) {
	IloInt nextnode = 0;
	cnode = getNodeId(nextnode)._id;
	depth = getDepth(nextnode);
	//myNodeData * a_new_node_data;
	//a_new_node_data = dynamic_cast <myNodeData *> (getNodeData(cnode));
	//IloInt dad = a_new_node_data->Id;
	//cout << "Father is: " << dad << endl;
	//cout << "We are at node: " << cnode << " and " << ti.count() << endl;
}


ILOMIPINFOCALLBACK7(Oinfocallback_test, vector<vector<size_t>>&, partitions, MyClock&, control_clock, size_t&, current_slot, size_t&, cnode, bool&, oofrootnode, solFeat&, org, size_t&, cont) {
	double best_integer;
	double best_bound;
	double explored_nodes;
	double rem_nodes;
	double GAP;
	control_clock.end = control_clock.sc.now();
	auto ti = static_cast<chrono::duration<double>>(control_clock.end - control_clock.start);
	if (cnode > 0 && oofrootnode == false) {
		//Do procedure to print infor when we went out of the root node
		cout << endl;
		cout << "EXPLORATION OF ROOT NODE FINISHED:";
		oofrootnode = true;
		best_integer = getIncumbentObjValue();
		best_bound = getBestObjValue();
		explored_nodes = getNnodes64();
		rem_nodes = getNremainingNodes64();
		GAP = getMIPRelativeGap();
		cout << "PartialReport_" << cont << ": "
			<< 8 << " " //8 means status running
			<< ti.count() << " "
			<< best_integer << " "
			<< best_bound << " "
			<< GAP << " "
			<< explored_nodes << " "
			<< rem_nodes << " "
			<< partitions.size() << " "
			<< org.optCuts << " "
			<< org.user_optCuts
			<< endl;
		cont++;
	}
	if (ti.count() > current_slot) {
		cout << "Testing info callback: time spent: " << ti.count() << endl;
		current_slot = current_slot + timeslot;
		best_integer = getIncumbentObjValue();
		best_bound = getBestObjValue();
		explored_nodes = getNnodes64();
		rem_nodes = getNremainingNodes64();
		GAP = getMIPRelativeGap();
		cout << "PartialReport_" << cont << ": "
			<< 8 << " "
			<< ti.count() << " "
			<< best_integer << " "
			<< best_bound << " "
			<< GAP << " "
			<< explored_nodes << " "
			<< rem_nodes << " "
			<< partitions.size() << " "
			<< org.optCuts << " "
			<< org.user_optCuts
			<< endl;
		cont++;
	}
}




OuterBendersSFLP::OuterBendersSFLP(string &inst_name, string &stoch_inst_name) : SFLP_GAPM(inst_name, stoch_inst_name) {
	cp.resize(0);
	current_theta.resize(nScenarios, 0.0);
	num_old_cuts = 0;
}

IloModel OuterBendersSFLP::MasterProblemCreation(const char &algo, bool removeobjs) {
	//Add variables
	AddVarsMaster(*this, algo);

	//Valid inequalities
	ValidInequalities(*this, algo);

	//Adding feasibility constraint
	FeasibilityConstraint(*this);

	return this->Mast_mod;
}

void OuterBendersSFLP::remove_inactive() {
	//remove inactive cuts from structure
	for (auto it = cp.begin(); it < cp.end(); it++) {
		if (it->active == false) {
			cp.erase(it);
		}
	}
}

void OuterBendersSFLP::MasterProblemModification(IloCplex *cplex_master, IloModel &master, const char &algo, bool removeobjs) {
	//Remove inactive cuts
	size_t current_num_cuts = cp.size();
	size_t to_remove = 0;
	if (removeobjs == true) {
		for (size_t p = 0; p < num_old_cuts; p++) {
			if (cp[p].active == false) {
				to_remove++;
				master.remove(cp[p].cut_pool_constr);
			}
		}
	}

	//Add new cuts
	if (removeobjs == false) {
		for (size_t p = num_old_cuts; p < current_num_cuts; p++) {
			master.add(cp[p].cut_pool_constr);
			num_old_cuts++;
		}
	}

	//remove inactive cuts from structure
	if (removeobjs == true) {
		cp.erase(
			remove_if(cp.begin(), cp.end(), [&](cut_pool const & c_p) {
			return c_p.active == false;
		}), cp.end());
		/*for (auto it = cp.begin(); it < cp.begin()+num_old_cuts; it++) {
			if (it->active == false) {
				cp.erase(it);
			}
		}*/
		num_old_cuts = num_old_cuts - to_remove;
	}
}

//Solving the master problem
void OuterBendersSFLP::MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL)
{
	//Setting up cplex
	//cplex_master->setParam(IloCplex::Param::Preprocessing::Reduce, CPX_PREREDUCE_NOPRIMALORDUAL);
	//cplex_master->setParam(IloCplex::Param::Preprocessing::Linear, IloFalse);
	cplex_master->setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
	//cplex_master->setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
	//cplex_master->setParam(IloCplex::Param::Threads, 1);
	cplex_master->setParam(IloCplex::Param::TimeLimit, TL);
	//cplex_master->setOut(SFLP_sp_cpx.getNullStream());

	//Solving
	cplex_master->extract(master);

	//Mip start
	IloNumArray xb(cplex_master->getEnv());
	for (size_t i = 0; i < nFacilities; i++) {
		xb.add(x_bar[i]);
	}
	cplex_master->addMIPStart(ent_sflp.x, xb);
	xb.end();

	cplex_master->exportModel("test.lp");

	auto control_clock = MyClock();
	size_t current_slot = 0;
	bool oorootnode = false;
	size_t report = 1;

	solFeat org;
	org.algo = 'o';
	org.feasCuts = 0;
	org.optCuts = 0;
	org.user_optCuts = 0;

	bool first_out = false;
	cplex_master->use(Onodecallback_test(this->Mast_Bend, org.cnode, org.depth));
	cplex_master->use(OLazyOptCuts(this->Mast_Bend, *this, org, first_out, org.algo));
	cplex_master->use(OUserOptCuts(this->Mast_Bend, *this, org));
	cplex_master->use(Oinfocallback_test(this->Mast_Bend, partition, control_clock, current_slot, org.cnode, oorootnode, org, report));

	//cplex_master->exportModel("mastergapm.lp");
	cplex_master->solve();

	// Solution recovery
	x_bar.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		x_bar[i] = cplex_master->getValue(ent_sflp.x[i]);
		//if (x_bar[i] > 0)
			//cout << "Facility " << (i + 1) << " will be open" << endl;
	}

	current_theta.resize(nScenarios);
	for (size_t s = 0; s < nScenarios; s++) {
		current_theta[s] = cplex_master->getValue(ent_sflp.theta[s]);
		//if (current_theta[s] > 0)
			//cout << "Theta value " << (s + 1) <<  endl;
	}

	LB = cplex_master->getObjValue();
	cout << "Current Lower Bound is " << LB << endl;

}

//New iteration cuts
void OuterBendersSFLP::new_cuts(vector<vector<double>> &stoch_agg, const vector<solution_sps>* sp_info_agg, bool &violated, bool benders) {
	vector<double> lhs;
	lhs.resize(partition.size());
	//For each element into a partition we check if there exist a violated cut
	for (size_t p = 0; p < partition.size(); p++) {
		//Compute the current rhs to check if the cut is violated
		double element_prob = 0.0;
		double theta_ac = 0.0;
		for (size_t s = 0; s < partition[p].size(); s++) {
			element_prob += probability[partition[p][s]];
			theta_ac += current_theta[partition[p][s]] * probability[partition[p][s]];
		}

		//constant of current first stage multiplied by duals of capacity constarints
		double const_fs = 0.0;
		for (size_t i = 0; i < nFacilities; i++) {
			const_fs += element_prob * facil_capacities[i] * x_bar[i] * sp_info_agg->at(p).lambda[nClients + i];
		}
		double rhs = theta_ac;
		//compute the constant part of the cut
		lhs[p] = 0.0;
		for (size_t j = 0; j < nClients; j++) {
			lhs[p] += element_prob * stoch_agg[p][j] * sp_info_agg->at(p).lambda[j];
		}
		lhs[p] += const_fs;

		if (smallerWithTolerance(rhs, lhs[p])) {
			violated = true;
			IloExpr new_cut(Mast_Bend);
			vector<double> x_coefs;
			vector<double> theta_coefs;
			for (size_t i = 0; i < nFacilities; i++) {
				x_coefs.push_back(element_prob * facil_capacities[i] * sp_info_agg->at(p).lambda[nClients + i]);
				new_cut += x_coefs.back() * ent_sflp.x[i];
			}
			for (size_t s = 0; s < partition[p].size(); s++) {
				theta_coefs.push_back(probability[partition[p][s]]);
				new_cut -= ent_sflp.theta[partition[p][s]] * theta_coefs.back();
			}
			double derecho = -1.0 * (lhs[p] - const_fs);
			//cout << new_cut << " <= " << derecho << endl;
			string cut_name = "test_" + to_string(p + 1);
			cp.emplace_back(&Mast_Bend, &new_cut, derecho, partition[p], cut_name, true, x_coefs, theta_coefs);
			new_cut.end();
		}
	}
}

//check active cuts
void OuterBendersSFLP::check_active_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool benders) {
	size_t current_num_cuts = cp.size();
	vector<double> lhs;
	lhs.resize(cp.size());
	for (size_t p = 0; p < current_num_cuts; p++) {
		//Compute the current rhs to check if the cut is violated
		double theta_ac = 0.0;
		lhs[p] = 0.0;
		for (size_t s = 0; s < cp[p].elem_cut.size(); s++) {
			theta_ac += current_theta[cp[p].elem_cut[s]] * probability[cp[p].elem_cut[s]];

			//constant of current first stage multiplied by duals of capacity constarints
			double const_fs = 0.0;
			for (size_t i = 0; i < nFacilities; i++) {
				const_fs += probability[cp[p].elem_cut[s]] * facil_capacities[i] * x_bar[i] * sp_info->at(cp[p].elem_cut[s]).lambda[nClients + i];
			}

			//compute the constant part of the cut
			for (size_t j = 0; j < nClients; j++) {
				lhs[p] += probability[cp[p].elem_cut[s]] * stoch[cp[p].elem_cut[s]][j] * sp_info->at(cp[p].elem_cut[s]).lambda[j];
			}
			lhs[p] += const_fs;
		}

		if (fabs(lhs[p] - theta_ac) > tolrelcuts) {
			cp[p].active = false;
		}
		else {
			cp[p].active = true;
		}
	}
	//size_t a = 1;
}

