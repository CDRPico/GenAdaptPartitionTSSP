// Created by CDRPico
// 20/06/2020 00:48

#include"../inc/BendersAPM_SFLP.h"
#include<chrono>

ILOLAZYCONSTRAINTCALLBACK3(LazyOptCuts,
    BendersSFLP&, Instance,
    solFeat&, org,
	bool&, first_out)
{
	
	//Callback execution time
	chrono::steady_clock sc_cut;
	auto start_cut = sc_cut.now();
	auto end_cut = sc_cut.now();
	Instance.sp_info.resize(Instance.nScenarios);
	Instance.stoch.resize(Instance.nScenarios);
    size_t scenarios = 0;
    if (org.algo == 's') {
        scenarios = 1;
    } else {
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
    Instance.stoch_agg.clear();
	Instance.stoch_agg.resize(Instance.partition.size());
	Instance.sp_info_agg.clear();
	Instance.sp_info_agg.resize(Instance.partition.size());
    Instance.lhs.clear();
    Instance.lhs.resize(Instance.partition.size());

    //Checking if any cut should added
    bool violated = false;
    double single_lhs = 0.0;
    IloExpr new_single_cut(getEnv());
    double const_part_single = 0.0;

    for (size_t p = 0; p < Instance.partition.size(); p++) {
		Instance.sp_info_agg[p].scen = Instance.partition[p];
        if (p == 0) {
			Instance.SPProblemModification_GRB(Instance.partition[p], true);
		}
		else {
			Instance.SPProblemModification_GRB(Instance.partition[p]);
		}
		Instance.SPProbleSolution_GRB(Instance.stoch_agg[p], &Instance.sp_info_agg[p], true);
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
            const_fs += element_prob * Instance.facil_capacities[i] * Instance.x_bar[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
        }
        double rhs = theta_ac;
        //compute the constant part of the cut
        Instance.lhs[p] = 0.0;
        for (size_t j = 0; j < Instance.nClients; j++) {
            Instance.lhs[p] += element_prob * Instance.stoch_agg[p][j] * Instance.sp_info_agg[p].lambda[j];
        }
        Instance.lhs[p] += const_fs;
        if (org.algo != 's') {            
            if(smallerWithTolerance(rhs, Instance.lhs[p])) {
                violated = true;
                IloExpr new_cut(getEnv());
                for (size_t i = 0; i < Instance.nFacilities; i++) {
                    new_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
                }
                for (size_t s = 0; s < Instance.partition[p].size(); s++) {
                    new_cut -= Instance.ent_sflp.theta[Instance.partition[p][s]] * Instance.probability[Instance.partition[p][s]];
                }
                double derecho = -1.0 * (Instance.lhs[p] - const_fs);
				//cout << new_cut << " <= " << derecho << endl;
                add(new_cut <= derecho).end();
				org.optCuts += 1;
                new_cut.end();
            }
        } else {
            // do single cut
            single_lhs += Instance.lhs[p];
            const_part_single +=  Instance.lhs[p] - const_fs;
            for (size_t i = 0; i < Instance.nFacilities; i++) {
                new_single_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
            }
        }
    }
    if (smallerWithTolerance(current_theta[0],single_lhs) && org.algo == 's'){
        violated = true;
        new_single_cut -= Instance.ent_sflp.theta[0];
		//cout << new_single_cut << " <= " << const_part_single << endl;
        add(new_single_cut <= -const_part_single).end();
		org.optCuts += 1;
    }

	/*if (org.cnode > 0 && first_out == false && org.algo == 'a') {
		first_out = true;
		Instance.partition.clear();
		Instance.partition.push_back(vector<size_t>(create_partition(Instance.nScenarios)));
	}*/

    if (violated == false && org.algo == 'a' && Instance.partition.size() < Instance.nScenarios) {
		//chrono::steady_clock sc;
		//auto start = sc.now();
		//auto end = sc.now();
        for (size_t s = 0; s < Instance.nScenarios; s++) {
            vector<size_t> el(1, s);
            Instance.sp_info[s].scen = el;
            if (s == 0) {
                Instance.SPProblemModification_GRB(Instance.sp_info[s].scen, true);
            }
            else {
                Instance.SPProblemModification_GRB(Instance.sp_info[s].scen);
            }
            double obj = 0.0;
            Instance.SPProbleSolution_GRB(Instance.stoch[s], &Instance.sp_info[s]);
        }
		//end = sc.now();
		//auto time_span = static_cast<chrono::duration<double>>(end - start);
		//double runtime = time_span.count();
		//cout << "it takes " << runtime << "to solve " << Instance.nScenarios << " subproblems" << endl;;
        disag_procedure to_dis;
		to_dis.disaggregation(Instance.sp_info, Instance.partition, Instance.nScenarios);
    }
	Instance.sp_info.clear();
	Instance.stoch.clear();
	Instance.sp_info_agg.clear();
	Instance.stoch_agg.clear();
	end_cut = sc_cut.now();
	auto time_span_cut = static_cast<chrono::duration<double>>(end_cut - start_cut);
	double runtime_cut = time_span_cut.count();
	cout << "it takes " << runtime_cut << "to execute the lazy callback " << endl;;
}

ILOUSERCUTCALLBACK2(UserOptCuts,
	BendersSFLP&, Instance,
	solFeat&, org)
{
	if (org.cnode == 0) {
		Instance.sp_info.resize(Instance.nScenarios);
		Instance.stoch.resize(Instance.nScenarios);
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
		Instance.stoch_agg.clear();
		Instance.stoch_agg.resize(Instance.partition.size());
		Instance.sp_info_agg.clear();
		Instance.sp_info_agg.resize(Instance.partition.size());
		Instance.lhs.clear();
		Instance.lhs.resize(Instance.partition.size());

		//Checking if any cut should added
		bool violated = false;
		double single_lhs = 0.0;
		IloExpr new_single_cut(getEnv());
		double const_part_single = 0.0;

		for (size_t p = 0; p < Instance.partition.size(); p++) {
			Instance.sp_info_agg[p].scen = Instance.partition[p];
			if (p == 0) {
				Instance.SPProblemModification_GRB(Instance.partition[p], true);
			}
			else {
				Instance.SPProblemModification_GRB(Instance.partition[p]);
			}
			Instance.SPProbleSolution_GRB(Instance.stoch_agg[p], &Instance.sp_info_agg[p], true);
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
				const_fs += element_prob * Instance.facil_capacities[i] * Instance.x_bar[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
			}
			double rhs = theta_ac;
			//compute the constant part of the cut
			Instance.lhs[p] = 0.0;
			for (size_t j = 0; j < Instance.nClients; j++) {
				Instance.lhs[p] += element_prob * Instance.stoch_agg[p][j] * Instance.sp_info_agg[p].lambda[j];
			}
			Instance.lhs[p] += const_fs;
			if (org.algo != 's') {
				if (smallerWithTolerance(rhs, Instance.lhs[p])) {
					violated = true;
					IloExpr new_cut(getEnv());
					for (size_t i = 0; i < Instance.nFacilities; i++) {
						new_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
					}
					for (size_t s = 0; s < Instance.partition[p].size(); s++) {
						new_cut -= Instance.ent_sflp.theta[Instance.partition[p][s]] * Instance.probability[Instance.partition[p][s]];
					}
					double derecho = -1.0 * (Instance.lhs[p] - const_fs);
					//cout << new_cut << " <= " << derecho << endl;
					add(new_cut <= derecho).end();
					new_cut.end();
					org.user_optCuts += 1;
				}
			}
			else {
				// do single cut
				single_lhs += Instance.lhs[p];
				const_part_single += Instance.lhs[p] - const_fs;
				for (size_t i = 0; i < Instance.nFacilities; i++) {
					new_single_cut += element_prob * Instance.facil_capacities[i] * Instance.ent_sflp.x[i] * Instance.sp_info_agg[p].lambda[Instance.nClients + i];
				}
			}
		}
		if (smallerWithTolerance(current_theta[0], single_lhs) && org.algo == 's') {
			violated = true;
			new_single_cut -= Instance.ent_sflp.theta[0];
			add(new_single_cut <= -const_part_single).end();
			org.user_optCuts += 1;
		}

		Instance.sp_info.clear();
		Instance.stoch.clear();
		Instance.sp_info_agg.clear();
		Instance.stoch_agg.clear();
	}	
}


ILONODECALLBACK2(nodecallback_test, size_t&, cnode, size_t&, depth) {
	IloInt nextnode = 0;
	cnode = getNodeId(nextnode)._id;
	depth = getDepth(nextnode);
	//myNodeData * a_new_node_data;
	//a_new_node_data = dynamic_cast <myNodeData *> (getNodeData(cnode));
	//IloInt dad = a_new_node_data->Id;
	//cout << "Father is: " << dad << endl;
	//cout << "We are at node: " << cnode << " and " << ti.count() << endl;
}


ILOMIPINFOCALLBACK7(infocallback_test, vector<vector<size_t>>&, partitions, MyClock&, control_clock, size_t&, current_slot, size_t&, cnode, bool&, oofrootnode, solFeat&, org, size_t&, cont) {
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

void BendersSFLP::CreateMaster(const char &algo) {
	size_t nscen = 0;
	if (algo == 's') {
		nscen = 1;
	}
	else {
		nscen = nScenarios;
	}

	vector<double> avg_scenarios = expected_demand(partition[0]);

	if (algo != 'a') {
		partition.clear();
		for (size_t s = 0; s < nScenarios; s++) {
			partition.push_back(vector<size_t>({ s }));
		}
	}
	
	//Declaring variables
    ent_sflp.x = IloNumVarArray(Mast_Bend, nFacilities);
	ent_sflp.objective = IloExpr(Mast_Bend);
    for (size_t i = 0; i < nFacilities; i++) {
        ent_sflp.x[i] = IloNumVar(Mast_Bend, 0.0, 1.0, ILOINT);
        Mast_mod.add(ent_sflp.x[i]);
		ent_sflp.objective += fixed_costs[i] * ent_sflp.x[i];
    }

    ent_sflp.theta = IloNumVarArray(Mast_Bend, nscen);
    for (size_t s = 0; s < nscen; s++) {
        ent_sflp.theta[s] = IloNumVar(Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
        Mast_mod.add(ent_sflp.theta[s]);
		if (nscen != 1) {
			ent_sflp.objective += probability[s] * ent_sflp.theta[s];
		}
    }
	if (algo == 's') {
		ent_sflp.objective += ent_sflp.theta[0];
	}


	//Add objective function
	Mast_mod.add(IloMinimize(Mast_Bend, ent_sflp.objective));

	//Valid inequalities
	ent_sflp.y = IloNumVarArray2(Mast_Bend, nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		ent_sflp.y[i] = IloNumVarArray(Mast_Bend, nClients);
		for (size_t j = 0; j < nClients; j++) {
			ent_sflp.y[i][j] = IloNumVar(Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
			Mast_mod.add(ent_sflp.y[i][j]);
		}
	}
	//Adding constraint to average demand
	IloRangeArray avg_dem_cons(Mast_Bend);
	for (size_t j = 0; j < nClients; j++) {
		//LHS
		IloExpr constr1(Mast_Bend);
		for (size_t i = 0; i < nFacilities; i++) {
			constr1 += ent_sflp.y[i][j];
		}
		//Adding the constraint
		avg_dem_cons.add(constr1 >= avg_scenarios[j]);
		//remove linear expression
		constr1.end();
	}
	Mast_mod.add(avg_dem_cons);

	IloRangeArray avg_cap_cons(Mast_Bend);
	for (size_t i = 0; i < nFacilities; i++) {
		//RHS
		IloExpr constr2(Mast_Bend);
		for (size_t j = 0; j < nClients; j++) {
			constr2 += ent_sflp.y[i][j];
		}
		//RHS
		//Adding the constraint
		constr2 -= ent_sflp.x[i] * facil_capacities[i];
		avg_cap_cons.add(constr2 <= 0);
		//Remove the linear expression
		constr2.end();
	}
	Mast_mod.add(avg_cap_cons);

	IloExpr Avg_cons(Mast_Bend);
	IloRangeArray AVG(Mast_Bend);

	for (size_t j = 0; j < nClients; j++) {
		for (size_t i = 0; i < nFacilities; i++) {
			Avg_cons -= dist_costs[i][j] * ent_sflp.y[i][j];
		}
	}
	if (algo == 's') {
		Avg_cons += ent_sflp.theta[0];
	}
	else {
		for (size_t s = 0; s < nScenarios; s++) {
			Avg_cons += probability[s] * ent_sflp.theta[s];
		}
	}
	AVG.add(Avg_cons >= 0);
	Avg_cons.end();
	Mast_mod.add(AVG);

    //Adding feasibility constraint
    ent_sflp.Feasibility = IloRangeArray(Mast_Bend);
    //Feasibility constraint
		//LHS
	IloExpr constr3(Mast_Bend);
	for (size_t i = 0; i < nFacilities; i++) {
		constr3 += facil_capacities[i] * ent_sflp.x[i];
	}
	//the constraint
	ent_sflp.Feasibility.add(constr3 >= max_dem);
	Mast_mod.add(ent_sflp.Feasibility);
	constr3.end();

	//subproblem creation
	SPProblemCreation_GRB();

	//Gathering solution info
	solFeat org;
	org.algo = algo;
	org.feasCuts = 0;
	org.optCuts = 0;
	org.user_optCuts = 0;

    /* Create cplex enviroment to solve the problem */
    IloCplex cplex(Mast_Bend);
    //Setting parameters to solve master problem
    cplex.setParam(IloCplex::Param::Preprocessing::Reduce, CPX_PREREDUCE_NOPRIMALORDUAL);
    cplex.setParam(IloCplex::Param::Preprocessing::Linear, IloFalse);
    cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::Param::MIP::Display, 4);
    cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
    cplex.setParam(IloCplex::Param::TimeLimit, timelimit);

    cplex.extract(Mast_mod);
	//cplex.exportModel("masterBenders.lp");

	auto control_clock = MyClock();
	size_t current_slot = 0;
	bool oorootnode = false;
	size_t report = 1;

	bool first_out = false;
	cplex.use(nodecallback_test(this->Mast_Bend, org.cnode, org.depth));
    cplex.use(LazyOptCuts(this->Mast_Bend, *this, org, first_out));
	cplex.use(UserOptCuts(this->Mast_Bend, *this, org));
	cplex.use(infocallback_test(this->Mast_Bend, partition, control_clock, current_slot, org.cnode, oorootnode, org, report));
	

	//Solve the model
	cplex.solve();

	//Recovering the solution
	RecoverSolution(cplex);

	//Printing the solution
	PrintSolution(org);
}

void BendersSFLP::RecoverSolution(IloCplex &master_cpx) {
	try {
		for (size_t i = 0; i < nFacilities; i++) {
			x_bar[i] = master_cpx.getValue(ent_sflp.x[i]);
			if (x_bar[i] > 0) {
				cout << " The facility: " << (i + 1) << " will be open" << endl;
			}
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

void BendersSFLP::PrintSolution(solFeat &org) {
	cout << "Final partition size is: " << partition.size() << endl;
	cout << "cap_sflp:"
		<< status << " "
		<< cpx_runtime << " "
		<< obj_fin << " "
		<< LB << " "
		<< GAP << " "
		<< exploredNodes << " "
		<< 0 << " "//0 nodes pending to be explored
		<< partition.size() << " "
		<< org.optCuts << " "
		<< org.user_optCuts
		<< endl;
}