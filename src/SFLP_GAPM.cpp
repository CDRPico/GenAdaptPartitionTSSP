// Created by CDRPico
// 11/06/2020 13:47

#include"../inc/SFLP_GAPM.h"
//#include"../inc/UsefulFunctions.h"

//Given the instance name and the stochastic name instances we read the data for the problem
SFLP_GAPM::SFLP_GAPM(string &inst_name, string &stoch_inst_name) {
	//Call auxiliary function which is templated to read data of any problem
	//using the same structure
	read_instance_data(*this, inst_name, stoch_inst_name);
	compute_varcosts();
	//initialize partition
	partition.push_back(vector<size_t>(create_partition(nScenarios)));
	//initialize values of obj sps
	obj_opt_sps.resize(nScenarios, 0.0);
	x_bar.resize(nFacilities, 0.0);
	//instatiate subproblem
	//SPProblemCreation();
}


vector<double> SFLP_GAPM::expected_demand(vector<size_t> &element) {
	vector<double> agg_demand_elem(nClients, 0.0);
	double partProb = 0.0;
	for (size_t s = 0; s < element.size(); s++) {
		partProb += probability[s];
	}
	//Computing average demands
	for (size_t j = 0; j < nClients; j++) {
		for (size_t s = 0; s < element.size(); s++) {
			agg_demand_elem[j] += stoch_param[j][element[s]] * probability[s] / partProb;
		}
	}
	return agg_demand_elem;
}


//Create and solve master problem given a certain partition
IloModel SFLP_GAPM::MasterProblemCreation(bool removeobjs)
{
	size_t total_elements = partition.size();
	//Master env
	IloEnv SFLP = IloEnv();

	//Master model
	IloModel master = IloModel(SFLP);
	try {

		// Creating first stage variables
		// Facilities to be open
		master_entities.x = IloNumVarArray(SFLP, nFacilities);
		//Populate OF
		master_entities.objective = IloExpr(SFLP);
		// Second stage decision variables (distribution of clients)
		master_entities.y = IloNumVarArray3(SFLP, nFacilities);
		for (size_t i = 0; i < nFacilities; i++) {
			master_entities.x[i] = IloNumVar(SFLP, 0.0, 1.0, ILOINT);
			master.add(master_entities.x[i]);
			master_entities.objective += fixed_costs[i] * master_entities.x[i];

			//Completing second stage variables
			master_entities.y[i] = IloNumVarArray2(SFLP, nClients);
			for (size_t j = 0; j < nClients; j++) {
				master_entities.y[i][j] = IloNumVarArray(SFLP, total_elements);
				for (size_t s = 0; s < total_elements; s++) {
					master_entities.y[i][j][s] = IloNumVar(SFLP, 0.0, DBL_MAX, ILOFLOAT);
					master.add(master_entities.y[i][j][s]);
					double elementProb = double(partition[s].size()) / double(nScenarios);
					master_entities.objective += (elementProb * dist_costs[i][j] * master_entities.y[i][j][s]);
				}
			}
		}

		//Add objective function
		master.add(IloMinimize(SFLP, master_entities.objective));

		//Create the arrays corresponding to the constraints that will be added
		master_entities.dem_satisfaction = IloRangeArray(SFLP);
		master_entities.linking_constraints = IloRangeArray(SFLP);
		//Feasibility at the end constraint 3
		master_entities.Feasibility = IloRangeArray(SFLP);

		//Demand constraints
		double max_demand = 0.0;
		for (size_t s = 0; s < total_elements; s++) {
			double scenario_tot_demand = 0.0;
			vector<double> elem_demand = expected_demand(partition[s]);
			for (size_t j = 0; j < nClients; j++) {
				//LHS
				IloExpr constr1(SFLP);
				for (size_t i = 0; i < nFacilities; i++) {
					constr1 += master_entities.y[i][j][s];
				}
				//Adding the constraint
				master_entities.dem_satisfaction.add(constr1 >= elem_demand[j]);
				//remove linear expression
				constr1.end();
				//total scenario demand
				scenario_tot_demand += elem_demand[j];
			}
			//update the maximum demand across scenarios
			if (scenario_tot_demand > max_demand) {
				max_demand = scenario_tot_demand;
			}
		}
		master.add(master_entities.dem_satisfaction);

		//linking constraints
		for (size_t s = 0; s < total_elements; s++) {
			for (size_t i = 0; i < nFacilities; i++) {
				//RHS
				IloExpr constr2(SFLP);
				for (size_t j = 0; j < nClients; j++) {
					constr2 += master_entities.y[i][j][s];
				}
				//RHS
				constr2 -= facil_capacities[i] * master_entities.x[i];
				//Adding the constraint
				master_entities.linking_constraints.add(constr2 <= 0);
				//Remove the linear expression
				constr2.end();
			}
		}
		master.add(master_entities.linking_constraints);

		//Feasibility constraint
		//LHS
		IloExpr constr3(SFLP);
		for (size_t i = 0; i < nFacilities; i++) {
			constr3 += facil_capacities[i] * master_entities.x[i];
		}
		//the constraint
		master_entities.Feasibility.add(constr3 >= max_demand);
		master.add(master_entities.Feasibility);
		constr3.end();

	}
	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception" << endl;
	}
	return master;
}

void SFLP_GAPM::MasterProblemSolution(IloModel &master, double &LB, const double &TL)
{
	//Setting up cplex
	IloEnv SFLP = master.getEnv();
	IloCplex cplex_master(SFLP);
	cplex_master.setParam(IloCplex::Param::Threads, 1);
	cplex_master.setParam(IloCplex::Param::TimeLimit, TL);
	cplex_master.setParam(IloCplex::Param::Benders::Strategy, 3);

	//Solving
	cplex_master.extract(master);
	cplex_master.solve();

	// Solution recovery
	x_bar.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		x_bar[i] = cplex_master.getValue(master_entities.x[i]);
		if (x_bar[i] > 0)
			cout << "Facility " << (i + 1) << " will be open" << endl;
	}

	LB = cplex_master.getObjValue();

}

void SFLP_GAPM::SPProblemCreation_CPX() {
	try {
		//Create variables
		//Demand distribution
		sp_entities_cpx.y = IloNumVarArray2(SFLP_sp_cpx, nFacilities);
		//Populate OF
		sp_entities_cpx.objective = IloExpr(SFLP_sp_cpx);
		for (size_t i = 0; i < nFacilities; i++) {
			sp_entities_cpx.y[i] = IloNumVarArray(SFLP_sp_cpx, nClients);
			for (size_t j = 0; j < nClients; j++) {
				sp_entities_cpx.y[i][j] = IloNumVar(SFLP_sp_cpx, 0.0, DBL_MAX, ILOFLOAT);
				subprob_cpx.add(sp_entities_cpx.y[i][j]);
				sp_entities_cpx.objective += dist_costs[i][j] * sp_entities_cpx.y[i][j];
			}
		}

		//Add objective function
		subprob_cpx.add(IloMinimize(SFLP_sp_cpx, sp_entities_cpx.objective));
		//Constraints
		//Demand satisfaction
		sp_entities_cpx.dem_constraints = IloRangeArray(SFLP_sp_cpx);
		// Facilities capacity
		sp_entities_cpx.cap_constraints = IloRangeArray(SFLP_sp_cpx);
		for (size_t j = 0; j < nClients; j++) {
			//LHS
			IloExpr constr1(SFLP_sp_cpx);
			for (size_t i = 0; i < nFacilities; i++) {
				constr1 += sp_entities_cpx.y[i][j];
			}
			//Adding the constraint
			sp_entities_cpx.dem_constraints.add(constr1 >= stoch_param[j][0]);
			//remove linear expression
			constr1.end();
		}
		subprob_cpx.add(sp_entities_cpx.dem_constraints);

		for (size_t i = 0; i < nFacilities; i++) {
			//RHS
			IloExpr constr2(SFLP_sp_cpx);
			for (size_t j = 0; j < nClients; j++) {
				constr2 += sp_entities_cpx.y[i][j];
			}
			//RHS
			//Adding the constraint
			sp_entities_cpx.cap_constraints.add(constr2 <= x_bar[i] * facil_capacities[i]);
			//Remove the linear expression
			constr2.end();
		}
		subprob_cpx.add(sp_entities_cpx.cap_constraints);

	}
	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception " << endl;
	}
}

//Modify the subproblem to solve a different scenario
void SFLP_GAPM::SPProblemModification_CPX(vector<size_t> &element, bool mod_x) {
	// s is used to modify RHS
	//modifying demands
	vector<double> elem_demand = expected_demand(element);
	for (size_t j = 0; j < nClients; j++) {
		sp_entities_cpx.dem_constraints[j].setLB(elem_demand[j]);
		//cout << sp_entities_cpx.dem_constraints[j] << endl << endl;
	}
	//Modifying composition of facilities built
	if (mod_x == true) {
		for (size_t i = 0; i < nFacilities; i++) {
			sp_entities_cpx.cap_constraints[i].setUB(x_bar[i] * facil_capacities[i]);
		}
	}
}

void SFLP_GAPM::SPProbleSolution_CPX(vector<double> &stoch, solution_sps *sp_info) {
	//Setting up cplex
	IloCplex cplex_sp(SFLP_sp_cpx);
	//cplex_sp.setParam(IloCplex::RootAlg, IloCplex::Primal);
	cplex_sp.setParam(IloCplex::PreInd, 0);
	cplex_sp.setOut(SFLP_sp_cpx.getNullStream());
	cplex_sp.extract(subprob_cpx);
	cplex_sp.exportModel("prueba_subproblem.lp");
	cplex_sp.solve();

	stoch.clear();
	sp_info->lambda.clear();
	if (cplex_sp.getStatus() == IloAlgorithm::Optimal) {
		cout << "duales ";
		for (size_t j = 0; j < nClients; j++) {
			vector<double> elem_demand = expected_demand(sp_info->scen);
			stoch.push_back(elem_demand[j]);
			sp_info->lambda.push_back(cplex_sp.getDual(sp_entities_cpx.dem_constraints[j]));
			//cout << "demanda " << j << " " << elem_demand[j];
			cout << sp_info->lambda.back() << " ";
		}
		cout << endl;
	}
	else {
		cerr << "Optimal solution of the problem was not found!" << endl;
	}
	sp_info->obj = cplex_sp.getObjValue();
}


/*double SFLP_GAPM::compute_UB(vector<solution_sps> &sp_info, double solution_sps::*obj) {
	double UB_1 = 0.0;
	for (size_t i = 0; i < nFacilities; i++) {
		UB_1 += x_bar[i] * fixed_costs[i];
	}

	for (solution_sps & c : sp_info) {
		for (size_t s = 0; s < nScenarios; s++) {
			UB_1 += probability[s] * c[s].obj;
		}
	}
	
	return UB_1;
}*/

void SFLP_GAPM::Labeling_y()
{
	Label_y.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		Label_y[i].resize(nClients);
		for (size_t j = 0; j < nClients; j++) {
			Label_y[i][j] = "y_" + to_string(i + 1) + "_" + to_string(j + 1);
		}
	}
}


void SFLP_GAPM::Label_demand_constr()
{
	Label_demconst.resize(nClients);
	for (size_t j = 0; j < nClients; j++) {
		Label_demconst[j] = "demconst_" + to_string(j + 1);
	}
}

void SFLP_GAPM::Label_capacity_constr()
{
	Label_capconst.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		Label_capconst[i] = "capconst_" + to_string(i + 1);
	}
}

void SFLP_GAPM::SPProblemCreation_GRB() {
	try {
		//Creating variables 
		Labeling_y();
		sp_entities_grb.y = new GRBVar *[nFacilities];
		for (size_t i = 0; i < nFacilities; i++) {
			sp_entities_grb.y[i] = new GRBVar[nClients];
			for (size_t j = 0; j < nClients; j++) {
				sp_entities_grb.y[i][j] = subprob_grb.addVar(0.0, GRB_INFINITY, dist_costs[i][j],
					GRB_CONTINUOUS, Label_y[i][j].c_str());
			}
		}

		//Creating the constraints
		//Demand constraints
		Label_demand_constr();
		for (size_t j = 0; j < nClients; j++) {
			//LHS
			GRBLinExpr constr1;
			for (size_t i = 0; i < nFacilities; i++) {
				constr1 = constr1 + sp_entities_grb.y[i][j];
			}
			//Adding the constraint
			subprob_grb.addConstr(constr1 >= stoch_param[j][0], Label_demconst[j]);
			//Clear expression
			constr1.clear();
		}

		//Capacity constraints
		Label_capacity_constr();
		for (size_t i = 0; i < nFacilities; i++) {
			//LHS
			GRBLinExpr constr2;
			for (size_t j = 0; j < nClients; j++) {
				constr2 = constr2 + sp_entities_grb.y[i][j];
			}
			//Adding the constraint
			subprob_grb.addConstr(constr2 <= facil_capacities[i] * x_bar[i], Label_capconst[i]);
			//Clear the expression
			constr2.clear();
		}

		subprob_grb.update();
	}
	catch (GRBException e) {
		cerr << "Error number: " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error GRB" << endl;
	}
}


void SFLP_GAPM::SPProblemModification_GRB(vector<size_t> &element, bool mod_x) {
	// s is used to modify RHS
	//modifying demands
	vector<double> elem_demand = expected_demand(element);
	for (size_t j = 0; j < nClients; j++) {
		GRBConstr dem_con = subprob_grb.getConstrByName(Label_demconst[j]);
		//Modify the value of RHS
		dem_con.set(GRB_DoubleAttr_RHS, elem_demand[j]);
	}
	//Modifying composition of facilities built
	if (mod_x == true) {
		for (size_t i = 0; i < nFacilities; i++) {
			GRBConstr cap_con = subprob_grb.getConstrByName(Label_capconst[i]);
			//Modify the RHS
			cap_con.set(GRB_DoubleAttr_RHS, facil_capacities[i] * x_bar[i]);
		}
	}
	//Updating the model
	subprob_grb.update();
}

//Solving the subproblem via Gurobi
void SFLP_GAPM::SPProbleSolution_GRB(vector<double> &stoch, solution_sps *sp_info) {
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
		for (size_t j = 0; j < nClients; j++) {
			vector<double> elem_demand = expected_demand(sp_info->scen);
			stoch.push_back(elem_demand[j]);
			sp_info->lambda.push_back(subprob_grb.getConstrByName(Label_demconst[j].c_str()).get(GRB_DoubleAttr_Pi));
		}
		sp_info->obj = subprob_grb.get(GRB_DoubleAttr_ObjVal);
	}
	//If infeasible, get extreme ray
	else if (sp_info->F == 3) {
		for (size_t j = 0; j < nClients; j++) {
			vector<double> elem_demand = expected_demand(sp_info->scen);
			stoch.push_back(elem_demand[j]);
			sp_info->lambda.push_back(-1.0 * subprob_grb.getConstrByName(Label_demconst[j].c_str()).get(GRB_DoubleAttr_FarkasDual));
		}
		sp_info->obj = DBL_MAX;
	}
}