// Created by CDRPico
// 09/07/2020 06:34

#include"../inc/SFCMFP_GAPM.h"
//#include"../inc/UsefulFunctions.h"

//Given the instance name and the stochastic name instances we read the data for the problem
SFCMFP_GAPM::SFCMFP_GAPM(string &inst_name, string &stoch_inst_name) {
	//Call auxiliary function which is templated to read data of any problem
	//using the same structure
	read_instance_data(*this, inst_name, stoch_inst_name);
	//initialize partition
	partition.push_back(vector<size_t>(create_partition(numScen)));
	//initialize values of obj sps
	obj_opt_sps.resize(numScen, 0.0);
	x_bar.resize(numArcs, 0.0);
	//Origin and destination of arcs
    orig();
    dest();
}

vector<double> SFCMFP_GAPM::expected_demand(const vector<size_t> &scenarios)
{
	size_t Comm = numComm;
	size_t partitionSize = scenarios.size();
	vector<double> PromDemanda(Comm);
	//>Computing average demand
	for (size_t k = 0; k < Comm; k++) {
		for (size_t a = 0; a < partitionSize; a++) {
			PromDemanda[k] += DemandS[scenarios[a] - 1][k];
		}
		PromDemanda[k] *= 1.0 / double(partitionSize);
	}
	return PromDemanda;
}


//Create and solve master problem given a certain partition
IloModel SFCMFP_GAPM::MasterProblemCreation(const char &algo, bool removeobjs)
{
	size_t total_elements = partition.size();
	//Master env
	IloEnv SFCMFP = IloEnv();

	//Master model
	IloModel master = IloModel(SFCMFP);
    try {

		// Creating first stage variables
		// Facilities to be open
		master_entities.x = IloNumVarArray(SFCMFP, numArcs);
		//Populate OF
		master_entities.objective = IloExpr(SFCMFP);
		// Second stage decision variables (distribution of clients)
		master_entities.y = IloNumVarArray3(SFCMFP, numArcs);
		for (size_t i = 0; i < numArcs; i++) {
            master_entities.x[i] = IloNumVar(SFCMFP, 0.0, 1.0, ILOINT);
			master.add(master_entities.x[i]);
			master_entities.objective += Cij[i] * master_entities.x[i];

			//Completing second stage variables
			master_entities.y[i] = IloNumVarArray2(SFCMFP, numComm);
			for (size_t k = 0; k < numComm; k++) {
				master_entities.y[i][k] = IloNumVarArray(SFCMFP, total_elements);
				for (size_t s = 0; s < total_elements; s++) {
					master_entities.y[i][k][s] = IloNumVar(SFCMFP, 0.0, DBL_MAX, ILOFLOAT);
					master.add(master_entities.y[i][k][s]);
					double elementProb = double(partition[s].size()) / double(numScen);
					master_entities.objective += (elementProb * Vij[i] * master_entities.y[i][k][s]);
				}
			}
		}

		//Add objective function
		master.add(IloMinimize(SFCMFP, master_entities.objective));

		//Create the array s corresponding to the constraints that will be added
		master_entities.flow_satisfaction = IloRangeArray(SFCMFP);
		master_entities.cap_constraints = IloRangeArray(SFCMFP);
		
        //Flow Constraints
		for (size_t s = 0; s < total_elements; s++) {
			size_t NumCon = 0;
			for (size_t k = 0; k < numComm; k++)
			{
				for (size_t i = 0; i < numNodes; i++)
				{
					vector<double> PromDemanda = expected_demand(partition[s]);
					// LHS
					IloExpr constr2(SFCMFP);
					// RHS
					double RHS_FlowConst;
					// Filling expressions
					for (size_t j = posor[i]; j < posor[i + 1]; j++) {
						string varName = "x_" + to_string(incoming[j]) + "_" + to_string(i + 1) + "_" + to_string(k);
						size_t arc;
						for (size_t tt = 0; tt < numArcs; tt++) {
							if (OArcs[tt] == incoming[j] && DArcs[tt] == i + 1) {
								arc = tt;
							}
						}
						constr2 = constr2 - master_entities.y[arc][k][s];
					}
					for (size_t j = posde[i]; j < posde[i + 1]; j++) {
						string varName = "x_" + to_string(i + 1) + "_" + to_string(outgoing[j]) + "_" + to_string(k);
						size_t arc;
						for (size_t tt = 0; tt < numArcs; tt++) {
							if (OArcs[tt] == i + 1 && DArcs[tt] == outgoing[j]) {
								arc = tt;
							}
						}
						constr2 = constr2 + master_entities.y[arc][k][s];
					}

					if (i == DemFrom[k] - 1) {
						RHS_FlowConst = PromDemanda[k];						
					}
					else if (i == DemTo[k] - 1) {
						RHS_FlowConst = -PromDemanda[k];
					}
					else {
						RHS_FlowConst = 0;
					}
					// Adding constraint to the object
					master_entities.flow_satisfaction.add(constr2 == RHS_FlowConst);
					// Closing Linear expression
					constr2.end();
				}
			}
		}
		//Adding whole set of flow constraints to the model
		master.add(master_entities.flow_satisfaction);

		//Capacity constraints
		for (size_t s = 0; s < total_elements; s++) {
			size_t NumCon = 0;
			for (size_t i = 0; i < numArcs; i++) {
				// LHS
				IloExpr constr3(SFCMFP);
				// Filling expressions
				for (size_t k = 0; k < numComm; k++) {
					constr3 += master_entities.y[i][k][s];
				}
				constr3 -= (CAPij[i] * master_entities.x[i]);
				// Adding constraint to the object IloRange
				master_entities.cap_constraints.add(constr3 <= 0);
				// Closing Linear Expression
				constr3.end();
			}
		}
		master.add(master_entities.cap_constraints);
	}
	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...) {
		cerr << "Other Exception" << endl;
	}
	return master;
}


void SFCMFP_GAPM::MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL)
{
	//Setting up cplex
	cplex_master->setParam(IloCplex::Param::Threads, 1);
	cplex_master->setParam(IloCplex::Param::TimeLimit, TL);
	cplex_master->setParam(IloCplex::Param::Benders::Strategy, 3);

	//Solving
	cplex_master->extract(master);
	//cplex_master->exportModel("mastergapm.lp");
	cplex_master->solve();

	// Solution recovery
	x_bar.resize(numArcs);
	for (size_t i = 0; i < numArcs; i++) {
		x_bar[i] = cplex_master->getValue(master_entities.x[i]);
		if (x_bar[i] > 0)
			cout << "Arc " << (i + 1) << " from " << OArcs[i] << " to " << DArcs[i] << " will be built" << endl;
	}

	LB = cplex_master->getObjValue();

}

// TODO Create functions of the subproblem