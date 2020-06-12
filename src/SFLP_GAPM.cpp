// Created by CDRPico
// 11/06/2020 13:47

#include"../inc/SFLP_GAPM.h"
#include"../inc/UsefulFunctions.h"

//Given the instance name and the stochastic name instances we read the data for the problem
SFLP_GAPM::SFLP_GAPM(string &inst_name, string &stoch_inst_name){
    //Call auxiliary function which is templated to read data of any problem
    //using the same structure
    read_instance_data(*this, inst_name, stoch_inst_name);

    //initialize partition
    partition.push_back(vector<size_t>(create_partition(nScenarios)));
    //initialize values of obj sps
    obj_opt_sps.resize(nScenarios, 0.0);
}


vector<double> SFLP_GAPM::expected_demand(vector<size_t> &element){
    vector<double> agg_demand_elem(nClients, 0.0);
    //Computing average demands
    for (size_t j = 0; j < nClients; j++){
        for (size_t s = 0; s < element.size(); s++){
            agg_demand_elem[j] += stoch_param[j][element[s]] * probability[s];
        }
    }
    return agg_demand_elem;
}


//Create and solve master problem given a certain partition
void SFLP_GAPM::MasterProblemCreation ()
{
    size_t total_elements = partition.size();
    try{

        // Creating first stage variables
        // Facilities to be open
        master_entities.x = IloNumVarArray(SFLP, nFacilities);
        //Populate OF
        master_entities.objective = IloExpr(SFLP);
        // Second stage decision variables (distribution of clients)
        master_entities.y = IloNumVarArray3 (SFLP, nFacilities);
        for (size_t i = 0; i < nFacilities; i++){
            master_entities.x[i] = IloNumVar (SFLP, 0.0, 1.0, ILOINT);
            master.add(master_entities.x[i]);
            master_entities.objective += fixed_costs[i]*master_entities.x[i];
            
            //Completing second stage variables
            master_entities.y[i] = IloNumVarArray2 (SFLP, nClients);
            for (size_t j = 0; j < nClients; j++){
                master_entities.y[i][j] = IloNumVarArray(SFLP, total_elements);
                for (size_t s = 0; s < total_elements; s++){
                    master_entities.y[i][j][s] = IloNumVar (SFLP, 0.0, DBL_MAX, ILOFLOAT);
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
        for (size_t s = 0; s < total_elements; s++){
            double scenario_tot_demand = 0.0;
            vector<double> elem_demand = expected_demand(partition[s]);
            for (size_t j = 0; j < nClients; j++){
                //LHS
                IloExpr constr1(SFLP);
                for (size_t i = 0; i < nFacilities; i++){
                    constr1 += master_entities.y[i][j][s];
                }
                //Adding the constraint
                master_entities.dem_satisfaction.add(constr1 >= elem_demand[j]);
                //remove linear expression
                constr1. end();
                //total scenario demand
                scenario_tot_demand += elem_demand[j];
            }
            //update the maximum demand across scenarios
            if (scenario_tot_demand > max_demand){
                max_demand = scenario_tot_demand;
            }
        }
        master.add(master_entities.dem_satisfaction);

        //linking constraints
        for (size_t s = 0; s < total_elements; s++){
            for (size_t i = 0; i < nFacilities; i++) {
                //RHS
                IloExpr constr2(SFLP);
                for(size_t j = 0; j < nClients; j++){
                    constr2 += master_entities.y[i][j][s];
                }
                //RHS
                constr2 -= facil_capacities[i]*master_entities.x[i];
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
        for (size_t i = 0; i < nFacilities; i++){
            constr3 += facil_capacities[i]*master_entities.x[i];
        }
        //the constraint
        master_entities.Feasibility.add(constr3 >= max_demand);
        master.add(master_entities.Feasibility);
        constr3.end();
        
    } catch (IloException& e) {
        cerr << "Concert Exception: " << e << endl;
    } catch (...) {
        cerr << "Other Exception" << endl;
    }
}

void SFLP_GAPM::MasterProblemSolution (double &LB, const double &TL)
{
    //Setting up cplex
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
			cout << "Facility " << (i+1) << " will be open" << endl;
	}

    LB = cplex_master.getObjValue();
        
}

void SFLP_GAPM::SPProblemCreation (){
    try {
        //Create variables
        //Demand distribution
        sp_entities.y = IloNumVarArray2(SFLP_sp, nFacilities);
        //Populate OF
        sp_entities.objective = IloExpr(SFLP_sp);
        for (size_t i = 0; i < nFacilities; i++){
            sp_entities.y[i] = IloNumVarArray(SFLP_sp, nClients);
            for (size_t j = 0; j < nClients; j++) {
                sp_entities.y[i][j] = IloNumVar(SFLP_sp, 0.0, DBL_MAX, ILOFLOAT);
                subprob.add(sp_entities.y[i][j]);
                sp_entities.objective += dist_costs[i][j] * sp_entities.y[i][j];
            }
        }

        //Add objective function
        subprob.add(IloMinimize(SFLP, master_entities.objective));

        //Constraints
        //Demand satisfaction
        sp_entities.dem_constraints = IloRangeArray(SFLP_sp);
        // Facilities capacity
        sp_entities.cap_constraints = IloRangeArray(SFLP_sp);

        for (size_t j = 0; j < nClients; j++){
            //LHS
            IloExpr constr1(SFLP);
            for (size_t i = 0; i < nFacilities; i++){
                constr1 += sp_entities.y[i][j];
            }
            //Adding the constraint
            sp_entities.dem_constraints.add(constr1 >= stoch_param[j][0]);
            //remove linear expression
            constr1. end();
        }
        subprob.add(sp_entities.dem_constraints);

        for (size_t i = 0; i < nFacilities; i++) {
            //RHS
            IloExpr constr2(SFLP);
            for(size_t j = 0; j < nClients; j++){
                constr2 += sp_entities.y[i][j];
            }
            //RHS
            //Adding the constraint
            sp_entities.cap_constraints.add(constr2 <= x_bar[i] * facil_capacities[i]);
            //Remove the linear expression
            constr2.end();
        }
        subprob.add(sp_entities.cap_constraints);

    } catch (IloException& e) {
        cerr << "Concert Exception: " << e << endl;
    } catch (...){
        cerr << "Other Exception " << endl;
    }
}

//Modify the subproblem to solve a different scenario
void SFLP_GAPM::SPProblemModification(vector<size_t> &element, const bool& mod_x) {
    // s is used to modify RHS
    //modifying demands
    for (size_t j = 0; j < nClients; j++){
        vector<double> elem_demand = expected_demand(element);
        sp_entities.dem_constraints[j].setLB(elem_demand[j]);
    }

    //Modifying composition of facilities built
    if (mod_x == true) {
        for (size_t i = 0; i < nFacilities; i++) {
            sp_entities.cap_constraints[i].setUB(x_bar[i] * facil_capacities[i]);
        }
    }
}

void SFLP_GAPM::SPProbleSolution(vector<double> &stoch, vector<double> &lambda, double &obj, vector<size_t> &element) {
    //Setting up cplex
    IloCplex cplex_sp(SFLP_sp);
    cplex_sp.setParam(IloCplex::RootAlg, IloCplex::Primal);
	cplex_sp.setParam(IloCplex::PreInd, 0);
	cplex_sp.solve();

    if (cplex_sp.getStatus() == IloAlgorithm::Optimal){
        for (size_t j = 0; j < nClients; j++) {
            vector<double> elem_demand = expected_demand(element);
            stoch.push_back(elem_demand[j]);
            lambda.push_back(cplex_sp.getDual(sp_entities.dem_constraints[j]));
        }
    } else {
        cerr << "Optimal solution of the problem was not found!" << endl;
    }
}


double SFLP_GAPM::compute_UB(){
    double UB = 0.0;
    for (size_t i = 0; i < nFacilities; i++){
        UB += x_bar[i]*fixed_costs[i];
    }

    for (size_t s = 0; s < nScenarios; s++){
        UB += probability[s] * obj_opt_sps[s];
    }
    return UB;
}