// Created by CDRPico
// 20/06/2020 00:48

#include"../inc/BendersAPM_SFLP.h"

ILOLAZYCONSTRAINTCALLBACK2(LazyOptCuts,
    BendersSFLP&, Instance,
    const string&, algo)
{
	Instance.sp_info.resize(Instance.nScenarios);
	Instance.stoch.resize(Instance.nScenarios);
    size_t scenarios = 0;
    if (algo == "s") {
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
			Instance.SPProblemModification_CPX(Instance.partition[p], true);
		}
		else {
			Instance.SPProblemModification_CPX(Instance.partition[p]);
		}
		Instance.SPProbleSolution_CPX(Instance.stoch_agg[p], &Instance.sp_info_agg[p], true);
        //Compute the current rhs to check if the cut is violated
        double element_prob = 0.0;
        double theta_ac = 0.0;
        for (size_t s = 0; s < Instance.partition[p].size(); s++) {
            element_prob += Instance.probability[Instance.partition[p][s]];
            theta_ac += current_theta[Instance.partition[p][s]] * Instance.probability[Instance.partition[p][s]];
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
        if (algo != "s") {            
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
                add(new_cut <= derecho).end();
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
    if (smallerWithTolerance(current_theta[0],single_lhs) && algo == "s"){
        violated = true;
        new_single_cut -= Instance.ent_sflp.theta[0];
        add(new_single_cut <= const_part_single).end();
    }

    if (violated == true && algo == "a") {
        for (size_t s = 0; s < Instance.nScenarios; s++) {
            vector<size_t> el(1, s);
            Instance.sp_info[s].scen = el;
            if (s == 0) {
                Instance.SPProblemModification_CPX(Instance.sp_info[s].scen, true);
            }
            else {
                Instance.SPProblemModification_CPX(Instance.sp_info[s].scen);
            }
            double obj = 0.0;
            Instance.SPProbleSolution_CPX(Instance.stoch[s], &Instance.sp_info[s]);
        }
        disag_procedure to_dis;
		to_dis.disaggregation(Instance.sp_info, Instance.partition, Instance.nScenarios);
    }
	Instance.sp_info.clear();
	Instance.stoch.clear();
	Instance.sp_info_agg.clear();
	Instance.stoch_agg.clear();
}

void BendersSFLP::CreateMaster() {
    //Declaring variables
    ent_sflp.x = IloNumVarArray(Mast_Bend, nFacilities);
	ent_sflp.objective = IloExpr(Mast_Bend);
    for (size_t i = 0; i < nFacilities; i++) {
        ent_sflp.x[i] = IloNumVar(Mast_Bend, 0.0, 1.0, ILOINT);
        Mast_mod.add(ent_sflp.x[i]);
		ent_sflp.objective += fixed_costs[i] * ent_sflp.x[i];
    }

    ent_sflp.theta = IloNumVarArray(Mast_Bend, nScenarios);
    for (size_t s = 0; s < nScenarios; s++) {
        ent_sflp.theta[s] = IloNumVar(Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
        Mast_mod.add(ent_sflp.theta[s]);
        ent_sflp.objective += ent_sflp.theta[s];
    }

    //Add objective function
	Mast_mod.add(IloMinimize(Mast_Bend, ent_sflp.objective));

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
	SPProblemCreation_CPX();

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

    string algo = "a";
    cplex.use(LazyOptCuts(this->Mast_Bend, *this, algo));

	//Solve the model
	cplex.solve();
}