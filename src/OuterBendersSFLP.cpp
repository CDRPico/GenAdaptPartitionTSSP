// Created by CDRPico
// 23/06/2020 01:14

#include"../inc/OuterBendersSFLP.h"

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
    for (size_t p = 0; p < num_old_cuts; p++) {
        if (cp[p].active == false) {
            master.remove(cp[p].cut_pool_constr);
        }
    }
    //Add new cuts
	for (size_t p = num_old_cuts; p < current_num_cuts; p++) {
        master.add(cp[p].cut_pool_constr);
    }
    //remove inactive cuts from structure
    for (auto it = cp.begin(); it < cp.begin()+num_old_cuts; it++) {
        if (it->active == false) {
            cp.erase(it);
        }
    }
    num_old_cuts = cp.size();
}

//Solving the master problem
void OuterBendersSFLP::MasterProblemSolution(IloCplex *cplex_master, IloModel &master, double &LB, const double &TL)
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
	x_bar.resize(nFacilities);
	for (size_t i = 0; i < nFacilities; i++) {
		x_bar[i] = cplex_master->getValue(master_entities.x[i]);
		if (x_bar[i] > 0)
			cout << "Facility " << (i + 1) << " will be open" << endl;
	}

	LB = cplex_master->getObjValue();

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

        if(smallerWithTolerance(rhs, lhs[p])) {
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
            string cut_name = "test_" + to_string(p+1);
            cp.emplace_back(Mast_Bend, new_cut, derecho, partition[p], cut_name, true, x_coefs, theta_coefs);
            new_cut.end();
        }
    }
}

//check active cuts
void OuterBendersSFLP::check_active_cuts(vector<vector<double>> &stoch, const vector<solution_sps>* sp_info, bool = false) {
    size_t current_num_cuts = cp.size();
    vector<double> lhs;
    lhs.resize(cp.size());
    for (size_t p = 0; p < current_num_cuts; p++) {
        //Compute the current rhs to check if the cut is violated
        double theta_ac = 0.0;
		for (size_t s = 0; s < cp[p].elem_cut.size() ; s++) {
			theta_ac += current_theta[cp[p].elem_cut[s]] * probability[cp[p].elem_cut[s]];

            //constant of current first stage multiplied by duals of capacity constarints
            double const_fs = 0.0;
            for (size_t i = 0; i < nFacilities; i++) {
                const_fs += probability[cp[p].elem_cut[s]] * facil_capacities[i] * x_bar[i] * sp_info->at(cp[p].elem_cut[s]).lambda[nClients + i];
            }

            //compute the constant part of the cut
            lhs[p] = 0.0;
            for (size_t j = 0; j < nClients; j++) {
                lhs[p] += probability[cp[p].elem_cut[s]] * stoch[cp[p].elem_cut[s]][j] * sp_info->at(cp[p].elem_cut[s]).lambda[j];
            }
            lhs[p] += const_fs;
		}
        
        if (fabs(lhs[p]-theta_ac) > tolabscuts) {
            cp[p].active = false;
        } else {
            cp[p].active = true;
        }        
    }
}