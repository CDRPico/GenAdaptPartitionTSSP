// Created by CDRPico
// 23/06/2020 01:14

#include"../inc/OuterBendersSFLP.h"

OuterBendersSFLP::OuterBendersSFLP(string &inst_name, string &stoch_inst_name) : SFLP_GAPM(inst_name, stoch_inst_name) {
	size_t a = 1;
}

void OuterBendersSFLP::MasterProblemCreation(const char &algo, bool removeobjs) {
    //Add variables
	AddVarsMaster(*this, algo);

    //Valid inequalities
	ValidInequalities(*this, algo);

    //Adding feasibility constraint
	FeasibilityConstraint(*this);
}

void OuterBendersSFLP::MasterProblemModification(const char &algo, bool removeobjs) {
    //Remove inactive cuts
    //Add new cuts
	size_t a = 1;
}

