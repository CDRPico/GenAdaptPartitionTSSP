// Created by CDRPico
// 23/06/2020 01:14

#pragma once

#ifndef OuterBendersGAPM
#define OuterBendersGAPM

#include"BendersAPM_SFLP.h"

struct cut_pool{
    //CPLEX cuts preserved from one iteration to the next one
    IloRangeArray cut_pool;
    //Scenarios for a certain cut
    vector<size_t> elem_cut;
    //cut name
    string cut_name;
    //is active
    bool active;
};

class OuterBendersSFLP : public SFLP_GAPM 
{
    //entities of the master problem
    entitiesBendSFLP ent_sflp;

    //struct to store pool of cuts
    vector<cut_pool> cp;

    //Instantiate using constructor of base class
    OuterBendersSFLP(string &inst_name, string &stoch_inst_name);

    //Create master problem given a certain partition
	void MasterProblemCreation(const char &algo, bool = false);

    //Modify master problem using pool cuts
    void MasterProblemModification(const char &algo, bool = false);

    //Check active cuts
    void check_active_cuts();

    //Create master problem given a certain partition
    //Override the basis class function
	IloCplex MasterProblemCreation(bool = false);

    //Master Environment
    IloEnv Mast_OutBend = IloEnv();

	//Master model
	IloModel Mast_OutMod = IloModel(Mast_OutBend);

    //Master cplex
    IloCplex Mast_cplex = IloCplex(Mast_OutBend);
};

#endif