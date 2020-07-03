// Created by CDRPico
// 11/06/2020 14:06

#include"../inc/UsefulFunctions.h"
#include"../inc/SFLP_GAPM.h"
#include"../inc/BendersAPM_SFLP.h"
#include"../inc/OuterBendersSFLP.h"

MyClock::MyClock(void) {
	start = sc.now();
	end = sc.now();
}


template <typename T>
void read_instance_data(T &ProblemInstance, string &inst_name, string &stoch_inst_name){
    const char* fln = inst_name.c_str();
    ifstream file(fln);

    //read the general data of the instance
    ProblemInstance.read_instance(file);

    const char* fln_stoch = stoch_inst_name.c_str();
    ifstream file_stoch(fln_stoch);

    //read the scenarios info
    ProblemInstance.read_stoch_data(file_stoch);
}
template void read_instance_data(SFLP_GAPM &ProblemInstance, string &inst_name, string &stoch_inst_name);

//suitable function to create the initial partition
vector<size_t> create_partition(const size_t &Scenarios) {
	vector<size_t> partition;
	for (size_t i = 0; i < Scenarios; i++) {
		partition.push_back(i);
	}
	return partition;
}


//taxicab norm
double taxicab(vector<double> &vect) {
	double taxicab = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		taxicab += fabs(vect[i]);
	}
	return taxicab;
}

//norm-2
double norm2(vector<double> &vect) {
	double norm = 0;
	for (size_t i = 0; i < vect.size(); i++) {
		norm += pow(vect[i],2);
	}
	return norm;
}

double max_demand(vector<vector<double>> &stoch_dem) {
	vector<double> dem(stoch_dem.size());
	double mdem = 0.0;
	for (size_t s = 0; s < stoch_dem[0].size(); s++){
		double scenario_tot_demand = 0.0;
		for (size_t j = 0; j < stoch_dem.size(); j++) {
			scenario_tot_demand += stoch_dem[j][s];
		}
		if (scenario_tot_demand > mdem) {
			mdem = scenario_tot_demand;
		}
	}
	return mdem;
}

bool smallerWithTolerance(double smaller, double larger)
{
	if (smaller < larger - tolabscuts)
	{
		if (larger - smaller > tolrelcuts * fabs(larger))
		{
			return true;
		}
	}
	return false;
}

//function designed to refine the current partition
void disag_procedure::disaggregation(vector<solution_sps> &sp_info, vector<vector<size_t>> &partition, const double &nScenarios) {
	//a new partition will be created and will further replace the previous one
	vector<vector<size_t>> new_partition = refine(sp_info, partition, nScenarios);
	partition.clear();
	partition = new_partition;
	new_partition.clear();
	cout << "the size of the partition now is " << partition.size() << endl;
	/*for (size_t i = 0; i < partition.size(); i++) {
		for (size_t j = 0; j < partition[i].size(); j++) {
			cout << partition[i][j] << ' ';
		}
		cout << endl;
		//cin.get();
	}*/
}

vector<vector<size_t>> disag_procedure::refine(vector<solution_sps> &sp_info, vector<vector<size_t>> &partition, const double &nScenarios) {
	size_t partition_size = partition.size();
	vector<vector<size_t>> outpart;
	for (size_t k = 0; k < partition_size; k++) {
		//Get the refination of element k into the partition
		vector<vector<size_t>> partial = refine_element(sp_info, partition[k], nScenarios);
		//Insert in outpart the vector obtained above
		size_t tamano = outpart.size();
		// if k = 0 outpart is empty and will be initialized with partial
		if (k == 0) {
			outpart = partial;
		}
		// otherwise the data from partial will be pushed back one scenario at a time
		else {
			for (size_t la = 0; la < partial.size(); la++) {
				for (size_t lb = 0; lb < partial[la].size(); lb++) {
					if (lb == 0) {
						outpart.push_back(vector<size_t>({ partial[la][lb] }));
					}
					else {
						outpart[la + tamano].push_back(partial[la][lb]);
					}
				}
			}
			partial.clear();
		}
	}
	return outpart;
}

//function to refine one element into the partition
vector<vector<size_t>> disag_procedure::refine_element(vector<solution_sps> &sp_info, vector<size_t> &element, const double &nScenarios) {
	//Element size
	size_t el_size = element.size();
	//Vector to be returned
	vector<vector<size_t>> new_element;
	vector<double> taxicab_computation;
	vector<double> normtwo_computation;
	//Iterate through the element to compare the duals of the subproblems
	for (size_t la = 0; la < el_size; la++) {
		taxicab_computation.push_back(taxicab(sp_info[la].lambda));
		normtwo_computation.push_back(norm2(sp_info[la].lambda));
		if (la == 0) {
			new_element.push_back(vector<size_t>({ element[la] }));
		}
		else {
			size_t tamano = new_element.size();
			//Create a new loop to compare current values regarding the previous already computed
			for (size_t lb = 0; lb < tamano; lb++) {
				if (sp_info[la].F == sp_info[new_element[lb].back()].F) {
					bool signal = true;
					if (fabs(taxicab_computation[la] - taxicab_computation[lb]) < tol_diff_disag &&
						fabs(normtwo_computation[la] - normtwo_computation[lb]) < tol_diff_disag)
					{
						signal = compare_duals(sp_info[element[la]], sp_info[element[lb]], element[la], element[lb]);
					}
					// If there is no difference between lambdas, this scenario can be added to the current new element position
					if (signal == false) {
						new_element[lb].push_back(element[la]);
						break;
					}
					//otherwise, continue checking until the end of the vector
					//if no coincidence was found, another element is added
					else {
						if (lb == (tamano - 1)) {
							new_element.push_back(vector<size_t>({ element[la] }));
						}
					}
				}
				else if (lb == (tamano - 1)) {
					new_element.push_back(vector<size_t>({ element[la] }));
				}
			}
		}
	}
	return new_element;
}

//function to campare scenarios pairwise (their duals)
bool disag_procedure::compare_duals(solution_sps &sp_info_s1, solution_sps &sp_info_s2, const size_t &s1, const size_t &s2) {
	bool signal = false;
	//Check equality between elements
	for (size_t j = 0; j < sp_info_s1.lambda.size(); j++) {
		if (fabs(sp_info_s1.lambda[j] - sp_info_s2.lambda[j]) > tol_diff_disag) {
			signal = true;
			break;
		}
	}
	return signal;
}

disag_procedure::~disag_procedure() {};


template<class T>
void AddVarsMaster(T &BendersProb, const char &algo) {
	size_t nscen = 0;
	if (algo == 's') {
		nscen = 1;
	}
	else {
		nscen = BendersProb.nScenarios;
	}

	BendersProb.avg_scenarios.resize(BendersProb.nClients);
	BendersProb.avg_scenarios = BendersProb.expected_demand(BendersProb.partition[0]);

	if (algo != 'a') {
		BendersProb.partition.clear();
		for (size_t s = 0; s < BendersProb.nScenarios; s++) {
			BendersProb.partition.push_back(vector<size_t>({ s }));
		}
	}

	//Declaring variables
	BendersProb.ent_sflp.x = IloNumVarArray(BendersProb.Mast_Bend, BendersProb.nFacilities);
	BendersProb.ent_sflp.objective = IloExpr(BendersProb.Mast_Bend);
	for (size_t i = 0; i < BendersProb.nFacilities; i++) {
		BendersProb.ent_sflp.x[i] = IloNumVar(BendersProb.Mast_Bend, 0.0, 1.0, ILOINT);
		BendersProb.Mast_mod.add(BendersProb.ent_sflp.x[i]);
		BendersProb.ent_sflp.objective += BendersProb.fixed_costs[i] * BendersProb.ent_sflp.x[i];
	}

	BendersProb.ent_sflp.theta = IloNumVarArray(BendersProb.Mast_Bend, nscen);
	for (size_t s = 0; s < nscen; s++) {
		BendersProb.ent_sflp.theta[s] = IloNumVar(BendersProb.Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
		BendersProb.Mast_mod.add(BendersProb.ent_sflp.theta[s]);
		if (nscen != 1) {
			BendersProb.ent_sflp.objective += BendersProb.probability[s] * BendersProb.ent_sflp.theta[s];
		}
	}
	if (algo == 's') {
		BendersProb.ent_sflp.objective += BendersProb.ent_sflp.theta[0];
	}


	//Add objective function
	BendersProb.Mast_mod.add(IloMinimize(BendersProb.Mast_Bend, BendersProb.ent_sflp.objective));
}
template void AddVarsMaster<BendersSFLP>(BendersSFLP &BendersProb, const char &algo);
template void AddVarsMaster<OuterBendersSFLP>(OuterBendersSFLP &BendersProb, const char &algo);

//Add constraints for the average scenario
template<class T>
void ValidInequalities(T &BendersProb, const char &algo) {
	BendersProb.ent_sflp.y = IloNumVarArray2(BendersProb.Mast_Bend, BendersProb.nFacilities);
	for (size_t i = 0; i < BendersProb.nFacilities; i++) {
		BendersProb.ent_sflp.y[i] = IloNumVarArray(BendersProb.Mast_Bend, BendersProb.nClients);
		for (size_t j = 0; j < BendersProb.nClients; j++) {
			BendersProb.ent_sflp.y[i][j] = IloNumVar(BendersProb.Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
			BendersProb.Mast_mod.add(BendersProb.ent_sflp.y[i][j]);
		}
	}
	//Adding constraint to average demand
	IloRangeArray avg_dem_cons(BendersProb.Mast_Bend);
	for (size_t j = 0; j < BendersProb.nClients; j++) {
		//LHS
		IloExpr constr1(BendersProb.Mast_Bend);
		for (size_t i = 0; i < BendersProb.nFacilities; i++) {
			constr1 += BendersProb.ent_sflp.y[i][j];
		}
		//Adding the constraint
		avg_dem_cons.add(constr1 >= BendersProb.avg_scenarios[j]);
		//remove linear expression
		constr1.end();
	}
	BendersProb.Mast_mod.add(avg_dem_cons);

	IloRangeArray avg_cap_cons(BendersProb.Mast_Bend);
	for (size_t i = 0; i < BendersProb.nFacilities; i++) {
		//RHS
		IloExpr constr2(BendersProb.Mast_Bend);
		for (size_t j = 0; j < BendersProb.nClients; j++) {
			constr2 += BendersProb.ent_sflp.y[i][j];
		}
		//RHS
		//Adding the constraint
		constr2 -= BendersProb.ent_sflp.x[i] * BendersProb.facil_capacities[i];
		avg_cap_cons.add(constr2 <= 0);
		//Remove the linear expression
		constr2.end();
	}
	BendersProb.Mast_mod.add(avg_cap_cons);

	IloExpr Avg_cons(BendersProb.Mast_Bend);
	IloRangeArray AVG(BendersProb.Mast_Bend);

	for (size_t j = 0; j < BendersProb.nClients; j++) {
		for (size_t i = 0; i < BendersProb.nFacilities; i++) {
			Avg_cons -= BendersProb.dist_costs[i][j] * BendersProb.ent_sflp.y[i][j];
		}
	}
	if (algo == 's') {
		Avg_cons += BendersProb.ent_sflp.theta[0];
	}
	else {
		for (size_t s = 0; s < BendersProb.nScenarios; s++) {
			Avg_cons += BendersProb.probability[s] * BendersProb.ent_sflp.theta[s];
		}
	}
	AVG.add(Avg_cons >= 0);
	Avg_cons.end();
	BendersProb.Mast_mod.add(AVG);
}
template void ValidInequalities<BendersSFLP>(BendersSFLP &BendersProb, const char &algo);
template void ValidInequalities<OuterBendersSFLP>(OuterBendersSFLP &BendersProb, const char &algo);


template<class T>
void FeasibilityConstraint(T &BendersProb) {
	BendersProb.ent_sflp.Feasibility = IloRangeArray(BendersProb.Mast_Bend);
	//Feasibility constraint
		//LHS
	IloExpr constr3(BendersProb.Mast_Bend);
	for (size_t i = 0; i < BendersProb.nFacilities; i++) {
		constr3 += BendersProb.facil_capacities[i] * BendersProb.ent_sflp.x[i];
	}
	//the constraint
	BendersProb.ent_sflp.Feasibility.add(constr3 >= BendersProb.max_dem);
	BendersProb.Mast_mod.add(BendersProb.ent_sflp.Feasibility);
	constr3.end();
}
template void FeasibilityConstraint<BendersSFLP>(BendersSFLP &BendersProb);
template void FeasibilityConstraint<OuterBendersSFLP>(OuterBendersSFLP &BendersProb);