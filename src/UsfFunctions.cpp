// Created by CDRPico
// 11/06/2020 14:06

#include"../inc/UsfFunctions.h"
#include"../inc/SFLP_GAPM.h"
#include"../inc/SFCMFP_GAPM.h"
#include"../inc/BendersAPM_SFLP.h"
#include"../inc/OuterBendersSFLP.h"
#include"../inc/BendersAPM_CP.h"
#include<algorithm>

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
template void read_instance_data(SFCMFP_GAPM &ProblemInstance, string &inst_name, string &stoch_inst_name);

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

bool compareSolutions(vector<double> x_prev, vector<double> cur_x) {
	bool r = true;
	for (size_t i = 0; i < x_prev.size(); i++) {
		if (fabs(x_prev[i] - cur_x[i]) > tolabscuts) {
			r = false;
			break;
		}
	}
	return r;
}

//function designed to refine the current partition
void disag_procedure::disaggregation(vector<solution_sps> &sp_info, vector<vector<size_t>> &partition, const double &nScenarios) {
	//a new partition will be created and will further replace the previous one
	/*partition.clear();
	partition.push_back(vector<size_t>(create_partition(nScenarios)));*/
	vector<vector<size_t>> new_partition;
	new_partition.reserve(nScenarios*2);
	
	new_partition = refine(sp_info, partition, nScenarios);
	partition.clear();
	partition = new_partition;
	new_partition.clear();
	/*cout << "the size of the partition now is " << partition.size() << endl;
	for (size_t i = 0; i < partition.size(); i++) {
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
	outpart.reserve(nScenarios*2);
	for (size_t s = 0; s < nScenarios; s++) {
		sp_info[s].taxicab_norm = taxicab(sp_info[s].lambda);
		sp_info[s].eucl_norm = norm2(sp_info[s].lambda);
	}
	
	for (size_t i = 0; i < partition_size; i++) {
		//Get the refination of element k into the partition
		vector<vector<size_t>> partial = refine_element(sp_info, partition[i], nScenarios);
		if (i == 0) {
			outpart = partial;
		}
		else {
			for (size_t j = 0; j < partial.size(); j++) {
				outpart.resize(outpart.size() + 1);
				outpart[outpart.size() - 1] = partial[j];
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
	// vector to copy the current scenarios into the element
	vector<solution_sps> elepart(el_size);
	for (size_t i = 0; i < el_size; i++) {
		//cout << element[i] << " - " << sp_info[element[i]].scen.size() << endl;
		elepart[i] = sp_info[element[i]];
	}
	//Vector to be returned
	vector<vector<size_t>> new_element;
	new_element.reserve(el_size*2);
	
	sort(elepart.begin(), elepart.end(), solution_sps::compsortDuals);
	
	vector<vector<size_t>> new_comp = check_equal_scen(elepart);
	
	new_element.resize(new_comp.size());
	for (size_t j = 0; j < new_element.size(); j++) {
		new_element[j].resize(new_comp[j][1]+1);
		//cout << elepart[new_comp[j][0]].scen.size() << endl;
		new_element[j][0] = elepart[new_comp[j][0]].scen[0];
		for (size_t s = 1; s < new_element[j].size(); s++) {
			new_element[j][s] = elepart[new_comp[j][0] + s].scen[0];
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

vector<double> disag_procedure::compute_part_prob(vector<vector<size_t>> &partition, const double &nScenarios) {
	vector<double> part_prob;
	part_prob.resize(partition.size());
	for (size_t p = 0; p < partition.size(); p++) {
		part_prob[p] = partition[p].size() / nScenarios;
	}
	return part_prob;
}

//Funcion para comparar y ordenar escenarios basado en la matriz de diferencia duales
bool solution_sps::compsortDuals(const solution_sps &s1, const solution_sps &s2) {
	return s1.lambda < s2.lambda;
}

// Method tocheck equality between pairwise scenarios once the dual class is sort
vector<vector<size_t>> disag_procedure::check_equal_scen(vector<solution_sps> &dif_dual){
	size_t scen = dif_dual.size();
	size_t comm = dif_dual[0].lambda.size();
	vector<double> cref_duales = dif_dual[0].lambda;
	size_t cref = 0;

	vector<vector<size_t>> r_part;
	r_part.reserve(scen*2);
	
	size_t cont = 0;
	for (size_t i = 1; i < scen; i++) {
		if (compareDuals_vector(dif_dual[cref], dif_dual[i])) {
			cont++;
		} else {
			r_part.resize(r_part.size()+1);
			r_part[r_part.size()-1] = {cref, cont};
			cref = i;
			cont = 0;
		}
	}

	r_part.resize(r_part.size()+1);
	r_part[r_part.size()-1] = {cref, cont};
	return r_part;
}

disag_procedure::~disag_procedure() {};

bool compareDuals_vector(const solution_sps &s1, const solution_sps &s2) {
	//Para recorrer nodos y commodities
	//Cantidad de nodos y commodities
	size_t Comm = s1.lambda.size();
	//Comparando duales
	if (fabs(s1.taxicab_norm - s2.taxicab_norm) > tol_diff_disag || fabs(s1.eucl_norm - s2.eucl_norm) > tol_diff_disag) {
		return false;
	} else {
		for (size_t k = 0; k < Comm; k++) {
			if (fabs(s1.lambda[k]- s2.lambda[k]) > tol_diff_disag) {
				return false;
			}
		}
	}
	return true;
}

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

	if (algo == 's' || algo == 'm') {
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
template void AddVarsMaster<BendersCP>(BendersCP &BendersProb, const char &algo);

template<class T>
void AddVarsMaster(T &BendersProb, const char &algo, const char &linear) {
	size_t nscen = 0;
	if (algo == 's' || algo == 'r') {
		nscen = 1;
	}
	else {
		nscen = BendersProb.nScenarios;
	}

	BendersProb.avg_scenarios.resize(BendersProb.nClients);
	BendersProb.avg_scenarios = BendersProb.expected_demand(BendersProb.partition[0]);

	if (algo == 's' || algo == 'm') {
		BendersProb.partition.clear();
		for (size_t s = 0; s < BendersProb.nScenarios; s++) {
			BendersProb.partition.push_back(vector<size_t>({ s }));
		}
	}

	//Declaring variables
	BendersProb.ent_sflp.x = IloNumVarArray(BendersProb.Mast_Bend, BendersProb.nFacilities);
	BendersProb.ent_sflp.objective = IloExpr(BendersProb.Mast_Bend);
	for (size_t i = 0; i < BendersProb.nFacilities; i++) {
		if (linear != 'l') {
			BendersProb.ent_sflp.x[i] = IloNumVar(BendersProb.Mast_Bend, 0.0, 1.0, ILOINT);
		}
		else {
			BendersProb.ent_sflp.x[i] = IloNumVar(BendersProb.Mast_Bend, 0.0, DBL_MAX, ILOFLOAT);
		}
		
		BendersProb.Mast_mod.add(BendersProb.ent_sflp.x[i]);
		BendersProb.ent_sflp.objective += BendersProb.fixed_costs[i] * BendersProb.ent_sflp.x[i];
	}

	BendersProb.ent_sflp.theta = IloNumVarArray(BendersProb.Mast_Bend, nscen);
	for (size_t s = 0; s < nscen; s++) {
		BendersProb.ent_sflp.theta[s] = IloNumVar(BendersProb.Mast_Bend, -DBL_MAX , DBL_MAX, ILOFLOAT);
		BendersProb.Mast_mod.add(BendersProb.ent_sflp.theta[s]);
		if (nscen != 1) {
			BendersProb.ent_sflp.objective += BendersProb.probability[s] * BendersProb.ent_sflp.theta[s];
		}
	}
	if (algo == 's' || algo == 'r') {
		BendersProb.ent_sflp.objective += BendersProb.ent_sflp.theta[0];
	}


	//Add objective function
	BendersProb.Mast_mod.add(IloMinimize(BendersProb.Mast_Bend, BendersProb.ent_sflp.objective));
}
template void AddVarsMaster<BendersCP>(BendersCP &BendersProb, const char &algo, const char &linear);

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


//Delete mant elements of a vector at once
void DeleteAll(vector<size_t>& data, const vector<size_t>& deleteIndices)
{
	vector<bool> markedElements(data.size(), false);
	vector<size_t> tempBuffer;
	tempBuffer.reserve(data.size() - deleteIndices.size());

	for (vector<size_t>::const_iterator itDel = deleteIndices.begin(); itDel != deleteIndices.end(); itDel++)
		markedElements[*itDel] = true;

	for (size_t i = 0; i < data.size(); i++)
	{
		if (!markedElements[i])
			tempBuffer.push_back(data[i]);
	}
	data = tempBuffer;
}


bool is_integer(float k)
{
	if (fabs(floor(k) - k) < 1e-6 || fabs(ceil(k) - k) < 1e-6) {
		return true;
	}
	else {
		return false;
	}
}

void remove_duplicates(vector<size_t> &v)
{
	vector<size_t>::iterator itr = v.begin();
	unordered_set<size_t> s;

	for (auto curr = v.begin(); curr != v.end(); ++curr) {
		if (s.insert(*curr).second)
			*itr++ = *curr;
	}

	v.erase(itr, v.end());
}