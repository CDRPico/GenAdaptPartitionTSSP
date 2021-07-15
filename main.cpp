
#include"inc/GAPM.h"
#include"inc/SFLP_GAPM.h"
#include"inc/BendersAPM_SFLP.h"
#include"inc/OuterBendersSFLP.h"
#include<fstream>
#include"inc/CapPlan_GAPM.h"
#include"SMPS_ElecPlan.h"
#include"inc/BendersAPM_CP.h"

using namespace std;

int main(int argc, char **argv) {
	//string file = argv[1];	
	//string file_stoch = argv[2];
	//const char algo = argv[3][0];

	//problem, 1 for Capacity planning
	//problem, 2 for facility location
	int problem = 1;

	if (problem == 1) {
		string file_ce_time = "cp_instances/time_nospace.mps";
		string file_ce_core = "cp_instances/core_nospace.mps";
		string file_ce_stoch = "cp_instances/stoch_nospace.mps";
		//size_t nnn = 1000;
		//testin.generate_instances(nnn);
		string testname = "capplan_small_1000.txt";
		const char algo = 'b';
		if (algo == 'n') {
			//Creating instance of the problem based on previous files
			CP_GAPM Instance(file_ce_core, file_ce_stoch, file_ce_time, testname);
			Instance.gapm_type = 's';
			//Instance.partition.push_back(vector<size_t>(create_partition(Instance.nScenarios)));
			/*Instance.partition.clear();
			Instance.partition.resize(Instance.nScenarios);
			for (size_t s = 0; s < Instance.nScenarios; s++) {
				Instance.partition.push_back({ s });
			}*/
			//Start the GAPM
			GAPM algor;
			algor.gapm_algorithm(Instance, algo, Instance.stoch_constr.size());
		}
		else if (algo == 'f') {
			CP_GAPM Instance(file_ce_core, file_ce_stoch, file_ce_time, testname);
			Instance.gapm_type = 's';
			Instance.FullCP();
		}
		else {
			//Initialize benders
			BendersCP Instance(file_ce_core, file_ce_stoch, file_ce_time, testname);
			Instance.gapm_type = 's';
			vector<vector<size_t>> part;
			//part.resize(Instance.nScenarios);
			part.push_back(vector<size_t>(create_partition(Instance.nScenarios)));
			vector<double> xb(Instance.nFacilities, 0.0);
			//Create and run procedure
			double optimal = Instance.runBenders(algo, part, xb);
			optimal = 0; // Agregado para evitar warning de que no se usaba
			//Instance.IterativeMaster(algo, file, file_stoch, part);
		}
	}
	else if (problem == 2) {
		string file = "cap111.txt";
		string file_stoch = "cap111_0.10_250.txt";
		const char algo = 'n';
		if (algo == 'n') {
			//Creating instance of the problem based on previous files
			SFLP_GAPM Instance(file, file_stoch);
			Instance.gapm_type = 's';
			//Start the GAPM
			GAPM algor;
			algor.gapm_algorithm(Instance, algo, Instance.nClients);
		}
		else {
			if (algo == 'o' || algo == 'l') {
				//Creating instance of the problem based on previous files
				OuterBendersSFLP Instance(file, file_stoch);
				//Start the GAPM
				GAPM algor;
				algor.gapm_algorithm(Instance, algo, Instance.nClients);
			}
			else {
				//Initialize benders
				BendersSFLP Instance(file, file_stoch);
				Instance.gapm_type = 's';
				vector<vector<size_t>> part;
				part.push_back(vector<size_t>(create_partition(Instance.nScenarios)));
				vector<double> xb(Instance.nFacilities, 0.0);
				//Create and run procedure
				double optimal = Instance.CreateMaster(algo, part, xb);
				optimal = 0; // Agregado para evitar warning de que no se usaba
				//Instance.IterativeMaster(algo, file, file_stoch, part);
			}
		}
	}
	//cin.get();
}