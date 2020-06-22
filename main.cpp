#include"inc/GAPM.h"
#include"inc/SFLP_GAPM.h"
#include"inc/BendersAPM_SFLP.h"
#include<fstream>
using namespace std;

int main(int argc, char **argv) {

	//string file = argv[1];
	string file = "cap111.txt";
	//string file_stoch = argv[2];
	string file_stoch = "cap111_50_50_0.30_250.txt";

	//Creating instance of the problem based on previous files
	//SFLP_GAPM Instance(file, file_stoch);
	BendersSFLP Instance(file, file_stoch);

	//Start the GAPM
	//GAPM algor;
	//algor.gapm_algorithm(Instance);
	Instance.CreateMaster();
	cin.get();
}