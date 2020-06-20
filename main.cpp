#include"inc/GAPM.h"
#include"inc/SFLP_GAPM.h"
#include<fstream>
using namespace std;

int main(int argc, char **argv) {

	string file = argv[1];
	string file_stoch = argv[2];

	//Creating instance of the problem based on previous files
	SFLP_GAPM Instance(file, file_stoch);

	//Start the GAPM
	GAPM algor;
	algor.gapm_algorithm(Instance);
}