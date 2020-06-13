#include"inc/GAPM.h"
#include"inc/SFLP_GAPM.h"
#include<fstream>
using namespace std;

int main(){

    string file = "cap111.txt";  
    string file_stoch = "cap111_50_50_0.10_250.txt";  
    
    //Creating instance of the problem based on previous files
    SFLP_GAPM Instance(file, file_stoch);

    //Start the GAPM
    GAPM algor;
    algor.gapm_algorithm(Instance);

}