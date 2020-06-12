// Created by CDRPico
// 11/06/2020 13:05

#include"../inc/GAPM.h"
#include"../inc/SFLP_GAPM.h"

//This method makes an iteration of the algorithm for a given Master and subproblem
template <typename T>
void GAPM::body_gapm(T &ProblemInstance){
    //Here we need to create the master problem
    // Notice, this need s to be created after each iteration because we change both
    // number of constraints and variables, cplex does not have some of these features
    ProblemInstance.MasterProblemCreation();

    //Solve the master model
    ProblemInstance.MasterProblemSolution(LB, timelimit);

    //Solve every subproblem
    for (size_t s = 0; s < ProblemInstance.nScenarios; s++){
        if (s == 0) {
            ProblemInstance.SPProblemModification(vector<site_t> a(1, s), true);
        } else {
            ProblemInstance.SPProblemModification(vector<site_t> a(1, s));
        }
        double obj = 0.0;
        ProblemInstance.SPProbleSolution(stoch[s], lambda[s], obj, vector<site_t> a(1, s));
        ProblemInstance.obj_opt_sps[s] = obj;
    }

    //Compute the upper bound
    ProblemInstance.compute_UB();

    //compute GAP
    //if (GAP > threshold){
        //run disaggregation procedure locally

        //run checking stopping criteria
        //here subproblems must be solved for each element into the partition
        //function of rpOrblemInstance
        //the stopping criteria should be checked locally
    //}


}
template void GAPM::body_gapm(SFLP_GAPM &ProblemInstance);


//function designed to refine the current partition
void GAPM::disaggregation(const double &nScenarios, vector<vector<double>> &lambda){
    //a new partition will be created and will further replace the previous one
    vector<vector<size_t>> new_partition = refine(nScenarios, lambda);
    partition.clear();
    partition = new_partition;
    new_partition.clear();
}

vector<vector<size_t>> GAPM::refine(const double &nScenarios, vector<vector<double>> &lambda) {
    size_t partition_size = partition.size();
    vector<vector<size_t>> outpart;
    for (size_t k = 0; k < partition_size; k++) {
        //Get the refination of element k into the partition
        vector<vector<size_t>> partial = refine_element(partition[k], nScenarios, lambda);
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
                    if (lb == 0){
                        outpart.push_back(vector<size_t>( { partial[la][lb] } ));
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
vector<vector<size_t>> GAPM::refine_element(vector<size_t> &element, const double &nScenarios, vector<vector<double>> &lambda) {
    //Element size
    size_t el_size = element.size();
    //Vector to be returned
    vector<vector<size_t>> new_element;
    vector<double> taxicab_computation;
    vector<double> normtwo_computation;
    //Iterate through the element to compare the duals of the subproblems
    for (size_t la = 0; la < el_size; la++) {
        taxicab_computation.push_back(taxicab(lambda[la]));
        normtwo_computation.push_back(norm2(lambda[la]));
        if (la == 0) {
            new_element.push_back(vector<size_t>( { element[la] } ));
        }
        else {
            size_t tamano = new_element.size();
            //Create a new loop to compare current values regarding the previous already computed
            for (size_t lb = 0; lb < tamano; lb++) {
                bool signal = true;
                if (fabs(taxicab_computation[la]-taxicab_computation[lb]) < tol_diff_disag && 
                    fabs(normtwo_computation[la]-normtwo_computation[lb]) < tol_diff_disag) 
                {
                    signal = compare_duals(lambda, element[la], element[lb]);
                }
                // If there is no difference between lambdas, this scenario can be added to the current new element position
                if (signal == false) {
                    new_element[lb].push_back(element[la]);
                    break;
                }
                //otherwise, continue checking until the end of the vector
                //if no coincidence was found, another element is added
                else {
                    if (lb == (tamano-1)) {
                        new_element.push_back(vector<size_t>( { element[la] } ));
                    }
                }
            }
        }
    }
    return new_element;
}

bool GAPM::compare_duals (const vector<vector<double>> &lambda, const size_t &s1, const size_t &s2) {
    bool signal = false;
    //Check equality between elements
    for (size_t j = 0; lambda[s1].size(); j++) {
        if ( fabs( lambda[s1][j] - lambda[s2][j] ) > tol_diff_disag ) {
            signal = true;
            break;
        }
    }
    return signal;
}