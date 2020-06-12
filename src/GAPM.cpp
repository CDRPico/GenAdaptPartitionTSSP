// Created by CDRPico
// 11/06/2020 13:05

#include"../inc/GAPM.h"
#include"../inc/SFLP_GAPM.h"
#include <chrono>

template<typename T>
void GAPM::gapm_algorithm(T &ProblemInstance) {
    //initialize size of stochastic parameters and duals
    lambda.resize(ProblemInstance.nScenarios);
    stoch.resize(ProblemInstance.nScenarios);
    partition = ProblemInstance.partition;
    //Creating the clock to control the execution time
    chrono::steady_clock sc;
	auto start = sc.now();
    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>>(end - start);
    execution_time = time_span.count();
    iterations = 0;
    GAP = 1;
    UB = DBL_MAX;
    LB = 0.0;
    termination = 9;
    //main loop of the algorithm
    while (GAP > GAP_threshold && execution_time < timelimit) {
        termination = gapm_algorithm(ProblemInstance);
        //update execution time
        end = sc.now();
        auto time_span = static_cast<chrono::duration<double>>(end - start);
        execution_time = time_span.count();
        if (termination == 2) {
            goto outwhile;
        }
        //restart vectors
        lambda.clear();
        stoch.clear();
        lambda.resize(ProblemInstance.nScenarios);
        stoch.resize(ProblemInstance.nScenarios);
    }
outwhile:
    /* DISPLAY */
	cout << "CRP: "
	<< termination << " "
	<< execution_time << " "
	<< UB << " "
	<< LB << " "
	<< (UB - LB) / (1e-10 + LB) << " "
	<< 0 << " "
	<< iterations << " "
	<< partitions.size()
	<< endl;
}
template void GAPM::gapm_algorithm(SFLP_GAPM &ProblemInstance);

//This method makes an iteration of the algorithm for a given Master and subproblem
template <typename T>
void GAPM::body_gapm(T &ProblemInstance){
    //Update number of iterations 
    iterations++;
    //Here we need to create the master problem
    // Notice, this need s to be created after each iteration because we change both
    // number of constraints and variables, cplex does not have some of these features
    ProblemInstance.MasterProblemCreation();

    //Solve the master model
    ProblemInstance.MasterProblemSolution(LB, timelimit - execution_time);

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

    //Update solution Gap
    compute_gap();

    double continue_algorithm = 0;
    if (GAP > GAP_threshold){
        //run disaggregation procedure locally
        disaggregation(ProblemInstance.nscenarios);

        for (size_t s = 0; s < partition.size(); s++) {
            if (s == 0) {
                ProblemInstance.SPProblemModification(partition[s], true);
            } else {
                ProblemInstance.SPProblemModification(partition[s]);
            }
            double obj = 0.0;
            ProblemInstance.SPProbleSolution(stoch_agg[s], lambda_agg[s], obj, vector<site_t> a(1, s));
        }

        //run checking stopping criteria
        bool optimality = stopping_criteria(ProblemInstance.nScenarios, ProblemInstance.probability);
        if (optimality == true) continue_algorithm = 2;
    } else {
        continue_algorithm = 1;
    }
}
template void GAPM::body_gapm(SFLP_GAPM &ProblemInstance);


//---------------- TODO include possibility of infeasible subproblems   -------------------//

//function designed to refine the current partition
void GAPM::disaggregation(const double &nScenarios){
    //a new partition will be created and will further replace the previous one
    vector<vector<size_t>> new_partition = refine(nScenarios);
    partition.clear();
    partition = new_partition;
    new_partition.clear();
    //reset stoch params agg and dual mult agg
    stoch_agg.clear();
    stoch_agg.resize(partition.size());
    lambda_agg.clear();
    lambda_agg.resize(partition.size());
}

vector<vector<size_t>> GAPM::refine(const double &nScenarios) {
    size_t partition_size = partition.size();
    vector<vector<size_t>> outpart;
    for (size_t k = 0; k < partition_size; k++) {
        //Get the refination of element k into the partition
        vector<vector<size_t>> partial = refine_element(partition[k], nScenarios);
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
vector<vector<size_t>> GAPM::refine_element(vector<size_t> &element, const double &nScenarios) {
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
                    signal = compare_duals(element[la], element[lb]);
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

//function to campare scenarios pairwise (their duals)
bool GAPM::compare_duals (const size_t &s1, const size_t &s2) {
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


void GAPM::compute_gap() {
    if (UB >= 1e300) {
		GAP = 1;
	}
	else {
		GAP = (UB - LB) / (1e-10 + LB);
	}
}

bool GAPM::stopping_criteria(const size_t &nScenarios, vector<double> &prob) {
    bool check_op = true;
    //aggregated expectation and atomized expectations for each element into the partition are computed
    for (size_t p = 0; p < partition.size(); p++){
        double atom_expectation = 0.0;
        double element_prob = partition[p].size()/nScenarios;
        double agg_expectation = 0.0;
        for (size_t j = 0; j < lambda_agg[p].size(); j++) {
            agg_expectation += (lambda_agg[p][j] * stoch_agg[p][j]);
            for (size_t s = 0; s < partition[p].size(); s++) {
                atom_expectation += prob[partition[p][s]] / element_prob * (lambda[partition[p][s]][j] * stoch[partition[p][s]][j]);
            }
        }

        //check equality of expectations, this is done for each element, when one element does not have the same value, procces is aborted
        if (fabs(agg_expectation - atom_expectation) > expectation_threshold) {
            check_op = false;
            break;
        }
    }
    return check_op;
}