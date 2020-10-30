// Created by CDRPico
// 09/07/2020 05:48

#pragma once

#ifndef INSTANCESFCMFP_H
#define INSTANCESFCMFP_H

using namespace std;
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>

class InstanceSFCMFP {
    public:
        
        // Sets size
        size_t numNodes = 0;
	    size_t numArcs = 0;
        size_t numComm = 0;
        size_t numScen = 0;

        // Deterministic data
        vector<size_t> OArcs;
        vector<size_t> DArcs;
        vector<size_t> DemFrom;
        vector<size_t> DemTo;
        vector<double> Vij, Cij, CAPij;
        vector<vector<double>> DemandS;
        vector<size_t> incoming;
        vector<size_t> posor;
        vector<size_t> index_incoming_arcs;
        vector<size_t> outgoing;
        vector<size_t> posde;
        vector<size_t> index_outgoing_arcs;

        //Stochastic entities
        vector<double> mean_demand;
        vector<vector<double>> stoch_param;
        vector<double> probability;
        
        //Solutions
        vector<double> x_bar;
        vector<vector<vector<double>>> y_bar;

        InstanceSFCMFP() = default;   

        //Given the file containing the info, we read the data
        void read_instance(ifstream &file);

        //This functions reads stochastic instance already generated (stoch demand)
        //The data is stored in stoch_param (it is used either to read or generate instance)
        void read_stoch_data(ifstream &file);
 
        //Incoming arcs (incoming vector) to each node (posor indicates position to read in incoming vector)
        void orig();

        //Outgoing arcs (outgoing vector) from each node (posde indicates position to read in outgoing vector)
        void dest(); 

};

#endif // INSTANCESFCMFP_H