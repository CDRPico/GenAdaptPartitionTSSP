// Created by CDRPico
// 09/07/2020 05:48

#include"../inc/InstancesSFCMFP.h"
#include"../inc/UsfFunctions.h"
#include"UsefulFunctions.h"

void InstanceSFCMFP::read_instance(ifstream &file) {
    string line;
	//First line is useless
	getline(file, line);
	//Main data of the problem, number of nodes, arcs and commodities
	getline(file, line);
	vector<string> row_values;
	split(line, row_values, ' ');
	//%zu specifies the format for reading size_t
	sscanf(row_values[0].c_str(), "%zu", &numNodes);
	sscanf(row_values[1].c_str(), "%zu", &numArcs);
	sscanf(row_values[2].c_str(), "%zu", &numComm);
	//Print the number of nodes, arcs and commodities to check if reading is correct
	cout << numNodes << numArcs << numComm;
	row_values.clear();

	//Initilize sizing of vector in the instances based on number of nodes, arcs and commodities
	//Origin node of every arc
	OArcs.resize(numArcs);
	//Destination node of every arc
	DArcs.resize(numArcs);
	//Variable routing cost for every arc
	Vij.resize(numArcs);
	//Fixed costs associated to open a certain arc
	Cij.resize(numArcs);
	//Arcs capacity
	CAPij.resize(numArcs);
	size_t carcs = 0;
	//Reading until there is no line containing arcs information
	while (carcs < numArcs)
	{
		getline(file, line);
		split(line, row_values, ' ');

		sscanf(row_values[0].c_str(), "%zu", &OArcs[carcs]);
		sscanf(row_values[1].c_str(), "%zu", &DArcs[carcs]);
		//atof convert strings to double
		Vij[carcs] = atof(row_values[2].c_str());
		CAPij[carcs] = atof(row_values[3].c_str());
		Cij[carcs] = atof(row_values[4].c_str());
		carcs++;
		row_values.clear();
	}

	//Resizing those vector which depend on the number of commodities
	//Origin nodes for each commodity
	DemFrom.resize(numComm);
	//Destination nodes for commodities
	DemTo.resize(numComm);
	size_t cdem = 0;
	while (cdem < numComm) {
		getline(file, line);
		split(line, row_values, ' ');

		sscanf(row_values[0].c_str(), "%zu", &DemFrom[cdem]);
		sscanf(row_values[1].c_str(), "%zu", &DemTo[cdem]);
		cdem++;
		row_values.clear();
	}
}

void InstanceSFCMFP::read_stoch_data(ifstream &fileScen) {
    string line;
    //Now the file which contains information related to scenarios is read
	getline(fileScen, line);
    vector<string> row_values;
	split(line, row_values, ' ');
	//First line specifies the total number of scenarios
	sscanf(row_values[0].c_str(), "%zu", &numScen);
	//Resize vectors which depend on the number of scenarios
	//Scenarios probability
	probability.resize(numScen);
	//Demands is a 2d array, scenarios are rows
	DemandS.resize(numScen);
	row_values.clear();

	size_t cscen = 0;
	while (cscen < numScen)
	{
		getline(fileScen, line);
		vector<string> row_values;
		split(line, row_values, ' ');
		//Demands columns are commodities
		DemandS[cscen].resize(numComm);
		probability[cscen] = atof(row_values[0].c_str());
		for (size_t k = 0; k < numComm; k++)
		{
			DemandS[cscen][k] = atof(row_values[k + 1].c_str());
		}
		cscen++;
		row_values.clear();
	}
}


//Incoming arcs (incoming vector) to each node (posor indicates position to read in incoming vector)
void InstanceSFCMFP::orig() {
	incoming.resize(numArcs);
	index_incoming_arcs.resize(numArcs);
	posor.resize(numNodes + 1);
	size_t count = 0;
	for (size_t i = 0; i < numNodes; i++) {
		posor[i] = count;
		for (size_t j = 0; j < numArcs; j++) {
			if (DArcs[j] == i + 1) {
				incoming[count] = OArcs[j];
				index_incoming_arcs[count] = j;
				count++;
			}
		}
	}
	posor[numNodes] = numArcs;
}

//Outgoing arcs (outgoing vector) from each node (posde indicates position to read in outgoing vector)
void InstanceSFCMFP::dest() {
	outgoing.resize(numArcs);
	index_outgoing_arcs.resize(numArcs);
	posde.resize(numNodes + 1);
	size_t count = 0;
	for (size_t i = 0; i < numNodes; i++) {
		posde[i] = count;
		for (size_t j = 0; j < numArcs; j++) {
			if (OArcs[j] == i + 1) {
				outgoing[count] = DArcs[j];
				index_outgoing_arcs[count] = j;
				count++;
			}
		}
	}
	posde[numNodes] = numArcs;
}
