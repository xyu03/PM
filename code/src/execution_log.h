/*
存储每次出现的改进解以及他们的解决方案
*/
#pragma once

#include <vector>
#include <ctime>

using namespace std;


struct ExecutionLog {

	ExecutionLog(Graph *g, double t, double tl) : timer(t), bestSolution(g), timeLimit(tl) {}

	typedef pair<EdgeCost, double> EntryType;

	inline void AddSolution(SteinerSolution &sol) {
		// Obey the time limit.
		if (timeLimit > 0 && GetTime()-begin_time > timeLimit)
			return;

		if (solCost.empty()) {
			solCost.push_back(EntryType(sol.GetCost(), GetTime()-begin_time));
			bestSolution.CopyFrom(&sol);
		}
		else {
			if (solCost.back().first > sol.GetCost() + EDGE_COST_PRECISION) {
				solCost.push_back(EntryType(sol.GetCost(), GetTime()-begin_time));
				bestSolution.CopyFrom(&sol);
			}
		}
	}

	// These are the cost and running times for each incumbent solution.
	vector<EntryType> solCost;
	double timer = 0;
	SteinerSolution bestSolution;
	const double timeLimit;

};