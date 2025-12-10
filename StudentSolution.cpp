#include "acequia_manager.h"
#include <algorithm>
#include <vector>

/*Instructions for this problem:

	The intend of this project is to simulate water management conservation principles in the context of New Mexico

	In this simulation, there exists several Regions (North, South, etc.). Each region class includes both:
	- a given water level
	- a given water need
	- a indicator boolean for if the region is flooded
	- an indicator boolean for if the region is in drought

	With each region, there are given waterSources provided (perhaps a .dat file list each waterSource to  a region) 
	and certain waterSources have canals connected to them to deliver water across regions.

	Given the current state of the region, students wil be asked to utlize the canals that connect regions to
	develop the logic and algorithm for finding a solution. The simulation has a fixed time



	The student solution will be evaluated on the criteria that each region meets the following:
	- a given region is NOT flooded
	- a given region is NOT in drought
	- the region waterNeed does not exceed the region waterLevel 
*/

/*This will be how the solveProblems function is set up. The student may enter their on  */
/*
void solveProblems(AcequiaManager& manager)
{
	//the student can call the members of the canals object such as name of canal. sourceRegion, and destinationRegion
	//This could be helpful in informing the students strategy to solve the problem
	auto canals = manager.getCanals();
	//students may call to get Region and WaterSource objects to inform decisions 


	while(!manager.isSolved && manager.hour!=manager.SimulationMax)
	{	
		//enter student code here


		manager.nexthour();
	}
}
*/


// Helper function to find canals by source and destination regions: George
std::vector<Canal*> findCanalsByRoute(const std::vector<Canal*>& canals, 
                                       const std::string& sourceName, 
                                       const std::string& destName) {
    std::vector<Canal*> matchingCanals;
    for(auto canal : canals) {
        if(canal->sourceRegion->name == sourceName && 
           canal->destinationRegion->name == destName) {
            matchingCanals.push_back(canal);
        }
    }
    return matchingCanals;
}

// Helper function to find all canals from a source region: George
std::vector<Canal*> findCanalsFromSource(const std::vector<Canal*>& canals, 
                                          const std::string& sourceName) {
    std::vector<Canal*> matchingCanals;
    for(auto canal : canals) {
        if(canal->sourceRegion->name == sourceName) {
            matchingCanals.push_back(canal);
        }
    }
    return matchingCanals;
}

// Helper function to find all canals to a destination region: George
std::vector<Canal*> findCanalsToDestination(const std::vector<Canal*>& canals, 
                                              const std::string& destName) {
    std::vector<Canal*> matchingCanals;
    for(auto canal : canals) {
        if(canal->destinationRegion->name == destName) {
            matchingCanals.push_back(canal);
        }
    }
    return matchingCanals;
}

// Helper function to calculate optimal flow rate based on deficit/surplus: George
double calculateFlowRate(double deficit, double surplus, double maxFlow = 1.0) {
    // Use the smaller of deficit or surplus to avoid over-transfer
    double targetAmount = (deficit < surplus) ? deficit : surplus;
    
    // Convert to flow rate (approximate: 1.0 flow rate = 3.6 units/hour)
    // Flow rate is in gal/sec, and updateWater uses 3600 seconds
    // So flowRate * 3600 / 1000 = flowRate * 3.6 units per hour
    double flowRate = targetAmount / 3.6;
    
    // For urgent cases (large deficits), use maximum flow rate
    if(deficit > 15.0) {
        flowRate = maxFlow; // Use maximum flow for urgent cases
    }
    // For moderate deficits, use high flow rate
    else if(deficit > 5.0) {
        flowRate = std::min(flowRate * 1.2, maxFlow * 0.8);
    }
    
    // Clamp between 0.3 and maxFlow (minimum 0.3 to ensure meaningful transfer)
    if(flowRate < 0.3) flowRate = 0.3;
    if(flowRate > maxFlow) flowRate = maxFlow;
    
    return flowRate;
}

// Main adaptive water management solution - balanced and deterministic: George
void solveProblems(AcequiaManager& manager)
{
    auto canals = manager.getCanals();
    auto regions = manager.getRegions();

    // Helper to find a region by name
    auto findRegion = [&](const std::string& name) -> Region* {
        for(auto r : regions){
            if(r->name == name) return r;
        }
        return nullptr;
    };

    // Cache common region pointers (may be nullptr if names differ)
    Region* north = findRegion("North");
    Region* south = findRegion("South");
    Region* east  = findRegion("East");

    auto clampFlow = [](double amt){
        double f = amt / 3.6;           // convert hourly amount to flowRate
        if(f < 0.15) f = 0.15;          // minimum meaningful flow
        if(f > 1.0)  f = 1.0;
        return f;
    };

    // Route from source to destination if helpful
    auto route = [&](Region* src, Region* dst){
        if(!src || !dst) return;
        if(src == dst) return;

        // Skip if destination already satisfied
        if(dst->waterLevel >= dst->waterNeed && !dst->isInDrought) return;

        // Compute surplus at source (above need + small buffer)
        double surplus = src->waterLevel - (src->waterNeed + 1.0);
        if(src->isFlooded){
            surplus = std::max(surplus, src->waterLevel - src->waterCapacity + 0.5);
        }
        if(src->isInDrought || surplus <= 0.0) return;

        // Compute deficit at destination (need + small buffer)
        double deficit = (dst->waterNeed + 1.0) - dst->waterLevel;
        if(dst->isInDrought) deficit += 2.0; // push harder for drought
        if(deficit <= 0.0) return;

        // Respect destination capacity
        double capMargin = dst->waterCapacity - dst->waterLevel - 0.5;
        if(capMargin <= 0.0) return;

        double amount = std::min({surplus, deficit, capMargin});
        if(amount <= 0.2) return; // too small to matter

        double flowRate = clampFlow(amount);

        auto connecting = findCanalsByRoute(canals, src->name, dst->name);
        for(auto canal : connecting){
            canal->setFlowRate(flowRate);
            canal->toggleOpen(true);
        }
    };

    while(!manager.isSolved && manager.hour != manager.SimulationMax)
    {
        // Reset all canals each hour to avoid stale settings: Jacob
        for(auto canal : canals){
            canal->toggleOpen(false);
            canal->setFlowRate(0.0);
        }

        // 1) Relieve any flooded regions first
        for(auto src : regions){
            if(src->isFlooded || src->waterLevel >= src->waterCapacity){
                for(auto dst : regions){
                    if(dst == src) continue;
                    if(dst->isFlooded) continue;
                    route(src, dst);
                }
            }
        }

        // 2) Use known network paths (A: N->S, B: S->E, C: N->E, D: E->N)
        route(north, south);
        route(south, east);
        route(north, east);
        route(east, north);

        // 3) General balancing across all pairs (best-effort)
        for(auto src : regions){
            for(auto dst : regions){
                if(src == dst) continue;
                route(src, dst);
            }
        }

        // Advance simulation time
        manager.nexthour();
    }
}


/*example 2*/
/*
void solveProblems(AcequiaManager& manager)
{
	auto canals = manager.getCanals();
	while(!manager.isSolved && manager.hour!=manager.SimulationMax)
	{
	//Students will implement this function to solve the probelms
	//Example: Adjust canal flow rates and directions
		if(manager.hour==1)
		{
			canals[0]->setFlowRate(0.1);
			canals[0]->toggleOpen(true);
			canals[1]->setFlowRate(0.2);
			canals[1]->toggleOpen(true);
		}
		else if(manager.hour==3)
		{
			canals[0]->toggleOpen(false);
			canals[1]->toggleOpen(false);
		}
	//student may add any necessary functions or check on the progress of each region as the simulation moves forward. 
	//The manager takes care of updating the waterLevels of each region and waterSource while the student is just expected
	//to solve how to address the state of each region

		
		manager.nexthour();
	}
}
*/


//In this solution, students can make functions that aid in identifying the best course of action for moving
//water resources. They can set conditions that then work on the canal vectors based on the information reported

//This can help in optimizing solutions for dynamic constraints such as weather (rain, or dried up waterSources) and
//make the solution to the problem more abstract, allowing the logic and algorithm to be the sole priority of the student
//while the computation is left for the Acequia Manager

//This would be the perfect opportunity to identify the tools learned from ECE 231L such as:
//data structures (stacks, queues, trees(?)), templates, vector class functions, etc... to aid in the algorithm solution

/*
int findCanal(std::vector<Canal *> canals, std::string region)
{
	int match;
	for(int i = 0; i< canals.size();i++)
	{
		if(canals[i]->sourceRegion->name == region)
		{
			match = i;
		}
	}
	return match;
}

void release(std::vector<Canal *> canals, std::string region)
{
	int match = findCanal(canals, region);
	canals[match]->toggleOpen(true);
	canals[match]->setFlowRate(1);
	return;
}

void close(std::vector<Canal *> canals, std::string region)
{
	int match = findCanal(canals, region);
	canals[match]->toggleOpen(false);
}


void solveProblems(AcequiaManager& manager)
{

	bool check = false;
	auto canals = manager.getCanals();
	auto regions = manager.getRegions();
	while(!manager.isSolved && manager.hour!=manager.SimulationMax)
	{
		
		if(manager.hour == 0)
		{
			for(int i = 0; i<canals.size(); i++)
			{
				canals[i]->toggleOpen(true);
				canals[i]->setFlowRate(1);
			}
		}

		for(int i =0 ; i<regions.size(); i++)
		{
			if(regions[i]->isFlooded == true)
			{
				//release water from that region by a canal
				release(canals, regions[i]->name);
			}
			else if(regions[i]->isInDrought = true)
			{
				//find canal to move water
				close
			}

			else if(regions[i]->isFlooded == true && regions[i]->isInDrought == true)
			{
				close(canals, regions[i]->name);
			}
		}
		
		manager.nexthour();
	}
}
*/
