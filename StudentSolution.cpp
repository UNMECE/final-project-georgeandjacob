#include "acequia_manager.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>

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


// Helper function to find canals by source and destination regions
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

// Helper function to find all canals from a source region
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

// Helper function to find all canals to a destination region
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

// Helper function to calculate optimal flow rate based on deficit/surplus
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

// Main adaptive water management solution - Maximum efficiency version
void solveProblems(AcequiaManager& manager)
{
    auto canals = manager.getCanals();
    auto regions = manager.getRegions();
    
    while(!manager.isSolved && manager.hour != manager.SimulationMax)
    {
        // Strategy: Maximize water flow to needy regions while preventing overflow/drought
        // Use all available canals aggressively but intelligently
        
        // Close problematic canals only
        for(auto canal : canals) {
            if(canal->isOpen) {
                // Close if destination is flooded and source is not
                if(canal->destinationRegion->isFlooded && !canal->sourceRegion->isFlooded) {
                    canal->toggleOpen(false);
                    continue;
                }
                
                // Close if source is in drought
                if(canal->sourceRegion->isInDrought) {
                    canal->toggleOpen(false);
                    continue;
                }
                
                // Close if source needs water much more than destination
                double sourceDeficit = canal->sourceRegion->waterNeed - canal->sourceRegion->waterLevel;
                double destDeficit = canal->destinationRegion->waterNeed - canal->destinationRegion->waterLevel;
                if(sourceDeficit > destDeficit + 15.0 && 
                   canal->sourceRegion->waterLevel < canal->sourceRegion->waterNeed - 5.0) {
                    canal->toggleOpen(false);
                    continue;
                }
            }
        }
        
        // Open canals for flooded regions first (highest priority)
        for(auto region : regions) {
            if(region->isFlooded || region->waterLevel >= region->waterCapacity) {
                // Find all needy destinations and route to them
                for(auto dest : regions) {
                    if(dest == region || dest->isFlooded) continue;
                    if(dest->waterLevel >= dest->waterNeed + 3.0) continue;
                    
                    auto connectingCanals = findCanalsByRoute(canals, region->name, dest->name);
                    for(auto canal : connectingCanals) {
                        if(!canal->isOpen) {
                            canal->setFlowRate(1.0);
                            canal->toggleOpen(true);
                        }
                    }
                }
            }
        }
        
        // Now route to all regions that need water
        for(auto needyRegion : regions) {
            if(needyRegion->isFlooded) continue;
            if(needyRegion->waterLevel >= needyRegion->waterNeed) continue;
            
            double needyDeficit = needyRegion->waterNeed - needyRegion->waterLevel;
            bool isCritical = needyRegion->isInDrought || 
                             needyRegion->waterLevel < 0.2 * needyRegion->waterCapacity;
            
            // Find all possible sources and rank them
            std::vector<std::pair<Region*, double>> sources;
            
            for(auto candidate : regions) {
                if(candidate == needyRegion) continue;
                
                double score = 0.0;
                bool canSupply = false;
                
                // Flooded regions - highest priority
                if(candidate->isFlooded || candidate->waterLevel >= candidate->waterCapacity) {
                    canSupply = true;
                    score = 100000.0;
                }
                // Well above need
                else if(candidate->waterLevel > candidate->waterNeed + 8.0) {
                    canSupply = true;
                    score = 10000.0 + (candidate->waterLevel - candidate->waterNeed);
                }
                // Above need
                else if(candidate->waterLevel > candidate->waterNeed + 3.0) {
                    canSupply = true;
                    score = 1000.0 + (candidate->waterLevel - candidate->waterNeed);
                }
                // At or just above need
                else if(candidate->waterLevel >= candidate->waterNeed) {
                    if(isCritical || needyDeficit > 20.0) {
                        canSupply = true;
                        score = 100.0 + (candidate->waterLevel - candidate->waterNeed);
                    }
                }
                // For critical needs, consider regions close to need
                else if(isCritical && needyDeficit > 25.0 &&
                        candidate->waterLevel > candidate->waterNeed - 8.0 &&
                        candidate->waterLevel > 0.35 * candidate->waterCapacity) {
                    canSupply = true;
                    score = 10.0;
                }
                
                // Never use sources in drought or very low
                if(candidate->isInDrought || candidate->waterLevel < 0.25 * candidate->waterCapacity) {
                    canSupply = false;
                }
                
                if(canSupply) {
                    sources.push_back({candidate, score});
                }
            }
            
            // Sort sources by score (best first)
            std::sort(sources.begin(), sources.end(),
                      [](const std::pair<Region*, double>& a, const std::pair<Region*, double>& b) {
                          return a.second > b.second;
                      });
            
            // Open canals from best sources
            for(auto& sourcePair : sources) {
                auto connectingCanals = findCanalsByRoute(canals, sourcePair.first->name, needyRegion->name);
                for(auto canal : connectingCanals) {
                    if(!canal->isOpen) {
                        canal->setFlowRate(1.0);
                        canal->toggleOpen(true);
                    }
                }
                // Use multiple sources if deficit is very large
                if(needyDeficit < 20.0) break; // Only use one source for smaller deficits
            }
        }
        
        // Advance to next hour
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