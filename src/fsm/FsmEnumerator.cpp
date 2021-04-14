#include "FsmEnumerator.h"

#include <iostream>

#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/FsmLabel.h"

FsmEnumerator::FsmEnumerator(int maxInput, int maxOutput, int maxState, const std::shared_ptr<FsmPresentationLayer>& presentationLayer, bool generateSmallerFsms) 
    : maxInput(maxInput), maxOutput(maxOutput), maxState(maxState), generateSmallerFsms(generateSmallerFsms), fsm("E",maxInput,maxOutput,std::vector<std::shared_ptr<FsmNode>>(),presentationLayer) 
{
    // initialize table
    for (int i = 0; i <= maxState; ++i) {
        std::vector<std::vector<int>> vi;
        for (int x = 0; x <= maxInput; ++x) {
            std::vector<int> vx;
            for (int y = 0; y <= maxOutput; ++y) {
                vx.push_back(-1);
            }
            vi.push_back(vx);
        }
        currentTable.push_back(vi);
    }
    currentTable[maxState][maxInput][maxOutput] = -2; // create a table "before" the actual first table
}


bool FsmEnumerator::updateTable() {

    std::cout << "update table:" << std::endl;    
    for (int i = 0; i <= maxState; ++i) {
        std::cout << i << ": [ ";
        for (int x = 0; x <= maxInput; ++x) {
            for (int y = 0; y <= maxOutput; ++y) {
                std::cout << " " << x << "/" << y << "->" << currentTable[i][x][y];
            }
        }
        std::cout << " ]" << std::endl;
    }


    int i;
    int x;
    int y;
    bool updatePossible = false;

    for (i = maxState; i >= 0; --i) {
        for (x = maxInput; x >= 0; --x) {
            for (y = maxOutput; y >= 0; --y) {
                if (currentTable[i][x][y] < maxState) {
                    currentTable[i][x][y] += 1; // just re-target the transition to the next state
                    updatePossible = true;
                    break;
                }
            }
            if (updatePossible) break;
        }
        if (updatePossible) break;
    }

    // return false of no transition can be re-targeted
    if (!updatePossible) return false;

    // otherwise re-target all entries following (i,x,y) to -1
    for (int di = i; di <= maxState; ++di) {
        for (int dx = (di==i) ? x : 0; dx <= maxInput; ++dx) {
            for (int dy = (di==i && dx==x) ? y+1 : 0; dy <= maxOutput; ++dy) {
                currentTable[di][dx][dy] = -1;
            }
        }
    }
    return true;
}

Fsm FsmEnumerator::generateFsmFromTable() {
    
    std::vector<std::shared_ptr<FsmNode>> nodes;
    for (int i = 0; i <= maxState; ++i) {
        nodes.push_back(std::make_shared<FsmNode>(i,fsm.getPresentationLayer()));
    }

    for (int i = 0; i <= maxState; ++i) {
        for (int x = 0; x <= maxInput; ++x) {
            for (int y = 0; y <= maxOutput; ++y) {
                int tgt = currentTable[i][x][y];
                if (tgt >= 0) {
                    auto label = std::make_shared<FsmLabel>(x,y,fsm.getPresentationLayer());
                    auto transition = std::make_shared<FsmTransition>(nodes[i], nodes[tgt], label);
                    nodes[i]->addTransition(transition);
                }
            }
        }
    }
    
    return Fsm(std::to_string(candidateNum),maxInput,maxOutput,nodes,fsm.getPresentationLayer());
}

bool FsmEnumerator::hasNext() {
    
    while (updateTable()) {
        fsm = generateFsmFromTable();

        // early exit if reachability is not to be checked
        if (generateSmallerFsms) {
            ++candidateNum;
            return true;
        }

        std::vector<std::shared_ptr<FsmNode>> unreachableNodes;
        fsm.removeUnreachableNodes(unreachableNodes);
        if (unreachableNodes.empty()) {
            ++candidateNum;
            return true;
        }
    }

    // no new Fsm could be generated
    return false;
}


Fsm FsmEnumerator::getNext() {
    return fsm;
}