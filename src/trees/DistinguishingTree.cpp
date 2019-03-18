/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#include <queue>
#include <unordered_map>
#include <functional>
#include "trees/DistinguishingTree.h"
#include "utils/Logger.hpp"

using namespace std;

DistinguishingTree::DistinguishingTree(const std::shared_ptr <Dfsm> &dfsm)
        : distinguishingSequence(vector<int>())
{
    //Logging
    LogCoordinator& logger = LogCoordinator::getStandardLogger();
    logger.setDefaultStream(cout);

    root = make_shared<DistinguishingTreeNode>();

    //A distinguishing tree can only be derived for truly completely specified, deterministic FSM
    if(!dfsm->isCompletelyDefined() || !dfsm->isDeterministic()) {
        return;
    }
    //Mapping from id to FsmNode to speed up the access
    idToFsmNode = vector<shared_ptr<FsmNode>>(dfsm->getMaxNodes(),nullptr);

    //Initialize the current uncertainty of the root node
    auto &initialCurrentUncertainty = root->getCurrentUncertainty();
    auto initialBlock = set<int>();

    for(auto node:dfsm->getNodes()) {
        initialBlock.insert(node->getId());
        idToFsmNode[node->getId()] = node;
    }

    initialCurrentUncertainty.insert(initialBlock);

    //Working list
    queue<shared_ptr<DistinguishingTreeNode>> workingList;
    workingList.push(root);

    // The hash table that is used to check if a current uncertainty already exists
    unordered_map< int, vector <multiset<set<int>>*> > uncTable;
    uncTable.insert({createHashForCurrentUncertainty(initialCurrentUncertainty),vector<multiset<set<int>>*> {&initialCurrentUncertainty}});

    //Algorithm to create the distinguishing tree and derive a distinguishing sequence for the dfsm
    while(!workingList.empty()) {
        LOG("VERBOSE_2") << " Next One please! " << std::endl;
        auto currentNode = workingList.front();
        workingList.pop();
        LOG("VERBOSE_2") << uncToString(currentNode->getCurrentUncertainty()) << std::endl;
        LOG("VERBOSE_2") << traceToString(currentNode->getInputTrace()) << std::endl;


        for(int x=0; x<=dfsm->getMaxInput(); ++x) {
            auto newNode = make_shared<DistinguishingTreeNode>();
            auto &nextCurrentUncertainty = newNode->getCurrentUncertainty();

            //create the next current Uncertainty from the current node with regard to x
            bool isTerminal = false;
            LOG("VERBOSE_2") << "compute next unc" << std::endl;
            computeNextCurrentUncertainty(currentNode->getCurrentUncertainty(),x,nextCurrentUncertainty,isTerminal);

            //the input 'x' wasnt valid for the current node, try the next one
            if(nextCurrentUncertainty.empty()) {
                continue;
            }

            //check whether the current uncertainty already exists earlier in the tree
            bool alreadyExists = false;
            int uncHash = createHashForCurrentUncertainty(nextCurrentUncertainty);
            auto uncIt = uncTable.find(uncHash);
            if(uncIt == uncTable.end()) {
                uncTable.insert({uncHash,vector<multiset<set<int>>*> {&nextCurrentUncertainty}});
            } else {
                for(auto& uncertainty:uncIt->second) {
                    if(nextCurrentUncertainty == *uncertainty) {
                        alreadyExists = true;
                        break;
                    }
                }
                if(alreadyExists) {
                    continue;
                } else {
                    uncIt->second.push_back(&nextCurrentUncertainty);
                }
            }

            auto &nextInputTrace = newNode->getInputTrace();
            nextInputTrace = currentNode->getInputTrace();
            nextInputTrace.push_back(x);

            /**
             * Check whether the current uncertainty consists of |number of dfsm states| singletons only.
             * In that case the distinguishing sequence has been found.
             */
            if(isTerminal) {
                distinguishingSequence = nextInputTrace;
                return;
            }

            auto newEdge = make_shared<DistinguishingTreeEdge>(x,newNode);
            currentNode->add(newEdge);

            workingList.push(newNode);

        }
    }

}

void DistinguishingTree::computeNextCurrentUncertainty(const multiset<set<int>> &currentUncertainty, int x,
                                                       multiset<set<int>> &nextCurrentUncertainty, bool& isTerminal) {
    //clear the 'nextCurrentUncertainty' for good measure. Might be unnecessary since it is always empty in the context the function is used.
    nextCurrentUncertainty.clear();

    isTerminal = true;

    for(auto& block:currentUncertainty) {
        //Partition the block of target states after reading x respecting the produced output
        unordered_map< int, shared_ptr<set<int>> > blockPartition;
        LOG("VERBOSE_2") << "starting new block" << std::endl;
        for(int id:block) {
            auto producedOutputs = vector<int>();
            /**
             * Apply 'x' to the node with id 'id'. note that there must exactly one output and one
             * target state produced, since 'dfsm' is completely specified and deterministic
             **/
            auto targetStates = idToFsmNode[id]->after(x,producedOutputs);
            int output = producedOutputs.front();
            LOG("VERBOSE_2") << id << " going to " << targetStates.front()->getId() << " with output " << output  << std::endl;
            auto bit = blockPartition.find(output);
            if(bit == blockPartition.end()) {
                auto partition = make_shared<set<int>>();
                blockPartition.insert({output,partition});
                partition->insert(targetStates.front()->getId());
            } else {
                auto it = bit->second->insert(targetStates.front()->getId());
                //the targetstate already exists inside that partition, which means that 'x' is not valid for 'currentUncertainty'
                if(!it.second) {
                    isTerminal = false;
                    return;
                }
            }
        }

        for(auto& p:blockPartition) {
            isTerminal &= (p.second->size() == 1);
            nextCurrentUncertainty.insert(*p.second);
        }

    }
}

int DistinguishingTree::createHashForCurrentUncertainty(const multiset<set<int>>& currentUncertainty) {
    string stringHash = "";
    hash<string> hash_fn;
    for(auto& block:currentUncertainty) {
        for(int id:block) {
            stringHash += to_string(id) + ",";
        }
        stringHash += "|";
    }
    return hash_fn(stringHash);
}

string DistinguishingTree::uncToString(const multiset<set<int>>& currentUncertainty) {
    string serialized = "{";
    for(auto& block:currentUncertainty) {
        serialized += " {";
        for(int id:block) {
            serialized += to_string(id) + ",";
        }
        serialized += "} ";
    }
    serialized += "}";

    return serialized;
}

string DistinguishingTree::traceToString(const vector<int> &trace) {
    string serialized = "";
    for(auto io:trace) {
        serialized += to_string(io) + ".";

    }
    return serialized;
}

const std::vector<int>& DistinguishingTree::getDistinguishingSequence() {
    return distinguishingSequence;
}
