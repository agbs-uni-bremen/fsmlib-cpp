/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

// TODO: remove superfluous includes
#include <chrono>
#include <algorithm>
#include <numeric>
#include <regex>
#include <math.h>
#include <cassert>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <functional>
#include <queue>

#include "fsm/StrongReductionTestSuiteGenerator.h"
#include "fsm/Fsm.h"
#include "fsm/FsmTransition.h"
#include "fsm/FsmLabel.h"
#include "fsm/FsmNode.h"
#include "fsm/OFSMTable.h"
#include "fsm/FsmVisitor.h"
#include "fsm/RDistinguishability.h"
#include "fsm/VPrimeLazy.h"
#include "fsm/IOTrace.h"
#include "sets/HittingSet.h"
#include "trees/AdaptiveTreeNode.h"
#include "trees/TreeNode.h"
#include "trees/TreeEdge.h"
#include "trees/IOListContainer.h"
#include "utils/Logger.hpp"
#include "utils/generic-equivalence-class-calculation.hpp"
#include "trees/TestSuite.h"
#include "trees/InputOutputTree.h"
#include "trees/InputTree.h"
#include "trees/OutputTree.h"
#include "trees/IOTreeContainer.h"
#include "fsm/IOTraceContainer.h"
#include "interface/FsmPresentationLayer.h"





void StrongReductionTestSuiteGenerator::calcDeterministicallyReachingSequences() {
    dReachingSequences = std::unordered_map<std::shared_ptr<FsmNode>, std::vector<int>>();

    // copy the fsm to perform modifications
    Fsm m(*fsm);
    vector<shared_ptr<FsmNode>> nodes = fsm->getNodes();
    int maxState = fsm->getMaxState();

    // add sink state 
    std::shared_ptr<FsmNode> sink = m.addNode("bot");

    // add transitions to the sink state for undefined inputs
    for ( auto node : m.getNodes()) {
        node->setColor(node->white); // set color to white for later search
        for (int x = 0; x <= m.getMaxInput(); ++x) {
            if (!node->hasTransition(x)) {
                shared_ptr<FsmLabel> label = make_shared<FsmLabel>(x, 0, fsm->getPresentationLayer());
                node->addTransition(make_shared<FsmTransition>(node, sink, label));
            }
        }

        // set all transition-outputs to 0
        for (auto transition : node->getTransitions()) {
            shared_ptr<FsmLabel> label = make_shared<FsmLabel>(transition->getLabel()->getInput(), 0, fsm->getPresentationLayer());
            transition->setLabel(label);
        }        
    }

    // make m observable, which effectively determinises the automaton obtained by dropping the outputs, as all outputs are 0
    Fsm mObs = m.transformToObservableFSM();

    // labels of states of mObs corresponding to singletons of states in fsm
    unordered_map<int, string> node2AutomatonLabel;
    for ( int n = 0; n <= maxState; ++n) {

        auto node = nodes[n];
        unordered_set<shared_ptr<FsmNode>> expectedLabel;
        expectedLabel.insert(node);
        string nodeName = Fsm::labelString(expectedLabel);
        node2AutomatonLabel[n] = nodeName;
    }
        
    // get deterministically reaching sequences via breadth first search,
    // where a sequence reaching {q} in mObs d-reaches q in fsm
    std::queue<std::shared_ptr<FsmNode>> todo;
    std::queue<std::vector<int>> prevSequences;
    todo.push(mObs.getInitialState());
    prevSequences.push(std::vector<int> ());
    while (!todo.empty()) {
        std::shared_ptr<FsmNode> curNode = todo.front();
        todo.pop();
        std::vector<int> dReachingSequence = prevSequences.front();
        prevSequences.pop();

        if (curNode->getColor() == curNode->black) {
            continue;
        }
        curNode->setColor(curNode->black);

        // check if current node is a singleton
        for ( int n = 0; n <= fsm->getMaxNodes(); ++n) {
            if (curNode->getName().compare(node2AutomatonLabel[n]) == 0) {
                dReachingSequences[nodes[n]] = dReachingSequence;
                break;
            } 
        }

        for (auto transition : curNode->getTransitions()) {
            todo.push(transition->getTarget());
            std::vector<int> nextSeq(dReachingSequence);
            nextSeq.push_back(transition->getLabel()->getInput());
            prevSequences.push(nextSeq);
        }
        
    }
}




void StrongReductionTestSuiteGenerator::calcRDistinguishingGraph() {

    rDistGraph = std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>>();

    vector<shared_ptr<FsmNode>> nodes = fsm->getNodes();
    int maxState = fsm->getMaxState();

    std::vector<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>> todo;
    std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>> emptyTargets = make_shared<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>();
    for (int i = 0; i < maxState; ++i) {
        auto ni = nodes[i];
        auto definedInputsI = ni->getDefinedInputs();

        for (int j = i+1; j <= maxState; ++j) {
            auto nj = nodes[j];

            auto pair = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(ni,nj);

            if (definedInputsI != nj->getDefinedInputs()) {
                // nodes are r(0)-distinguishable if their sets of defined inputs differ
                // thus no input needs to be applied (marked as -1)
                auto value = std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>(-1,emptyTargets);
                rDistGraph[pair] = value;

                // also add reversed pair
                auto pairR = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(nj,ni);
                rDistGraph[pairR] = value;
            } else {
                todo.push_back(pair);
            }
        }
    }

    bool reachedFixpoint = false;

    while (!reachedFixpoint) {
        // to be set to false again if result is updated in this iteration
        reachedFixpoint = true;

        std::vector<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>> todoNext;

        for (auto nodePair : todo) {            

            if (rDistGraph.count(nodePair) == 1) {
                // skip pairs that have already been added
                continue;
            }

            auto n1 = nodePair.first;
            auto n2 = nodePair.second;
            
            for (int x : n1->getDefinedInputs()) {
                std::unordered_set<int> outputs1;
                for (auto output : n1->getPossibleOutputs(x)) {
                    outputs1.insert(output->get().front());
                }
                
                std::unordered_set<int> sharedOutputs;
                for (auto output : n2->getPossibleOutputs(x)) {
                    int y = output->get().front();
                    if (outputs1.find(y) != outputs1.end()) {
                        sharedOutputs.insert(y);
                    }
                }

                std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>> targets = make_shared<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>();

                bool distinguishesAllSharedOutputs = true;
                for (int y : sharedOutputs) {

                    // as x is defined in both nodes (as the are in todo and hence not r(0)-distinguishable)
                    // and y is a shared output, both nodes must exhibit a transition for x/y
                    std::shared_ptr<FsmNode> t1;
                    for (shared_ptr<FsmTransition> trans : n1->getTransitions()) {
                        if (trans->getLabel()->getInput() == x && trans->getLabel()->getOutput() == y) {
                            t1 = trans->getTarget();
                        }
                    }
                    std::shared_ptr<FsmNode> t2;
                    for (shared_ptr<FsmTransition> trans : n2->getTransitions()) {
                        if (trans->getLabel()->getInput() == x && trans->getLabel()->getOutput() == y) {
                            t2 = trans->getTarget();
                        }
                    }

                    auto targetPair = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(t1,t2); 
                    targets->insert(targetPair);

                    if (rDistGraph.count(targetPair) == 0) {
                        distinguishesAllSharedOutputs = false;
                        break;
                    }
                }

                if (distinguishesAllSharedOutputs) {
                    // the pair can be distinguished by applying x as first input
                    rDistGraph[nodePair] = std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>(x, targets);
                    // also add reversed pair
                    auto pairR = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(n2,n1);
                    rDistGraph[pairR] = std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>(x, targets);

                    // result has been updated
                    reachedFixpoint = false;
                    break;
                } else {
                    todoNext.push_back(nodePair);
                }
            }
        }
        todo = todoNext;
    }
}



void StrongReductionTestSuiteGenerator::calcRDistinguishingTrees() {
    calcRDistinguishingGraph();
    rDistTrees = std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::shared_ptr<InputTree>>();
    
    std::function<std::shared_ptr<InputTree>(std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>)> helper = [this,&helper](std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>> pair)->std::shared_ptr<InputTree> {
        auto res = make_shared<InputTree>(fsm->getPresentationLayer());
        auto entry = rDistGraph.find(pair);
        int x = entry->second.first;
        auto targets = entry->second.second;

        if (x == -1) {
            return res;
        }

        if (targets->empty()) {
            res->addToRoot(std::vector<int> {x});
        }

        for (auto targetPair : *targets) {
            auto targetTree = helper(targetPair);
            for (auto trace : targetTree->getInputTraces()) {
                trace.prepend(std::vector<int> {x});   
                res->addToRoot(trace.get());
            }
        }

        return res;
    };

    for (auto entry : rDistGraph) {
        rDistTrees[entry.first] = helper(entry.first);
    }
}

void StrongReductionTestSuiteGenerator::calcMaximalRDistinguishableSets() {
    maximalRDistSets = std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>>();

    for (auto node : fsm->getNodes()) {
        std::unordered_set<std::shared_ptr<FsmNode>> set;
        set.insert(node);

        bool fixpointReached = false;
        while (!fixpointReached) {
            fixpointReached = true;
            for (auto otherNode : fsm->getNodes()) {
                if (set.count(otherNode) != 0) continue;

                bool isRdFromAllCurrentNodes = true;
                for (auto nodeToBeRD : set) {
                    auto pair = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(nodeToBeRD,otherNode);
                    if (rDistTrees.count(pair) == 0) {
                        isRdFromAllCurrentNodes = false;
                        break;
                    }
                }
                if (isRdFromAllCurrentNodes) {
                    set.insert(otherNode);
                    fixpointReached = false;
                }            
            }
        }

        maximalRDistSets.insert(set);
    }
}

void StrongReductionTestSuiteGenerator::calcAllMaximalRDistinguishableSets() {
    maximalRDistSets = std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>>();

    vector<shared_ptr<FsmNode>> nodes = fsm->getNodes();

    std::function<std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>>(unsigned int)> helper = [this,nodes,&helper](unsigned int n)->std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> {
        if (n == 0) {
            std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> result;
            std::unordered_set<std::shared_ptr<FsmNode>> emptySet;
            std::unordered_set<std::shared_ptr<FsmNode>> singleton;
            singleton.insert(nodes[n]);
            result.insert(singleton);
            result.insert(emptySet);
            return result;
        }

        auto node = nodes[n];
        auto sets = helper(n-1);

        auto result = std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>>(sets);
        for (auto set : sets) {
            bool isRdFromAllCurrentNodes = true;
            for (auto otherNode : set) {
                auto pair = std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>(node,otherNode);
                if (rDistTrees.count(pair) == 0) {
                    isRdFromAllCurrentNodes = false;
                    break;
                }
                           
            }
            if (isRdFromAllCurrentNodes) {
                auto setWithNode = std::unordered_set<std::shared_ptr<FsmNode>>(set);
                setWithNode.insert(node);
                result.insert(setWithNode);
            } 
        }
        return result;
    };

    auto candidates = helper(nodes.size()-1);

    // TODO: inefficient
    for (auto candidate : candidates) {
        bool isMaximal = true;
        for (auto otherCandidate : candidates) {
            if (candidate == otherCandidate) continue;
            
            bool isProperSubset = true;
            for (auto elem : candidate) {
                if (otherCandidate.count(elem) == 0) {
                    isProperSubset = false;
                    break;
                }
            }
            if (isProperSubset) {
                isMaximal = false;
                break;
            }
        }
        if (isMaximal) {
            maximalRDistSets.insert(candidate);
        }
    }
}


StrongReductionTestSuiteGenerator::StrongReductionTestSuiteGenerator(const std::shared_ptr<Fsm> fsm, bool calculateAllMaximalRDistinguishableSets) : fsm(fsm)
{
    calcDeterministicallyReachingSequences();
    calcRDistinguishingTrees();
    if (calculateAllMaximalRDistinguishableSets) {
        calcAllMaximalRDistinguishableSets();
    } else {
        calcMaximalRDistinguishableSets();
    }
}


std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::shared_ptr<InputTree>> StrongReductionTestSuiteGenerator::getRDistinguishingTrees() const {
    return rDistTrees;
}

std::unordered_map<std::shared_ptr<FsmNode>, std::vector<int>> StrongReductionTestSuiteGenerator::getDeterministicallyReachingSequences() const {
    return dReachingSequences;
}

std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>> StrongReductionTestSuiteGenerator::getRDistGraph() const {
    return rDistGraph;
}

std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> StrongReductionTestSuiteGenerator::getMaximalRDistinguishableSets() {
    return maximalRDistSets;
}

std::vector<std::pair<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>,int>> StrongReductionTestSuiteGenerator::getTerminationTuples(int m) {
    std::vector<std::pair<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>,int>> result;

    for (auto rdSet : maximalRDistSets) {
        int dr = 0;
        for (auto node : rdSet) {
            if (dReachingSequences.count(node) != 0) {
                ++ dr;
            }
        }
        result.push_back(std::make_pair(make_shared<std::unordered_set<std::shared_ptr<FsmNode>>>(rdSet),m-dr+1));
    }

    return result;
}


std::vector<std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>> StrongReductionTestSuiteGenerator::calcTraversalSet(std::shared_ptr<FsmNode> node, int m) {
    std::vector<std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>> result;

    std::vector<std::pair<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>,int>> terminationTuples = getTerminationTuples(m);
    std::vector<int> initialMissingVisits;
    for (unsigned int i = 0; i < terminationTuples.size(); ++i) {
        initialMissingVisits.push_back(terminationTuples[i].second);
    }

    std::function<std::vector<std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>>(std::shared_ptr<FsmNode>, std::vector<int>)> helper = [this,terminationTuples,&helper](std::shared_ptr<FsmNode> node, std::vector<int> missingVisits)->std::vector<std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>> {

        std::vector<std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>> result;
        std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>> terminatingSets;
        for (unsigned int i = 0; i < missingVisits.size(); ++i) {
            if (missingVisits[i] <= 0) {
                terminatingSets.push_back(terminationTuples[i].first);
            }
        }

        // terminate if all required visits have been performed for some set
        if (!terminatingSets.empty()) {
            IOTrace emptyTrace = IOTrace(InputTrace(fsm->getPresentationLayer()), OutputTrace(fsm->getPresentationLayer()), node);
            std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>> singleEntry = std::make_pair(emptyTrace,terminatingSets);
            result.push_back(singleEntry);
            return result;
        }

        for (auto transition : node->getTransitions()) {            
            // update visits based on current target
            std::vector<int> nextMissingVisits;
            for (unsigned int i = 0; i < missingVisits.size(); ++i) {
                if (terminationTuples[i].first->count(transition->getTarget()) == 0) {
                    nextMissingVisits.push_back(missingVisits[i]);
                } else {
                    nextMissingVisits.push_back(missingVisits[i]-1);
                }
            }

            auto targetResult = helper(transition->getTarget(), nextMissingVisits);
            for (auto entry : targetResult) {
                // prepend results for the transition target by the transitions IO
                IOTrace curTrace = IOTrace(entry.first);
                IOTrace transitionTrace = IOTrace(transition->getLabel()->getInput(), transition->getLabel()->getOutput(), fsm->getPresentationLayer());
                curTrace.prepend(transitionTrace);
                result.push_back(std::make_pair(curTrace,entry.second));
            }
        }
        return result;    
    };

    return helper(node,initialMissingVisits);
}


std::shared_ptr<InputTree> StrongReductionTestSuiteGenerator::augmentToRDistSet(std::shared_ptr<FsmNode> n1, std::shared_ptr<FsmNode> n2, std::shared_ptr<InputTree> currentlyAppliedSequences) {
    std::shared_ptr<InputTree> result = make_shared<InputTree>(fsm->getPresentationLayer());

    // do nothing if the currently applied sequences are already sufficient to r-distinguish n1 and n2
    if (n1->idRDistinguishedBy(n2,currentlyAppliedSequences)) {
        return result;
    } 

    // TODO: currently uses no heuristic, just applied the pre-calculated r-dist set for n1 and n2
    return rDistTrees[std::make_pair(n1,n2)];
}


InputTree StrongReductionTestSuiteGenerator::initialTestSuite(int m) {
    InputTree result(fsm->getPresentationLayer());

    LOG("VERBOSE_2") << "initialTestSuite" << endl;

    for (auto drEntry : dReachingSequences) {
        auto node = drEntry.first;
        auto v = drEntry.second;
        result.addToRoot(v);

        auto vTrace = InputTrace(v,fsm->getPresentationLayer());
        auto travSet = calcTraversalSet(node,m);

        LOG("VERBOSE_2") << "\tnode " << node->getId() << " with d-r sequence " << vTrace << endl;


        IOListContainer cont(fsm->getPresentationLayer());
        for (auto trEntry : travSet) {
            cont.add(trEntry.first.getInputTrace());
            LOG("VERBOSE_2") << "\t\tadding " << trEntry.first << endl;
        }
        result.addAfter(vTrace,cont);

        LOG("VERBOSE_2") << "\tintermediate result: " << endl << result << endl;        
    }
    
    LOG("VERBOSE_2") << "finished initialising test suite" << endl;

    return result;
}

void StrongReductionTestSuiteGenerator::updateTestSuite(const std::shared_ptr<FsmNode> node, const std::pair<IOTrace, std::vector<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>& nextElementOfD, InputTree& currentTestSuite) {
    // TODO: add heuristic to choose an S, currently the first one is always chosen
    auto trace = nextElementOfD.first;
    auto rdSet = *nextElementOfD.second.cbegin();
    
    auto vTrace = InputTrace(dReachingSequences[node],fsm->getPresentationLayer());

    LOG("VERBOSE_2") << "update for node " << node->getId() << " (d-reached by " << vTrace << "), traversal-trace " << trace << " and set { ";
    for (auto n : *rdSet) { LOG("VERBOSE_2") << n->getId() << " "; }
    LOG("VERBOSE_2") << "}" << endl;

    // TODO: maybe pre-calculate
    std::unordered_set<std::shared_ptr<FsmNode>> drrdNodes; // nodes in the rdSet that are also d-reachable
    LOG("VERBOSE_2") << "\tdrrdNodes are { "; 
    for (auto drEntry : dReachingSequences) {  
        if (rdSet->count(drEntry.first) != 0) {
            LOG("VERBOSE_2") << drEntry.first->getId() << " ";
            drrdNodes.insert(drEntry.first);
        }
    }
    LOG("VERBOSE_2") << "}" << endl;

    // check all pairs of prefixes of trace and the statecover
    // first consider all prefixes of trace against its proper prefixes and against the state cover
    // note: prefixes of trace need to be applied AFTER the d-reaching sequence of node
    auto prefixes = trace.getPrefixes();
    prefixes.push_back(IOTrace(fsm->getPresentationLayer())); // add empty sequence to the prefix list
    for (unsigned int i = 0; i < prefixes.size(); ++i) {
        auto trace1 = prefixes[i];  
        InputTrace preTrace1(trace1.getInputTrace());
        preTrace1.prepend(vTrace);      

        // TODO: pre-compute targets
        auto target1 = *node->after(trace1).cbegin(); // assumes that the trace is in the language of node
        LOG("VERBOSE_2") << "\tcheck trace1 " << trace1 << " reaching " << target1->getId() << " (pre-trace: " << preTrace1 << ")" << endl;
        if (rdSet->count(target1) == 0) {
            LOG("VERBOSE_2") << "\t\treaches non-rd target" << endl;   
            continue; // skip transitions that reach states not in the rdSet 
        }; 

        for (unsigned int j = i+1; j < prefixes.size(); ++j) {
            auto trace2 = prefixes[j];
            auto target2 = *node->after(trace2).cbegin();
            InputTrace preTrace2(trace2.getInputTrace());
            preTrace2.prepend(vTrace);
            LOG("VERBOSE_2") << "\tagainst trace2 " << trace2 << " reaching " << target2->getId() << " (pre-trace: " << preTrace2 << ")" << endl;
            if (target1 == target2 || rdSet->count(target2) == 0) {
                LOG("VERBOSE_2") << "\t\treaches same target as trace1 or a non-rd target" << endl;
                continue;
            }

            auto sharedTraceExtensions = currentTestSuite.sharedExtensions(preTrace1, preTrace2);
            LOG("VERBOSE_2") << "\t\tshared extensions in current test suite:" << endl;
            LOG("VERBOSE_2") << *sharedTraceExtensions;

            auto testSuiteExtension = augmentToRDistSet(target1, target2, sharedTraceExtensions);
            LOG("VERBOSE_2") << "\t\textension to be applied after both traces:" << endl;
            LOG("VERBOSE_2") << *testSuiteExtension;

            IOListContainer cont(fsm->getPresentationLayer());
            for (auto ext : testSuiteExtension->getInputTraces()) {
                cont.add(ext);
            }
            currentTestSuite.addAfter(preTrace1,cont);
            currentTestSuite.addAfter(preTrace2,cont);

            LOG("VERBOSE_2") << "\tintermediate result: " << endl << currentTestSuite << endl;
        }

        for (auto drrdNode : drrdNodes) {
            auto drTrace = InputTrace(dReachingSequences[drrdNode],fsm->getPresentationLayer());
            LOG("VERBOSE_2") << "\tagainst d-r r-d node " << drrdNode->getId() << " with d-r sequence " << drTrace << endl;

            if (target1 == drrdNode) {
                LOG("VERBOSE_2") << "\t\tis target of trace1" << endl;
                continue;
            }            

            auto sharedTraceExtensions = currentTestSuite.sharedExtensions(preTrace1, drTrace);
            LOG("VERBOSE_2") << "\t\tshared extensions in current test suite:" << endl;
            LOG("VERBOSE_2") << *sharedTraceExtensions;

            auto testSuiteExtension = augmentToRDistSet(target1, drrdNode, sharedTraceExtensions);
            LOG("VERBOSE_2") << "\t\textension to be applied after both traces:" << endl;
            LOG("VERBOSE_2") << *testSuiteExtension;

            IOListContainer cont(fsm->getPresentationLayer());
            for (auto ext : testSuiteExtension->getInputTraces()) {
                cont.add(ext);
            }
            currentTestSuite.addAfter(preTrace1,cont);
            currentTestSuite.addAfter(drTrace,cont);

            LOG("VERBOSE_2") << "\tintermediate result: " << endl << currentTestSuite << endl;
        }
    }

    // second consider pairs of sequences in the state cover
    // translating drrdNodes to a vector to avoid mirrored comparisons
    std::vector<std::shared_ptr<FsmNode>> drrdNodesVector;
    for (auto drrdNode : drrdNodes) {
        drrdNodesVector.push_back(drrdNode);
    }
    for (unsigned int i = 0; i < drrdNodesVector.size()-1; ++i) {
        auto drrdNode1 = drrdNodesVector[i];
        auto drTrace1 = InputTrace(dReachingSequences[drrdNode1],fsm->getPresentationLayer());
        LOG("VERBOSE_2") << "\tcheck d-r node " << drrdNode1->getId() << " with d-r sequence " << drTrace1 << endl;

        for (unsigned int j = i+1; j < drrdNodesVector.size(); ++j) {
            auto drrdNode2 = drrdNodesVector[j];
            auto drTrace2 = InputTrace(dReachingSequences[drrdNode2],fsm->getPresentationLayer());
            LOG("VERBOSE_2") << "\tagainst d-r node " << drrdNode2->getId() << " with d-r sequence " << drTrace2 << endl;

            if (drrdNode1 == drrdNode2) {
                LOG("VERBOSE_2") << "\t\tis same as drrdNode1" << endl;
                continue;
            }            

            auto sharedTraceExtensions = currentTestSuite.sharedExtensions(drTrace1, drTrace2);
            LOG("VERBOSE_2") << "\t\tshared extensions in current test suite:" << endl;
            LOG("VERBOSE_2") << *sharedTraceExtensions;

            auto testSuiteExtension = augmentToRDistSet(drrdNode1, drrdNode2, sharedTraceExtensions);
            LOG("VERBOSE_2") << "\t\textension to be applied after both traces:" << endl;
            LOG("VERBOSE_2") << *testSuiteExtension;

            IOListContainer cont(fsm->getPresentationLayer());
            for (auto ext : testSuiteExtension->getInputTraces()) {
                cont.add(ext);
            }
            currentTestSuite.addAfter(drTrace1,cont);
            currentTestSuite.addAfter(drTrace2,cont);

            LOG("VERBOSE_2") << "\tintermediate result: " << endl << currentTestSuite << endl;
        }
    }
}


InputTree StrongReductionTestSuiteGenerator::generateTestSuite(int m) {
    InputTree ts = initialTestSuite(m);        

    // for all d-reachable states s ...
    for (auto drEntry : dReachingSequences) {
        auto node = drEntry.first;
        auto travSet = calcTraversalSet(node,m);

        // ... and all (s,trace,rdSets) in Tr(s,m) ...
        for (auto trEntry : travSet) {

            // ... update the test suite via on-the-fly extension with r-distinguishing sets
            updateTestSuite(node, trEntry, ts);
        }
    }
    
    return ts;
}