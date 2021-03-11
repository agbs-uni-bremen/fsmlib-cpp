/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_STRONGSEMIREDUCTIONTESTSUITEGENERATOR_H_
#define FSM_FSM_STRONGSEMIREDUCTIONTESTSUITEGENERATOR_H_

// TODO: remove superfluous includes
#include <memory>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"

class Fsm;
class FsmLabel;
class FsmNode;
class FsmTransition;
class IOTrace;
class OutputTrace;
class FsmPresentationLayer;
class Tree;
class InputTree;
class IOTreeContainer;
class IOListContainer;
class InputTrace;
class IOTraceContainer;

class StrongSemiReductionTestSuiteGenerator
{
protected:
    
    const std::shared_ptr<Fsm> fsm;

    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>> rDistGraph;
    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::shared_ptr<InputTree>> rDistTrees;
    std::unordered_map<std::shared_ptr<FsmNode>, std::vector<int>> dReachingSequences;
    std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> maximalRDistSets;

    /**
     * Calculates a mapping from pairs of states to inputs and target sets that 
     * represents the graph such that any contained pair of states (s1,s2) that
     * maps to (x, tgts) can be r-distinguished by applying x and after x recursively
     * applying the inputs obtained by applying the mapping to the targets.
     * 
     * If x is -1, then the states are already r(0)-distinguishably (they differ in
     * their defined inputs) and no inputs needs to be applied.
     */ 
    void calcRDistinguishingGraph();

    

    /**
     * Creates a mapping that maps each pair of r-distinguishable states to a
     * set of input sequences that r-distinguishes them.
     */ 
    void calcRDistinguishingTrees();

    

    /**
     * Calculates a vector of states and input sequences (q,xs) is contained
     * in the result if and only if q is a state of this FSM and is deterministically
     * reachable (for possibly partial FSMs) via xs.
     * Assumes that "this" is observable.
     */
    void calcDeterministicallyReachingSequences();

    
    /**
     * Calculates a set of sets of states of fsm that are pairwise r-distinguishable
     * such that every state of fsm is contained in at least one of the resulting sets.
     * 
     * NOTE: does not solve the maximal clique problem and hence may not return all
     *       maximal sets of states of fsm that are r-distinguishable.
     */
    void calcMaximalRDistinguishableSets();

    /**
     * Calculates the set of all maximal sets of states of fsm that are pairwise r-distinguishable.
     * 
     * NOTE: performs a brute force calculation that should be replaced by
     *       a more sophisticated algorithm to calculate maximal pairwise r-distinguishable sets.
     */
    void calcAllMaximalRDistinguishableSets();

    
public:
    
    StrongSemiReductionTestSuiteGenerator(const Fsm& other, bool calculateAllMaximalRDistinguishableSets = false);

    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::shared_ptr<InputTree>> getRDistinguishingTrees() const;
    std::unordered_map<std::shared_ptr<FsmNode>, std::vector<int>> getDeterministicallyReachingSequences() const;
    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>> getRDistGraph() const;
    std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> getMaximalRDistinguishableSets();
     
};





#endif  
