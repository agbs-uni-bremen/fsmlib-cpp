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
class OutputTree;
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
    // NOTE: as we only consider input sequences and all targets of a d-reaching sequences for state q are q, we do not explicitly store V'
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
     * in the result only if q is a state of this FSM and is deterministically
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
    
    StrongSemiReductionTestSuiteGenerator(const std::shared_ptr<Fsm> other, bool calculateAllMaximalRDistinguishableSets = false);

    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::shared_ptr<InputTree>> getRDistinguishingTrees() const;
    std::unordered_map<std::shared_ptr<FsmNode>, std::vector<int>> getDeterministicallyReachingSequences() const;
    std::unordered_map<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>, std::pair<int,std::shared_ptr<std::unordered_set<std::pair<std::shared_ptr<FsmNode>,std::shared_ptr<FsmNode>>>>>> getRDistGraph() const;
    std::unordered_set<std::unordered_set<std::shared_ptr<FsmNode>>> getMaximalRDistinguishableSets();


    /**
     * Calculate all tuples (S,rep) where S is a maximal pairwise r-distinguishable set of states of fsm,
     * and rep = m - dr + 1, where dr is the number of d-reachable states in S.
     */
    std::vector<std::pair<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>,int>> getTerminationTuples(int m);

    /**
     * Generate the traversal set T(s,m).
     */
    std::vector<std::pair<IOTrace, std::unordered_set<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>> calcTraversalSet(std::shared_ptr<FsmNode> node, int m);

    /**
     * Create a set W of input sequences that r-distinguishes the given states.
     * If the inputs currently applied after both states are already r-distinguishing, then the returned set is empty.
     * 
     * TODO: use heuristic to choose an augmenting set that is small, rather than just using the pre-calculated r-distinguishing set
     */
    std::shared_ptr<InputTree> augmentToRDistSet(std::shared_ptr<FsmNode> n1, std::shared_ptr<FsmNode> n2, std::shared_ptr<InputTree> currentlyAppliedSequences);

    /**
     * Compute the initial test suite by extending for each d-reachable state s of fsm its d-reaching sequence with all inputs in Tr(s,m).
     */ 
    InputTree initialTestSuite(int m);

    /**
     * Updates the test suite for a given terminated pair.
     * 
     * Assumes that the testSuite already contains the result of initialTestSuite(m).
     */
    void updateTestSuite(const std::shared_ptr<FsmNode> node, const std::pair<IOTrace, std::unordered_set<std::shared_ptr<std::unordered_set<std::shared_ptr<FsmNode>>>>>& nextElementOfD, InputTree& currentTestSuite);

    /**
     * Generate an m-complete test suite for fsm.
     */ 
    InputTree generateTestSuite(int m);
};





#endif  
