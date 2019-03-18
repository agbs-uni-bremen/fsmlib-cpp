/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#ifndef FSM_TREES_DISTINGUISHINGTREE_H_
#define FSM_TREES_DISTINGUISHINGTREE_H_

#include <vector>
#include <set>
#include <string>
#include "trees/DistinguishingTreeNode.h"
#include "fsm/Dfsm.h"

using namespace std;

/**
 * The distinguishing tree is a datastructure to compute a distinguishing sequence for
 * a dfsm, that distinguishes every state from every other state. For more information the
 * book Model-Based Testing of Reactive Systems: Advanced Lectures by Broy, Manfred (Chapter 2) is
 * recommended.
 */
class DistinguishingTree final
{
private:
    /**
    * The root of this tree
    */
    shared_ptr<DistinguishingTreeNode> root;

    /**
    * The distinguishing sequence associated with this distinguishing tree, if one exists
    */
    vector<int> distinguishingSequence;

    /**
    * A mapping of the FsmNode ids to the actual FsmNode object. It is used to help speed up the access
    */
    vector<shared_ptr<FsmNode>> idToFsmNode;

    /**
    * Creates a hash value for a current uncertainty multiset.
    */
    int createHashForCurrentUncertainty(const multiset<set<int>>& currentUncertainty);

    /**
    * Create the next current uncertainty based on the given current uncertainty and input.
    * If the input is not valid for 'currentUncertainty' , 'nextCurrentUncertainty' is empty afterwards.
    * @param currentUncertainty the 'current' current uncertainty
    * @param x the next input symbol to apply
    * @param nextCurrentUncertainty reference to the next current uncertainty to be updated
    * @param isTerminal flag that is set to true, if the next input symbol is valid for 'currentUncertainty' and 'nextCurrentUncertainty' consists of singletons only
    */
    void computeNextCurrentUncertainty(const multiset<set<int>>& currentUncertainty,int x,multiset<set<int>>& nextCurrentUncertainty, bool& isTerminal);

    string uncToString(const multiset<set<int>>& currentUncertainty);

    string traceToString(const vector<int>& trace);

public:

    /**
     * Returns the distinguishing sequence associated with this distinguishing tree, if one exists.
     *
     * @return the distinguishing sequence if one exists, otherwise an empty sequence
     */
    const std::vector<int>& getDistinguishingSequence();

    /**
	* Creates a new distinguishing tree from a DFSM, if it is completely specified. A distinguishing sequence
    * is derived from it, if one exists.
	* @param dfsm  DFSM, for which the distinguishing tree is being created
	*/
    DistinguishingTree(const std::shared_ptr<Dfsm>& dfsm);
};

#endif //FSM_TREES_DISTINGUISHINGTREE_H_
