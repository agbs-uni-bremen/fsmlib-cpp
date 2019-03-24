/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_TREES_SPLITTINGTREE_H_
#define FSM_TREES_SPLITTINGTREE_H_

#include <vector>
#include <memory>
#include "fsm/Dfsm.h"
#include "trees/SplittingTreeNode.h"

using namespace std;

/**
 * The splitting tree is a datastructure to compute an adaptive distinguishing sequence for
 * a dfsm, that distinguishes every state from every other state. For more information the
 * book Model-Based Testing of Reactive Systems: Advanced Lectures by Broy, Manfred (Chapter 2) is
 * recommended.
 */
class SplittingTree final
{
private:
    /**
    * The root of this tree
    */
    shared_ptr<SplittingTreeNode> root;

    /**
    * The adaptive distinguishing sequence associated with this splitting tree, if one exists.
    */
    shared_ptr<InputOutputTree> adaptiveDistinguishingSequence;


public:

    /**
	* Creates a new distinguishing tree from a DFSM, if it is completely specified. A distinguishing sequence
    * is derived from it, if one exists.
	* @param dfsm  DFSM, for which the distinguishing tree is being created
	*/
    SplittingTree(const shared_ptr<Dfsm>& dfsm);
};
#endif //FSM_TREES_SPLITTINGTREE_H_
