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

    /**
    * A mapping of the FsmNode ids to the actual FsmNode object. It is used to help speed up the access
    */
    vector<shared_ptr<FsmNode>> idToFsmNode;

    /**
    * checks if an input is a-valid for the underlying block of the node
    * @param blockNode the node, whose underlying block is to be checked against x for a-validity
    * @param x the input that is checked for a-validity against `blockNode`
    * @param blockPartition it is the partition of the block induced by the outputs of `x`, while reading `x` in the states of `blockNode`, if `x` is a-valid for `blockNode`
    * @param blockToTarget mapping of the block in `blockNode` to the target states block
    * @param isValid flag that is set, if the input `x` is valid at all for the underlying block of `blockNode`
    * @return true, if `x` is a-valid, otherwise false
    */
    bool checkAValid(shared_ptr<SplittingTreeNode> &blockNode, const int x,
                     unordered_map<int, shared_ptr<set<int>>> &blockPartition,
                     shared_ptr<unordered_map<int, int>> &blockToTarget, bool &isValid);

public:

    /**
	* Creates a new distinguishing tree from a DFSM, if it is completely specified. A distinguishing sequence
    * is derived from it, if one exists.
	* @param dfsm  DFSM, for which the distinguishing tree is being created
	*/
    SplittingTree(const shared_ptr<Dfsm>& dfsm);
};
#endif //FSM_TREES_SPLITTINGTREE_H_
