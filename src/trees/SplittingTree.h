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
#include "AdaptiveTreeNode.h"

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
    * Reference to the dfsm, the splitting tree is build for
    */
    shared_ptr<Dfsm> dfsm;

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
                     map<int, shared_ptr<set<int>>> &blockPartition,
                     shared_ptr<unordered_map<int, int>> &blockToTarget, bool &isValid);

    /**
    * find the corresponding splitting tree node to an adaptive tree node, that contains all
    * the current set of the adaptive tree node
    * @param rootNode the root node of the splitting tree to search in
    * @param adaptiveTreeNode the adaptive tree node, for which the corresponding splitting tree node is searched
    * @return the corresponding splitting tree node of `adaptiveTreeNode`
    */
    shared_ptr<SplittingTreeNode> findDeepestNode(shared_ptr<SplittingTreeNode> &rootNode,shared_ptr<AdaptiveTreeNode> &adaptiveTreeNode);


public:

    /**
	* Creates a new instance of a distinguishing tree from a DFSM, if it is completely specified. A distinguishing sequence
    * is derived from it, if one exists.
	* @param dfsm  DFSM, for which the distinguishing tree has to be created
	*/
    SplittingTree(const shared_ptr<Dfsm>& dfsm);

    /**
    * An auxilliary function to compose auxilliary dfsm state mappings, in course of the creation of the splitting tree
    * @param firstMap the first mapping
    * @param secondMap the second mapping, which is composed to `firstMap`
    * @return the composition mapping of `firstMap` and `secondMap`
    */
    static shared_ptr<unordered_map<int, int>> composeBlockToTarget(shared_ptr<unordered_map<int, int>>& firstMap,shared_ptr<unordered_map<int, int>>& secondMap);

    /**
    * Gets the adaptive distinguishing sequence, derived from the splitting tree. if none exists the smart pointer to it is empty.
    * @return the adaptive distinguishing sequence
    */
    shared_ptr<InputOutputTree>& getAdaptiveDistinguishingSequence();

    /**
    * builds the splitting tree and tries to derive an adaptive distinguishing sequence if one exists
    */
    void build();

    /**
    * gets the root node of this splitting tree
    * @return the root node of this splitting tree
    */
    const shared_ptr<SplittingTreeNode>& getRoot() const;

    /**
    * Prints the splitting tree in .dot format into the given filename (in the current path)
    * @param fname the filename (without extension) to print the tree into
    */
    void toDot(const string& fname);

    friend ostream & operator<<(ostream & out, const SplittingTree & splittingTree);
};
#endif //FSM_TREES_SPLITTINGTREE_H_
