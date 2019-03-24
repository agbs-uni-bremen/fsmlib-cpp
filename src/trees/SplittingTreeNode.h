/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_TREES_SPLITTINGTREENODE_H_
#define FSM_TREES_SPLITTINGTREENODE_H_

#include <memory>
#include <vector>
#include <set>

#include "trees/SplittingTreeEdge.h"

using namespace std;

class SplittingTreeEdge;

class SplittingTreeNode: public enable_shared_from_this<SplittingTreeNode>  {
protected:

    /**
    * The block of the partition of states of a dfsm, that is induced by this subtree.
    */
    set<int> block;

    /**
    * The input sequence that is associated with this node.
    */
    vector<int> trace;

    /**
    * An auxilliary mapping of the states in 'block' to the target states, that are reached after reading 'trace'
    */
    vector<int> blockToTarget;

    /**
    * The parent of this node. It is nullptr if this is a root node.
    */
    weak_ptr<SplittingTreeNode> parent;

    /**
	* A list containing the children of this node
	*/
    vector<shared_ptr<SplittingTreeEdge>> children;

public:

    /**
    * Create a new distinguishing tree node
    */
    SplittingTreeNode();

    /**
    * Add an edge to this nodes children
    * @param edge The edge to be added
    */
    void add(const shared_ptr<SplittingTreeEdge>& edge);

    /**
    * Set the parent of this node
    * @param parent the new parent of this node
    */
    void setParent(const weak_ptr<SplittingTreeNode>& parent);

    /**
    * Returns the associated block of this node.
    */
    set<int>& getBlock();


};

#endif //FSM_TREES_SPLITTINGTREENODE_H_
