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
#include <unordered_map>

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
    shared_ptr<unordered_map< int, int>> blockToTarget;

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
    * Create a new distinguishing tree node
    * @param block the block of a partition, that is associated with this subtree
    */
    SplittingTreeNode(const set<int>& block);

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

    /**
    * Sets the input trace associated with this node
    */
    void setTrace(const vector<int>& trace);

    /**
    * Sets the auxilliary mapping of states of `block`
    */
    void setBlockToTarget(const shared_ptr<unordered_map<int,int>>& blockToTarget);


};

#endif //FSM_TREES_SPLITTINGTREENODE_H_
