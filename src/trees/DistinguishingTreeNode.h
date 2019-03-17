/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#ifndef FSM_TREES_DISTINGUISHINGTREENODE_H_
#define FSM_TREES_DISTINGUISHINGTREENODE_H_

#include <vector>
#include <set>
#include "trees/DistinguishingTreeEdge.h"
#include "trees/DistinguishingTreeNode.h"

using namespace std;

class DistinguishingTreeEdge;

class DistinguishingTreeNode: public enable_shared_from_this<DistinguishingTreeNode> {
protected:

    /**
    * The parent of this node. It is nullptr if this is a root node.
    */
    weak_ptr<DistinguishingTreeNode> parent;

    /**
    * The input sequence that is labelling the path from the root to this node. For
    * performance reasons the input sequence is redundantly stored at the node.
    */
    vector<int> inputTrace;

    /**
    * The current state uncertainty multiset, that states, in which state the automaton might be after
    * execution of the associated input sequence and the initial state being unknown. For each possible
    * output sequence exists a Set of IDs of the DFSM states, that are the possible current states for this.
    * If there are exactly |Number of DFSM states|  singleton sets, the input trace associated with this
    * node is a distinguishing sequence.
    */
    multiset<set<int>> currentUncertainty;

    /**
	* A list containing the children of this node
	*/
    vector<shared_ptr<DistinguishingTreeEdge>> children;

public:

    /**
    * Create a new distinguishing tree node
    */
    DistinguishingTreeNode();

    /**
    * Add an edge to this nodes children
    * @param edge The edge to be added
    */
    void add(const shared_ptr<DistinguishingTreeEdge>& edge);

    /**
    * Set the parent of this node
    * @param parent the new parent of this node
    */
    void setParent(const weak_ptr<DistinguishingTreeNode>& parent);

    /**
    * Gets the input trace associated with this node.
    * @return the input trace of this node
    */
    vector<int>& getInputTrace();

    /**
    * Gets the current uncertainty of this node
    * @return the current uncertainty of this node
    */
    multiset<set<int>>& getCurrentUncertainty();

};

#endif //FSM_TREES_DISTINGUISHINGTREENODE_H_
