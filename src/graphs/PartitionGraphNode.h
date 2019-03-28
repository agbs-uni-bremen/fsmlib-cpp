/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_GRAPHS_PARTITIONGRAPHNODE_H_
#define FSM_GRAPHS_PARTITIONGRAPHNODE_H_

#include <memory>
#include "graphs/Node.h"
#include "trees/SplittingTreeNode.h"

class PartitionGraphNode: public Node {
protected:

    /**
    * The block of a partition this node represents in form of an underlying splitting tree node, this node manages.
    */
    shared_ptr<SplittingTreeNode> block;

public:

    /**
    * Create new partitiongraph node.
    * @param id the id of the new node
    * @param block the block of a partition the new node represents
    */
    PartitionGraphNode(const int id,const shared_ptr<SplittingTreeNode>& block);

    /**
    * Gets the block this node represents
    * @return the block this node represents
    */
    shared_ptr<SplittingTreeNode>& getBlock();

};

#endif //FSM_GRAPHS_PARTITIONGRAPHNODE_H_
