/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "PartitionGraphNode.h"

PartitionGraphNode::PartitionGraphNode(const int id, const shared_ptr<SplittingTreeNode> &block)
    : block(block), Node(id)
{

}

shared_ptr<SplittingTreeNode>& PartitionGraphNode::getBlock(){
    return block;
}