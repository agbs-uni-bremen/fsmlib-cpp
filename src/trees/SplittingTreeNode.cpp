/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "SplittingTreeNode.h"

SplittingTreeNode::SplittingTreeNode()
    :parent(weak_ptr<SplittingTreeNode>())
{

}

SplittingTreeNode::SplittingTreeNode(const set<int>& block)
    : block(block)
{

}
void SplittingTreeNode::setParent(const weak_ptr<SplittingTreeNode>& parent) {
    this->parent = parent;
}

void SplittingTreeNode::add(const shared_ptr<SplittingTreeEdge> &edge) {
    edge->getTarget()->setParent(shared_from_this());
    children.push_back(edge);
}

set<int> &SplittingTreeNode::getBlock() {
    return block;
}

void SplittingTreeNode::setTrace(const vector<int> &trace) {
    this->trace = trace;
}
