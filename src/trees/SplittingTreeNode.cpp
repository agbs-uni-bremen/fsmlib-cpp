/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "SplittingTreeNode.h"

SplittingTreeNode::SplittingTreeNode()
    :parent(weak_ptr<SplittingTreeNode>()),isAValid(false),isBValid(false)
{

}

SplittingTreeNode::SplittingTreeNode(const set<int>& block)
    : block(block),parent(weak_ptr<SplittingTreeNode>()),isAValid(false),isBValid(false)
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

vector<int> &SplittingTreeNode::getTrace() {
    return trace;
}

bool &SplittingTreeNode::getIsAValid() {
    return isAValid;
}

bool &SplittingTreeNode::getIsBValid() {
    return isBValid;
}

shared_ptr<unordered_map<int, int>> &SplittingTreeNode::getBlockToTarget() {
    return blockToTarget;
}

void SplittingTreeNode::setBlockToTarget(const shared_ptr<unordered_map<int, int>> &blockToTarget) {
    this->blockToTarget = blockToTarget;
}

weak_ptr<SplittingTreeNode> &SplittingTreeNode::getParent() {
    return parent;
}

vector<shared_ptr<SplittingTreeEdge>> &SplittingTreeNode::getChildren() {
    return children;
}

