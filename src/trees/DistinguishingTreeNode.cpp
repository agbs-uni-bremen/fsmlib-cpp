/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "trees/DistinguishingTreeNode.h"

using namespace std;

DistinguishingTreeNode::DistinguishingTreeNode()
    : parent(weak_ptr<DistinguishingTreeNode>()),inputTrace(vector<int>()),currentUncertainty(multiset<set<int>>()),
    children(vector<shared_ptr<DistinguishingTreeEdge>>())
{

}
void DistinguishingTreeNode::setParent(const weak_ptr<DistinguishingTreeNode>& parent) {
    this->parent = parent;
}

multiset<set<int>>& DistinguishingTreeNode::getCurrentUncertainty() {
    return currentUncertainty;
}

vector<int>& DistinguishingTreeNode::getInputTrace() {
    return inputTrace;
}

void DistinguishingTreeNode::add(const shared_ptr<DistinguishingTreeEdge> &edge) {
    edge->getTarget()->setParent(shared_from_this());
    children.push_back(edge);
}
