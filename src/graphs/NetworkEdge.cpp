/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "NetworkEdge.h"

NetworkEdge::NetworkEdge(const vector<int> &trace, const weak_ptr<Node> &source, const weak_ptr<Node> &target,
                         int capacity,int cost)
                         :Edge(trace,source,target),capacity(capacity),flow(0),isDs(false),isAlpha(false),
                         isReverse(false),referenceEdge(weak_ptr<NetworkEdge>())
{
    this->cost = cost;
}

bool NetworkEdge::getIsAlpha() {
    return isAlpha;
}

void NetworkEdge::setIsAlpha(bool isAlpha) {
    this->isAlpha = isAlpha;
}

bool NetworkEdge::getIsDs() {
    return isDs;
}

void NetworkEdge::setIsDs(bool isDs) {
    this->isDs = isDs;
}

int NetworkEdge::getCapacity() const {
    return capacity;
}

void NetworkEdge::setCapacity(int capacity) {
    NetworkEdge::capacity = capacity;
}

unsigned int NetworkEdge::getFlow() const {
    return flow;
}

void NetworkEdge::setFlow(unsigned int flow) {
    NetworkEdge::flow = flow;
}

bool NetworkEdge::getIsReverse() const {
    return isReverse;
}

void NetworkEdge::setIsReverse(bool isReverse) {
    NetworkEdge::isReverse = isReverse;
}

const weak_ptr<NetworkEdge> &NetworkEdge::getReferenceEdge() const {
    return referenceEdge;
}

void NetworkEdge::setReferenceEdge(const weak_ptr<NetworkEdge> &reverseEdge) {
    NetworkEdge::referenceEdge = reverseEdge;
}
