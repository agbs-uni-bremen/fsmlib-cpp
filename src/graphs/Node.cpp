/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "graphs/Node.h"

using namespace std;

Node::Node(const int id)
    : id(id)
{

}

void Node::addInEdge(const shared_ptr<Edge> &inEdge) {
    inEdges.push_back(inEdge);
}

void Node::addEdge(const shared_ptr<Edge> &edge) {
    edges.push_back(edge);
}

int Node::getId() {
    return id;
}