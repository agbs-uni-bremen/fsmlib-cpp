/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "graphs/Node.h"
#include "Node.h"


using namespace std;

Node::Node(const int id)
    : id(id),visited(false)
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

void Node::setVisited(bool visited) {
    this->visited = visited;
}

bool Node::isVisited() {
    return visited;
}

vector<shared_ptr<Edge>> &Node::getEdges() {
    return edges;
}

const vector<weak_ptr<Edge>> &Node::getInEdges() const {
    return inEdges;
}
