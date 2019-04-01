/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/Edge.h"
#include "Edge.h"

#include <cassert>

Edge::Edge(const vector<int> &trace, const weak_ptr<Node> &source, const weak_ptr<Node> &target)
    : trace(trace), source(source), target(target)
{
    //causes bad_weak_ptr exception if shared_from_this is called for an instance, that has no shared_ptr pointing to it yet
    //assert(!source.expired());
   // assert(!target.expired());

    //source.lock()->addEdge(shared_from_this());
    //target.lock()->addInEdge(shared_from_this());
}

weak_ptr<Node> Edge::getSource() {
    return source;
}

weak_ptr<Node> Edge::getTarget() {
    return target;
}

vector<int> &Edge::getTrace() {
    return trace;
}
