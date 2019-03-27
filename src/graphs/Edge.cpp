/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/Edge.h"

Edge::Edge(const vector<int> &trace, const weak_ptr<Node> &source, const weak_ptr<Node> &target)
    : trace(trace), source(source), target(target)
{
    if(!source.expired())
        source.lock()->addEdge(shared_from_this());

    if(!target.expired())
        target.lock()->addInEdge(shared_from_this());
}