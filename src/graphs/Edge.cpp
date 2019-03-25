/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/Edge.h"

Edge::Edge(const vector<int> &trace, const shared_ptr<Node> &target)
    : trace(trace),target(target)
{
    target->addInEdge(shared_from_this());
}