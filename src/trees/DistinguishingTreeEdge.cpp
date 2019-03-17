/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "trees/DistinguishingTreeEdge.h"

using namespace std;

DistinguishingTreeEdge::DistinguishingTreeEdge(const int input, const shared_ptr<DistinguishingTreeNode> &target)
    : input(input),target(target)
{

}

shared_ptr<DistinguishingTreeNode> DistinguishingTreeEdge::getTarget() {
    return target;
}