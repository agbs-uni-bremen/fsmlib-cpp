/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "SplittingTreeEdge.h"

using namespace std;

SplittingTreeEdge::SplittingTreeEdge(const int output, const shared_ptr<SplittingTreeNode> &target)
        : output(output),target(target)
{

}

shared_ptr<SplittingTreeNode> SplittingTreeEdge::getTarget() {
    return target;
}

int SplittingTreeEdge::getOutput() {
    return output;
}
