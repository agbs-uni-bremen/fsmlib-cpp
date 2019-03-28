/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/PartitionGraphEdge.h"
#include "PartitionGraphEdge.h"


PartitionGraphEdge::PartitionGraphEdge(const vector<int> &trace, const weak_ptr<Node> &source,
                                       const weak_ptr<Node> &target,
                                       shared_ptr<unordered_map<int, int>> &blockToTarget)
                                       :blockToTarget(blockToTarget),Edge(trace,source,target)
{

}

shared_ptr<unordered_map<int, int>> &PartitionGraphEdge::getBlockToTarget() {
    return blockToTarget;
}
