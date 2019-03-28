/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#include "graphs/PartitionGraph.h"
#include "graphs/PartitionGraphNode.h"
#include "graphs/PartitionGraphEdge.h"
#include <queue>
#include <cassert>
#include <utility>

PartitionGraph::PartitionGraph()
    : Graph(vector<shared_ptr<Node>>())
{

}

void PartitionGraph::addEdge(shared_ptr<SplittingTreeNode> &nodeA, shared_ptr<SplittingTreeNode> &nodeB, int x,shared_ptr<unordered_map< int, int>>& blockToTarget) {

    shared_ptr<PartitionGraphNode> partitionNodeA,
                                partitionNodeB;

    for(auto& node:nodes){
        auto castedNode = static_pointer_cast<PartitionGraphNode>(node);
        if(castedNode->getBlock() == nodeA) {
            partitionNodeA = castedNode;
        } else if(castedNode->getBlock() == nodeB) {
            partitionNodeB = castedNode;
        }
        if(partitionNodeA && partitionNodeB) break;
    }

    if(!partitionNodeA) {
        //create new partition graph node; dont care for id
        partitionNodeA = make_shared<PartitionGraphNode>(0,nodeA);
        nodes.push_back(partitionNodeA);
    }

    if(!partitionNodeB) {
        if(nodeA == nodeB) {
            partitionNodeB = partitionNodeA;
        } else {
            //create new partition graph node; dont care for id
            partitionNodeB = make_shared<PartitionGraphNode>(0,nodeB);
            nodes.push_back(partitionNodeB);
        }
    }


    auto edge = make_shared<PartitionGraphEdge>(vector<int>{x},partitionNodeA,partitionNodeB,blockToTarget);

}

shared_ptr<SplittingTreeNode> PartitionGraph::findPathToAOrBValidNode(shared_ptr<SplittingTreeNode> &sNode, vector<int>& cTrace,shared_ptr<unordered_map< int, int>>& blockToTarget) {

    shared_ptr<PartitionGraphNode> partitionNode;

    for(auto& n:nodes) {
        auto castedNode = static_pointer_cast<PartitionGraphNode>(n);
        castedNode->setVisited(false);
        if(castedNode->getBlock() == sNode)
            partitionNode = castedNode;
    }
    if(!partitionNode)
        return shared_ptr<SplittingTreeNode>();

    shared_ptr<SplittingTreeNode> targetNode;
    //length of current path P to targetNode, such that |label(P).(targetNode.trace)| is minimal
    int minLength = 0;

    queue< shared_ptr< tuple< vector<int>,shared_ptr<PartitionGraphNode>,shared_ptr<unordered_map< int, int>> > > > workingQueue;
    workingQueue.push(make_shared<tuple<vector<int>,shared_ptr<PartitionGraphNode>,shared_ptr<unordered_map< int, int>>>>(vector<int>(),partitionNode, blockToTarget));
    partitionNode->setVisited(true);

    while(!workingQueue.empty()) {
        auto currentNode = workingQueue.front();
        workingQueue.pop();

        for(auto& edge: get<1>(*currentNode)->getEdges()) {
            assert(edge->getTarget().lock());
            auto castedEdge = static_pointer_cast<PartitionGraphEdge>(edge);

            auto currentTargetNode = static_pointer_cast<PartitionGraphNode>(edge->getTarget().lock());
            auto& currentTargetBlock = currentTargetNode->getBlock();
            if(!currentTargetNode->isVisited()) {
                auto currentTrace = get<0>(*currentNode);
                currentTrace.push_back(edge->getTrace().front());

                auto currentBlockToTarget = get<2>(*currentNode);
                if(currentBlockToTarget->size() == 0) {
                    currentBlockToTarget = castedEdge->getBlockToTarget();
                } else {
                    for(auto& idToTarget:*currentBlockToTarget) {

                    }
                }

                if(currentTargetBlock->getIsAValid() || currentTargetBlock->getIsBValid()) {
                    if(!targetNode || (currentTrace.size() + currentTargetBlock->getTrace().size() < minLength)) {
                        targetNode = currentTargetBlock;
                        cTrace = currentTrace;
                        minLength = cTrace.size() + currentTargetBlock->getTrace().size();
                    }
                } else {
                    currentTargetNode->setVisited(true);
                    //workingQueue.push(make_shared<tuple<vector<int>,shared_ptr<PartitionGraphNode>>>(currentTrace,currentTargetNode));
                }
            }
        }
    }

    return targetNode;
}
