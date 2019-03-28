/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#include "graphs/PartitionGraph.h"
#include "graphs/PartitionGraphNode.h"
#include <queue>
#include <cassert>
#include <utility>

PartitionGraph::PartitionGraph()
    : Graph(vector<shared_ptr<Node>>())
{

}

void PartitionGraph::addEdge(shared_ptr<SplittingTreeNode> &nodeA, shared_ptr<SplittingTreeNode> &nodeB, int x) {

    shared_ptr<PartitionGraphNode> partitionNodeA,
                                partitionNodeB;

    for(auto& node:nodes){
        auto castedNode = dynamic_pointer_cast<PartitionGraphNode>(node);
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


    auto edge = make_shared<Edge>(vector<int>{x},partitionNodeA,partitionNodeB);

}

shared_ptr<SplittingTreeNode> PartitionGraph::findPathToAOrBValidNode(shared_ptr<SplittingTreeNode> &node, vector<int>& cTrace) {

    shared_ptr<PartitionGraphNode> partitionNode;

    for(auto& node:nodes) {
        auto castedNode = dynamic_pointer_cast<PartitionGraphNode>(node);
        castedNode->setVisited(false);
        if(castedNode->getBlock() == node)
            partitionNode = castedNode;
    }
    if(!partitionNode)
        return shared_ptr<SplittingTreeNode>();

    shared_ptr<SplittingTreeNode> targetNode;
    //length of current path P to targetNode, such that |label(P).(targetNode.trace)| is minimal
    int minLength = 0;

    queue< shared_ptr< pair< vector<int>,shared_ptr<PartitionGraphNode> > > > workingQueue;
    workingQueue.push(make_shared<pair<vector<int>,shared_ptr<PartitionGraphNode>>>(vector<int>(),partitionNode));
    partitionNode->setVisited(true);

    while(!workingQueue.empty()) {
        auto currentNode = workingQueue.front();
        workingQueue.pop();

        for(auto& edge:currentNode->second->getEdges()) {
            assert(edge->getTarget().lock());

            auto currentTargetNode = dynamic_pointer_cast<PartitionGraphNode>(edge->getTarget().lock());
            auto& currentTargetBlock = currentTargetNode->getBlock();
            if(!currentTargetNode->isVisited()) {
                auto tempTrace = currentNode->first;
                tempTrace.push_back(edge->getTrace().front());
                if(currentTargetBlock->getIsAValid() || currentTargetBlock->getIsBValid()) {
                    if(!targetNode || (tempTrace.size() + currentTargetBlock->getTrace().size() < minLength)) {
                        targetNode = currentTargetBlock;
                        cTrace = tempTrace;
                        minLength = cTrace.size() + currentTargetBlock->getTrace().size();
                    }
                } else {
                    currentTargetNode->setVisited(true);
                    workingQueue.push(make_shared<pair<vector<int>,shared_ptr<PartitionGraphNode>>>(tempTrace,currentTargetNode));
                }
            }
        }
    }

}
