/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */


#include "graphs/PartitionGraph.h"
#include "graphs/PartitionGraphNode.h"

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

        }
        //create new partition graph node; dont care for id
        partitionNodeB = make_shared<PartitionGraphNode>(0,nodeB);

    }


    auto edge = make_shared<Edge>(vector<int>{x},partitionNodeA,partitionNodeB);

}
