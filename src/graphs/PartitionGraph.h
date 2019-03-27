/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_GRAPHS_PARTITIONGRAPH_H_
#define FSM_GRAPHS_PARTITIONGRAPH_H_

#include <vector>
#include <memory>
#include "fsm/Dfsm.h"
#include "graphs/Graph.h"
#include "trees/SplittingTreeNode.h"

using namespace std;

class PartitionGraph:public Graph {
protected:

public:

    /**
    * Creates a new partitiongraph
    */
    PartitionGraph();

    /**
    * Creates nodes for partition blocks represented by `nodeA` and `nodeB` respectively, if they do not exist, and adds an edge from
    * `nodeA` to `nodeB` with label `x`
    * @param nodeA the corresponding splitting tree node that contains the block, which represents the node from which the new edge starts
    * @param nodeB the corresponding splitting tree node that contains the block, which represents the node at which the new edge ends
    * @param x the label of the new edge
    */
    void addEdge(shared_ptr<SplittingTreeNode>& nodeA, shared_ptr<SplittingTreeNode>& nodeB,int x);
};


#endif //FSM_GRAPHS_PARTITIONGRAPH_H_
