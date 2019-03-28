/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_GRAPHS_PARTITIONGRAPHEDGE_H
#define FSM_GRAPHS_PARTITIONGRAPHEDGE_H

#include "graphs/Edge.h"
#include <unordered_map>

class PartitionGraphEdge:public Edge {
protected:

    /**
    * An auxilliary mapping of a subset of fsm states to the target states, that are reached after reading 'trace'
    */
    shared_ptr<unordered_map< int, int>> blockToTarget;
public:

    /**
    * Creates a new edge
    * @param trace the trace associated with this edge
    * @param source the source node of this edge
    * @param target the target node of this edge
    * @param blockToTarget the auxilliary mapping of the states in the underlying block of source to the states of target
    */
    PartitionGraphEdge(const vector<int>& trace, const weak_ptr<Node>& source, const weak_ptr<Node>& target,shared_ptr<unordered_map< int, int>>& blockToTarget);


    /**
    * Gets the auxilliary mapping of this edge
    * @return the auxilliary mapping of this edge
    */
    shared_ptr<unordered_map< int, int>>& getBlockToTarget();
};

#endif //FSM_GRAPHS_PARTITIONGRAPHEDGE_H
