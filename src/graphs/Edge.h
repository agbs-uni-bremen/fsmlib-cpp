/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_GRAPHS_EDGE_H_
#define FSM_GRAPHS_EDGE_H_

#include <vector>
#include <memory>
#include <utility>
#include "graphs/Node.h"

using namespace std;


class Edge {
protected:

    /**
    The input trace of this edge
    */
    vector<int> trace;

    /**
    The target node of this edge
    */
    shared_ptr<Node> target;

public:

    /**
    * Creates a new edge
    * @param trace the trace associated with this edge
    * @param target the target node of this edge
    */
    Edge(const vector<int>& trace, const shared_ptr<Node>& target);

};
#endif //FSM_GRAPHS_EDGE_H_
