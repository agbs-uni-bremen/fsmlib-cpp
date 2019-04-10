/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_GRAPHS_GRAPH_H_
#define FSM_GRAPHS_GRAPH_H_

#include <vector>
#include <deque>
#include <list>
#include "graphs/Node.h"

using namespace std;

class Graph {
protected:

    /**
    * The nodes of this graph
    */
    vector<shared_ptr<Node>> nodes;

    /**
    * Check whether the ids of the nodes are
    * 1. unique
    * 2. respecting the node set size 'n' and dont exceed that boundary
    * 3. matching their id to the indices of #nodes
    *
    * so that we have guarenteed a bijective indexing function from {0,...,n-1} to the set of nodes
    */
    bool validateNodeIds() const;

public:

    /**
    * Creates a new graph
    * @param nodes the nodes of this graph
    *
    * @note it is assumed that the node ids are valid according to the requirements in Graph::validateNodeIds()
     *      otherwise the result of all algorithms are undefined.
    */
    Graph(const vector<shared_ptr<Node>>& nodes);

    /**
     * Compute the shortest path with respect to edge costs
     * from a source node to a target node in the graph by
     * using the *Bellman-Ford* algorithm, if one exists.
     * The algorithm is able to detect negative cost cycles, in which case
     * there is no shortest path either.
     * @param source the source node of the shortest path to search for
     * @param target the target node of the shortest path to search for
     * @result a path from \a source to \a target, if one exists, otherwise an empty path
     *
     */
    shared_ptr<deque<shared_ptr<Edge>>> shortestPathByBellmanFord(const shared_ptr<Node>& source,
            const shared_ptr<Node>& target);

    /**
     * create an Euler Tour from the graph, provided that it is symmetric (all nodes possess same in- and outdegree),
     * otherwise an empty Euler Tour is delivered
     * @return euler tour for the graph, if it is symmetric, or an empty tour otherwise
     *
     * @note the algorithm of *Hierholz* is used for the extraction of the Euler Tour as described in
     *  > Srivastava, Sudhir:
     *   > Study of Different algorithm in Euler Graph.
     *   > In: International Journal on Recent and Innovation Trends in
     *   > Computing and Communication 4 (2016), Dez, Nr. 12, 309–311.
     *   > http://www.ijritcc.org/download/browse/Volume_4_Issues/December_16_Volume_4_Issue_12/1501310493_29-07-2017.pdf. – ISSN 2321–8169
     */
    shared_ptr<list<shared_ptr<Edge>>> generateEulerTour();

    /**
    * Prints the graph in .dot format into the given filename (in the current path)
    * @param fname the filename (without extension) to print the graph into
    */
    void toDot(const string& fname);

    friend ostream & operator<<(ostream & out, const Graph & graph);
};

#endif //FSM_GRAPHS_GRAPH_H_
