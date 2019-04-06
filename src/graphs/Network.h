/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_GRAPHS_NETWORK_H_
#define FSM_GRAPHS_NETWORK_H_

#include "graphs/Graph.h"

using namespace std;

/**
 * This class represents a network graph datastructure as described in <br>
 * Gibbons, Alan: Algorithmic Graph Theory. Cambridge University Press, 1985 .<br>
 * It is mainly used by the D-Method of Hierons et. al. Dfsm::hieronsDMethodOnMinimisedDfsm(bool)
 */
class Network:public Graph {
protected:

    /**
     * The id of the Node inside the #nodes list, that represents the source node of the network.
     * Note that #nodes indices are supposed to respect the node ids. the source node should not have
     * any ingoing edges.
     */
    int sourceNodeId;

    /**
     * The id of the Node inside the #nodes list, that represents the source node of the network.
     * Note that #nodes indices are supposed to respect the node ids. the sink node should not have
     * any outgoing edges.
     */
    int sinkNodeId;

    /**
     * the residual network of this network, if #this is a residual network itself, the pointer is empty.
     */
    shared_ptr<Network> residualNetwork;



public:

    /**
    * Creates a new network graph
    * @param nodes the nodes of this network
    * @param sourceNodeId the id of the source node of this network
    * @param sinkNodeId the id of the sink node of this network
    *
    * @note it is assumed that the ids of elements of \p nodes respect the boundaries of the \p nodes list size
    *       and are unique according to the requirements in Graph::validateNodeIds(). An explicit check is
    *       not performed in the constructor, and has to happen beforehand. furthermore the index of an element
    *       in \p nodes matches its id
    *
    * @note it is assumed that the elements of \p nodes all do use NetworkEdge as edges. Otherwise the behavior
    *       of the algorithms used on this datastructure are undefined.
    *
    * @note the specified source and sink nodes should balance the network graph, according to
    *       the description in https://www.topcoder.com/community/competitive-programming/tutorials/minimum-cost-flow-part-one-key-concepts/.
    *      the source node should not have any ingoing edges and the sink node no outgoing edges accordingly
    *
    */
    Network(const vector<shared_ptr<Node>>& nodes,int sourceNodeId, int sinkNodeId);

    /**
     * constructor for the creation of an residual network, whose field #residualNetwork is empty. it is called
     * by the public constructor to initialize the #residualNetwork field.
     * @param network the network to create the residual network for
     *
     * @note it is assumed that the ids of elements of \p nodes respect the boundaries of the \p nodes list size
     *       and are unique according to the requirements in Graph::validateNodeIds(). An explicit check is
     *       not performed in the constructor, and has to happen beforehand. furthermore the index of an element
     *       in \p nodes matches its id
     *
     * @note it is assumed that the elements of \p nodes all do use NetworkEdge as edges. Otherwise the behavior
     *       of the algorithms used on this datastructure are undefined.
     *
     * @note the specified source and sink nodes should balance the network graph, according to
     *       the description in https://www.topcoder.com/community/competitive-programming/tutorials/minimum-cost-flow-part-one-key-concepts/.
     *      the source node should not have any ingoing edges and the sink node no outgoing edges accordingly
     *
     */
    Network(Network* network);

    /**
     * Calculates a maximal flow with minimal costs(length of the edge traces), for a network, if the following preconditions are met:
     * 1. the network does not contain negative cost cycles
     * 2. the capacities of the outgoing source edges and ingoing sink edges are maxed out by the resulting flow
     *
     * otherwise the resulting flow might not be minimal in the costs
     *
     * The algorithm used is the *Successive Shortest Path* Algorithm by *Edmonds and Karp* described in
     *
     *   > Edmonds, Jack ; Karp, Richard M.:
     *   > Theoretical Improvements in Al-
     *   > gorithmic Efficiency for Network Flow Problems.
     *   > In: J. ACM 19 (1972),
     *   > April, Nr. 2, 248–264. http://dx.doi.org/10.1145/321694.321699. – DOI
     *   > 10.1145/321694.321699. – ISSN 0004–5411
     *
     * A description of the algorithm from a practical viewpoint can be found at the website
     * https://www.topcoder.com/community/competitive-programming/tutorials/minimum-cost-flow-part-two-algorithms/
     */
    void calculateMinimumCostMaximumFlow();

    /**
   * Prints the network in .dot format into the given filename (in the current path)
   * @param fname the filename (without extension) to print the graph into
   */
    void toDot(const string& fname);

    friend ostream & operator<<(ostream & out, const Network & network);
};


#endif //FSM_GRAPHS_NETWORK_H_
