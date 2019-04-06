/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_GRAPHS_NETWORKEDGE_H
#define FSM_GRAPHS_NETWORKEDGE_H

#include "graphs/Edge.h"

/**
* A directed edge of a network.
*/
class NetworkEdge: public Edge {
protected:

    /**
     * The capacity of the network edge. A negative value is interpreted as infinity.
     */
    int capacity;

    /**
     * the flow value of the network edge. it should not exceed the capacity and conserve the flow, which means
     * that the net flow of ingoing and outgoing edges of a node (except the source and sink) should be the same
     */
    unsigned int flow;

    /**
     * auxilliary flag that is set if this edge represents an optimized alpha sequence
     */
    bool isAlpha;

    /**
     * auxilliary flag that is set if this edge represents a trace that is a preset distinguishing sequence
     * or an input portion of a path of an adaptive distinguishing sequence
     */
    bool isDs;

    /**
     * flag that is set if this edge represents a reverse edge.
     */
    bool isReverse;

    /**
     * the reverse edge, if one exists. for the reverse edge it is the pointer to the original edge and vice versa
     */
    weak_ptr<NetworkEdge> reverseEdge;

public:

    /**
    * Creates a new network edge
    * @param trace the trace associated with this edge
    * @param source the source node of this edge
    * @param target the target node of this edge
    * @param capacity the capacity of this network edge
    * @param cost the cost of this network edge
    */
    NetworkEdge(const vector<int>& trace, const weak_ptr<Node>& source, const weak_ptr<Node>& target,int capacity,int cost);

    /**
     * gets the #isAlpha flag
     * @return the #isAlpha flag
     */
    bool getIsAlpha();

    /**
     * sets the #isAlpha flag
     * @param isAlpha new value for #isAlpha
     */
    void setIsAlpha(bool isAlpha);

    /**
     * gets the #isDs flag
     * @return the #isDs flag
     */
    bool getIsDs();

    /**
     * sets the #isDs flag
     * @param isDs the new value for #isDs
     */
    void setIsDs(bool isDs);

    int getCapacity() const;

    void setCapacity(int capacity);

    void setFlow(unsigned int flow);

    unsigned int getFlow() const;

    bool getIsReverse() const;

    void setIsReverse(bool isReverse);

    const weak_ptr<NetworkEdge> &getReverseEdge() const;

    void setReverseEdge(const weak_ptr<NetworkEdge> &reverseEdge);
};
#endif // FSM_GRAPHS_NETWORKEDGE_H
