/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_TREES_SPLITTINGTREEEDGE_H_
#define FSM_TREES_SPLITTINGTREEEDGE_H_

#include <memory>
#include "trees/SplittingTreeNode.h"

using namespace std;

class SplittingTreeNode;

class SplittingTreeEdge
{
protected:
    /**
    The output of this edge
    */
    int output;

    /**
    The target of this edge
    */
    shared_ptr<SplittingTreeNode> target;

public:
    /**
    * Create a new splitting tree edge
    * @param output the output of this edge
    * @param target the target node of this edge
    */
    SplittingTreeEdge(const int output,const shared_ptr<SplittingTreeNode>& target);

    /**
    * Gets the target node of this edge
    * @return the target node if this edge
    */
    shared_ptr<SplittingTreeNode> getTarget();

    /**
    * Gets the output of this edge
    * @return the output of this edge
    */
    int getOutput();

};


#endif //FSM_TREES_SPLITTINGTREEEDGE_H_
