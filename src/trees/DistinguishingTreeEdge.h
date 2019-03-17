/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FSM_TREES_DISTINGUISHINGTREEEDGE_H_
#define FSM_TREES_DISTINGUISHINGTREEEDGE_H_

#include <vector>
#include <memory>
#include "trees/DistinguishingTreeNode.h"

using namespace std;

class DistinguishingTreeNode;

class DistinguishingTreeEdge
{
protected:
    /**
    The input of this edge
    */
    int input;

    /**
    The target of this edge
    */
    shared_ptr<DistinguishingTreeNode> target;

public:
    /**
    * Create a new distinguishing tree edge
    * @param input the
    */
    DistinguishingTreeEdge(const int input,const shared_ptr<DistinguishingTreeNode>& target);

    /**
    * Gets the target node of this edge
    * @return the target node if this edge
    */
    shared_ptr<DistinguishingTreeNode> getTarget();

};

#endif //FSM_TREES_DISTINGUISHINGTREEEDGE_H_
