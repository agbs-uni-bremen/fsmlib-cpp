/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "SplittingTree.h"

SplittingTree::SplittingTree(const shared_ptr<Dfsm> &dfsm)
    :adaptiveDistinguishingSequence(shared_ptr<InputOutputTree>())
{

    root = make_shared<SplittingTreeNode>();

    //A distinguishing tree can only be derived for truly completely specified, deterministic FSM
    if(!dfsm->isCompletelyDefined() || !dfsm->isDeterministic()) {
        return;
    }

    //Initialize the root node
    auto rootBlock = root->getBlock();

    for(auto &node:dfsm->getNodes()) {
        rootBlock.insert(node->getId());
    }


}
