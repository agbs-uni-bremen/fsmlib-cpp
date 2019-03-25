/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "trees/SplittingTree.h"
#include <queue>

SplittingTree::SplittingTree(const shared_ptr<Dfsm> &dfsm)
    :adaptiveDistinguishingSequence(shared_ptr<InputOutputTree>())
{

    root = make_shared<SplittingTreeNode>();

    //A distinguishing tree can only be derived for truly completely specified, deterministic FSM
    if(!dfsm->isCompletelyDefined() || !dfsm->isDeterministic()) {
        return;
    }

    //Mapping from id to FsmNode to speed up the access
    idToFsmNode = vector<shared_ptr<FsmNode>>(dfsm->getMaxNodes(),nullptr);

    //Initialize the root node
    auto rootBlock = root->getBlock();

    for(auto &node:dfsm->getNodes()) {
        rootBlock.insert(node->getId());
        idToFsmNode[node->getId()] = node;
    }

    vector<shared_ptr<SplittingTreeNode>> leaves;
    vector<shared_ptr<SplittingTreeNode>> incompleteLeaves;
    leaves.push_back(root);
    incompleteLeaves.push_back(root);

    bool isDiscretePartition = false;

    while(!isDiscretePartition) {

        for(auto& currentLeaf:incompleteLeaves) {
            for(int x=0; x<=dfsm->getMaxInput(); ++x) {


            }
        }

    }



}

bool SplittingTree::isAValid(shared_ptr<SplittingTreeNode> blockNode, const int x,
                             vector<set<int>> &partitionedBlock)
{
    unordered_map< int, shared_ptr<set<int>> > blockPartition;
    for(auto id: blockNode->getBlock()) {
        auto producedOutputs = vector<int>();
        /**
         * Apply 'x' to the node with id 'id'. note that there must exactly one output and one
         * target state produced, since 'dfsm' is completely specified and deterministic
         **/
        auto targetStates = idToFsmNode[id]->after(x,producedOutputs);
        int output = producedOutputs.front();
        auto bit = blockPartition.find(output);
    }

}