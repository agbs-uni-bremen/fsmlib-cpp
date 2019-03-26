/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "trees/SplittingTree.h"
#include <queue>
#include <pair>

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
        vector<shared_ptr<SplittingTreeNode>> newLeaves;
        vector<shared_ptr<SplittingTreeNode>> newIncompleteLeaves;
        vector<shared_ptr<SplittingTreeNode>> cValidNodes;

        for(auto& currentLeaf:incompleteLeaves) {
            unordered_map< int, shared_ptr<unordered_map< int, int>> > validInputs;
            bool isAValid = false;
            bool isBValid = false;
            //check if there is an a-valid input for currentLeaf
            for(int x=0; x<=dfsm->getMaxInput(); ++x) {
                unordered_map< int, shared_ptr<set<int>> > blockPartition;
                shared_ptr<unordered_map< int, int>> blockToTarget = make_shared<unordered_map<int,int>>();
                bool isValid = false;
                isAValid = checkAValid(currentLeaf,x,blockPartition,blockToTarget,isValid);
                if(isAValid) {
                    currentLeaf->setTrace(vector<int> {x});
                    currentLeaf->setBlockToTarget(blockToTarget);
                    for(auto& partitionBlock:blockPartition) {
                        auto newNode = make_shared<SplittingTreeNode>(*partitionBlock.second);
                        auto newEdge = make_shared<SplittingTreeEdge>(partitionBlock.first,newNode);
                        currentLeaf->add(newEdge);

                        newLeaves.push_back(newNode);
                        if(partitionBlock.second->size() > 1) {
                            newIncompleteLeaves.push_back(newNode);
                        }
                    }
                    break;
                } else if(isValid) {
                    validInputs.insert({x,blockToTarget});
                }
            }
            if(isAValid) continue;
            //if there are no valid inputs there is no ads at this point
            if(validInputs.size() == 0) return;

            //check if there is a b-valid input for currentLeaf with respect to the partition represented by `leaves`
            pair<int, shared_ptr<unordered_map<int,int>>> bValidInput;
            for(auto& validInput:validInputs) {
                for(auto& leaf:leaves) {
                    auto block = leaf->getBlock();
                    bool containsTarget = false,
                         missedOneTarget = false;
                    for(auto& idToTarget:*validInput.second) {
                        auto it = block.find(idToTarget.second);
                        containsTarget |= it != block.end();
                        missedOneTarget |= it == block.end();
                        if(containsTarget && missedOneTarget) {
                            isBValid = true;
                            bValidInput = validInput;
                            break;
                        }
                    }
                    if(isBValid) break;

                }
                if(isBValid) break;
            }

        }

    }



}

bool SplittingTree::checkAValid(shared_ptr<SplittingTreeNode> &blockNode, const int x,
                                unordered_map<int, shared_ptr<set<int>>> &blockPartition,
                                shared_ptr<unordered_map<int, int>> &blockToTarget,
                                bool &isValid)
{
    isValid = true;
    blockPartition.clear();
    blockToTarget->clear();
    unordered_map< int, shared_ptr<set<int>> > targetPartition;

    for(auto id: blockNode->getBlock()) {
        auto producedOutputs = vector<int>();
        /**
         * Apply 'x' to the node with id 'id'. note that there must exactly one output and one
         * target state produced, since 'dfsm' is completely specified and deterministic
         **/
        auto targetStates = idToFsmNode[id]->after(x,producedOutputs);
        int target = targetStates.front()->getId();
        int output = producedOutputs.front();
        auto tpit = targetPartition.find(output);
        auto bpit = blockPartition.find(output);
        if(tpit == targetPartition.end()) {
            auto targetPartitionBlock = make_shared<set<int>>();
            targetPartition.insert({output,targetPartitionBlock});
            targetPartitionBlock->insert(target);

            auto partitionBlock = make_shared<set<int>>();
            blockPartition.insert({output,partitionBlock});
            partitionBlock->insert(id);

        } else {
            auto it = tpit->second->insert(target);
            //the targetstate already exists inside that partition, which means that 'x' is not valid for 'blockNode'
            if(!it.second) {
                isValid = false;
                return false;
            }
            bpit->second->insert(id);
        }
        blockToTarget->insert({id,target});
    }
    return blockPartition.size() > 1;

}