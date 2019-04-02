/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "trees/SplittingTree.h"
#include "graphs/PartitionGraph.h"
#include "trees/SplittingTree.h"
#include "trees/InputOutputTree.h"
#include "SplittingTree.h"


#include <queue>
#include <cassert>
#include <utility>

SplittingTree::SplittingTree(const shared_ptr<Dfsm> &dfsm)
    :root(make_shared<SplittingTreeNode>(0)),adaptiveDistinguishingSequence(shared_ptr<InputOutputTree>()),
    dfsm(dfsm),idToFsmNode(vector<shared_ptr<FsmNode>>(dfsm->getMaxNodes(),nullptr))
{

}

bool SplittingTree::checkAValid(shared_ptr<SplittingTreeNode> &blockNode, const int x,
                                map<int, shared_ptr<set<int>>> &blockPartition,
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

shared_ptr<unordered_map<int, int>> SplittingTree::composeBlockToTarget(shared_ptr<unordered_map<int, int>> &firstMap,
                                                                        shared_ptr<unordered_map<int, int>> &secondMap)
{
    auto nextBlockToTarget = make_shared<unordered_map< int, int>>();
    for(auto& idToTarget:*firstMap) {
        auto it = secondMap->find(idToTarget.second);
        assert(it != secondMap->end());
        nextBlockToTarget->insert({idToTarget.first,it->second});
    }
    return nextBlockToTarget;
}

shared_ptr<SplittingTreeNode> SplittingTree::findDeepestNode(shared_ptr<SplittingTreeNode> &rootNode,
                                                              shared_ptr<AdaptiveTreeNode> &adaptiveTreeNode)
{
    shared_ptr<SplittingTreeNode> currentNode = rootNode;
    while(true) {

        bool foundDeeperNode = false;
        for(auto& edge:currentNode->getChildren()) {
            auto targetNode = edge->getTarget();
            auto& targetNodeBlock = targetNode->getBlockToTarget();
            bool containsAllTarget = true;
            for(auto& idToTarget:*adaptiveTreeNode->getInitialToCurrentSet()) {
                auto it = targetNodeBlock->find(idToTarget.second);
                containsAllTarget &= it != targetNodeBlock->end();
            }
            if(containsAllTarget) {
                currentNode = targetNode;
                foundDeeperNode = true;
                break;
            }
        }
        if(foundDeeperNode) continue;

        break;
    }

    return currentNode;
}

shared_ptr<InputOutputTree> &SplittingTree::getAdaptiveDistinguishingSequence() {
    return adaptiveDistinguishingSequence;
}

void SplittingTree::build() {

    //A splitting tree can only be derived for truly completely specified, deterministic FSM
    if(!dfsm->isCompletelyDefined() || !dfsm->isDeterministic()) {
        return;
    }

    //Initialize the root node of the splitting tree and the adaptive distinguishing sequence
    auto adaptiveRootNode = make_shared<AdaptiveTreeNode>();
    auto& rootBlock = root->getBlock();
    auto& rootInitialToCurrentSet = adaptiveRootNode->getInitialToCurrentSet();
    rootInitialToCurrentSet = make_shared<unordered_map< int, int>>();

    for(auto &node:dfsm->getNodes()) {
        rootBlock.insert(node->getId());
        rootInitialToCurrentSet->insert({node->getId(),node->getId()});
        idToFsmNode[node->getId()] = node;
    }

    vector<shared_ptr<SplittingTreeNode>> leaves;
    vector<shared_ptr<SplittingTreeNode>> largestLeaves;
    leaves.push_back(root);
    largestLeaves.push_back(root);

    int nodeIncrementId = 1;

    unsigned long maxBlockSize = rootBlock.size();
    //create the splitting tree
    while(!largestLeaves.empty()) {
        vector<shared_ptr<SplittingTreeNode>> newLeaves;
        vector<shared_ptr<SplittingTreeNode>> cValidNodes;
        auto implicationGraph = make_shared<PartitionGraph>();
        unsigned long maxLeafBlockSize = 0;

        for(auto& currentLeaf:largestLeaves) {
            map< int, shared_ptr<unordered_map< int, int>> > validInputs;
            bool& isAValid = currentLeaf->getIsAValid();
            bool& isBValid = currentLeaf->getIsBValid();
            //check if there is an a-valid input for currentLeaf
            for(int x=0; x<=dfsm->getMaxInput(); ++x) {
                map< int, shared_ptr<set<int>> > blockPartition;
                shared_ptr<unordered_map< int, int>> blockToTarget = make_shared<unordered_map<int,int>>();
                bool isValid = false;
                isAValid = checkAValid(currentLeaf,x,blockPartition,blockToTarget,isValid);
                if(isAValid) {
                    currentLeaf->setTrace(vector<int> {x});
                    currentLeaf->setBlockToTarget(blockToTarget);
                    for(auto& partitionBlock:blockPartition) {
                        auto newNode = make_shared<SplittingTreeNode>(nodeIncrementId++,*partitionBlock.second);
                        auto newEdge = make_shared<SplittingTreeEdge>(partitionBlock.first,newNode);
                        currentLeaf->add(newEdge);

                        newLeaves.push_back(newNode);
                        maxLeafBlockSize = partitionBlock.second->size() > maxLeafBlockSize?partitionBlock.second->size():maxLeafBlockSize;
                    }
                    break;
                } else if(isValid) {
                    validInputs.insert({x,blockToTarget});
                }
            }
            if(isAValid) continue;
            //if there are no valid inputs there is no ads at this point
            if(validInputs.empty()) return;

            //check if there is a b-valid input for currentLeaf with respect to the partition represented by `leaves`
            pair<int, shared_ptr<unordered_map<int,int>>> bValidInput;
            shared_ptr<SplittingTreeNode> bValidTargetNode;
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
                            bValidTargetNode = leaf;
                            break;
                        }
                    }
                    if(isBValid) break;
                    if(containsTarget) {
                        implicationGraph->addEdge(currentLeaf,leaf,validInput.first,validInput.second);
                        break;
                    }
                }
                if(isBValid) break;
            }
            if(isBValid) {
                //Traverse up the tree until u get a node, that includes all target states of `bValidInput`
                do {
                    assert(!bValidTargetNode->getParent().expired());
                    bValidTargetNode = bValidTargetNode->getParent().lock();

                    bool missedOne = false;
                    auto block = bValidTargetNode->getBlock();
                    for(auto& idToTarget:*bValidInput.second) {
                        auto it = block.find(idToTarget.second);
                        missedOne |= it == block.end();
                        if(missedOne) break;
                    }
                    if(missedOne) continue;
                    break;
                } while(true);

                //Set the trace of the current leaf to x.(bValidTargetNode->trace)
                auto& inputTrace = currentLeaf->getTrace();
                inputTrace.push_back(bValidInput.first);
                for(int input:bValidTargetNode->getTrace()) {
                    inputTrace.push_back(input);
                }
                auto& currentBlockToTarget = currentLeaf->getBlockToTarget();
                currentBlockToTarget = composeBlockToTarget(bValidInput.second,bValidTargetNode->getBlockToTarget());

                for(auto& child:bValidTargetNode->getChildren()) {
                    auto newNode = make_shared<SplittingTreeNode>(nodeIncrementId++);
                    auto& newBlock = newNode->getBlock();
                    auto& childBlock = child->getTarget()->getBlock();
                    for(auto& idToTarget:*bValidInput.second) {
                        auto it = childBlock.find(idToTarget.second);
                        if(it != childBlock.end())
                            newBlock.insert(idToTarget.first);
                    }
                    auto newEdge = make_shared<SplittingTreeEdge>(child->getOutput(),newNode);
                    currentLeaf->add(newEdge);

                    newLeaves.push_back(newNode);
                    maxLeafBlockSize = newBlock.size() > maxLeafBlockSize?newBlock.size():maxLeafBlockSize;
                }

            } else {
                cValidNodes.push_back(currentLeaf);
            }
        }
        for(auto& currentLeaf:cValidNodes) {
            auto& cTrace = currentLeaf->getTrace();
            auto& currentBlockToTarget = currentLeaf->getBlockToTarget();
            auto cValidTargetNode = implicationGraph->findPathToAOrBValidNode(currentLeaf,cTrace,currentBlockToTarget);
            if(cValidTargetNode) {
                for(int input:cValidTargetNode->getTrace()) {
                    cTrace.push_back(input);
                }
                for(auto& child:cValidTargetNode->getChildren()) {
                    auto newNode = make_shared<SplittingTreeNode>(nodeIncrementId++);
                    auto& newBlock = newNode->getBlock();
                    auto& childBlock = child->getTarget()->getBlock();
                    for(auto& idToTarget:*currentBlockToTarget) {
                        auto it = childBlock.find(idToTarget.second);
                        if(it != childBlock.end())
                            newBlock.insert(idToTarget.first);
                    }
                    auto newEdge = make_shared<SplittingTreeEdge>(child->getOutput(),newNode);
                    currentLeaf->add(newEdge);

                    newLeaves.push_back(newNode);
                    maxLeafBlockSize = newBlock.size() > maxLeafBlockSize?newBlock.size():maxLeafBlockSize;
                }
                currentBlockToTarget = composeBlockToTarget(currentBlockToTarget,cValidTargetNode->getBlockToTarget());
            } else {
                //There exists no ads
                return;
            }
        }
        for(auto& leaf:leaves) {
            unsigned long blockSize = leaf->getBlock().size();
            if(blockSize < maxBlockSize) {
                newLeaves.push_back(leaf);
                maxLeafBlockSize = blockSize > maxLeafBlockSize?blockSize:maxLeafBlockSize;
            }
        }
        if(maxLeafBlockSize == 1) break;
        maxBlockSize = maxLeafBlockSize;
        largestLeaves.clear();
        leaves.clear();
        for(auto& newLeaf:newLeaves) {
            if(newLeaf->getBlock().size() == maxBlockSize) {
                largestLeaves.push_back(newLeaf);
            }
            leaves.push_back(newLeaf);
        }
    }

    queue<shared_ptr<AdaptiveTreeNode>> indistinctLeaves;
    indistinctLeaves.push(adaptiveRootNode);

    //generate the adaptive distinguishing sequence from the splitting tree
    while(!indistinctLeaves.empty()) {
        auto currentLeaf = indistinctLeaves.front();
        indistinctLeaves.pop();

        //find the deepest node in the splitting tree that contains every state of the current set of `currentLeaf`
        auto v = findDeepestNode(root,currentLeaf);
        auto& trace = v->getTrace();
        //start the single path chain
        shared_ptr<AdaptiveTreeNode> u = currentLeaf;
        //representive state of the current set of `u`, get the first one
        int q = u->getInitialToCurrentSet()->begin()->second;
        for(int i=0;i<trace.size()-1;++i) {
            u->setInput(trace[i]);
            vector<int> outputs;
            auto targetStates = idToFsmNode[q]->after(trace[i],outputs);
            q = targetStates.front()->getId();
            int output = outputs.front();
            auto w = make_shared<AdaptiveTreeNode>();
            auto edge = make_shared<TreeEdge>(output,w);
            u->add(edge);
            u = w;
        }
        u->setInput(trace[trace.size()-1]);
        for(auto& edge:v->getChildren()) {
            auto childNode = edge->getTarget();
            auto& childNodeBlock = childNode->getBlock();
            auto& vGetBlockToTarget = v->getBlockToTarget();

            auto w = make_shared<AdaptiveTreeNode>();
            auto& initialToCurrentSet = w->getInitialToCurrentSet();
            initialToCurrentSet = make_shared<unordered_map< int, int>>();

            bool containsTarget = false;
            for(auto& idToTarget:*currentLeaf->getInitialToCurrentSet()) {
                auto it = childNodeBlock.find(idToTarget.second);
                if(it != childNodeBlock.end()) {
                    containsTarget = true;
                    auto vit = vGetBlockToTarget->find(idToTarget.second);
                    initialToCurrentSet->insert({idToTarget.first,vit->second});
                }
            }
            if(!containsTarget)
                continue;

            auto adaptiveEdge = make_shared<TreeEdge>(edge->getOutput(),w);
            u->add(adaptiveEdge);
            if(initialToCurrentSet->size() > 1)
                indistinctLeaves.push(w);
        }
    }

    adaptiveDistinguishingSequence = make_shared<InputOutputTree>(adaptiveRootNode,dfsm->getPresentationLayer());
}

ostream & operator<<(ostream & out, const SplittingTree & splittingTree) {
    out << "digraph g {" << endl << endl;

    out << "node [shape = ellipse]" << endl;

    queue<shared_ptr<SplittingTreeNode>> workingQueue;
    workingQueue.push(splittingTree.getRoot());

    while(!workingQueue.empty()) {
        auto node = workingQueue.front();
        workingQueue.pop();

        string traceString = "";
        for(int input:node->getTrace()) {
            traceString += to_string(input) + ".";
        }
        traceString = traceString.substr(0,traceString.length()-1);

        string blockString = "(";
        for(int state:node->getBlock()) {
            blockString += to_string(state) + ",";
        }
        if(blockString.size() > 1) {
            blockString = blockString.substr(0,blockString.length()-1);
        }
        blockString += ")";

        out << node->getId() << "[label=\"" << blockString << "\\n" << traceString << "\"];" <<  endl;

        for(auto& edge:node->getChildren()) {
            auto targetNode = edge->getTarget();
            workingQueue.push(targetNode);
        }
    }

    workingQueue.push(splittingTree.getRoot());

    while(!workingQueue.empty()) {
        auto node = workingQueue.front();
        workingQueue.pop();

        for(auto& edge:node->getChildren()) {
            auto targetNode = edge->getTarget();
            out << node->getId() << " -> " << targetNode->getId() << "[label=\"" << edge->getOutput() << "\"];"
                << "  //" << node->getId() << " -> " << targetNode->getId() << endl;
            workingQueue.push(targetNode);
        }
    }

    out << endl << "}" << endl;
    return out;
}

const shared_ptr<SplittingTreeNode>& SplittingTree::getRoot() const {
    return root;
}

void SplittingTree::toDot(const string &fname) {
    ofstream out(fname + ".dot");
    out << *this;
    out.close();
}

