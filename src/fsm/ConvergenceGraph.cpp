#include "fsm/ConvergenceGraph.h"
#include "fsm/Dfsm.h"
#include "fsm/FsmNode.h"
#include "trees/Tree.h"
#include "trees/TreeNode.h"

#include <memory>
#include <vector>


ConvergenceGraph::ConvergenceGraph(const Dfsm& dfsm, const std::shared_ptr<Tree> testSuite, const std::shared_ptr<ConvergenceNode> root)
 : dfsm(dfsm), testSuite(testSuite), root(root) 
{

}

std::shared_ptr<ConvergenceNode> ConvergenceGraph::after(const std::vector<int>& trace) {
    std::shared_ptr<ConvergenceNode> node = root;
    for (int x : trace) {
        node = node->nextForInput[x];
    }
    return node;
}

void ConvergenceGraph::replaceAfter(const std::vector<int>& trace, int input, std::shared_ptr<ConvergenceNode> newNode) {
    std::shared_ptr<ConvergenceNode> node = after(trace);
    node->nextForInput[input] = newNode;
}


void ConvergenceGraph::add(const std::vector<int>& trace) 
{
    std::shared_ptr<ConvergenceNode> node = root;
    std::shared_ptr<FsmNode> state = dfsm.getInitialState();
    std::shared_ptr<TreeNode> treeNode = testSuite->getRoot();
    for (int x : trace) {

        // the dfsm is assumed to be complete and any trace
        // is assumed to contain only inputs accepted by the dfsm
        // - thus exactly one state is reached via input x
        state = *state->after({x}).begin();

        // it is also assumed that trace is contained in the test suite
        treeNode = treeNode->after(x);

        if (!node->nextForInput[x]) {
            node->nextForInput[x] = std::make_shared<ConvergenceNode>(state,treeNode,dfsm.getMaxInput()+1);
        } else {
            // if a convergent node already exists, then the current prefix of
            // the trace converges with that node and should be added 
            node->nextForInput[x]->convergentNodes.insert(treeNode);
        }

        node = node->nextForInput[x];
    }
}


void ConvergenceGraph::mergeInto(const std::shared_ptr<ConvergenceNode> node1, const std::shared_ptr<ConvergenceNode> node2) {
    for (auto conv : node1->convergentNodes) {
        node2->convergentNodes.insert(conv);
    }
    
    for (int x = 0; x <= dfsm.getMaxInput(); ++x) {
        mergeInto(node1->nextForInput[x], node2->nextForInput[x]);
    }
}


void ConvergenceGraph::merge(const std::vector<int>& trace1, int input, const std::vector<int>& trace2) 
{
    std::shared_ptr<ConvergenceNode> node  = after(trace1);
    std::shared_ptr<ConvergenceNode> nodeL = node->nextForInput[input];
    std::shared_ptr<ConvergenceNode> nodeR = after(trace2);

    node->nextForInput[input] = nodeR;    
    mergeInto(nodeL,nodeR);
}


std::unordered_set<std::shared_ptr<TreeNode>> ConvergenceGraph::getConvergentTraces(const std::vector<int>& trace) 
{
    return after(trace)->convergentNodes;
}

bool ConvergenceGraph::hasLeaf(const std::vector<int>& trace) 
{
    for (auto treeNode : getConvergentTraces(trace)) {
        if (treeNode->isLeaf())
            return true;
    }
    return false;
}
