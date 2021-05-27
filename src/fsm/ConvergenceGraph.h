/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_CONVERGENCEGRAPH_H_
#define FSM_FSM_CONVERGENCEGRAPH_H_

#include <memory>
#include <vector>
#include <unordered_set>

#include <trees/TreeNode.h>
#include <fsm/FsmNode.h>
#include <fsm/Dfsm.h>


// Type of nodes in the convergence graph
struct ConvergenceNode {
    // traces converging in this node
    std::unordered_set<std::shared_ptr<TreeNode>> convergentNodes;
    // nextForInput[x] is the convergence node reached by input x 
    std::vector<std::shared_ptr<ConvergenceNode>> nextForInput;
    // state of the underlying DFSM corresponding to this node
    std::shared_ptr<FsmNode> state;

    bool isLeaf = true;

    ConvergenceNode(const std::shared_ptr<FsmNode>& node,const std::shared_ptr<TreeNode> convergentNode, size_t numInputs) : state(node) {
        convergentNodes.insert(convergentNode);
        nextForInput.resize(numInputs);
    }
};

class ConvergenceGraph
{
protected:
    const Dfsm& dfsm;
    const std::shared_ptr<Tree> testSuite;
    std::shared_ptr<ConvergenceNode> root;

    // assumes that trace is aready defined in the graph
    std::shared_ptr<ConvergenceNode> after(const InputTrace& trace);

    // assumes that trace followed by input is aready defined in the graph
    void replaceAfter(const InputTrace& trace, int input, std::shared_ptr<ConvergenceNode> newNode);

    void mergeInto(const std::shared_ptr<ConvergenceNode> node1, const std::shared_ptr<ConvergenceNode> node2);

public:
	ConvergenceGraph(const Dfsm& dfsm, const std::shared_ptr<Tree> testSuite, const std::shared_ptr<ConvergenceNode> root);

    void add(const InputTrace& trace);

    void merge(const InputTrace& trace1, int input, const InputTrace& trace2);

    std::unordered_set<std::shared_ptr<TreeNode>> getConvergentTraces(const InputTrace& trace);
	
};
#endif //FSM_FSM_CONVERGENCEGRAPH_H_
