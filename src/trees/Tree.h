/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_TREE_H_
#define FSM_TREES_TREE_H_

#include <algorithm>
#include <memory>
#include <vector>

#include "fsm/InputTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"
#include "trees/TreeEdge.h"
#include "trees/TreeNode.h"

class Tree
{
protected:
	/**
	The root of this tree
	*/
	std::shared_ptr<TreeNode> root;

	/**
	The list of the leaves of this tree (empty, unless you call calcleaves)
	*/
	std::vector<std::shared_ptr<TreeNode>> leaves;

	/**
	The presentation layer used by this tree
	*/
	const std::shared_ptr<FsmPresentationLayer> presentationLayer;

	/**
	Calculate the leaves of the all tree, calling calcLeaves on the root of the tree
	*/
	void calcLeaves();

	//TODO
	void remove(const std::shared_ptr<TreeNode> thisNode, const std::shared_ptr<TreeNode> otherNode);

	/**
	Print every childran of this tree to a dot format into a standard output stream
	@param out The standard output stream to use
	@param top The root of the tree
	@param idNode The id of this node, used to differenciate node in dot format
	*/
	void printChildren(std::ostream & out, const std::shared_ptr<TreeNode> top, const std::shared_ptr<int> idNode) const;
public:
	/**
	Create a new tree, with a root and a presenation layer
	@param root  root of the tree
	@param presentationLayer The presentation layer to use
	*/
	Tree(const std::shared_ptr<TreeNode> root,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Calculate the leaves, then give the leaves back
	@return The leaves of this tree
	*/
	std::vector<std::shared_ptr<TreeNode>> getLeaves();

	/**
	Getter for the root of this tree
	@return The root of this tree
	*/
	std::shared_ptr<TreeNode> getRoot() const;

	//TODO
	IOListContainer getIOLists();

	/**
	Special remove operation.
	@param otherTree For all edges in otherTree that correspond to
	an edge in this tree, the corresponding source
	node and target node in this tree are marked as deleted.
	*/
	void remove(const std::shared_ptr<Tree> otherTree);

	/**
	Output this tree to a dot format, into a standard output stream
	@param out The standard output stream to use
	*/
	void toDot(std::ostream & out);

	/**
	Get the test cases of this tree
	@return the test cases
	*/
	IOListContainer getTestCases();

	IOListContainer getDeterministicTestCases();

	/**
	Append a list of input traces to EVERY node of the input tree.
	Do not create redundant input sequences that are already contained
	(possibly as a prefix) in the existing tree.
	*/
	void add(const IOListContainer & tcl);

	/**
	 * Insert a list of input traces at the root of the input tree.
	 * Do not create redundant input sequences that are already contained
	 * (possibly as a prefix) in the existing tree.
	 */
	void addToRoot(const IOListContainer & tcl);
    
    /**
     * Add a single input trace represented as 
     * vector of int to the root of the tree
     */
    void addToRoot(const std::vector<int> & lst);

	/**
	Construct the union of this Tree and otherTree by adding
	every maximal input trace of otherTree to this inputTree.
	*/
	void unionTree(const std::shared_ptr<Tree> otherTree);

	//TODO
	void addAfter(const InputTrace & tr, const IOListContainer & cnt);
    
    /** Return number of nodes in the tree */
    size_t size();
    
};
#endif //FSM_TREES_TREE_H_
