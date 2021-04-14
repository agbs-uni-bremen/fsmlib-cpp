/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_TREE_H_
#define FSM_TREES_TREE_H_

#include <memory>
#include <vector>

#include "cloneable/ICloneable.h"

class InputTrace;
class FsmPresentationLayer;
class IOListContainer;
class TreeEdge;
class TreeNode;
class SegmentedTrace;

class Tree: public ICloneable, public std::enable_shared_from_this<Tree>
{
protected:
    Tree(const Tree* other);
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
	 * Calculate the leaves of the tree, calling calcLeaves on the root of the tree
	 */
	void calcLeaves();
	std::vector<std::shared_ptr<TreeNode const>> calcLeaves() const;

	//TODO
    void remove(const std::shared_ptr<TreeNode>& thisNode,
                const std::shared_ptr<TreeNode>& otherNode);

	/**
	 * Print all children of this tree to a dot format into a standard output stream
	 * @param out The standard output stream to use
	 * @param top The root of the tree
	 * @param idNode The id of this node, used to differenciate node in dot format
	*/
    void printChildren(std::ostream & out, const std::shared_ptr<TreeNode>& top, const std::shared_ptr<int>& idNode) const;

    /**
     *  @return true if one of the input traces is prefix of the other one, false otherwise.
     */
    bool inPrefixRelation(std::vector<int> aPath, std::vector<int> bPath);
public:
	/**
	Create a new tree, with a root and a presentation layer
	@param root  root of the tree
	@param presentationLayer The presentation layer to use
	*/
    Tree(const std::shared_ptr<TreeNode>& root,
         const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

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

	/**
     * Get vector of all I/O lists in the tree.
     * Each list is represented as a vector.
     * Lists which are prefixes of other lists are NOT
     * represented as separate vectors.
     *
     * @return vector of I/O lists, each list again represented
     *         as a vector.
     */
	IOListContainer getIOLists() const;
    
    /**
     * Get vector of all I/O lists in the tree, including all prefixes.
     * For example, the returned container also contains the empty list.
     * Each list is represented as a vector. 
     *
     * @return vector of I/O lists, each list again represented
     *         as a vector.
     */
    IOListContainer getIOListsWithPrefixes();

	/**
	Special remove operation.
	@param otherTree For all edges in otherTree that correspond to
	an edge in this tree, the corresponding source
	node and target node in this tree are marked as deleted.
	*/
    void remove(const std::shared_ptr<Tree>& otherTree);

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

    /**
     * Calculates all test cases that are being represented by this tree.
     * Prefixes are not being removed, since this is needed for deterministic
     * state cover. Every d-reachable state has to be d-reached (not just visited)
     * by some input sequence.
     * @return Test cases in form of input sequences.
     */
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
    void unionTree(const std::shared_ptr<Tree>& otherTree);

	//TODO
	void addAfter(const InputTrace & tr, const IOListContainer & cnt);

    /**
     * Determines wether this tree's root node has an outgoing edge with the given io.
     * @param y The given io
     * @return {@code true}, if there is an outgoing edge with the given io,
     * {@code false} otherwise.
     */
    bool isDefined(int y) const;
    
    /** Return number of nodes in the tree */
    size_t size();

    virtual Tree* _clone() const;
    std::shared_ptr<Tree> Clone() const;

    /**
     *  Construct a pseudo-intersection tree of this Tree and b.
     *
     *  Appending one of the pathes of the resulting tree to the roots of both
     *  trees does not result in more breadth (more test cases).
     *  Therefore it is useful to search in the resulting for distinguishing
     *  traces for the two root nodes.
     *
     *  @param b For every path of one of the two trees (this and b) that is
     *           a prefix of a path of the other tree we add the longer path to
     *           the resulting tree.
     *  @return Tree
     */
    std::shared_ptr<Tree> getPrefixRelationTree(const std::shared_ptr<Tree> &b);
    
    /**
     * create a deep copy of a subtree that is reached by alpha
     * @param alpha InputTrace that leads to the root of the new subtree
     * @return Tree subtree with after-alpha as the new root node
     *         or null if no tree node could be found after applying alpha
     */
    std::shared_ptr<Tree> getSubTree(const std::shared_ptr<InputTrace>& alpha);
    
    /**
     *  return the TreeNode where the subtree after input trace alpha starts
     */
    std::shared_ptr<TreeNode> getSubTree(std::shared_ptr< std::vector<int> > alpha);

    /**
     *  return the TreeNode where the subtree after input trace alpha starts
     */
    std::shared_ptr<TreeNode> getSubTree(const std::vector<int>& alpha);

    /**
     * Construct the intersection of this tree with another tree.
     * 
     * The resulting tree contains a sequence xs if and only if both
     * input trees contain sequences that xs is a prefix of, while this
     * property holds for no extension of xs.
     */
    std::shared_ptr<Tree> getIntersectionTree(const std::shared_ptr<Tree> &b);
    
    
    /**
     *   Check the effect of adding a trace at the tree root,
     *   using addToRoot(),
     *   without actually changing the tree.
     *   @param alpha The trace to be checked w.r.t. insertion into the tree
     *
     *   @return 0 If the alpha is already contained in the tree,
     *             so it's not necessary to add it at all
     *   @return 1 If adding alpha will just extend one path of the tree at
     *             one of its leaves.
     *   @return 2 If adding alpha will create a new branch in the tree,
     *             leading to an additional test case.
     */
    int tentativeAddToRoot(const std::vector<int>& alpha);
    int tentativeAddToRoot(SegmentedTrace& alpha) const;

    friend std::ostream & operator<<(std::ostream & out, Tree & t);

    std::string str();
};
#endif //FSM_TREES_TREE_H_
