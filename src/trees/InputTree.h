/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_INPUTTREE_H_
#define FSM_TREES_INPUTTREE_H_

#include <vector>
#include <memory>

#include "trees/Tree.h"

class IOTrace;
class InputTrace;
class FsmPresentationLayer;

class InputTree : public Tree
{
private:
    
protected:
    /**
	 * Create a new InputTree by copying a given tree
	 * @param other The tree to copy
	*/
    InputTree(const InputTree* other);
    
    /**
	 * Print this one and every child (to a dot format)
	 * @param out The standard output stream to use
	 * @param top The current "root"
	 * @param idNode The current id of the node, incremented for EVERY node,
     * used to differentiate them in the dot file
	 * @param idInput The current id of the input trace, incremented each time you go deeper in the tree
	 */
	void printChildrenInput(std::ostream& out,
                             const std::shared_ptr<TreeNode>& top,
                             const std::shared_ptr<int>& idNode,
                             const int idInput) const;			   

public:
    /**
	 * Create a new InputTree rooted at a given node
	 * @param root The root of the output tree
	 * @param presentationLayer The presentation layer to use
	*/
    InputTree(const std::shared_ptr<TreeNode>& root,
               const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	 * Create a new empty InputTree representing the singleton set containing the empty sequence
	 * @param presentationLayer The presentation layer to use
	*/
    InputTree(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    
    /**
     * Returns a vector of all input traces stored in this tree
     * @return The vector of input traces stored in this tree
     */
    std::vector<InputTrace> getInputTraces() const;

	/**
	 * Check whether this InputTree instance is a superset of the output traces
	 * contained in the other instance ot.
	 * @param ot The other InputTree instance
	 * @return false if the trees have been produced by different input traces,
	 * @return false if they are associated with the same input traces,
	 *         but ot contains an output trace that does not exist in this InputTree instance,
	 * @return true otherwise.
     *
     * @note This operations re-calculates the leaves of this and of ot,
     *        so it changes the internal state of both objects.
	*/
    bool contains(const InputTree& ot) const;

	/**
	 * Store the InputTree to a standard output file in dot format
	 * @param out The standard output file to use
	 */
	void toDot(std::ostream& out) const;

	/**
	 * Store the InputTree to a standard output file
	 * @param out The standard output file to use
	 */
	void store(std::ofstream& file);

	/**
	 * Retrieve the subtree of this tree for a given edge from the root.
	 * @param input The input to get the subtree for.
	 * @return The subtree of the root reached by a edge labelled with the input,
	 *         if such a edge exists.
	 * 		   Otherwise the empty tree is returned.
	 */
	std::shared_ptr<InputTree> getSubtreeForInput(int input);

	/**
	 * @return The inputs of all edges from the root.  
	 */
	std::vector<int> getInputsAtRoot() const;

	/**
 	 * Get all maximal sequences applied in this tree after both sequences.
	 */ 
	std::shared_ptr<InputTree> sharedExtensions(const InputTrace& t1, const InputTrace& t2);

	/**
	 * Get the number of maximal sequences contained in this tree.
	 */
	unsigned int getNumberOfSequences();

	/**
	 * Get the combined length of all maximal sequences contained in this tree.
	 */
	unsigned int getTotalLengthOfSequences();
    
    virtual InputTree* _clone() const;
    std::shared_ptr<InputTree> Clone() const;

	/**
	 * Input the InputTree to a standard output stream
	 * @param out The standard output stream to use
	 * @param ot The InputTree to print
	 * @return The standard output stream used, to allow user to cascade <<
	 */
	friend std::ostream& operator<<(std::ostream& out, InputTree& ot);

	/**
	 * Check this InputTree instance and the instance ot for equality
	 * @param inputTree1 The first InputTree instance
	 * @param inputTree2 The other InputTree instance
	 * @return      false if they do not contain the exact same set of 
     *              input sequences
	 * @return      true otherwise.
     *
     * @note Checking for equality (or inequality) has a side effect on the input trees involved:
     *       Their leaves are calculated again.
	 */
	friend bool operator==(InputTree const &inputTree1, InputTree const &inputTree2);
    
    /** complementary operator to == */
    friend bool operator!=(InputTree const &inputTree1, InputTree const &inputTree2);
};

bool operator==(InputTree const &inputTree1, InputTree const &inputTree2);

/** complementary operator to == */
bool operator!=(InputTree const &inputTree1, InputTree const &inputTree2);
#endif //FSM_TREES_INPUTTREE_H_
