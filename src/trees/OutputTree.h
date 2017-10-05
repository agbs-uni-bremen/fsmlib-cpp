/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_OUTPUTTREE_H_
#define FSM_TREES_OUTPUTTREE_H_

#include <vector>

#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/Tree.h"

class OutputTree : public Tree
{
private:
	/**
	The inputTrace of this output tree (one input trace, maybe more than one output possible)
	*/
	InputTrace inputTrace;

	/**
	Print this one and every child (to a dot format)
	@param out The standard output stream to use
	@param top The current "root"
	@param idNode The current id of the node, incremented for EVERY node, used to differenciate them in the dot file
	@param idInput The current id of the input trace, incremented each time you go deeper in the tree
	@return The standard output stream used, to allow user to cascade <<
	*/
	void printChildrenOutput(std::ostream & out, const std::shared_ptr<TreeNode> top, const std::shared_ptr<int> idNode, const int idInput) const;//TODO NOT PRESENT IN JAVA
protected:
    OutputTree(const OutputTree* other);
public:
	/**
	Create a new OutputTree
	@param root The root of the output tree
	@param inputTrace The inputTrace of this output tree (one input trace, maybe more than one output possible)
	@param presentationLayer The presentation layer to use
	*/
	OutputTree(const std::shared_ptr<TreeNode> root, const InputTrace & inputTrace, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Check whether this OutputTree instance is a superset of the output traces
	contained the other instance ot.
	@param ot The other OutputTree instance
	@return false if the trees have been produced by different input traces,
	false if they are associated with the same input traces,
	but ot contains an output trace that does not exist in this OutputTree instance,
	true otherwise.
	*/
	bool contains(const OutputTree & ot) const;

    std::vector<IOTrace> getOutputsIntersection(OutputTree & ot);

	/**
	Store the OutputTree to a standard output file in dot format
	@param out The standard output file to use
	*/
	void toDot(std::ostream & out) const;//TODO NOT PRESENT IN JAVA

	/**
	Store the OutputTree to a standard output file
	@param out The standard output file to use
	*/
	void store(std::ofstream & file);//TODO NOT PRESENT IN JAVA

    
    void toIOTrace(std::vector<IOTrace>& iotrVec);

    virtual OutputTree* _clone() const;
    std::shared_ptr<OutputTree> Clone() const;

	/**
	Output the OutputTree to a standard output stream
	@param out The standard output stream to use
	@param ot The OutputTree to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, OutputTree & ot);

	/**
	Check this OutputTree instance and the instance ot for equality
	@param outputTree1 The first OutputTree instance
	@param outputTree2 The other OutputTree instance
	@return		false if the trees have been produced by different input traces,
	false if they are associated with the same input traces,
	but their is no one-to-one correspondence between the edges
	and their labels,
	true otherwise.
	*/
	friend bool operator==(OutputTree const & outputTree1, OutputTree const & outputTree2);
    
    /** complementary operator to == */
    friend bool operator!=(OutputTree const & outputTree1, OutputTree const & outputTree2);
};
#endif //FSM_TREES_OUTPUTTREE_H_
