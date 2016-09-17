/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_TREEEDGE_H_
#define FSM_TREES_TREEEDGE_H_

#include <memory>

class TreeNode;

class TreeEdge
{
private:
	/**
	The input or output of this tree edge
	*/
	int io;

	/**
	The target of this tree edge
	*/
	std::shared_ptr<TreeNode> target;
public:
	/**
	Create a new tree edge
	\param io The input or output of this tree edge
	\param target The target of this tree edge
	*/
	TreeEdge(const int io, const std::shared_ptr<TreeNode> target);

	/**
	Getter for the input or ouput
	\return The input or output of this tree edge
	*/
	int getIO() const;

	/**
	Getter for the target
	\return The target of this tree edge
	*/
	std::shared_ptr<TreeNode> getTarget() const;
};
#endif //FSM_TREES_TREEEDGE_H_