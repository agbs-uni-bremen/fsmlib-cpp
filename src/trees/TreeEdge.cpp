/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/TreeEdge.h"
#include "trees/TreeNode.h"

TreeEdge::TreeEdge(const int io, const std::shared_ptr<TreeNode> target)
	: io(io), target(target)
{

}

int TreeEdge::getIO() const
{
	return io;
}

std::shared_ptr<TreeNode> TreeEdge::getTarget() const
{
	return target;
}
