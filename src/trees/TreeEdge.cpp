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

TreeEdge::TreeEdge(std::shared_ptr<TreeEdge const> other)
{
    if (other->target != nullptr)
    {
        target = std::shared_ptr<TreeNode>(other->target->clone());
    }
    io = other->io;

}

int TreeEdge::getIO() const
{
	return io;
}

std::shared_ptr<TreeNode> TreeEdge::getTarget() const
{
	return target;
}

TreeEdge* TreeEdge::clone() const
{
    return new TreeEdge( *this );
}
