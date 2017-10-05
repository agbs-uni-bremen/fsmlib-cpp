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

TreeEdge::TreeEdge(const TreeEdge* other)
{
    if (other->target != nullptr)
    {
        target = other->target->Clone();
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

TreeEdge* TreeEdge::_clone() const
{
    return new TreeEdge( this );
}

std::shared_ptr<TreeEdge> TreeEdge::Clone() const
{
    return std::shared_ptr<TreeEdge>(_clone());
}

bool operator==(TreeEdge const & edge1, TreeEdge const & edge2)
{
    return edge1.io == edge2.io && *edge1.target == *edge2.target;
}

bool operator!=(TreeEdge const & edge1, TreeEdge const & edge2)
{
    return !(edge1 == edge2);
}
