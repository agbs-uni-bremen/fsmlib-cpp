/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/TreeNode.h"

TreeNode::TreeNode()
	: parent(std::weak_ptr<TreeNode>()), children(std::make_shared<std::vector<std::shared_ptr<TreeEdge>>>()), deleted(false)
{

}

TreeNode::TreeNode(const TreeNode* other):
    parent(std::weak_ptr<TreeNode>()), children(std::make_shared<std::vector<std::shared_ptr<TreeEdge>>>())
{
    for (std::shared_ptr<TreeEdge> child: *other->children)
    {
        std::shared_ptr<TreeEdge> childCopy = child->Clone();
        children->push_back(childCopy);
    }
    deleted = other->deleted;
}

void TreeNode::setParent(const std::weak_ptr<TreeNode> pparent)
{
	parent = pparent;
}

std::weak_ptr<TreeNode> TreeNode::getParent() const
{
	return parent;
}

void TreeNode::deleteNode()
{
	deleted = true;

	if (!isLeaf())
	{
		return;
	}

	std::shared_ptr<TreeNode> c = shared_from_this();
	std::shared_ptr<TreeNode> t = parent.lock();

	while (t != nullptr)
	{
		t->remove(c);
		if (!t->isLeaf())
		{
			break;
		}
		if (!t->isDeleted())
		{
			break;
		}
		c = t;
		t = t->getParent().lock();
	}
}

bool TreeNode::isDeleted() const
{
	return deleted;
}

std::shared_ptr<std::vector<std::shared_ptr<TreeEdge>>> TreeNode::getChildren() const
{
	return children;
}

void TreeNode::remove(const std::shared_ptr<TreeNode> node)
{
	for (std::shared_ptr<TreeEdge> e : *children)
	{
		if (e->getTarget() == node)
		{
			children->erase(std::find(children->begin(), children->end(), e));
			break;
		}
	}
}

void TreeNode::calcLeaves(std::vector<std::shared_ptr<TreeNode>>& leaves)
{
	if (isLeaf())
	{
		leaves.push_back(shared_from_this());
	}
	else
	{
		for (std::shared_ptr<TreeEdge> e : *children)
		{
			e->getTarget()->calcLeaves(leaves);
		}
	}
}

void TreeNode::add(const std::shared_ptr<TreeEdge> edge)
{
	edge->getTarget()->setParent(shared_from_this());
	children->push_back(edge);
}

bool TreeNode::isLeaf() const
{
	return children->empty();
}

int TreeNode::getIO(const std::shared_ptr<TreeNode> node) const
{
	for (std::shared_ptr<TreeEdge> e : *children)
	{
		if (e->getTarget() == node)
		{
			return e->getIO();
		}
	}
	exit(EXIT_FAILURE);
}

std::shared_ptr<TreeEdge> TreeNode::hasEdge(const std::shared_ptr<TreeEdge> edge) const
{
	for (std::shared_ptr<TreeEdge> g : *children)
	{
		if (g->getIO() == edge->getIO())
		{
			return g;
		}
	}
	return nullptr;
}

std::vector<int> TreeNode::getPath()
{
	std::vector<int> path;

	std::shared_ptr<TreeNode> m = shared_from_this();
	std::shared_ptr<TreeNode> n = parent.lock();

	while (n != nullptr)
	{
		path.insert(path.begin(), n->getIO(m));
		m = n;
		n = n->getParent().lock();
	}

	return path;
}

bool TreeNode::superTreeOf(const std::shared_ptr<TreeNode> otherNode) const
{
	if (children->size() < otherNode->children->size())
	{
		return false;
	}

	for (std::shared_ptr<TreeEdge> eOther : *otherNode->children)
	{
		int y = eOther->getIO();
		bool yFound = false;

		for (std::shared_ptr<TreeEdge> eMine : *children)
		{
			if (y == eMine->getIO())
			{
				if (!eMine->getTarget()->superTreeOf(eOther->getTarget()))
				{
					return false;
				}
				yFound = true;
				break;
			}
		}

		/*If this node does not have an outgoing edge labelled with y, the nodes differ.*/
		if (!yFound)
		{
			return false;
		}
	}
	return true;
}

bool operator==(TreeNode const & treeNode1, TreeNode const & treeNode2)
{
	if (treeNode1.children->size() != treeNode2.children->size())
	{
		return false;
	}

	if (treeNode1.deleted != treeNode2.deleted)
	{
		return false;
	}

	/*Now compare the child nodes linked by edges with the same output label.
	Since we are only dealing with observable FSMs, the output label
	uniquely determines the edge and the target node: all outputs
	have been generated from the SAME input.*/
	for (std::shared_ptr<TreeEdge> e : *treeNode1.children)
	{
		int y = e->getIO();
		bool yFound = false;

		for (std::shared_ptr<TreeEdge> eOther : *treeNode2.children)
		{
			if (y == eOther->getIO())
			{
				if (!(*e->getTarget() == *eOther->getTarget()))
				{
					return false;
				}
				yFound = true;
				break;
			}
		}

		/*If otherNode does not have an outgoing edge labelled with y, the nodes differ.*/
		if (!yFound)
		{
			return false;
		}

	}
	return true;
}

bool operator!=(TreeNode const & treeNode1, TreeNode const & treeNode2)
{
    return !(treeNode1 == treeNode2);
}

std::shared_ptr<TreeNode> TreeNode::add(const int x)
{
	for (std::shared_ptr<TreeEdge> e : *getChildren())
	{
		if (e->getIO() == x)
		{
			return e->getTarget();
		}
	}

	std::shared_ptr<TreeNode> tgt = std::make_shared<TreeNode>();
	add(std::make_shared<TreeEdge>(x, tgt));
	return tgt;
}

void TreeNode::add(std::vector<int>::const_iterator lstIte, const std::vector<int>::const_iterator end)
{
	/*There may be no next list element, when this method is called*/
	if (lstIte == end)
	{
		return;
	}

	/*Which input is represented by the list iterator?*/
	int x = *lstIte++;

	for (std::shared_ptr<TreeEdge> e : *getChildren())
	{
		/*Is there already an edge labelled with this input?*/
		if (e->getIO() == x)
		{
			/*We do not need to extend the tree, but follow the existing edge*/
			std::shared_ptr<TreeNode> nTgt = e->getTarget();
			nTgt->add(lstIte, end);
			return;
		}
	}

	/*No edge labelled with x exists for this node.
	Therefore one has to be created*/
	std::shared_ptr<TreeNode> newNode = std::make_shared<TreeNode>();
	newNode->setParent(shared_from_this());
	getChildren()->push_back(std::make_shared<TreeEdge>(x, newNode));
	newNode->add(lstIte, end);
}

void TreeNode::add(const IOListContainer & tcl)
{
	/*First delegate the work to the children*/
	for (std::shared_ptr<TreeEdge> e : *getChildren())
	{
		std::shared_ptr<TreeNode> nTgt = e->getTarget();
		nTgt->add(tcl);
	}

	/*Now append each input sequence in tcl to this node,
	using the special strategy of the add(lstIte) operation*/
	for (std::vector<int>& lst : *tcl.getIOLists())
	{
		add(lst.cbegin(), lst.cend());
	}
}

void TreeNode::addToThisNode(const IOListContainer& tcl)
{
	/*Append each input sequence in tcl to this node,
	using the special strategy of the add(lstIte) operation*/
	for (std::vector<int>& lst : *tcl.getIOLists())
	{
		add(lst.cbegin(), lst.cend());
	}
}

void TreeNode::addToThisNode(const std::vector<int> &lst)
{
    add(lst.cbegin(), lst.cend());
}

std::shared_ptr<TreeNode> TreeNode::after(std::vector<int>::const_iterator lstIte, const std::vector<int>::const_iterator end)
{
	if (lstIte != end)
	{
		int x = *lstIte ++;

		for (std::shared_ptr<TreeEdge> e : *getChildren())
		{
			if (e->getIO() == x)
			{
				return e->getTarget()->after(lstIte, end);
			}
		}

		/*Could not find an edge labelled by x*/
		return nullptr;
	}

	/*This is the last node reached by the input trace
	represented by iterator trIte. We have processed the
	whole input trace and can therefore return this
	as the final node.*/
	return shared_from_this();
}

void TreeNode::calcSize(size_t& theSize) {
    
    theSize++;
    
    for ( auto t : *children) {
        t->getTarget()->calcSize(theSize);
    }
    
}

TreeNode* TreeNode::clone() const
{
    return new TreeNode( this );
}

std::shared_ptr<TreeNode> TreeNode::Clone() const
{
    std::shared_ptr<TreeNode> copy = std::shared_ptr<TreeNode>(clone());
    for (auto child: *(copy->children))
    {
        child->getTarget()->setParent(copy);
    }
    return copy;
}












