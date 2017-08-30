/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/Tree.h"

void Tree::calcLeaves()
{
	leaves.clear();
	root->calcLeaves(leaves);
}

void Tree::remove(const std::shared_ptr<TreeNode> thisNode, const std::shared_ptr<TreeNode> otherNode)
{
	thisNode->deleteNode();

	for (std::shared_ptr<TreeEdge> e : *thisNode->getChildren())
	{
		std::shared_ptr<TreeEdge> eOther = otherNode->hasEdge(e);
		if (eOther != nullptr)
		{
			remove(e->getTarget(), eOther->getTarget());
		}
	}
}

void Tree::printChildren(std::ostream & out, const std::shared_ptr<TreeNode> top, const std::shared_ptr<int> idNode) const
{
	int idNodeBase = *idNode;
	for (std::shared_ptr<TreeEdge> edge : *top->getChildren())
	{
		out << idNodeBase << " -> " << ++ *idNode << "[label = \"" << edge->getIO() << "\" ];" << std::endl;
		printChildren(out, edge->getTarget(), idNode);
	}
}

Tree::Tree(const std::shared_ptr<TreeNode> root, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: root(root), presentationLayer(presentationLayer)
{

}

Tree::Tree(const Tree* other): presentationLayer(other->presentationLayer)
{
    if (other->root != nullptr)
    {
        root = other->root->Clone();
    }
    if (other->leaves.size() > 0)
    {
        calcLeaves();
    }
}

std::vector<std::shared_ptr<TreeNode>> Tree::getLeaves()
{
	calcLeaves();
	return leaves;
}

std::shared_ptr<TreeNode> Tree::getRoot() const
{
	return root;
}

IOListContainer Tree::getIOLists()
{
	std::shared_ptr<std::vector<std::vector<int>>> ioll = std::make_shared<std::vector<std::vector<int>>>();
	calcLeaves();

	for (std::shared_ptr<TreeNode> n : leaves)
	{
		ioll->push_back(n->getPath());
	}

	return IOListContainer(ioll, presentationLayer);
}

void Tree::remove(const std::shared_ptr<Tree> otherTree)
{
	remove(getRoot(), otherTree->getRoot());
}

void Tree::toDot(std::ostream & out)
{
	out << "digraph Tree {" << std::endl;
	out << "\trankdir=TB;" << std::endl;//Top -> Bottom, to create a vertical graph
	out << "\tnode [shape = circle];" << std::endl;
	std::shared_ptr<int> id = std::make_shared<int>(0);
	printChildren(out, root, id);
	out << "}";
}

IOListContainer Tree::getTestCases()
{
	return getIOLists();
}

IOListContainer Tree::getDeterministicTestCases()
{
	std::shared_ptr<std::vector<std::vector<int>>> ioll = std::make_shared<std::vector<std::vector<int>>>();
	std::shared_ptr<TreeNode> currentNode = root;
	ioll->push_back({-1});
	std::shared_ptr<std::vector<std::shared_ptr<TreeEdge>>> edges = currentNode->getChildren();

	for (size_t i=0; i < edges->size(); ++i)
	{
		std::shared_ptr<TreeEdge> edge = (*edges)[i];
		std::shared_ptr<TreeNode> node = edge->getTarget();
		ioll->push_back(node->getPath());
		std::shared_ptr<std::vector<std::shared_ptr<TreeEdge>>> newEdges = node->getChildren();
		std::copy (newEdges->begin(), newEdges->end(), std::back_inserter(*edges));
	}
	return IOListContainer(ioll, presentationLayer);
}

void Tree::add(const IOListContainer & tcl)
{
	std::shared_ptr<TreeNode> r = getRoot();
	r->add(tcl);
}

void Tree::addToRoot(const IOListContainer & tcl)
{
	std::shared_ptr<TreeNode> r = getRoot();
	r->addToThisNode(tcl);
}

void Tree::addToRoot(const std::vector<int> &lst)
{
    std::shared_ptr<TreeNode> r = getRoot();
    r->addToThisNode(lst);
}

void Tree::unionTree(const std::shared_ptr<Tree> otherTree)
{
	addToRoot(otherTree->getIOLists());
}

void Tree::addAfter(const InputTrace & tr, const IOListContainer & cnt)
{
	std::shared_ptr<TreeNode> r = getRoot();
	std::shared_ptr<TreeNode> n = r->after(tr.cbegin(), tr.cend());

	if (n == nullptr)
	{
		return;
	}
	n->addToThisNode(cnt);
}


size_t Tree::size() {
    
    size_t theSize = 0;
    
    std::shared_ptr<TreeNode> r = getRoot();
    r->calcSize(theSize);
    
    
    return theSize;
}

Tree* Tree::clone() const
{
    return new Tree( this );
}

std::shared_ptr<Tree> Tree::Clone() const
{
    return std::shared_ptr<Tree>(clone());
}
