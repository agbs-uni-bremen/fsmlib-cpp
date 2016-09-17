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