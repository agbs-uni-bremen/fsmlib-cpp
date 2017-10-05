/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/Tree.h"

using namespace std;

void Tree::calcLeaves()
{
	leaves.clear();
	root->calcLeaves(leaves);
}

void Tree::remove(const shared_ptr<TreeNode> thisNode, const shared_ptr<TreeNode> otherNode)
{
	thisNode->deleteNode();

	for (shared_ptr<TreeEdge> e : *thisNode->getChildren())
	{
		shared_ptr<TreeEdge> eOther = otherNode->hasEdge(e);
		if (eOther != nullptr)
		{
			remove(e->getTarget(), eOther->getTarget());
		}
	}
}

void Tree::printChildren(ostream & out, const shared_ptr<TreeNode> top, const shared_ptr<int> idNode) const
{
	int idNodeBase = *idNode;
	for (shared_ptr<TreeEdge> edge : *top->getChildren())
	{
		out << idNodeBase << " -> " << ++ *idNode << "[label = \"" << edge->getIO() << "\" ];" << endl;
		printChildren(out, edge->getTarget(), idNode);
    }
}

bool Tree::inPrefixRelation(std::vector<int> aPath, std::vector<int> bPath)
{
    if (aPath.size() == 0 || bPath.size() == 0)
        return false;
    for (unsigned i = 0; i<aPath.size() && i < bPath.size(); i++)
    {
        if (aPath[i] != bPath[i])
        {
            return false;
        }
    }
    return true;
}

Tree::Tree(const shared_ptr<TreeNode> root, const shared_ptr<FsmPresentationLayer> presentationLayer)
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

shared_ptr<TreeNode> Tree::getRoot() const
{
	return root;
}

std::shared_ptr<Tree> Tree::getSubTree(const shared_ptr<InputTrace> alpha)
{
    shared_ptr<TreeNode> afterAlpha = getRoot()->after(alpha->cbegin(), alpha->cend());
    shared_ptr<TreeNode> cpyNode = afterAlpha->clone();
    return make_shared<Tree>(cpyNode, presentationLayer);
}

IOListContainer Tree::getIOLists()
{
	shared_ptr<vector<vector<int>>> ioll = make_shared<vector<vector<int>>>();
	calcLeaves();

	for (shared_ptr<TreeNode> n : leaves)
	{
		ioll->push_back(n->getPath());
	}

	return IOListContainer(ioll, presentationLayer);
}



IOListContainer Tree::getIOListsWithPrefixes()
{
    shared_ptr<vector<vector<int>>> ioll = make_shared<vector<vector<int>>>();
    
    // Create empty I/O-list as vector
    vector<int> thisVec; 
    
    // Perform in-order traversal of the tree
    // and create all I/O-lists.
    root->traverse(thisVec,ioll);
    
    return IOListContainer(ioll, presentationLayer);
}

void Tree::remove(const shared_ptr<Tree> otherTree)
{
	remove(getRoot(), otherTree->getRoot());
}

void Tree::toDot(ostream & out)
{
	out << "digraph Tree {" << endl;
	out << "\trankdir=TB;" << endl;//Top -> Bottom, to create a vertical graph
	out << "\tnode [shape = circle];" << endl;
	shared_ptr<int> id = make_shared<int>(0);
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
	shared_ptr<TreeNode> r = getRoot();
	r->add(tcl);
}

void Tree::addToRoot(const IOListContainer & tcl)
{
	shared_ptr<TreeNode> r = getRoot();
	r->addToThisNode(tcl);
}

void Tree::addToRoot(const vector<int> &lst)
{
    shared_ptr<TreeNode> r = getRoot();
    r->addToThisNode(lst);
}

void Tree::unionTree(const shared_ptr<Tree> otherTree)
{
	addToRoot(otherTree->getIOLists());
}

void Tree::addAfter(const InputTrace & tr, const IOListContainer & cnt)
{
	shared_ptr<TreeNode> r = getRoot();
	shared_ptr<TreeNode> n = r->after(tr.cbegin(), tr.cend());

	if (n == nullptr)
	{
		return;
	}
	n->addToThisNode(cnt);
}

bool Tree::isDefined(int y) const
{
    return root->isDefined(y);
}

size_t Tree::size() {
    
    size_t theSize = 0;
    
    shared_ptr<TreeNode> r = getRoot();
    r->calcSize(theSize);
    
    
    return theSize;
}

Tree* Tree::_clone() const
{
    return new Tree( this );
}

std::shared_ptr<Tree> Tree::Clone() const
{
    return std::shared_ptr<Tree>(_clone());
}

std::shared_ptr<Tree> Tree::getPrefixRelationTree(const std::shared_ptr<Tree> & b)
{
    IOListContainer aIOlst = getIOLists();
    IOListContainer bIOlst = b->getIOLists();

    std::shared_ptr<std::vector<std::vector<int>>> aPrefixes = aIOlst.getIOLists();
    std::shared_ptr<std::vector<std::vector<int>>> bPrefixes = bIOlst.getIOLists();

    shared_ptr<TreeNode> r = make_shared<TreeNode>();
    shared_ptr<Tree> tree = make_shared<Tree>(r, presentationLayer);

    if (aPrefixes->at(0).size() == 0 && bPrefixes->at(0).size() == 0)
    {
        return tree;
    }

    if (aPrefixes->at(0).size() == 0)
    {
        return b;
    }
    if (bPrefixes->at(0).size() == 0)
    {
        return shared_from_this();
    }


    for (auto aPrefix : *aPrefixes)
    {
        for (auto bPrefix : *bPrefixes)
        {
            if (inPrefixRelation(aPrefix, bPrefix))
            {
                r->addToThisNode(aPrefix);
                r->addToThisNode(bPrefix);
            }
        }
    }
    return tree;

}
