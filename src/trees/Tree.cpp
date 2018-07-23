/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmLabel.h"
#include "trees/Tree.h"
#include "logging/easylogging++.h"

using namespace std;

void Tree::calcLeaves()
{
	leaves.clear();
	root->calcLeaves(leaves);
}

vector<shared_ptr<TreeNode const>> Tree::calcLeaves() const
{
    vector<shared_ptr<TreeNode const>> result;
	root->calcLeaves(result);
    return result;
}

void Tree::remove(const shared_ptr<TreeNode>& thisNode, const shared_ptr<TreeNode>& otherNode)
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

void Tree::printChildren(ostream & out, const shared_ptr<TreeNode>& top, const shared_ptr<int>& idNode) const
{
	int idNodeBase = *idNode;
	for (shared_ptr<TreeEdge> edge : *top->getChildren())
	{
		out << idNodeBase << " -> " << ++ *idNode << "[label = \"" << edge->getIO() << "\" ];" << endl;
		printChildren(out, edge->getTarget(), idNode);
    }
}

bool Tree::inPrefixRelation(vector<int> aPath, vector<int> bPath)
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

Tree::Tree(const shared_ptr<TreeNode>& root, const shared_ptr<FsmPresentationLayer>& presentationLayer)
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

shared_ptr<Tree> Tree::getSubTree(const shared_ptr<InputTrace>& alpha)
{
    shared_ptr<TreeNode> afterAlpha = root->after(alpha->cbegin(), alpha->cend());
    shared_ptr<TreeNode> cpyNode = afterAlpha->clone();
    return make_shared<Tree>(cpyNode, presentationLayer);
}

shared_ptr<TreeNode> Tree::getSubTree(shared_ptr< vector<int> > alpha) {
    
    return root->after(alpha->begin(),alpha->end());
}

IOListContainer Tree::getIOLists() const
{
	shared_ptr<vector<vector<int>>> ioll = make_shared<vector<vector<int>>>();

	for (shared_ptr<TreeNode const> n : calcLeaves())
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

void Tree::remove(const shared_ptr<Tree>& otherTree)
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
    TIMED_FUNC(timerObj);
	std::shared_ptr<std::vector<std::vector<int>>> ioll = std::make_shared<std::vector<std::vector<int>>>();
	std::shared_ptr<TreeNode> currentNode = root;
    ioll->push_back({FsmLabel::EPSILON});
	std::shared_ptr<std::vector<std::shared_ptr<TreeEdge>>> edges = currentNode->getChildren();

	for (size_t i=0; i < edges->size(); ++i)
	{
		std::shared_ptr<TreeEdge> edge = (*edges)[i];
		std::shared_ptr<TreeNode> node = edge->getTarget();
        const vector<int>& path = node->getPath();
        if (!IOListContainer::contains(ioll, path))
        {
            ioll->push_back(path);
        }
		std::shared_ptr<std::vector<std::shared_ptr<TreeEdge>>> newEdges = node->getChildren();
		std::copy (newEdges->begin(), newEdges->end(), std::back_inserter(*edges));
	}
	return IOListContainer(ioll, presentationLayer);
}

void Tree::add(const IOListContainer & tcl)
{
	root->add(tcl);
}

void Tree::addToRoot(const IOListContainer & tcl)
{
	root->addToThisNode(tcl);
}

void Tree::addToRoot(const vector<int> &lst)
{
    root->addToThisNode(lst);
}

void Tree::unionTree(const shared_ptr<Tree>& otherTree)
{
	addToRoot(otherTree->getIOLists());
}

void Tree::addAfter(const InputTrace & tr, const IOListContainer & cnt)
{
	shared_ptr<TreeNode> n = root->after(tr.cbegin(), tr.cend());

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
    root->calcSize(theSize);
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

shared_ptr<Tree> Tree::getPrefixRelationTree(const shared_ptr<Tree> & b)
{
    IOListContainer aIOlst = getIOLists();
    IOListContainer bIOlst = b->getIOLists();

    shared_ptr<vector<vector<int>>> aPrefixes = aIOlst.getIOLists();
    shared_ptr<vector<vector<int>>> bPrefixes = bIOlst.getIOLists();

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

int Tree::tentativeAddToRoot(const std::vector<int>& alpha) {
    return root->tentativeAddToThisNode(alpha.cbegin(), alpha.cend());
}





int Tree::tentativeAddToRoot(SegmentedTrace& alpha) const {
    
    int r;
    shared_ptr<TreeNode> n = root;
    
    for ( size_t i = 0; i < alpha.size(); i++ ) {
        shared_ptr<TraceSegment> seg = alpha.getSegments().at(i);
        r = n->tentativeAddToThisNode(seg->get()->cbegin(), seg->get()->cend(),n);
        if ( r > 0 ) return r;
    }
    
    return 0;

}

ostream & operator<<(ostream & out, Tree & t)
{
    out << t.str();
    return out;
}

string Tree::str()
{
    string str;
    calcLeaves();

    for (size_t k = 0; k < leaves.size(); ++k)
    {
        shared_ptr<TreeNode> leave = leaves.at(k);
        shared_ptr<TreeNode> l = static_pointer_cast<TreeNode>(leave);
        vector<int> path = l->getPath();
        for (size_t i = 0; i < path.size(); ++i)
        {
            if ( i > 0 )
            {
                str += ".";
            }

            str += "(" + presentationLayer->getInId(path.at(i)) + ")";
        }
        if (k != leaves.size() - 1)
        {
            str += ", ";
        }
    }
    return str;
}
