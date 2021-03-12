/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "trees/TreeNode.h"
#include "trees/TreeEdge.h"
#include "trees/IOListContainer.h"
#include <deque>
#include <algorithm>

using namespace std;

TreeNode::TreeNode()
: parent(weak_ptr<TreeNode>()), children(make_shared<vector<shared_ptr<TreeEdge>>>()), deleted(false)
{
    
}

std::shared_ptr<TreeNode> TreeNode::clone() const
{
    shared_ptr<TreeNode> clone = make_shared<TreeNode>();
    clone->getChildren()->reserve(children->size());
    
    for (const auto& c : *children)
    {
        std::shared_ptr<TreeEdge> childClone = c->clone();
        childClone->getTarget()->setParent(clone);
        clone->getChildren()->push_back(childClone);
    }
    return clone;
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

void TreeNode::setParent(const weak_ptr<TreeNode>& parent)
{
    this->parent = parent;
}

weak_ptr<TreeNode> TreeNode::getParent() const
{
    return parent;
}

void TreeNode::deleteSingleNode()
{
    deleted = true;
    
    if (!isLeaf())
    {
        return;
    }
    
    shared_ptr<TreeNode> c = shared_from_this();
    shared_ptr<TreeNode> t = parent.lock();
    
    if (t != nullptr)
    {
        t->remove(c);
    }
}

void TreeNode::deleteNode()
{
    deleted = true;
    
    if (!isLeaf())
    {
        return;
    }
    
    shared_ptr<TreeNode> c = shared_from_this();
    shared_ptr<TreeNode> t = parent.lock();
    
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

shared_ptr<vector<shared_ptr<TreeEdge>>> TreeNode::getChildren() const
{
    return children;
}

void TreeNode::remove(const shared_ptr<TreeNode>& node)
{
    for (shared_ptr<TreeEdge> e : *children)
    {
        if (e->getTarget() == node)
        {
            children->erase(find(children->begin(), children->end(), e));
            break;
        }
    }
}

void TreeNode::calcLeaves(vector<shared_ptr<TreeNode>>& leaves)
{
    if (isLeaf())
    {
        leaves.push_back(shared_from_this());
    }
    else
    {
        for (shared_ptr<TreeEdge> e : *children)
        {
            e->getTarget()->calcLeaves(leaves);
        }
    }
}

void TreeNode::calcLeaves(vector<shared_ptr<TreeNode const>>& leaves) const
{
    if (isLeaf())
    {
        leaves.push_back(shared_from_this());
    }
    else
    {
        for (shared_ptr<TreeEdge> e : *children)
        {
            e->getTarget()->calcLeaves(leaves);
        }
    }
}

void TreeNode::add(const shared_ptr<TreeEdge>& edge)
{
    edge->getTarget()->setParent(shared_from_this());
    children->push_back(edge);
}

bool TreeNode::isLeaf() const
{
    return children->empty();
}

int TreeNode::getIO(const shared_ptr<TreeNode const>& node) const
{
    for (shared_ptr<TreeEdge> e : *children)
    {
        if (e->getTarget() == node)
        {
            return e->getIO();
        }
    }
    exit(EXIT_FAILURE);
}

shared_ptr<TreeEdge> TreeNode::hasEdge(const shared_ptr<TreeEdge>& edge) const
{
    for (shared_ptr<TreeEdge> g : *children)
    {
        if (g->getIO() == edge->getIO())
        {
            return g;
        }
    }
    return nullptr;
}

vector<int> TreeNode::getPath() const
{
    //As all path elements are inserted at the front, a deque is better to use
    //than a vector
    deque<int> path;
    
    shared_ptr<TreeNode const> m = shared_from_this();
    shared_ptr<TreeNode const> n = parent.lock();
    
    while (n != nullptr)
    {
        path.insert(path.begin(), n->getIO(m));
        m = n;
        n = n->getParent().lock();
    }
    //Copy deque to vector. As this is a sequence of contiguous memory access
    //operations, this should be quite fast.
    vector<int> result(path.begin(), path.end());
    
    return result;
}

bool TreeNode::superTreeOf(const shared_ptr<TreeNode>& otherNode) const
{
    
    if (children->size() < otherNode->children->size())
    {
        return false;
    }
    
    for (shared_ptr<TreeEdge> eOther : *otherNode->children)
    {
        int y = eOther->getIO();
        bool yFound = false;
        
        for (shared_ptr<TreeEdge> eMine : *children)
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
    for (shared_ptr<TreeEdge> e : *treeNode1.children)
    {
        int y = e->getIO();
        bool yFound = false;
        
        for (shared_ptr<TreeEdge> eOther : *treeNode2.children)
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

shared_ptr<TreeNode> TreeNode::add(const int x)
{
    for (shared_ptr<TreeEdge> e : *getChildren())
    {
        if (e->getIO() == x)
        {
            return e->getTarget();
        }
    }
    
    shared_ptr<TreeNode> tgt = make_shared<TreeNode>();
    add(make_shared<TreeEdge>(x, tgt));
    return tgt;
}

void TreeNode::add(vector<int>::const_iterator lstIte, const vector<int>::const_iterator end)
{
    /*There may be no next list element, when this method is called*/
    if (lstIte == end)
    {
        return;
    }
    
    /*Which input is represented by the list iterator?*/
    int x = *lstIte++;
    
    for (shared_ptr<TreeEdge> e : *getChildren())
    {
        /*Is there already an edge labelled with this input?*/
        if (e->getIO() == x)
        {
            /*We do not need to extend the tree, but follow the existing edge*/
            shared_ptr<TreeNode> nTgt = e->getTarget();
            nTgt->add(lstIte, end);
            return;
        }
    }
    
    /*No edge labelled with x exists for this node.
     Therefore one has to be created*/
    shared_ptr<TreeNode> newNode = make_shared<TreeNode>();
    newNode->setParent(shared_from_this());
    getChildren()->push_back(make_shared<TreeEdge>(x, newNode));
    newNode->add(lstIte, end);
}

void TreeNode::add(const IOListContainer & tcl)
{
    /*First delegate the work to the children*/
    for (shared_ptr<TreeEdge> e : *getChildren())
    {
        shared_ptr<TreeNode> nTgt = e->getTarget();
        nTgt->add(tcl);
    }
    
    /*Now append each input sequence in tcl to this node,
     using the special strategy of the add(lstIte) operation*/
    for (vector<int>& lst : *tcl.getIOLists())
    {
        add(lst.cbegin(), lst.cend());
    }
}


int TreeNode::tentativeAddToThisNode(vector<int>::const_iterator start,
                                     vector<int>::const_iterator stop) {
    
    // If we have reached the end, the trace is fully contained in the tree
    if ( start == stop ) return 0;
    
    // If this is a node without children (i.e., a leaf),
    // we just have to extend the tree, without creating a new branch.
    if ( children->empty() ) return 1;
    
    // Now we have to check whether an existing edge
    // is labelled with inout *start
    int x = *start;
    for ( auto e : *children ) {
		if (e->getIO() == x) {
			shared_ptr<TreeNode> next = e->getTarget();
			return next->tentativeAddToThisNode(++start, stop);
		}
    }
    
    // Adding this trace requires a new branch in the tree,
    // this means, an additional test case.
    return 2;
}

int TreeNode::tentativeAddToThisNode(vector<int>::const_iterator start,
                                     vector<int>::const_iterator stop,
                                     std::shared_ptr<TreeNode>& n) {
    
    n = shared_from_this();
    
    // If we have reached the end, the trace is fully contained in the tree
    if ( start == stop ) return 0;
    
    // If this is a node without children (i.e., a leaf),
    // we just have to extend the tree, without creating a new branch.
    if ( children->empty() ) return 1;
    
    // Now we have to check whether an existing edge
    // is labelled with inout *start
    int x = *start;
    for ( auto e : *children ) {
        if ( e->getIO() == x ) {
            shared_ptr<TreeNode> next = e->getTarget();
            return next->tentativeAddToThisNode(++start,stop,n);
        }
    }
    
    // Adding this trace requires a new branch in the tree,
    // this means, an additional test case.
    return 2;
}



void TreeNode::addToThisNode(const IOListContainer& tcl)
{
    /*Append each input sequence in tcl to this node,
     using the special strategy of the add(lstIte) operation*/
    for (vector<int>& lst : *tcl.getIOLists())
    {
        add(lst.cbegin(), lst.cend());
    }
}

void TreeNode::addToThisNode(const vector<int> &lst)
{
    add(lst.cbegin(), lst.cend());
}

shared_ptr<TreeNode> TreeNode::after(vector<int>::const_iterator lstIte, const vector<int>::const_iterator end)
{
    if (lstIte != end)
    {
        int x = *lstIte ++;
        
        for (shared_ptr<TreeEdge> e : *getChildren())
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

std::shared_ptr<TreeNode> TreeNode::after(const int y) const
{
    for (std::shared_ptr<TreeEdge> edge : *children)
    {
        if (edge->getIO() == y)
        {
            return edge->getTarget();
        }
    }
    return nullptr;
}

bool TreeNode::isDefined(int y) const
{
    for (std::shared_ptr<TreeEdge> edge : *children)
    {
        if (edge->getIO() == y)
        {
            return true;
        }
    }
    return false;
}

void TreeNode::calcSize(size_t& theSize) {
    
    theSize++;
    
    for ( auto t : *children) {
        t->getTarget()->calcSize(theSize);
    }
    
}

TreeNode* TreeNode::_clone() const
{
    return new TreeNode( this );
}

std::shared_ptr<TreeNode> TreeNode::Clone() const
{
    std::shared_ptr<TreeNode> copy = std::shared_ptr<TreeNode>(_clone());
    for (auto child: *(copy->children))
    {
        child->getTarget()->setParent(copy);
    }
    return copy;
}

void TreeNode::traverse(vector<int>& v,
                        shared_ptr<vector<vector<int>>> ioll) {
    
    // traverse all edges to child nodes
    for ( auto e : *children ) {
        
        int io = e->getIO();
        shared_ptr<TreeNode> n = e->getTarget();
        v.push_back(io);
        
        n->traverse(v,ioll);
        
        // Pop the last element in v
        v.pop_back();
        
    }
    
    // add v to vector of I/O-lists
    ioll->push_back(v);
    
    
}



std::shared_ptr<TreeNode> TreeNode::getIntersectionNode(const std::shared_ptr<TreeNode> &otherNode) {
    std::shared_ptr<TreeNode> root = std::make_shared<TreeNode>();

    for (auto thisEdge : *children) {

        int thisIO = thisEdge->getIO();
        for (auto otherEdge : *otherNode->children) {
            if (thisIO == otherEdge->getIO()) {
                std::shared_ptr<TreeNode> target = thisEdge->getTarget()->getIntersectionNode(otherEdge->getTarget());
                std::shared_ptr<TreeEdge> edge = std::make_shared<TreeEdge>(thisIO,target);
                root->add(edge);
            }
        }
    }
    return root;
}







