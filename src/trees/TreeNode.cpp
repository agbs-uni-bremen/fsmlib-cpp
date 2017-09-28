/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "trees/TreeNode.h"

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

void TreeNode::setParent(const weak_ptr<TreeNode> pparent)
{
    parent = pparent;
}

weak_ptr<TreeNode> TreeNode::getParent() const
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

void TreeNode::remove(const shared_ptr<TreeNode> node)
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

void TreeNode::add(const shared_ptr<TreeEdge> edge)
{
    edge->getTarget()->setParent(shared_from_this());
    children->push_back(edge);
}

bool TreeNode::isLeaf() const
{
    return children->empty();
}

int TreeNode::getIO(const shared_ptr<TreeNode> node) const
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

shared_ptr<TreeEdge> TreeNode::hasEdge(const shared_ptr<TreeEdge> edge) const
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

vector<int> TreeNode::getPath()
{
    vector<int> path;
    
    shared_ptr<TreeNode> m = shared_from_this();
    shared_ptr<TreeNode> n = parent.lock();
    
    while (n != nullptr)
    {
        path.insert(path.begin(), n->getIO(m));
        m = n;
        n = n->getParent().lock();
    }
    
    return path;
}

bool TreeNode::superTreeOf(const shared_ptr<TreeNode> otherNode) const
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

void TreeNode::calcSize(size_t& theSize) {
    
    theSize++;
    
    for ( auto t : *children) {
        t->getTarget()->calcSize(theSize);
    }
    
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














