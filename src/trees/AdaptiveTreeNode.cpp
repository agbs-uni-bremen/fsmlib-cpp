#include "AdaptiveTreeNode.h"

using namespace std;

AdaptiveTreeNode::AdaptiveTreeNode(int input):
    input(input)
{

}

AdaptiveTreeNode::AdaptiveTreeNode()
{

}

AdaptiveTreeNode::AdaptiveTreeNode(const AdaptiveTreeNode* other):
    TreeNode (other)
{
    input = other->input;
}

int AdaptiveTreeNode::getInput()
{
    return input;
}


std::vector<int> AdaptiveTreeNode::getInputPath()
{
    vector<int> path;

    shared_ptr<AdaptiveTreeNode> m = static_pointer_cast<AdaptiveTreeNode>(shared_from_this());
    shared_ptr<AdaptiveTreeNode> n = static_pointer_cast<AdaptiveTreeNode>(parent.lock());

    while (n != nullptr)
    {
        path.insert(path.begin(), n->getInput());
        m = n;
        n = static_pointer_cast<AdaptiveTreeNode>(n->getParent().lock());
    }

    return path;
}

std::vector<int> AdaptiveTreeNode::getOutputPath()
{
    vector<int> path;

    shared_ptr<AdaptiveTreeNode> m = static_pointer_cast<AdaptiveTreeNode>(shared_from_this());
    shared_ptr<AdaptiveTreeNode> n = static_pointer_cast<AdaptiveTreeNode>(parent.lock());

    while (n != nullptr)
    {
        path.insert(path.begin(), n->getIO(m));
        m = n;
        n = static_pointer_cast<AdaptiveTreeNode>(n->getParent().lock());
    }

    return path;
}

bool AdaptiveTreeNode::superTreeOf(const std::shared_ptr<AdaptiveTreeNode> otherNode) const
{
    return input == otherNode->getInput() && TreeNode::superTreeOf(otherNode);
}

AdaptiveTreeNode* AdaptiveTreeNode::clone() const
{
    return new AdaptiveTreeNode( this );
}

std::shared_ptr<AdaptiveTreeNode> AdaptiveTreeNode::Clone() const
{
    std::shared_ptr<AdaptiveTreeNode> copy = std::shared_ptr<AdaptiveTreeNode>(clone());
    for (auto child: *(copy->children))
    {
        child->getTarget()->setParent(copy);
    }
    return copy;
}

bool operator==(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2)
{
    if (node1.input != node2.input)
    {
        return false;
    }
    if (node1.parent.expired()  || node2.parent.expired())
    {
        if (!node1.parent.expired() || !node2.parent.expired())
        {
            return false;
        }
    }
    if (*node1.parent.lock() != *node2.parent.lock())
    {
        return false;
    }
    if (node1.children == nullptr || node2.children == nullptr)
    {
        if (node1.children != nullptr || node2.children != nullptr)
        {
            return false;
        }
    }
    if (node1.children->size() != node2.children->size())
    {
        return false;
    }
    for (size_t i = 0; i < node1.children->size(); ++i)
    {
        if (*node1.children->at(i) != *node2.children->at(i))
        {
            return false;
        }
    }
    if (node1.deleted != node2.deleted)
    {
        return false;
    }
    return true;
}

bool operator!=(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2)
{
    return !(node1 == node2);
}
