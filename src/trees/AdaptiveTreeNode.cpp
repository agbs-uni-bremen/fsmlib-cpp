#include "AdaptiveTreeNode.h"

using namespace std;

AdaptiveTreeNode::AdaptiveTreeNode(int input):
    input(input)
{

}

AdaptiveTreeNode::AdaptiveTreeNode()
{

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

