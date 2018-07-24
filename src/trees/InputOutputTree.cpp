#include "InputOutputTree.h"
#include "AdaptiveTreeNode.h"

using namespace std;

InputOutputTree::InputOutputTree(const std::shared_ptr<AdaptiveTreeNode>& root,
                                 const std::shared_ptr<FsmPresentationLayer>& presentationLayer)
    : OutputTree(root->Clone(),
                 InputTrace(root->getInput(), presentationLayer),
                 presentationLayer)
{

}

InputOutputTree::InputOutputTree(const InputOutputTree* other):
    OutputTree(other)
{

}

bool InputOutputTree::isEmpty() const
{
    return root == nullptr || root->isLeaf();
}

bool InputOutputTree::contains(const InputOutputTree& other) const
{
    if (isEmpty() || other.isEmpty())
    {
        if (!isEmpty() || !other.isEmpty())
        {
            return false;
        }
    }

    if (!OutputTree::contains(other))
    {
        return false;
    }

    shared_ptr<AdaptiveTreeNode> thisRoot = static_pointer_cast<AdaptiveTreeNode>(getRoot());
    shared_ptr<AdaptiveTreeNode> otherRoot = static_pointer_cast<AdaptiveTreeNode>(other.getRoot());
    return thisRoot->superTreeOf(otherRoot);
}

IOListContainer InputOutputTree::getInputLists()
{
    std::shared_ptr<std::vector<std::vector<int>>> ioll = std::make_shared<std::vector<std::vector<int>>>();
    calcLeaves();

    for (std::shared_ptr<TreeNode> n : leaves)
    {
        std::shared_ptr<AdaptiveTreeNode> an = static_pointer_cast<AdaptiveTreeNode>(n);
        ioll->push_back(an->getInputPath());
    }

    return IOListContainer(ioll, presentationLayer);
}

IOListContainer InputOutputTree::getOutputLists()
{
    std::shared_ptr<std::vector<std::vector<int>>> ioll = std::make_shared<std::vector<std::vector<int>>>();
    calcLeaves();

    for (std::shared_ptr<TreeNode> n : leaves)
    {
        std::shared_ptr<AdaptiveTreeNode> an = static_pointer_cast<AdaptiveTreeNode>(n);
        ioll->push_back(an->getOutputPath());
    }

    return IOListContainer(ioll, presentationLayer);
}

InputOutputTree* InputOutputTree::_clone() const
{
    return new InputOutputTree( this );
}

std::shared_ptr<InputOutputTree> InputOutputTree::Clone() const
{
    return std::shared_ptr<InputOutputTree>(_clone());
}

ostream & operator<<(ostream & out, InputOutputTree & ot)
{
    out << ot.str();
    return out;
}

string InputOutputTree::str()
{
    string str;
    calcLeaves();

    //for (shared_ptr<TreeNode> leave : ot.leaves)
    for (size_t k = 0; k < leaves.size(); ++k)
    {
        shared_ptr<TreeNode> leave = leaves.at(k);
        shared_ptr<AdaptiveTreeNode> l = static_pointer_cast<AdaptiveTreeNode>(leave);
        vector<int> inputPath = l->getInputPath();
        vector<int> outputPath = l->getPath();
        for (size_t i = 0; i < inputPath.size(); ++i)
        {
            if ( i > 0 )
            {
                str += ".";
            }

            str += "(" + presentationLayer->getInId(inputPath.at(i)) + "/" + presentationLayer->getOutId(outputPath.at(i)) + ")";
        }
        if (k != leaves.size() - 1)
        {
            str += ", ";
        }
    }
    return str;
}
