#include "InputOutputTree.h"
#include "AdaptiveTreeNode.h"

using namespace std;

InputOutputTree::InputOutputTree(const std::shared_ptr<AdaptiveTreeNode> root, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
    : OutputTree(root, InputTrace({root->getInput()}, presentationLayer), presentationLayer)
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

InputOutputTree* InputOutputTree::clone() const
{
    return new InputOutputTree( this );
}

std::shared_ptr<InputOutputTree> InputOutputTree::Clone() const
{
    return std::shared_ptr<InputOutputTree>(clone());
}

ostream & operator<<(ostream & out, InputOutputTree & ot)
{
    ot.calcLeaves();
    vector<shared_ptr<TreeNode>> leaves = ot.getLeaves();

    for (shared_ptr<TreeNode> leave : leaves)
    {
        shared_ptr<AdaptiveTreeNode> l = static_pointer_cast<AdaptiveTreeNode>(leave);
        vector<int> inputPath = l->getInputPath();
        vector<int> outputPath = l->getPath();
        for (size_t i = 0; i < inputPath.size(); ++i)
        {
            if ( i > 0 )
            {
                out << ".";
            }

            out << "(" << ot.presentationLayer->getInId(inputPath.at(i)) << "/" << ot.presentationLayer->getOutId(outputPath.at(i)) << ")";
        }
        out << "\n";

    }
    return out;
}
