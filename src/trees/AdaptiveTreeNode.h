#ifndef ADAPTIVETREENODE_H
#define ADAPTIVETREENODE_H

#include "trees/TreeNode.h"


class AdaptiveTreeNode : public TreeNode
{
private:
    int input;
protected:
    AdaptiveTreeNode(const AdaptiveTreeNode* other);
public:
    AdaptiveTreeNode(int input);
    AdaptiveTreeNode();
    int getInput();
    std::vector<int> getInputPath();

    virtual AdaptiveTreeNode* clone() const;
    std::shared_ptr<AdaptiveTreeNode> Clone() const;
};

#endif // ADAPTIVETREENODE_H
