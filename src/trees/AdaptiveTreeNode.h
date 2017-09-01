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
    std::vector<int> getOutputPath();

    virtual AdaptiveTreeNode* clone() const;
    std::shared_ptr<AdaptiveTreeNode> Clone() const;

    friend bool operator==(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2);
    friend bool operator!=(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2);
};

#endif // ADAPTIVETREENODE_H
