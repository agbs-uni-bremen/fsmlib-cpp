#ifndef ADAPTIVETREENODE_H
#define ADAPTIVETREENODE_H

#include "trees/TreeNode.h"


class AdaptiveTreeNode : public TreeNode
{
private:
    int input;
public:
    AdaptiveTreeNode(int input);
    AdaptiveTreeNode();
    int getInput();
    std::vector<int> getInputPath();
};

#endif // ADAPTIVETREENODE_H
