#ifndef INPUTOUTPUTTREE_H
#define INPUTOUTPUTTREE_H

#include "trees/OutputTree.h"

class AdaptiveTreeNode;


class InputOutputTree : public OutputTree
{
public:
    InputOutputTree(const std::shared_ptr<AdaptiveTreeNode> root, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    friend std::ostream & operator<<(std::ostream & out, InputOutputTree & ot);
};

#endif // INPUTOUTPUTTREE_H
