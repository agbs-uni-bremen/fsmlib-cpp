#ifndef INPUTOUTPUTTREE_H
#define INPUTOUTPUTTREE_H

#include "trees/OutputTree.h"

class AdaptiveTreeNode;


class InputOutputTree : public OutputTree
{
protected:
    InputOutputTree(const InputOutputTree* other);
public:
    InputOutputTree(const std::shared_ptr<AdaptiveTreeNode> root, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    virtual InputOutputTree* clone() const;

    std::shared_ptr<InputOutputTree> Clone() const;

    friend std::ostream & operator<<(std::ostream & out, InputOutputTree & ot);
};

#endif // INPUTOUTPUTTREE_H
