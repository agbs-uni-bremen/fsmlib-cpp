#ifndef INPUTOUTPUTTREE_H
#define INPUTOUTPUTTREE_H

#include "trees/OutputTree.h"

class AdaptiveTreeNode;

/**
 * This class represents an adaptive test case. Every node has an associated input,
 * the edges represent outputs
 */
class InputOutputTree : public OutputTree
{
protected:
    /**
     * Constructs a new InputOutputTree by cloning a given InputOutputTree.
     * @param other The given tree.
     */
    InputOutputTree(const InputOutputTree* other);
public:
    /**
     * Constructs a new InputOutputTree by cloning the given root.
     * @param root The given root.
     * @param presentationLayer The presentation layer.
     */
    InputOutputTree(const std::shared_ptr<AdaptiveTreeNode> root, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    /**
     * Determines wether the tree is empty (the tree has no root or the root is a leaf).
     * @return `true`, if the tree is empty, `false`, otherwise.
     */
    bool isEmpty() const;

    /**
     * Determines wether this tree is the super tree of a given tree.
     * @param other The given tree.
     * @return `true`, if this tree is a super tree of the given tree,
     * `false`, otherwise.
     */
    bool contains(const InputOutputTree& other) const;

    virtual InputOutputTree* _clone() const;

    /**
     * Calculates all input traces, that are being represented by the tree.
     * @return All input traces.
     */
    IOListContainer getInputLists();

    /**
     * Calculates all output traces, that are being represented by the tree.
     * @return All output traces.
     */
    IOListContainer getOutputLists();    

    std::shared_ptr<InputOutputTree> Clone() const;

    friend std::ostream & operator<<(std::ostream & out, InputOutputTree & ot);
    std::string str();

};

#endif // INPUTOUTPUTTREE_H
