#ifndef ADAPTIVETREENODE_H
#define ADAPTIVETREENODE_H

#include "trees/TreeNode.h"
#include <unordered_map>

/**
 * Represents a tree node in an adaptive test case.
 * The difference to a "normal" tree node is that the adaptive tree node
 * contains an associated input.
 */
class AdaptiveTreeNode : public TreeNode
{
private:
    /**
     * The associated input.
     */
    int input;

    /**
    * An auxilliary mapping of the possible initial states, that can be assumed at this point in the adaptive experiment,
    * this node represents a subtree of, to the possible current states. This datastructure is used in the course of the creation
    * of an adaptive distinguishing sequence, which is represented by an InputOutputTree
    */
    std::shared_ptr<std::unordered_map<int, int>> initialToCurrentSet;
protected:
    /**
     * Constructs a new adaptive tree node by cloning a given one.
     * @param other The given tree node.
     */
    AdaptiveTreeNode(const AdaptiveTreeNode* other);
public:

    /**
     * Constructs a new adaptive tree node with the given associated input.
     * @param input The given input.
     */
    AdaptiveTreeNode(int input);

    /**
     * Constructs a new adaptive tree ndoe with no specified input.
     */
    AdaptiveTreeNode();

    /**
     * Returns the associated input.
     * @return The associated input.
     */
    int getInput();

    /**
     * Sets the associated input.
     * @param input The associated input.
     */
    void setInput(int input);

    /**
     * Calculates the input sequence that reaches this node.
     * @return The input sequence that reaches this node.
     */
    std::vector<int> getInputPath();

    /**
     * Calculates the output sequence that reaches this node.
     * @return The output sequence that reaches this node.
     */
    std::vector<int> getOutputPath();

    /**
    * Gets the auxilliary mapping of possible initial states to possible current states
    * @return the auxilliary mapping of possible initial states to possible current states
    */
    std::shared_ptr<std::unordered_map<int, int>>& getInitialToCurrentSet();

    /**
     * Determines wether this tree (node) is a super tree (node) of a given tree (node).
     * @param otherNode The given tree node
     * @return `true`, if this tree (node) is a super tree of the given tree (node),
     * `false`, otherwise.
     */
    bool superTreeOf(const std::shared_ptr<AdaptiveTreeNode>& otherNode) const;

    virtual AdaptiveTreeNode* _clone() const;
    std::shared_ptr<AdaptiveTreeNode> Clone() const;

    friend bool operator==(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2);
    friend bool operator!=(AdaptiveTreeNode const & node1, AdaptiveTreeNode const & node2);
};

#endif // ADAPTIVETREENODE_H
