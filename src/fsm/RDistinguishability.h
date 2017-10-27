#ifndef RDISTINGUISHABILITY_H
#define RDISTINGUISHABILITY_H

#include <vector>
#include <memory>
#include <map>

class FsmNode;
class InputOutputTree;
class FsmPresentationLayer;

/**
 * This class serves as data structure that is being used to calculate the
 * r-distinguishability for FSM nodes.
 */
class RDistinguishability
{
private:
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
    bool hBeenCalculated;
public:
    RDistinguishability(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

protected:
    /**
     * This maps holds for every `i` the states that are r(i)-distinguishable from the state
     * that corresponds to the `RDistinguishability` instance.
     */
    std::map<size_t, std::vector<std::shared_ptr<FsmNode>>> rDistinguishableWith;
    /**
     * This maps holds for every `i` the states that are not r(i)-distinguishable from the state
     * that corresponds to the `RDistinguishability` instance.
     */
    std::map<size_t, std::vector<std::shared_ptr<FsmNode>>> notRDistinguishableWith;

    /**
     * This map holds for every state an adaptive input sequence, that r-distinguishes
     * that state from the state that corresponds to the `RDistinguishability` instance.
     */
    std::map<std::shared_ptr<FsmNode>, std::shared_ptr<InputOutputTree>> adaptiveIOSequences;
public:
    /**
     * Removes a given state from the set that holds the states that are not r(i)-distinguishable
     * for a given `i`.
     * @param i The given non-r(i)-distinguishability index
     * @param node The given node to be removed
     * @return An iterator, pointing to the removed node, if the node has been found in the
     * non-r(i)-distinguisability set, an iterator, pointing to the end of the non-r(i)-distinguisability set,
     * if the node has not beend found.
     */
    std::vector<std::shared_ptr<FsmNode>>::iterator removeNotRDistinguishable(size_t i, std::shared_ptr<FsmNode> node);

    void initRDistinguishable(size_t i);

    /**
     * Adds a node to the set that holds the states that are r(i)-distinguishable
     * for a given `i`.
     * @param i The given r(i)-distinguishability index
     * @param node The given node to be added
     */
    void addRDistinguishable(size_t i, std::shared_ptr<FsmNode> node);

    /**
     * Adds a node to the set that holds the states that are not r(i)-distinguishable
     * for a given `i`.
     * @param i The given non-r(i)-distinguishability index
     * @param node The given node to be added
     */
    void addNotRDistinguishable(size_t i, std::shared_ptr<FsmNode> node);

    /**
     * Sets an empty set as non-r(i)-distinguishability set for a given index `i`.
     * @param i The given non-r(i)-distinguishability index
     */
    void addNotRDistinguishable(size_t i);

    /**
     * Sets the r-distinguishing adaptive input sequence for a given state.
     * @param otherNode The goven state
     * @param tree The given adaptive input sequence
     */
    void addAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode, std::shared_ptr<InputOutputTree> tree);

    /**
     * Returns all states that are r(i)-distinguishable for a given index `i` from the state
     * that corresponds to the `RDistinguishability` instance.
     *
     * @param i The given r(i)-distinguishability index
     * @return All r(i)-distinguishable states
     */
    std::vector<std::shared_ptr<FsmNode>> getRDistinguishableWith(size_t i);

    /**
     * Returns all states that are r-distinguishable from the state that corresponds to the
     * `RDistinguishability` instance.
     * @return All r-distinguishable states
     */
    std::vector<std::shared_ptr<FsmNode>> getRDistinguishableWith();

    /**
     * Returns all states that are non-r(i)-distinguishable for a given index `i` from the state
     * that corresponds to the `RDistinguishability` instance.
     *
     * @param i The given non-r(i)-distinguishability index
     * @return All non-r(i)-distinguishable states
     */
    std::vector<std::shared_ptr<FsmNode>> getNotRDistinguishableWith(size_t i);

    /**
     * Determines wether a given state is r(i)-distinguishable from the state that
     * corresponds to the `RDistinguishability` instance for a given index `i`.
     * @param i The given r(i)-distinguishability index
     * @param node The given state
     * @return `true`, if the given state is r(i)-distinguishable from the state that
     * corresponds to the `RDistinguishability` instance, `false`, otherwise
     */
    bool isRDistinguishableWith(size_t i, std::shared_ptr<FsmNode> node);

    /**
     * Determines wether a given state is r-distinguishable from the state that
     * corresponds to the `RDistinguishability` instance.
     * @param node The given state
     * @return `true`, if the given state is r-distinguishable from the state that
     * corresponds to the `RDistinguishability` instance, `false`, otherwise
     */
    bool isRDistinguishableWith(std::shared_ptr<FsmNode> node);

    /**
     * Determines wether every state of a given set of states is r-distinguishable from
     * the state that corresponds to the `RDistinguishability` instance.
     * @param node The given set of states
     * @return `true`, if every state from the given set of states is r-distinguishable
     * from the state that corresponds to the `RDistinguishability` instance, `false`, otherwise
     */
    bool isRDistinguishableWith(std::vector<std::shared_ptr<FsmNode>> nodes);
    /**
     * Determines if the state that corresponds to the `RDistinguishability` instance is
     * r-distinguishable.
     * @return `true`, if the state is r-distinguishable, `false`, otherwise.
     */
    bool isNotRDistinguishable();

    /**
     * Returns an adaptive input sequence, that r-distinguishes a given state from the state
     * that corresponds to the `RDistinguishability` instance.
     * @param otherNode The given state
     * @return An empty sequence, if the two states are not r-distinguishable, an adaptive
     * input sequence that r-distinguishes them, otherwise.
     */
    std::shared_ptr<InputOutputTree> getAdaptiveIOSequence(std::shared_ptr<FsmNode> otherNode);

    /**
     * Sets the r(i)-distinguishability and the non-r(i)-distinguishability sets for a given
     * index `i` to the same values as the r(i-1)-distinguishability and the
     * non-r(i-1)-distinguishability.
     * @param i The given r(i)-distinguishability index
     */
    void inheritDistinguishability(size_t i);

    /**
     * States wether the r-distinguishability has been calculated or not.
     * @return `true`, if the r-distinguishability has been calculated, `false`,
     * otherwise.
     */
    bool hasBeenCalculated() const;

    /**
     * Sets the flag that states wether the r-distinguishability has been calculated or not.
     * @param hasBeen The new flag value
     */
    void hasBeenCalculated(bool hasBeen);

};

#endif // RDISTINGUISHABILITY_H
