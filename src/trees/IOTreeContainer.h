#ifndef IOTREECONTAINER_H
#define IOTREECONTAINER_H

#include <memory>
#include <vector>

#include "interface/FsmPresentationLayer.h"
#include "trees/InputOutputTree.h"

/**
 * This class is a container for InputOutputTrees.
 */
class IOTreeContainer
{
private:
    /**
     * The list containing the InputOutputTrees
     */
    std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> list;
    const std::shared_ptr<FsmPresentationLayer> presentationLayer;

    /**
     * Removes all trees from the container, that are prefixes of the given tree.
     * @param tree The given tree
     */
    void removePrefixes(std::shared_ptr<InputOutputTree> tree);
public:
    /**
     * Constructs an empty container, with no elements.
     * @param presentationLayer The presentation layer
     */
    IOTreeContainer(const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    /**
     * Constructs a container that acquires the elements of the given list.
     * @param list The given list
     * @param presentationLayer The presentation layer
     */
    IOTreeContainer(const std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> list, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    /**
     * Returns the underlying list.
     * @return The lsit containing the elements
     */
    std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> getList() const;

    /**
     * Calculates all possible input traces from all tree elements from this container.
     * Prefixes and duplicates of resulting elements are being removed.
     * @return A set of input traces
     */
    IOListContainer toIOList() const;

    /**
     * Add a new tee to the IOTreeContainer if it doesn't exist already
     * @param tree The tree that will be added.
     */
    void addUnique(std::shared_ptr<InputOutputTree> tree);

    /**
     * Add a new tee to the IOTreeContainer if it doesn't exist already.
     * Additionally, remove all real prefixes, regardless of whether
     * the {@code tree} has been added.
     * @param tree The tree that will be added.
     */
    void addUniqueRemovePrefixes(std::shared_ptr<InputOutputTree> tree);

    /**
     * Returns the number of elements in the container.
     * @return The number of elements in thee container
     */
    size_t size() const;

    /**
     * Output the IOTreeContainer to a standard output stream
     * @param out The standard output stream to use
     * @param ot The IOTreeContainer to print
     * @return The standard output stream used, to allow user to cascade <<
     */
    friend std::ostream & operator<<(std::ostream & out, const IOTreeContainer & ot);
};

#endif // IOTREECONTAINER_H
