#ifndef IOTREECONTAINER_H
#define IOTREECONTAINER_H

#include <memory>
#include <vector>

#include "interface/FsmPresentationLayer.h"
#include "trees/InputOutputTree.h"

class IOTreeContainer
{
private:
    /**
     * The list containing the InputOutputTrees
     */
    std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> list;
    const std::shared_ptr<FsmPresentationLayer> presentationLayer;
    void removePrefixes(std::shared_ptr<InputOutputTree> tree);
public:
    IOTreeContainer(const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    IOTreeContainer(const std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> list, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    std::shared_ptr<std::vector<std::shared_ptr<InputOutputTree>>> getList() const;

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
