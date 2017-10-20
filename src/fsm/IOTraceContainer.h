#ifndef IOTRACECONTAINER_H
#define IOTRACECONTAINER_H

#include "fsm/IOTrace.h"
#include "trees/OutputTree.h"

class IOTraceContainer
{
private:
    std::shared_ptr<std::vector<IOTrace>> list;
    void removeRealPrefixes(const IOTrace & trc);
public:
    IOTraceContainer();
    IOTraceContainer(std::shared_ptr<std::vector<IOTrace>>& list);
    IOTraceContainer(std::shared_ptr<IOTrace> trace);
    std::shared_ptr<std::vector<IOTrace>> getList() const;
    /**
     * Adds the given trace to the container, only if the container
     * does not already contain a trace with the given inputs and outputs.
     * @param trc The given trace
     */
    void addUnique(IOTrace& trc);

    /**
     * Adds the given trace to the container, only if the container
     * does not already contain a trace with the given inputs and outputs.
     * Additionally all real prefixes of `trc` are being removed from the container.
     * @param trc The given trace.
     */
    void addUniqueRemovePrefixes(const IOTrace& trc);

    /**
     * Adds every trace from the given container to the container, only if the container
     * does not already contain a trace with the corresponding inputs and outputs.
     * Additionally all real prefixes of every trace from `cont` are being removed
     * from the container.
     * @param cont The given container.
     */
    void addUniqueRemovePrefixes(const IOTraceContainer& cont);

    /**
     * Adds the given trace to the container.
     * @param trc The given trace
     */
    void add(IOTrace& trc);

    /**
     * Every trace from the given container that is not already being held by
     * this container gets added to this container.
     * @param trc The given container
     */
    void addUnique(IOTraceContainer& container);

    /**
     * Adds every trace that is being represented by the given tree to the container,
     * only if the container dies not already contain a trace with the given inputs
     * and outputs.
     * @param tree The given tree.
     */
    void addUnique(OutputTree& tree);

    /**
     * Adds every trace from the given container to this container.
     * @param container The given container.
     */
    void add(IOTraceContainer& container);

    /**
     * Adds every trace that is being represented by the given tree to this container.
     * @param container The given tree.
     */
    void add(OutputTree& tree);

    /**
     * Checks if the container contains the given trace .
     * @param trace The given trace
     * @return {@code true}, if the container contains the given trace,
     * {@code false}, otherwise.
     */
    bool contains(IOTrace& trace) const;

    /**
     * Concatenates a given trace with each element of this container.
     * This modifies the container.
     * @param trace The given trace.
     */
    void concatenate(IOTrace& trace);

    /**
     * Concatenates each element of a given trace container with each
     * element of this container. This container will be modified.
     * @param container The given container of traces.
     */
    void concatenate(IOTraceContainer& container);

    /**
     * Concatenates the given input trace and the given output trace as an `IOTrace`
     * to the front of each trace in this container.
     * @param inputTrace The given input trace
     * @param outputTrace The given output trace
     */
    void concatenateToFront(InputTrace& inputTrace, OutputTrace outputTrace);

    /**
      Concatenates the given input/output trace to the front of each trace in this container.
     * @param iOTrace The given input/output trace.
     */
    void concatenateToFront(IOTrace& iOTrace);

    /**
     * Clears this container.
     */
    void clear();

    /**
     * Removes all occurrences of the given trace from this container.
     * @param trace The given trace
     */
    void remove (IOTrace& trace);

    /**
     * Removes all occurrebces of the traces in the given container from this container.
     * @param container The given container of traces that will be removed.
     */
    void remove (IOTraceContainer& container);

    /**
     * Returns all output traces.
     * @return All output traces
     */
    std::vector<OutputTrace> getOutputTraces() const;

    /**
     * Returns the size of the conteiner.
     * @return The size of the container
     */
    size_t size() const { return list->size(); }

    /**
     * Determines wether the container is empty.
     * @return `true`, if the container is empty, `false`, otherwise.
     */
    bool isEmpty() const {return size() == 0; }

    /**
     * Output the IOTraceContainer to a standard output stream
     * @param out The standard output stream to use
     * @param ot The IOTraceContainer to print
     * @return The standard output stream used, to allow user to cascade <<
     */
    friend std::ostream & operator<<(std::ostream & out, const IOTraceContainer & iot);
};

#endif // IOTRACECONTAINER_H
