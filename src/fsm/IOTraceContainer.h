#ifndef IOTRACECONTAINER_H
#define IOTRACECONTAINER_H

#include "fsm/IOTrace.h"

class IOTraceContainer
{
private:
    std::shared_ptr<std::vector<IOTrace>> list;
    const std::shared_ptr<FsmPresentationLayer> presentationLayer;
public:
    IOTraceContainer(const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    IOTraceContainer(std::shared_ptr<std::vector<IOTrace>>& list, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    IOTraceContainer(std::shared_ptr<IOTrace> trace, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    std::shared_ptr<std::vector<IOTrace>> getList() const;
    /**
     * Adds the given trace to the container, only if the container
     * does not already contain a trace with the given inputs and outputs.
     * @param trc The given trace
     */
    void addUnique(IOTrace& trc);

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
     * Adds every trace from the given container to this container.
     * @param container The given container.
     */
    void add(IOTraceContainer& container);

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
     * Returns the size of the conteiner.
     * @return The size of the container
     */
    size_t size() const { return list->size(); }

    /**
     * Output the IOTraceContainer to a standard output stream
     * @param out The standard output stream to use
     * @param ot The IOTraceContainer to print
     * @return The standard output stream used, to allow user to cascade <<
     */
    friend std::ostream & operator<<(std::ostream & out, const IOTraceContainer & iot);
};

#endif // IOTRACECONTAINER_H
