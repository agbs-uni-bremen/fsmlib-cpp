/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_IOLISTCONTAINER_H_
#define FSM_TREES_IOLISTCONTAINER_H_

#include <memory>
#include <vector>

#include "fsm/Trace.h"
#include "interface/FsmPresentationLayer.h"

class IOListContainer
{
private:
    /**
     * The input list of this test cases
     */
    std::shared_ptr<std::vector<std::vector<int>>> iolLst;
    
    /**
     * The presentation layer used by this test cases
     */
    const std::shared_ptr<FsmPresentationLayer> presentationLayer;
    
    /**
     * Check whether or not it is the last list
     * @param maxInput The maximum input, to know the last one
     * @param lst The list
     * @return true if it is the last one, false otherwise
     */
    bool isLastLst(const int maxInput, const std::vector<int>& lst) const;
    
    /**
     * Return the next input trace of the same length as lst, in the order induced
     * by the input values 0..maxInput. For example, if maxInput = 1 and lst.size() = 2,
     * then the ordered enumeration of traces is 0.0, 0.1, 1.0, 1.1.
     * @param maxInput maximal input value to be used as trace element
     * @param lst preceding list, where the successor is to be constructed
     * @return Empty list, if no successor exists,
     * the successor list otherwise.
     */
    std::vector<int> nextLst(const int maxInput, const std::vector<int>& lst) const;
public:
    /**
     * Create a new IOListContainer (test cases)
     * @param iolLst The list of input traces
     * @param presentationLayer The presentation layer to use
     */
    IOListContainer(const std::shared_ptr<std::vector<std::vector<int>>> iolLst, const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     * Create an IOListContainer with input traces from length minLength
     * up to length maxLength.
     * For each length, all sequences with arbitrary inputs in range 0..maxInput
     * are created.
     * @param maxInput maximal input value to be created in an input trace.
     * @param minLength minimal length of the input traces to be created.
     * @param maxLength maximal length of a trace to be created.
     * @param presentationLayer The presentation layer to use
     */
    IOListContainer(const int maxInput,
                    const int minLength,
                    const int maxLength,
                    const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    IOListContainer(const std::shared_ptr<FsmPresentationLayer>
                        presentationLayer);
    
    /**
     * Getter for the input list
     * @return The input list
     */
    std::shared_ptr<std::vector<std::vector<int>>> getIOLists() const;
    
    /**
     * Add a new trace to the IOListContainer
     * @param trc The trace to add
     */
    void add(const Trace & trc);

    /**
     * Getter for the size of the IOListContainer
     * @return The size of the IOListContainer
     */
    int size() const;
    
    /**
     * Output the IOListContainer to a standard output stream
     * @param out The standard output stream to use
     * @param ot The IOListContainer to print
     * @return The standard output stream used, to allow user to cascade <<
     */
    friend std::ostream & operator<<(std::ostream & out, const IOListContainer & ot);
};
#endif //FSM_TREES_IOLISTCONTAINER_H_
