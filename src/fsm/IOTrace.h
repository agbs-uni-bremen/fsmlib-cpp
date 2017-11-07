/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_IOTRACE_H_
#define FSM_FSM_IOTRACE_H_

#include <iostream>
#include <vector>

#include "fsm/InputTrace.h"
#include "fsm/OutputTrace.h"

class FsmNode;

class IOTrace
{
private:
	/**
	The input trace contained into the iotrace
	*/
	InputTrace inputTrace;

	/**
	The output trace contained into the iotrace
	*/
	OutputTrace outputTrace;

    /**
     * The node that can be reached via this trace
     */
    std::shared_ptr<FsmNode> targetNode;
public:
	/**
	Create an iotrace from one inut trace and one output trace
	@param i The input trace contained into the iotrace
	@param o The output trace contained into the iotrace
	*/
    IOTrace(const InputTrace & i, const OutputTrace & o, std::shared_ptr<FsmNode> targetNode = nullptr);
    IOTrace(const int i, const int o, std::shared_ptr<FsmPresentationLayer> pl);
    IOTrace(const int i, const int o, std::shared_ptr<FsmNode> targetNode, std::shared_ptr<FsmPresentationLayer> pl);
    IOTrace(const IOTrace & ioTrace);
    IOTrace(const IOTrace & ioTrace, int n, std::shared_ptr<FsmNode> targetNode = nullptr);
    IOTrace(std::shared_ptr<FsmPresentationLayer> pl);

	/**
	Getter for the input trace
	@return The input trace contained into the iotrace
	*/
	InputTrace getInputTrace() const;

	/**
	Getter for the output trace
	@return The output trace contained into the iotrace
	*/
	OutputTrace getOutputTrace() const;

    std::vector<IOTrace> getPrefixes(bool proper = false) const;

    void setTargetNode(std::shared_ptr<FsmNode> target)
    {
        targetNode = target;
    }

    std::shared_ptr<FsmNode> getTargetNode()
    {
        return targetNode;
    }

    size_t size() const;

    /**
     * Appends the given trace to this trace.
     * @param other The given trace
     */
    void append(IOTrace& other);
    /**
     * Prepends the given trace to this trace.
     * @param other The given trace
     */
    void prepend(IOTrace& other);

    void append(int input, int output);

    /**
     * Determines if the given trace is a prefix of this trace.
     * @param other The given trace
     * @return {@code true}, if the given trace is a prefix of
     * this trace, {@code false}, otherwise.
     */
    bool isPrefix(const IOTrace& other) const;

    /**
     * Determines if the given trace is a suffix of this trace.
     * @param other The given trace
     * @return `true`, if the given trace is a suffix of
     * this trace, `false`,, otherwise.
     */
    bool isSuffix(const IOTrace& other) const;
    /**
     * Determines if this trace is a prefix of the given trace.
     * @param other The given trace
     * @return {@code true}, if this trace is a prefix of
     * the given trace, {@code false}, otherwise.
     */
    bool isPrefixOf(const IOTrace& other) const;

    static std::shared_ptr<IOTrace> getEmptyTrace(std::shared_ptr<FsmPresentationLayer> pl);

	/**
	Output the IOTrace to a standard output stream
	@param out The standard output stream to use
	@param trace The IOTrace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const IOTrace & trace);
    friend bool operator==(IOTrace const & iOTrace1, IOTrace const & iOTrace2);
    IOTrace& operator= (IOTrace&& other);
    std::string toRttString() const;
};
#endif //FSM_FSM_IOTRACE_H_
