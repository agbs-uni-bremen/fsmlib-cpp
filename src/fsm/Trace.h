/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_TRACE_H_
#define FSM_FSM_TRACE_H_

#include <memory>
#include <vector>

#include "interface/FsmPresentationLayer.h"

class Trace
{
protected:
	/**
	The trace itself, represented by a list of int
	*/
	std::vector<int> trace;

	/**
	The presentation layer used by the trace
	*/
    std::shared_ptr<FsmPresentationLayer> presentationLayer;
public:
	/**
	Create an empty trace, with only one presentation layer
	@param presentationLayer The presentation layer used by the trace
	*/
    Trace(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	Create a trace
	@param trace The trace itself, represented by a list of int
	@param presentationLayer The presentation layer used by the trace
	*/
	Trace(const std::vector<int>& trace,
          const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

    Trace(const Trace& other);
	
	/**
	 * Add an element, at the end of the trace
	 */
	void add(const int e);
    
    /** 
     * Append a vector of int to the trace
     */
    void append(const std::vector<int>& traceToAppend);
    /**
     * Prepends the given vector of int to this trace.
     * @param traceToPrepend The given vector of int
     */
    void prepend(const std::vector<int>& traceToPrepend);

    /**
     * Appends the given trace to this trace.
     * @param traceToPrepend The given trace
     */
    void append(const Trace& traceToAppend);
    /**
     * Prepends the given trace to this trace.
     * @param traceToPrepend The given trace
     */
    void prepend(const Trace& traceToPrepend);

    /**
     * Determines if the given trace is a prefix of this trace.
     * @param other The given trace
     * @return {@code true}, if the given trace is a prefix of
     * this trace, {@code false}, otherwise.
     */
    bool isPrefix(const Trace& other) const;

    /**
     * Determines if the given trace is a suffix of this trace.
     * @param other The given trace
     * @return `true`, if the given trace is a suffix of
     * this trace, `false`, otherwise.
     */
    bool isSuffix(const Trace& other) const;
    /**
     * Determines if this trace is a prefix of the given trace.
     * @param other The given trace
     * @return {@code true}, if this trace is a prefix of
     * the given trace, {@code false}, otherwise.
     */
    bool isPrefixOf(const Trace& other) const;

    void removeElements(int n);

	/**
	Getter for the trace itself
	@return The trace itself, represented by a list of int
	*/
	std::vector<int> get() const;

	/**
	Getter for an iterator of the trace, pointing at the beginning
	@return The iterator
	*/
	std::vector<int>::const_iterator cbegin() const;

	/**
	Getter for an iterator of the trace, pointing at the end
	@return The iterator
	*/
	std::vector<int>::const_iterator cend() const;
    
    const std::shared_ptr<FsmPresentationLayer> getPresentationLayer() const { return presentationLayer; }

    std::vector<Trace> getPrefixes() const;

    static bool contains(const std::vector<std::shared_ptr<Trace>>& list, const Trace& trace);
    bool contains(const int io) const;

	/**
	Check wheter or not, the 2 trace are the same
	@param trace1 The first trace
	@param trace2 The second trace
	@return True if they are the same, false otherwise
	*/
	friend bool operator==(Trace const & trace1, Trace const & trace2);

    /**
     * Check whether or not the list of integers in trace1
     * is the same as the list of integers specified by trace2
     * @param trace1 The trace
     * @param trace2 A vector of integers, representing a trace
     * @return True if they are the same, false otherwise
     */
    friend bool operator==(Trace const & trace1, std::vector<int> const & trace2);

	/**
	Output the Trace to a standard output stream
	@param out The standard output stream to use
	@param trace The Trace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const Trace & trace);

    Trace& operator=(Trace&& other);
};
#endif //FSM_FSM_TRACE_H_
