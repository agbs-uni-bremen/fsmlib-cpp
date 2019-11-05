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
	std::shared_ptr<FsmPresentationLayer const> presentationLayer;

    Trace(const std::vector<int>::const_iterator& begin,
          const std::vector<int>::const_iterator& end,
          const std::shared_ptr<FsmPresentationLayer const>& presentationLayer);

public:
	/**
	Create an empty trace, with only one presentation layer
	@param presentationLayer The presentation layer used by the trace
	*/
    Trace(const std::shared_ptr<FsmPresentationLayer const>& presentationLayer);

	/**
	Create a trace
	@param trace The trace itself, represented by a list of int
	@param presentationLayer The presentation layer used by the trace
	*/
	Trace(const std::vector<int>& trace,
          const std::shared_ptr<FsmPresentationLayer const>& presentationLayer);
	
    Trace(const Trace& other);
    /**
     * Creates a trace based on the given trace `other` and skips the first `n` symbols.
     * @param other The given trace the new trace will be based on.
     * @param n The number of symbols to be skipped.
     * @param defaultToEmpty If `true` and `n` is larger than the `other` trace's size, a trace containing ε
     * will be created. If `false` and `n` is larger than the `other` trace's size, a trace not containing
     * any symbol will be created.
     */
    Trace(const Trace& other, size_t n, bool defaultToEmpty = false);

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

    bool isEmptyTrace() const;

    /**
     * Determines if the given trace is a prefix of this trace.
     * @param other The given trace
     * @return {@code true}, if the given trace is a prefix of
     * this trace, {@code false}, otherwise.
     */
    bool isPrefix(const Trace& other, bool proper = false, bool allowEmpty = true) const;

    Trace removeEpsilon() const;
    Trace removeLeadingEpsilons() const;

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

    /**
     * Creates a trace based on this trace and skips the first `n` symbols.
     * @param n The number of symbols to be skipped.
     * @param defaultToEmpty If `true` and `n` is larger than the `other` trace's size, a trace containing ε
     * will be created. If `false` and `n` is larger than the `other` trace's size, a trace not containing
     * any symbol will be created.
     */
    std::shared_ptr<Trace> getSuffix(size_t n, bool defaultToEmpty = false) const;

    /**
     * Creates a trace based on this trace and takes only the first `n` symbols into account.
     * @param n The number of symbols to be taken into account.
     * @param defaultToEmpty If `true` and `n` is larger than the `other` trace's size, a trace containing ε
     * will be created. If `false` and `n` is larger than the `other` trace's size, a trace not containing
     * any symbol will be created.
     */
    std::shared_ptr<const Trace> getPrefix(size_t n, bool defaultToEmpty = false) const;

    Trace getSuffix(const Trace& prefix) const;

    void removeElements(int n);

	/**
	Getter for the trace itself
	@return The trace itself, represented by a list of int
	*/
	std::vector<int> get() const;

    size_t size() const;

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
    
    std::shared_ptr<FsmPresentationLayer const> getPresentationLayer() const { return presentationLayer; }

    std::vector<Trace> getPrefixes(bool proper = false) const;

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
     * Define a total order on traces (used in sorting algorithms).
     * @param trace to compare against
     * @return True if this trace is defined to be less than the other trace, false otherwise
     */
    bool operator<(Trace const &other) const;

	/**
	Output the Trace to a standard output stream
	@param out The standard output stream to use
	@param trace The Trace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const Trace & trace);

    Trace& operator=(Trace& other);
    Trace& operator=(Trace&& other);

};

namespace std {
    template <> struct hash<Trace>
    {
        size_t operator()(const Trace& trace) const
        {
            std::hash<int> hasher;
            size_t seed = 0;
            for (int i : trace.get()) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }
            return seed;
        }
        size_t operator()(const shared_ptr<Trace>& trace) const
        {
            std::hash<int> hasher;
            size_t seed = 0;
            for (int i : trace->get()) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }
            return seed;
        }
    };
    template <> struct equal_to<Trace>
    {
        bool operator() (const Trace& a, const Trace& b) const 
        {
            return a == b;
        }

        bool operator() (const shared_ptr<const Trace>& a, const shared_ptr<const Trace>& b) const
        {
            if (a == b)
            {
                return true;
            }
            return *a == *b;
        }
    };
}

#endif //FSM_FSM_TRACE_H_
