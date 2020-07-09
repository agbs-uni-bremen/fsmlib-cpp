/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_IOTRACE_H_
#define FSM_FSM_IOTRACE_H_

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
    std::weak_ptr<FsmNode> targetNode;
public:
	/**
	Create an iotrace from one inut trace and one output trace
	@param i The input trace contained into the iotrace
	@param o The output trace contained into the iotrace
	*/
    IOTrace(const InputTrace & i, const OutputTrace & o, std::shared_ptr<FsmNode> targetNode = nullptr);
    IOTrace(const int i, const int o, std::shared_ptr<FsmPresentationLayer const> pl);
    IOTrace(const Trace& i, const Trace& o);
    IOTrace(const int i, const int o, std::shared_ptr<FsmNode> targetNode, std::shared_ptr<FsmPresentationLayer> pl);
    IOTrace(const IOTrace & ioTrace);
    IOTrace(const IOTrace& ioTrace, const IOTrace& append, bool prepend = false);
    IOTrace(const IOTrace & ioTrace, int n, std::shared_ptr<FsmNode> targetNode = nullptr);
    IOTrace(std::shared_ptr<FsmPresentationLayer> pl);
    IOTrace(const IOTrace& other, size_t n, bool defaultToEmpty = false);

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

    std::shared_ptr<FsmNode> getTargetNode() const
    {
        return targetNode.lock();
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

    bool isEmptyTrace() const;

    /**
     * Determines if the given trace is a prefix of this trace.
     * @param other The given trace
     * @return {@code true}, if the given trace is a prefix of
     * this trace, {@code false}, otherwise.
     */
    bool isPrefix(const IOTrace& other, bool proper = false, bool allowEmpty = true) const;

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

    /**
     * Returns the remaining suffix of this trace after removing its prefix.
     * @param prefix The given prefix
     * @return  The remaining suffix.
     */
    IOTrace getSuffix(const IOTrace& prefix) const;

    static std::shared_ptr<IOTrace> getEmptyTrace(std::shared_ptr<FsmPresentationLayer> pl);

    IOTrace removeEpsilon() const { return IOTrace(inputTrace.removeEpsilon(), outputTrace.removeEpsilon(), targetNode.lock()); }
    IOTrace removeLeadingEpsilons() const { return IOTrace(inputTrace.removeLeadingEpsilons(), outputTrace.removeLeadingEpsilons(), targetNode.lock()); }

	/**
	Output the IOTrace to a standard output stream
	@param out The standard output stream to use
	@param trace The IOTrace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const IOTrace & trace);
    friend bool operator==(IOTrace const & iOTrace1, IOTrace const & iOTrace2);
    friend bool operator!=(IOTrace const & iOTrace1, IOTrace const & iOTrace2)
    {
        return !(iOTrace1 == iOTrace2);
    }
    friend bool operator<=(IOTrace const & iOTrace1, IOTrace const & iOTrace2);
    IOTrace& operator= (IOTrace& other);
    IOTrace& operator= (IOTrace&& other);
    std::string toRttString() const;
    
    /**
     * Check equality of IOTraces
     */
    friend bool operator==(IOTrace const& trc1, IOTrace const& trc2);
    
    /**
     * Define a total order on IOTraces (used in sorting algorithms).
     * @param trace to compare against
     * @return True if this trace is defined to be less than the other trace, false otherwise
     */
    bool operator<(IOTrace const &other) const;
};

namespace std {
  template <> struct hash<IOTrace>
  {
    size_t operator()(const IOTrace& trace) const
    {

        size_t seed = 0;
        std::hash<size_t> hasher;

        size_t iH = hash<InputTrace>()(trace.getInputTrace());
        size_t oH = hash<OutputTrace>()(trace.getOutputTrace());

        seed ^= hasher(iH) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(oH) + 0x9e3779b9 + (seed<<6) + (seed>>2);

        return seed;
    }
  };
}

#endif //FSM_FSM_IOTRACE_H_
