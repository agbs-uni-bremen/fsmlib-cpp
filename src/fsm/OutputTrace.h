/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_OUTPUTTRACE_H_
#define FSM_FSM_OUTPUTTRACE_H_

#include <vector>
#include <memory>

#include "fsm/Trace.h"
#include "interface/FsmPresentationLayer.h"

class OutputTrace : public Trace
{
public:
	/**
	Create an empty output trace, with only one presentation layer
	@param presentationLayer The presentation layer used by the trace
	*/
    OutputTrace(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	Create an output trace
	@param trace The trace itself, represented by a list of int
	@param presentationLayer The presentation layer used by the trace
	*/
    OutputTrace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer const>& presentationLayer);
    OutputTrace(const OutputTrace& other, size_t n, bool defaultToEmpty = false);
    OutputTrace(const OutputTrace& other);
    OutputTrace(const Trace& other);

    /**
     * Determines wether a given list of output traces contains a given output trace.
     * @param list The given list of output traces.
     * @param trace The given trace.
     * @return `true`, if `list` contains an element equal to `trace`, `false`, otherwise.
     */
    static bool contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const OutputTrace& trace);

    /**
     * Determines wether a given list of output traces contains an output trace that
     * contains the given (single) output.
     * @param list The given list of output traces.
     * @param output The given output.
     * @return `true`, if `list` contains an element that contains `output`,
     * `false`, otherwise.
     */
    static bool contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const int output);

    /**
     * Determines wether this output trace contains a given (single) output.
     * @param output The given output.
     * @return `true`, if this output trace contains `output`, `false`, otherwise.
     */
    bool contains(const int output) const;

    OutputTrace removeEpsilon() const { return OutputTrace(static_cast<Trace>(*this).removeEpsilon()); }
    OutputTrace removeLeadingEpsilons() const { return OutputTrace(static_cast<Trace>(*this).removeLeadingEpsilons()); }

	/**
	Output the OutputTrace to a standard output stream
	@param out The standard output stream to use
	@param trace The OutputTrace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
    friend std::ostream & operator<<(std::ostream & out, const OutputTrace & trace);

    OutputTrace& operator=(OutputTrace& other);
    OutputTrace& operator=(OutputTrace&& other);
};

namespace std {
  template <> struct hash<OutputTrace>
  {
    size_t operator()(const OutputTrace& trace) const
    {
        const Trace& t = static_cast<Trace>(trace);
        return std::hash<Trace>()(t);
    }
  };
}

#endif //FSM_FSM_OUTPUTTRACE_H_
