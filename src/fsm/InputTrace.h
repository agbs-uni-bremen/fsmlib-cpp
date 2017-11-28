/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_INPUTTRACE_H_
#define FSM_FSM_INPUTTRACE_H_

#include <iostream>
#include <vector>

#include "fsm/Trace.h"
#include "interface/FsmPresentationLayer.h"

class InputTrace : public Trace
{
public:
	/**
	Create an empty input trace, with only one presentation layer
	@param presentationLayer The presentation layer used by the trace
	*/
    InputTrace(const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

	/**
	Create an input trace
	@param trace The trace itself, represented by a list of int
	@param presentationLayer The presentation layer used by the trace
	*/
	InputTrace(const std::vector<int>& trace,
               const std::shared_ptr<FsmPresentationLayer>& presentationLayer);
    InputTrace(int x, const std::shared_ptr<FsmPresentationLayer>& presentationLayer);

    InputTrace(const InputTrace& other, size_t n, bool defaultToEmpty = false);
    InputTrace(const InputTrace& other);
    InputTrace(const Trace& other);

    static bool contains(const std::vector<std::shared_ptr<InputTrace>>& list, const InputTrace& trace);
    bool isEmptyTrace() const;

    InputTrace removeEpsilon() const { return InputTrace(static_cast<Trace>(*this).removeEpsilon()); }
    InputTrace removeLeadingEpsilons() const { return InputTrace(static_cast<Trace>(*this).removeLeadingEpsilons()); }

	/**
	Output the InputTrace to a standard output stream
	@param out The standard output stream to use
	@param trace The InputTrace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out,
                                     const InputTrace & trace);
    InputTrace& operator=(InputTrace& other);
    InputTrace& operator=(InputTrace&& other);

};

namespace std {
  template <> struct hash<InputTrace>
  {
    size_t operator()(const InputTrace& trace) const
    {
        const Trace& t = static_cast<Trace>(trace);
        return std::hash<Trace>()(t);
    }
  };
}

#endif //FSM_FSM_INPUTTRACE_H_
