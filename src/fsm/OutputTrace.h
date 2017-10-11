/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_OUTPUTTRACE_H_
#define FSM_FSM_OUTPUTTRACE_H_

#include <iostream>
#include <vector>

#include "fsm/Trace.h"
#include "interface/FsmPresentationLayer.h"

class OutputTrace : public Trace
{
public:
	/**
	Create an empty output trace, with only one presentation layer
	@param presentationLayer The presentation layer used by the trace
	*/
	OutputTrace(const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Create an output trace
	@param trace The trace itself, represented by a list of int
	@param presentationLayer The presentation layer used by the trace
	*/
    OutputTrace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer> presentationLayer);

    static bool contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const OutputTrace& trace);
    static bool contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const int output);
    bool contains(const int output) const;

	/**
	Output the OutputTrace to a standard output stream
	@param out The standard output stream to use
	@param trace The OutputTrace to print
	@return The standard output stream used, to allow user to cascade <<
	*/
    friend std::ostream & operator<<(std::ostream & out, const OutputTrace & trace);
};
#endif //FSM_FSM_OUTPUTTRACE_H_
