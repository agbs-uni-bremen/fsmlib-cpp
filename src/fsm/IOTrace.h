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
public:
	/**
	Create an iotrace from one inut trace and one output trace
	@param i The input trace contained into the iotrace
	@param o The output trace contained into the iotrace
	*/
    IOTrace(const InputTrace & i, const OutputTrace & o);
    IOTrace(const int i, const int o, std::shared_ptr<FsmPresentationLayer> pl);
    IOTrace(const IOTrace & ioTrace);

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

    std::vector<IOTrace> getPrefixes() const;

    size_t size() const;

    void append(IOTrace& other);

    void append(int input, int output);

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
