/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/IOTrace.h"

IOTrace::IOTrace(const InputTrace & i, const OutputTrace & o)
	: inputTrace(i), outputTrace(o)
{

}

InputTrace IOTrace::getInputTrace() const
{
	return inputTrace;
}

OutputTrace IOTrace::getOutputTrace() const
{
	return outputTrace;
}

std::ostream & operator<<(std::ostream & out, const IOTrace & trace)
{
	out << trace.inputTrace << "/" << trace.outputTrace;
	return out;
}
