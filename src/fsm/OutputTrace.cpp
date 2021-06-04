/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include <ostream>

#include "interface/FsmPresentationLayer.h"
#include "fsm/OutputTrace.h"

OutputTrace::OutputTrace(const std::shared_ptr<FsmPresentationLayer>& presentationLayer)
	: Trace(presentationLayer)
{

}

OutputTrace::OutputTrace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer const>& presentationLayer)
	: Trace(trace, presentationLayer)
{

}

OutputTrace::OutputTrace(const OutputTrace& other, size_t n, bool defaultToEmpty):
    Trace(other, n, defaultToEmpty)
{

}

OutputTrace::OutputTrace(const OutputTrace& other):
    Trace(other)
{

}

OutputTrace::OutputTrace(const Trace& other):
    Trace(other)
{

}

bool OutputTrace::contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const OutputTrace& trace)
{
    for (std::shared_ptr<OutputTrace> t : list)
    {
        if (*t == trace)
        {
            return true;
        }
    }
    return false;
}

bool OutputTrace::contains(const std::vector<std::shared_ptr<OutputTrace>>& list, const int output)
{
    for (std::shared_ptr<OutputTrace> t : list)
    {
        if (t->contains(output))
        {
            return true;
        }
    }
    return false;
}


bool OutputTrace::contains(const int output) const
{
    return this->Trace::contains(output);
}


OutputTrace& OutputTrace::operator=(OutputTrace& other)
{
    if (this != &other)
    {
        trace = other.trace;
    }
    return *this;
}

OutputTrace& OutputTrace::operator=(OutputTrace&& other)
{
    if (this != &other)
    {
        Trace::operator=(std::move(other));
    }
    return *this;
}

std::ostream & operator<<(std::ostream & out, const OutputTrace & trace)
{
	for (auto it = trace.cbegin(); it != trace.cend(); ++ it)
	{
		if (it != trace.cbegin())
		{
			out << ".";
		}
        if (*it == -1)
        {
            out << "ε";
        }
        else
        {
            out << trace.presentationLayer->getOutId(*it);
        }
	}
	return out;
}
