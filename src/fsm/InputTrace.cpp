/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/InputTrace.h"
#include "fsm/FsmLabel.h"

InputTrace::InputTrace(const std::shared_ptr<FsmPresentationLayer>& presentationLayer)
	: Trace(presentationLayer)
{

}

InputTrace::InputTrace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer const>& presentationLayer)
	: Trace(trace, presentationLayer)
{

}

InputTrace::InputTrace(int x,
           const std::shared_ptr<FsmPresentationLayer>& presentationLayer)
    :Trace(std::vector<int>({x}), presentationLayer)
{

}

InputTrace::InputTrace(const InputTrace& other, size_t n, bool defaultToEmpty):
    Trace(other, n, defaultToEmpty)
{

}

InputTrace::InputTrace(const InputTrace& other):
    Trace(other)
{

}

InputTrace::InputTrace(const Trace& other):
    Trace(other)
{

}

bool InputTrace::contains(const InputTraceSet& list,
                          const std::shared_ptr<InputTrace>& trace)
{
    return list.find(trace) != list.end();
}

bool InputTrace::isEmptyTrace() const
{
    return trace.size() == 1 && trace.at(0) == FsmLabel::EPSILON;
}

std::ostream & operator<<(std::ostream & out, const InputTrace & trace)
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
            out << trace.presentationLayer->getInId(*it);
        }
	}
	return out;
}

InputTrace& InputTrace::operator=(InputTrace& other)
{
    if (this != &other)
    {
        trace = other.trace;
    }
    return *this;
}

InputTrace& InputTrace::operator=(InputTrace&& other)
{
    if (this != &other)
    {
        Trace::operator=(std::move(other));
    }
    return *this;
}
