/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/Trace.h"

Trace::Trace(const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: presentationLayer(presentationLayer)
{

}

Trace::Trace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: trace(trace), presentationLayer(presentationLayer)
{

}

Trace::Trace(const Trace& other):
    trace(other.trace), presentationLayer(other.presentationLayer)
{

}

void Trace::add(const int e)
{
	trace.push_back(e);
}

void Trace::append(const std::vector<int>& traceToAppend) {
    trace.reserve(trace.size() + traceToAppend.size());
    for ( size_t i = 0; i < traceToAppend.size(); i++ ) {
        trace.push_back(traceToAppend.at(i));
    }
    
}

void Trace::append(const Trace& traceToAppend) {
    append(traceToAppend.get());
}

std::vector<int> Trace::get() const
{
	return trace;
}

std::vector<int>::const_iterator Trace::cbegin() const
{
	return trace.cbegin();
}

std::vector<int>::const_iterator Trace::cend() const
{
	return trace.cend();
}

std::vector<Trace> Trace::getPrefixes() const
{
     std::vector<Trace> result;
     if (trace.size() > 1)
     {
         for (size_t i = 1; i < trace.size(); ++i)
         {
             Trace prefix = Trace(std::vector<int>(trace.begin(), trace.end() - static_cast<std::vector<int>::difference_type>(i)), presentationLayer);
            result.push_back(prefix);
         }
     }
     return result;
}

bool operator==(Trace const & trace1, Trace const & trace2)
{
	if (trace1.get().size() != trace2.get().size())
	{
		return false;
	}

	for (unsigned int i = 0; i < trace1.get().size(); ++ i)
	{
		if (trace1.get().at(i) != trace2.get().at(i))
		{
			return false;
		}
	}
	return true;
}

bool operator==(Trace const & trace1, std::vector<int> const & trace2)
{
    if (trace1.get().size() != trace2.size())
    {
        return false;
    }

    for (unsigned int i = 0; i < trace1.get().size(); ++ i)
    {
        if (trace1.get()[i] != trace2[i])
        {
            return false;
        }
    }
    return true;
}

std::ostream & operator<<(std::ostream & out, const Trace & trace)
{
	for (auto it = trace.cbegin(); it != trace.cend(); ++ it)
	{
		if (it != trace.cbegin())
		{
			out << ".";
		}
		out << *it;
	}
	return out;
}

Trace& Trace::operator=(Trace&& other)
{
    if (this != &other)
    {
        trace = std::move(other.trace);
        presentationLayer = std::move(other.presentationLayer);
    }
    return *this;
}
