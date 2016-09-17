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

void Trace::add(const int e)
{
	trace.push_back(e);
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
