/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/OutputTrace.h"

OutputTrace::OutputTrace(const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Trace(presentationLayer)
{

}

OutputTrace::OutputTrace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: Trace(trace, presentationLayer)
{

}

std::ostream & operator<<(std::ostream & out, const OutputTrace & trace)
{
	for (auto it = trace.cbegin(); it != trace.cend(); ++ it)
	{
		if (it != trace.cbegin())
		{
			out << ".";
		}
		out << trace.presentationLayer->getOutId(*it);
	}
	return out;
}