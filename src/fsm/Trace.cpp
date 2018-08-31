/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include <sstream>

#include "fsm/Trace.h"
#include "fsm/FsmLabel.h"
#include "utils/Logger.hpp"

Trace::Trace(const std::shared_ptr<FsmPresentationLayer const>& presentationLayer)
	: presentationLayer(presentationLayer)
{

}

Trace::Trace(const std::vector<int>& trace, const std::shared_ptr<FsmPresentationLayer const>& presentationLayer)
	: trace(trace), presentationLayer(presentationLayer)
{

}

Trace::Trace(const Trace& other):
    trace(other.trace), presentationLayer(other.presentationLayer)
{

}

Trace::Trace(const std::vector<int>::const_iterator& begin,
             const std::vector<int>::const_iterator& end,
             const std::shared_ptr<FsmPresentationLayer const>& presentationLayer)
    : presentationLayer(presentationLayer)
{
    trace = std::vector<int>(begin, end);
}

Trace::Trace(const Trace& other, size_t n, bool defaultToEmpty):
    presentationLayer(other.presentationLayer)
{
    std::vector<int> otherTrace = other.get();
    if (otherTrace.size() == 0)
    {
        if (defaultToEmpty)
        {
            trace = {FsmLabel::EPSILON};
        }
    }
    else
    {
        if (n >= otherTrace.size())
        {
            if (defaultToEmpty)
            {
                trace = {FsmLabel::EPSILON};
            }
        }
        else
        {
            trace = std::vector<int>(otherTrace.begin() + static_cast<std::vector<int>::difference_type>(n), otherTrace.end());
        }
    }
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

void Trace::prepend(const std::vector<int>& traceToPrepend) {
    trace.insert(trace.begin(), traceToPrepend.begin(), traceToPrepend.end());
}

void Trace::append(const Trace& traceToAppend) {
    append(traceToAppend.get());
}

void Trace::prepend(const Trace& traceToPrepend) {
    prepend(traceToPrepend.get());
}

Trace Trace::removeEpsilon() const
{
    Trace result = Trace(presentationLayer);
    for (int symbol : trace)
    {
        if (symbol != FsmLabel::EPSILON)
        {
            result.add(symbol);
        }
    }
    return result;
}

Trace Trace::removeLeadingEpsilons() const
{
    Trace result = Trace(presentationLayer);
    bool foundSymbol = false;
    for (int symbol : trace)
    {
        if (foundSymbol || symbol != FsmLabel::EPSILON)
        {
            result.add(symbol);
            foundSymbol = true;
        }
    }
    return result;
}

bool Trace::isEmptyTrace() const
{
    for (int symbol : trace)
    {
        if (symbol != FsmLabel::EPSILON)
        {
            return false;
        }
    }
    return true;
}

bool Trace::isPrefix(const Trace& other, bool proper, bool allowEmpty) const
{
    if (allowEmpty && other.isEmptyTrace())
    {
        return true;
    }

    if (!allowEmpty && other.isEmptyTrace())
    {
        return false;
    }

    const Trace& thisCopy = removeEpsilon();
    const Trace& otherCopy = other.removeEpsilon();

    if (proper && thisCopy == otherCopy)
    {
        return false;
    }

    if (otherCopy.get().size() > thisCopy.get().size())
    {
        return false;
    }
    for (size_t i = 0; i < otherCopy.get().size(); ++i)
    {
        if (otherCopy.get().at(i) != thisCopy.get().at(i))
        {
            return false;
        }
    }
    return true;
}

bool Trace::isSuffix(const Trace& other) const
{
    if (other.get().size() > trace.size())
    {
        return false;
    }
    const int thisMax = static_cast<int>(get().size()) - 1;
    const int otherMax = static_cast<int>(other.get().size()) - 1;
    for (int i = 0; i <= otherMax ; ++i)
    {
        if (other.get().at(static_cast<size_t>(otherMax - i)) != trace.at(static_cast<size_t>(thisMax - i)))
        {
            return false;
        }
    }
    return true;
}

bool Trace::isPrefixOf(const Trace& other) const
{
 return other.isPrefix(*this);
}

std::shared_ptr<Trace> Trace::getSuffix(size_t n, bool defaultToEmpty) const
{
    const std::vector<int>::const_iterator& begin = trace.begin() + static_cast<std::vector<int>::difference_type>(n);
    const std::vector<int>::const_iterator& end = trace.end();
    Trace t = Trace(begin, end, presentationLayer);
    std::shared_ptr<Trace> trace = std::make_shared<Trace>(t);
    if (defaultToEmpty && trace->size() == 0)
    {
        trace->add(FsmLabel::EPSILON);
    }
    return trace;
}

std::shared_ptr<const Trace> Trace::getPrefix(size_t n, bool defaultToEmpty) const
{
    const std::vector<int>::const_iterator& begin = trace.begin();
    const std::vector<int>::const_iterator& end = trace.begin() + static_cast<std::vector<int>::difference_type>(n);
    Trace t = Trace(begin, end, presentationLayer);
    std::shared_ptr<Trace> trace = std::make_shared<Trace>(t);
    if (defaultToEmpty && trace->size() == 0)
    {
        trace->add(FsmLabel::EPSILON);
    }
    return trace;
}

Trace Trace::getSuffix(const Trace& prefix) const
{
    if (!isPrefix(prefix))
    {
        std::stringstream ss;
        ss << "The given prefix is not a prefix of this trace.";
        LOG("FATAL") << ss.str() << std::endl;
        throw ss.str();
    }

    const Trace& thisCopy = removeEpsilon();
    const Trace& prefixCopy = prefix.removeEpsilon();

    std::vector<int> newTrace(thisCopy.cbegin() + static_cast<std::vector<int>::difference_type>(prefixCopy.get().size()), thisCopy.cend());
    return Trace(newTrace, presentationLayer);
}

void Trace::removeElements(int n)
{
    if (n > 0)
    {
        trace.erase(trace.begin(), trace.begin() + n);
    }
    else if (n < 0)
    {
        trace.erase(trace.end() + n, trace.end());
    }
}

std::vector<int> Trace::get() const
{
	return trace;
}

size_t Trace::size() const
{
    Trace tmp = this->removeEpsilon();
    return tmp.trace.size();
}

std::vector<int>::const_iterator Trace::cbegin() const
{
	return trace.cbegin();
}

std::vector<int>::const_iterator Trace::cend() const
{
	return trace.cend();
}

std::vector<Trace> Trace::getPrefixes(bool proper) const
{
     std::vector<Trace> result;
     size_t prefixIndex = proper ? 1 : 0;
     if (trace.size() > prefixIndex)
     {
         for (size_t i = prefixIndex; i < trace.size(); ++i)
         {
             Trace prefix = Trace(std::vector<int>(trace.begin(), trace.end() - static_cast<std::vector<int>::difference_type>(i)), presentationLayer);
            result.push_back(prefix);
         }
     }
     return result;
}

bool Trace::contains(const std::vector<std::shared_ptr<Trace>>& list, const Trace& trace)
{
    for (std::shared_ptr<Trace> t : list)
    {
        if (*t == trace)
        {
            return true;
        }
    }
    return false;
}

bool Trace::contains(const int io) const
{
    for (int i : trace)
    {
        if (i == io)
        {
            return true;
        }
    }
    return false;
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

bool Trace::operator<(Trace const &other) const {
    return trace < other.trace;
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

Trace& Trace::operator=(Trace& other)
{
    if (this != &other)
    {
        trace = other.trace;
        presentationLayer = other.presentationLayer;
    }
    return *this;
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
