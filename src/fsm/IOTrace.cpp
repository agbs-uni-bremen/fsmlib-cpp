/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/IOTrace.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmLabel.h"
#include "interface/FsmPresentationLayer.h"

#include <iostream>
#include <sstream>

using namespace std;

IOTrace::IOTrace(const InputTrace & i, const OutputTrace & o, std::shared_ptr<FsmNode> targetNode)
    : inputTrace(i), outputTrace(o), targetNode(targetNode)
{
    if (i.get().size() != o.get().size())
    {
        cerr << "Input trace length and output trace length differ." << endl;
        exit(EXIT_FAILURE);
    }
}

IOTrace::IOTrace(const int i, const int o, shared_ptr<FsmPresentationLayer const> pl):
    inputTrace({i}, pl), outputTrace({o}, pl)
{

}

IOTrace::IOTrace(const Trace& i, const Trace& o):
    inputTrace(i), outputTrace(o)
{

}

IOTrace::IOTrace(const int i, const int o, std::shared_ptr<FsmNode> targetNode, shared_ptr<FsmPresentationLayer> pl):
    IOTrace(i, o, pl)
{
    this->targetNode = targetNode;
}

IOTrace::IOTrace(const IOTrace & ioTrace):
inputTrace(ioTrace.inputTrace), outputTrace(ioTrace.outputTrace), targetNode(ioTrace.targetNode)
{

}

IOTrace::IOTrace(const IOTrace& ioTrace, const IOTrace& append, bool prepend):
    inputTrace(ioTrace.inputTrace), outputTrace(ioTrace.outputTrace), targetNode(ioTrace.targetNode)
{
    if (prepend)
    {
        inputTrace.prepend(append.getInputTrace());
        outputTrace.prepend(append.getOutputTrace());
    }
    else
    {
        inputTrace.append(append.getInputTrace());
        outputTrace.append(append.getOutputTrace());
    }
}

IOTrace::IOTrace(const IOTrace & ioTrace, int n, std::shared_ptr<FsmNode> targetNode):
    inputTrace(ioTrace.getInputTrace()), outputTrace(ioTrace.getOutputTrace()), targetNode(targetNode)
{
    inputTrace.removeElements(n);
    outputTrace.removeElements(n);
}

IOTrace::IOTrace(std::shared_ptr<FsmPresentationLayer> pl):
    inputTrace(vector<int>(), pl), outputTrace(vector<int>(), pl)
{

}

IOTrace::IOTrace(const IOTrace& other, size_t n, bool defaultToEmpty):
    inputTrace(other.inputTrace, n, defaultToEmpty), outputTrace(other.outputTrace, n, defaultToEmpty)
{
    if ( inputTrace.size() != outputTrace.size() ) {
        cerr << "IOTrace(): Non-matching sizes of input and output traces."
        << endl << "In : " << inputTrace << endl << "Out: " << outputTrace <<
        endl;
    }
}

InputTrace IOTrace::getInputTrace() const
{
	return inputTrace;
}

OutputTrace IOTrace::getOutputTrace() const
{
	return outputTrace;
}

vector<IOTrace> IOTrace::getPrefixes(bool proper) const
{
    vector<IOTrace> result;
    vector<int> inputRaw = inputTrace.get();
    vector<int> outputRaw = outputTrace.get();
    if (inputRaw.size() != outputRaw.size())
    {
        cerr << "Input trace and output trace differ in size." << endl;
        exit(EXIT_FAILURE);
    }

    size_t prefixIndex = proper ? 1 : 0;
    if (inputRaw.size() > prefixIndex) {
        vector<Trace> inPre = inputTrace.getPrefixes(proper);
        vector<Trace> outPre = outputTrace.getPrefixes(proper);
        for (size_t i = 0; i < inPre.size(); ++i)
        {
            InputTrace ip = InputTrace(inPre.at(i).get(), inPre.at(i).getPresentationLayer());
            OutputTrace op = OutputTrace(outPre.at(i).get(), outPre.at(i).getPresentationLayer());
            IOTrace prefix = IOTrace(ip, op);
            result.push_back(prefix);
        }

    }
    return result;
}


size_t IOTrace::size() const
{
    return inputTrace.get().size();
}

void IOTrace::append(IOTrace& other)
{
    inputTrace.append(other.getInputTrace());
    outputTrace.append(other.getOutputTrace());
}

void IOTrace::prepend(IOTrace& other)
{
    inputTrace.prepend(other.getInputTrace());
    outputTrace.prepend(other.getOutputTrace());
}

void IOTrace::append(int input, int output)
{
    inputTrace.add(input);
    outputTrace.add(output);
}

bool IOTrace::isEmptyTrace() const
{
    return inputTrace.isEmptyTrace() && outputTrace.isEmptyTrace();
}

bool IOTrace::isPrefix(const IOTrace& other, bool proper, bool allowEmpty) const
{
    return inputTrace.isPrefix(other.inputTrace, proper, allowEmpty) &&
            outputTrace.isPrefix(other.outputTrace, proper, allowEmpty);
}

bool IOTrace::isSuffix(const IOTrace& other) const
{
    return inputTrace.isSuffix(other.inputTrace) &&
            outputTrace.isSuffix(other.outputTrace);
}

bool IOTrace::isPrefixOf(const IOTrace& other) const
{
    return inputTrace.isPrefixOf(other.inputTrace) &&
            outputTrace.isPrefixOf(other.outputTrace);
}

IOTrace IOTrace::getSuffix(const IOTrace& prefix) const
{
    if (!isPrefix(prefix))
    {
        stringstream ss;
ss << "The given prefix is not a prefix of this trace.";
std::cerr << ss.str();
throw ss.str();
    }
    const Trace& in = inputTrace.getSuffix(prefix.getInputTrace());
    const Trace& out = outputTrace.getSuffix(prefix.getOutputTrace());

    if (in.get().size() == 0 && out.get().size() == 0)
    {
        return IOTrace(FsmLabel::EPSILON, FsmLabel::EPSILON, inputTrace.getPresentationLayer());
    }
    else{
        return IOTrace(in, out);
    }
}

shared_ptr<IOTrace> IOTrace::getEmptyTrace(shared_ptr<FsmPresentationLayer> pl)
{
    InputTrace i(FsmLabel::EPSILON, pl);
    OutputTrace o({FsmLabel::EPSILON}, pl);
    return make_shared<IOTrace>(i, o);
}

ostream & operator<<(ostream & out, const IOTrace & trace)
{
    out << trace.inputTrace << "/" << trace.outputTrace;
    if (!trace.targetNode.expired())
    {
        out << " -> " << trace.targetNode.lock()->getName();
    }
	return out;
}

bool operator==(IOTrace const & iOTrace1, IOTrace const & iOTrace2)
{
    return iOTrace1.getInputTrace() == iOTrace2.getInputTrace() && iOTrace1.getOutputTrace() == iOTrace2.getOutputTrace();
}

bool operator<=(IOTrace const & trace1, IOTrace const & trace2)
{
    if (trace1.size() < trace2.size())
    {
        return true;
    }
    else if (trace1.size() > trace2.size())
    {
        return false;
    }
    else
    {
        vector<int> comp1 = trace1.getInputTrace().get();
        vector<int> comp2 = trace2.getInputTrace().get();
        for (size_t i = 0; i < comp1.size(); ++i)
        {
            if (comp1.at(i) < comp2.at(i))
            {
                return true;
            }
            else if (comp1.at(i) > comp2.at(i))
            {
                return false;
            }
        }
        comp1 = trace1.getOutputTrace().get();
        comp2 = trace2.getOutputTrace().get();
        for (size_t i = 0; i < comp1.size(); ++i)
        {
            if (comp1.at(i) < comp2.at(i))
            {
                return true;
            }
            else if (comp1.at(i) > comp2.at(i))
            {
                return false;
            }
        }
    }
    return false;
}

IOTrace& IOTrace::operator= (IOTrace& other)
{
    if (this != &other)
    {
        inputTrace = other.inputTrace;
        outputTrace = other.outputTrace;
        targetNode = other.targetNode;
    }
    return *this;
}

IOTrace& IOTrace::operator= (IOTrace&& other)
{
    if (this != &other)
    {
        inputTrace = std::move(other.inputTrace);
        outputTrace = std::move(other.outputTrace);
        targetNode = other.targetNode;
    }
    return *this;
}

string IOTrace::toRttString() const {
    
    string s;
    float ts = 0.0;
    
    vector<int> inputs = inputTrace.get();
    vector<int> outputs = outputTrace.get();
    std::shared_ptr<FsmPresentationLayer const> pl = inputTrace.getPresentationLayer();
    
    for ( size_t i = 0; i < inputs.size(); i++ ) {
        ostringstream ossIn;
        ostringstream ossOut;

        string xStr = pl->getInId(inputs[i]);
        string yStr = pl->getOutId(outputs[i]);
        
        ossIn << ts << ";" << xStr << ";0" << endl;
        s += ossIn.str();
        
        ts += 0.100;
        ossOut << ts << ";" << yStr  << ";0" << endl;
        s += ossOut.str();
        
    }
    
    return s;
}

bool IOTrace::operator<(IOTrace const &other) const {
    if(inputTrace < other.inputTrace) {
        return true;
    }
    if(other.inputTrace < inputTrace) {
        return false;
    }
    return outputTrace < other.outputTrace;
}

