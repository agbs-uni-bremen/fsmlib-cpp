/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/IOTrace.h"

using namespace std;

IOTrace::IOTrace(const InputTrace & i, const OutputTrace & o)
	: inputTrace(i), outputTrace(o)
{
    if (i.get().size() != o.get().size())
    {
        cerr << "Input trace length and output trace length differ." << endl;
        exit(EXIT_FAILURE);
    }
}

IOTrace::IOTrace(const int i, const int o, shared_ptr<FsmPresentationLayer> pl):
    inputTrace({i}, pl), outputTrace({o}, pl)
{

}

IOTrace::IOTrace(const IOTrace & ioTrace):
inputTrace(ioTrace.inputTrace), outputTrace(ioTrace.outputTrace)
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

vector<IOTrace> IOTrace::getPrefixes() const
{
    vector<IOTrace> result;
    vector<int> inputRaw = inputTrace.get();
    vector<int> outputRaw = outputTrace.get();
    if (inputRaw.size() != outputRaw.size())
    {
        cerr << "Input trace and output trace differ in size." << endl;
        exit(EXIT_FAILURE);
    }

    //vector<InputTrace> inPrefixes = inputTrace.getPrefixes();

    if (inputRaw.size() > 1) {
        vector<Trace> inPre = inputTrace.getPrefixes();
        vector<Trace> outPre = outputTrace.getPrefixes();
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

void IOTrace::append(int input, int output)
{
    inputTrace.add(input);
    outputTrace.add(output);
}

shared_ptr<IOTrace> IOTrace::getEmptyTrace(shared_ptr<FsmPresentationLayer> pl)
{
    InputTrace i = InputTrace({-1}, pl);
    OutputTrace o = OutputTrace({-1}, pl);
    return make_shared<IOTrace>(i, o);
}

ostream & operator<<(ostream & out, const IOTrace & trace)
{
	out << trace.inputTrace << "/" << trace.outputTrace;
	return out;
}

bool operator==(IOTrace const & iOTrace1, IOTrace const & iOTrace2)
{
    return iOTrace1.getInputTrace() == iOTrace2.getInputTrace() && iOTrace1.getOutputTrace() == iOTrace2.getOutputTrace();
}

IOTrace& IOTrace::operator= (IOTrace&& other)
{
    if (this != &other)
    {
        inputTrace = std::move(other.inputTrace);
        outputTrace = std::move(other.outputTrace);
    }
    return *this;
}

string IOTrace::toRttString() const {
    
    string s;
    float ts = 0.0;
    
    vector<int> inputs = inputTrace.get();
    vector<int> outputs = outputTrace.get();
    const std::shared_ptr<FsmPresentationLayer> pl = inputTrace.getPresentationLayer();
    
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


