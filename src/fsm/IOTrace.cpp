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

}

InputTrace IOTrace::getInputTrace() const
{
	return inputTrace;
}

OutputTrace IOTrace::getOutputTrace() const
{
	return outputTrace;
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


