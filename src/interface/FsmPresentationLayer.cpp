/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "interface/FsmPresentationLayer.h"
#include <algorithm>

const std::string FsmPresentationLayer::FAIL_STATE_NAME = "Fail";
const std::string FsmPresentationLayer::FAIL_OUTPUT = "fail";

FsmPresentationLayer::FsmPresentationLayer()
{

}

FsmPresentationLayer::FsmPresentationLayer(const FsmPresentationLayer& pl):
    in2String(pl.in2String), out2String(pl.out2String), state2String(pl.state2String)
{

}

FsmPresentationLayer::FsmPresentationLayer(const std::vector<std::string>& in2String, const std::vector<std::string>& out2String, const std::vector<std::string>& state2String)
	: in2String(in2String), out2String(out2String), state2String(state2String)
{

}

FsmPresentationLayer::FsmPresentationLayer(const std::string& inputs, const std::string& outputs, const std::string& states)
{
	std::string line;

	std::ifstream inputsFile(inputs);
	while (getline(inputsFile, line))
	{
		in2String.push_back(line);
	}
	inputsFile.close();


	std::ifstream outputsFile(outputs);
	while (getline(outputsFile, line))
	{
		out2String.push_back(line);
	}
	outputsFile.close();

	std::ifstream statesFile(states);
	while (getline(statesFile, line))
	{
		state2String.push_back(line);
	}
}

std::shared_ptr<FsmPresentationLayer> FsmPresentationLayer::mergeAlphabets(
        std::shared_ptr<const FsmPresentationLayer> plA,
        std::shared_ptr<const FsmPresentationLayer> plB)
{
    std::shared_ptr<FsmPresentationLayer> result = std::make_shared<FsmPresentationLayer>();
    size_t minNumberInputs = std::min(plA->in2String.size(), plB->in2String.size());
    const std::vector<std::string>& biggerInputs = plA->in2String.size() > plB->in2String.size() ? plA->in2String : plB->in2String;

    size_t minNumberOutputs = std::min(plA->out2String.size(), plB->out2String.size());
    const std::vector<std::string>& biggerOutputs = plA->out2String.size() > plB->out2String.size() ? plA->out2String : plB->out2String;

    for (size_t i = 0; i < minNumberInputs; ++i)
    {
        if (plA->in2String.at(i) != plB->in2String.at(i))
        {
            std::cerr << "Could not merge presentation layers. Input alphabets differ." << std::endl;
            exit(EXIT_FAILURE);
        }
        result->in2String.push_back(plA->in2String.at(i));
    }
    for (size_t i = minNumberInputs; i < biggerInputs.size(); ++i)
    {
        result->in2String.push_back(biggerInputs.at(i));
    }

    for (size_t i = 0; i < minNumberOutputs; ++i)
    {
        if (plA->out2String.at(i) != plB->out2String.at(i))
        {
            std::cerr << "Could not merge presentation layers. Output alphabets differ." << std::endl;
            exit(EXIT_FAILURE);
        }
        result->out2String.push_back(plA->out2String.at(i));
    }
    for (size_t i = minNumberOutputs; i < biggerOutputs.size(); ++i)
    {
        result->out2String.push_back(biggerOutputs.at(i));
    }
    return result;
}

void FsmPresentationLayer::setState2String(std::vector<std::string> state2String)
{
    this->state2String = state2String;
}

void FsmPresentationLayer::addState2String(std::string name)
{
    state2String.push_back(name);
}

void FsmPresentationLayer::removeState2String(const int index)
{
    if (index >= 0 && state2String.size() > static_cast<size_t>(index))
    {
        state2String.erase(state2String.begin() + index);
    }
}

int FsmPresentationLayer::addOut2String(std::string name)
{
    out2String.push_back(name);
    return static_cast<int>(out2String.size() - 1);
}

int FsmPresentationLayer::addIn2String(std::string name)
{
    in2String.push_back(name);
    return static_cast<int>(in2String.size() - 1);
}

void FsmPresentationLayer::truncateState2String(const int index)
{
    if (state2String.size() > static_cast<size_t>(index))
    {
        state2String.erase(state2String.begin() + index, state2String.end());
    }
}

void FsmPresentationLayer::truncateIn2String(const int index)
{
    if (in2String.size() > static_cast<size_t>(index))
    {
        in2String.erase(in2String.begin() + index, in2String.end());
    }
}

void FsmPresentationLayer::truncateOut2String(const int index)
{
    if (out2String.size() > static_cast<size_t>(index))
    {
        out2String.erase(out2String.begin() + index, out2String.end());
    }
}

std::string FsmPresentationLayer::getInId(const unsigned int id) const
{
	if (id >= in2String.size())
	{
		return std::to_string(id);
	}
	return in2String.at(id);
}

std::string FsmPresentationLayer::getOutId(const unsigned int id) const
{
	if (id >= out2String.size())
	{
		return std::to_string(id);
	}
	return out2String.at(id);
}

std::string FsmPresentationLayer::getStateId(const unsigned int id, const std::string & prefix) const
{
	if (id >= state2String.size())
	{
		if (prefix.empty())
		{
			return std::to_string(id);
		}
		return prefix + std::to_string(id);
	}
	return state2String.at(id);
}

void FsmPresentationLayer::dumpIn(std::ostream & out) const
{
	for (unsigned int i = 0; i < in2String.size(); ++ i)
	{
		if (i != 0)
		{
			out << std::endl;
		}
		out << in2String.at(i);
	}
}

void FsmPresentationLayer::dumpOut(std::ostream & out) const
{
	for (unsigned int i = 0; i < out2String.size(); ++ i)
	{
		if (i != 0)
		{
			out << std::endl;
		}
		out << out2String.at(i);
	}
}

void FsmPresentationLayer::dumpState(std::ostream & out) const
{
	for (unsigned int i = 0; i < state2String.size(); ++ i)
	{
		if (i != 0)
		{
			out << std::endl;
		}
		out << state2String.at(i);
	}
}

bool FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer> otherPresentationLayer)
{
	if (in2String.size() != otherPresentationLayer->in2String.size())
	{
		return false;
	}

	if (out2String.size() != otherPresentationLayer->out2String.size())
	{
		return false;
	}

	for (unsigned int i = 0; i < in2String.size(); ++ i)
	{
		if (in2String.at(i) != otherPresentationLayer->in2String.at(i))
		{
			return false;
		}
	}

	for (unsigned int i = 0; i < out2String.size(); ++ i)
	{
		if (out2String.at(i) != otherPresentationLayer->out2String.at(i))
		{
			return false;
		}
	}
	return true;
}


int FsmPresentationLayer::in2Num(const std::string& name) {
    
    for ( size_t i = 0; i < in2String.size(); i++ ) {
        if ( in2String[i] == name ) return (int)i;
    }
    
    return -1;
    
}

int FsmPresentationLayer::out2Num(const std::string& name) {
    
    for ( size_t i = 0; i < out2String.size(); i++ ) {
        if ( out2String[i] == name ) return (int)i;
    }
    
    return -1;
}

int FsmPresentationLayer::state2Num(const std::string& name) {
    
    for ( size_t i = 0; i < state2String.size(); i++ ) {
        if ( state2String[i] == name ) return (int)i;
    }
    
    return -1;
}

FsmPresentationLayer& FsmPresentationLayer::operator=(FsmPresentationLayer& other)
{
    if (this != &other)
    {
        in2String = other.in2String;
        out2String = other.out2String;
        state2String = other.state2String;
    }
    return *this;
}












