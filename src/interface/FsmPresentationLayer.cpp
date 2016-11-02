/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "interface/FsmPresentationLayer.h"

FsmPresentationLayer::FsmPresentationLayer()
{

}

FsmPresentationLayer::FsmPresentationLayer(const std::vector<std::string>& in2String, const std::vector<std::string>& out2String, const std::vector<std::string>& state2String)
	: in2String(in2String), out2String(out2String), state2String(state2String)
{

}

FsmPresentationLayer::FsmPresentationLayer(const std::string& inputs, const std::string& outputs, const std::string& states)
{
	/*Quite different from the Java version. It is because the way of reading file is very different*/
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
