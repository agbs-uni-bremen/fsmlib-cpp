/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmLabel.h"

FsmLabel::FsmLabel(const int input, const int output, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: input(input), output(output), presentationLayer(presentationLayer)
{

}

int FsmLabel::getInput() const
{
	return input;
}

int FsmLabel::getOutput() const
{
	return output;
}

bool operator==(FsmLabel const & label1, FsmLabel const & label2)
{
	return label1.getInput() == label2.getInput() && label1.getOutput() == label2.getOutput();
}

bool operator<(FsmLabel const & label1, FsmLabel const & label2)
{
	if (label1.getInput() < label2.getInput())
	{
		return true;
	}
	else if (label1.getInput() == label2.getInput() && label1.getOutput() < label2.getOutput())
	{
		return true;
	}
	return false;
}

std::ostream & operator<<(std::ostream & out, const FsmLabel & label)
{
	out << label.presentationLayer->getInId(label.input) << "/" << label.presentationLayer->getOutId(label.output);
	return out;
}