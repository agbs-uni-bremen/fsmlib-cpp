/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmLabel.h"

const int FsmLabel::EPSILON = -1;
const int FsmLabel::ERROR_OUTPUT = -2;
const int FsmLabel::UNDEFINED_OUTPUT = -1000;

FsmLabel::FsmLabel(const int input,
                   const int output,
                   const std::shared_ptr<FsmPresentationLayer>& presentationLayer)
	: input(input), output(output), presentationLayer(presentationLayer)
{

}

FsmLabel::FsmLabel(const FsmLabel& other) :
   input(other.input), output(other.output), presentationLayer(other.presentationLayer)
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


void FsmLabel::accept(FsmVisitor &v) {
    
    v.visit(*this);
}
