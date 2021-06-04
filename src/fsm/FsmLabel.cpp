/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmLabel.h"
#include "fsm/FsmVisitor.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"

const int FsmLabel::EPSILON = -1;

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

std::shared_ptr<IOTrace> FsmLabel::toIOTrace() const
{
    return std::make_shared<IOTrace>(input, output, presentationLayer);
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
