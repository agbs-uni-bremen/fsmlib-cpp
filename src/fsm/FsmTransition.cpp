/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmTransition.h"
#include "fsm/FsmNode.h"

FsmTransition::FsmTransition(const std::shared_ptr<FsmNode> source, const std::shared_ptr<FsmNode> target, const FsmLabel & label)
	: source(source), target(target), label(label)
{

}

std::shared_ptr<FsmNode> FsmTransition::getSource() const
{
	return source.lock();
}

std::shared_ptr<FsmNode> FsmTransition::getTarget() const
{
	return target.lock();
}

FsmLabel FsmTransition::getLabel() const
{
	return label;
}

std::ostream & operator<<(std::ostream & out, const FsmTransition & transition)
{
	out << transition.getSource()->getId() << " -> " << transition.getTarget()->getId() << "[label=\" " << transition.label << "   \"];";
	return out;
}
