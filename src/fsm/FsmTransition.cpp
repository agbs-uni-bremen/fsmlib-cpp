/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmTransition.h"
#include "fsm/FsmNode.h"

using namespace std;

FsmTransition::FsmTransition(const shared_ptr<FsmNode>  source,
                             const shared_ptr<FsmNode>  target,
                             const shared_ptr<FsmLabel> label)
	: source(source), target(target), label(label)
{

}

shared_ptr<FsmNode> FsmTransition::getSource()
{
	return source;
}

void FsmTransition::setSource(shared_ptr<FsmNode> src) {
    source = src;
}

shared_ptr<FsmNode> FsmTransition::getTarget()
{
	return target;
}

void FsmTransition::setTarget(shared_ptr<FsmNode> tgt) {
    target = tgt;
}

void FsmTransition::setLabel(std::shared_ptr<FsmLabel> lbl) {
    label = lbl;
}


shared_ptr<FsmLabel> FsmTransition::getLabel()
{
	return label;
}

ostream & operator<<(ostream& out, FsmTransition& transition)
{
	out << transition.getSource()->getId() << " -> " << transition.getTarget()->getId() << "[label=\" " << *transition.label << "   \"];";
	return out;
}

