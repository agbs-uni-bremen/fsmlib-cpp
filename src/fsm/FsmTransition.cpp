/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/FsmTransition.h"
#include "fsm/FsmNode.h"

#include <iostream>
#include <sstream>

using namespace std;

FsmTransition::FsmTransition(const shared_ptr<FsmNode>& source,
                             const shared_ptr<FsmNode>& target,
                             const shared_ptr<FsmLabel>& label)
	: source(source), target(target), label(label)
{
    
    if ( source == nullptr ) {
        cerr << "ERROR: Constructor FsmTransition() called with null pointer as source node" << endl;
    }
    
    if ( target == nullptr ) {
        cerr << "ERROR: Constructor FsmTransition() called with null pointer as target node" << endl;
    }
    
    if ( label == nullptr ) {
        cerr << "ERROR: Constructor FsmTransition() called with null pointer as label" << endl;
    }

}

shared_ptr<FsmNode> FsmTransition::getSource()
{
    return source.lock();
}

void FsmTransition::setSource(shared_ptr<FsmNode> src) {
    source = src;
}

shared_ptr<FsmNode> FsmTransition::getTarget()
{
    return target.lock();
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
    out << transition.str();
    return out;
}

string FsmTransition::str()
{
    stringstream out;
    if (auto src = getSource())
    {
        if (auto tar = getTarget())
        {
            out << src->getId() << " -> " << tar->getId() << "[label=\"" << *label << "\"];"
                << "  //" << src->getName() << " -> " << tar->getName();
        }
        else
        {
            out << src->getId() << " -> " << "[label=\"" << *label << "\"];"
                << "  //" << src->getName() << " -> ";
        }
    }
    else if (auto tar = getTarget())
    {
        out << " -> " << tar->getId() << "[label=\"" << *label << "\"];"
            << "  //" << " -> " << tar->getName();
    }
    else
    {
        out << " -> " << "[label=\"" << *label << "\"];"
            << "  //" << " -> ";
    }
    return out.str();
}

void FsmTransition::accept(FsmVisitor &v) {
    
    v.visit(*this);
    label->accept(v);
    //target->accept(v);
    
    
}
