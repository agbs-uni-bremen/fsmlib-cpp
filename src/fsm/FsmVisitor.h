/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#ifndef FsmVisitor_h
#define FsmVisitor_h

#include "fsm/Fsm.h"

class FsmVisitor {
    
protected:
    
    FsmVisitor() { };
    
public:
    
    ~FsmVisitor() { } ;
    
    virtual void visit(Fsm& f) = 0;
    
    
};


#endif /* FsmVisitor_h */
