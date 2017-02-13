//
//  FsmSimVisitor.hpp
//  fsm
//
//  Created by Jan Peleska on 2017-02-12.
//
//

#ifndef FsmSimVisitor_hpp
#define FsmSimVisitor_hpp

#include <stdio.h>
#include "fsm/FsmVisitor.h"

class FsmSimVisitor : public FsmVisitor {
    
private:
    bool finalRun;
    
public:
    
    FsmSimVisitor() { finalRun = false; }
    
    virtual void visit(Fsm& f);
    virtual void visit(FsmNode& f);
    virtual void visit(FsmTransition& f);
    virtual void visit(FsmLabel& f);
    
    void setFinalRun(bool b) { finalRun = b; }
    
};

#endif /* FsmSimVisitor_hpp */
