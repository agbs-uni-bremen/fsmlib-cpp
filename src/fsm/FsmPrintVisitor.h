//
//  FsmPrintVisitor.hpp
//  fsm
//
//  Created by Jan Peleska on 2017-02-12.
//
//

#ifndef FsmPrintVisitor_hpp
#define FsmPrintVisitor_hpp

#include "fsm/FsmVisitor.h"

class FsmPrintVisitor : public FsmVisitor {
    
public:
    
    virtual void visit(Fsm& /*f*/);
    virtual void visit(FsmNode& /*f*/);
    virtual void visit(FsmTransition& /*f*/);
    virtual void visit(FsmLabel& /*f*/);
    
};

#endif /* FsmPrintVisitor_hpp */
