//
//  FsmOraVisitor.hpp
//  fsm
//
//  Created by Jan Peleska on 2017-02-12.
//
//

#ifndef FsmOraVisitor_hpp
#define FsmOraVisitor_hpp

#include <stdio.h>
#include "fsm/FsmVisitor.h"

class FsmOraVisitor : public FsmVisitor {
    
private:
    bool finalRun;
    
public:
    
    FsmOraVisitor() { finalRun = false; }
    
    virtual void visit(Fsm& f);
    virtual void visit(FsmNode& f);
    virtual void visit(FsmTransition& f);
    virtual void visit(FsmLabel& f);
    
    void setFinalRun(bool b) { finalRun = b; }
    
};

#endif /* FsmOraVisitor_hpp */
