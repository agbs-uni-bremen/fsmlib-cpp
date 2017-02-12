//
//  FsmPrintVisitor.hpp
//  fsm
//
//  Created by Jan Peleska on 2017-02-12.
//
//

#ifndef FsmPrintVisitor_hpp
#define FsmPrintVisitor_hpp

#include <stdio.h>
#include "fsm/FsmVisitor.h"

class FsmPrintVisitor : public FsmVisitor {
    
    
    virtual void visit(Fsm& f);
    
    
    
};

#endif /* FsmPrintVisitor_hpp */
