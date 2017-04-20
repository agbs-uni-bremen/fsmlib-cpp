//
//  FsmOraVisitor.cpp
//  fsm
//
//  Created by Jan Peleska on 2017-02-12.
//
//
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"
#include "fsm/FsmLabel.h"
#include "fsm/FsmOraVisitor.h"


using namespace std;


void FsmOraVisitor::visit(Fsm& f) {
    
    if ( finalRun ) {
        cout << endl
        << "default: break;" << endl
        << "} // switch (rttIOPost->event )  " << endl
        << "default: break;" << endl
        << "} // switch (thisState) " << endl
        << "@rttYield();" << endl
        << "}" << endl << "}";
    }
    else {
        cout << "@func void ora_" << f.getName() << "() {" << endl << endl;
        cout << "int thisState = " << f.getInitStateIdx() << ";" << endl;
        cout << endl << "while ( @rttIsRunning ) {" << endl << endl;
        cout << "rttIOPost->fsmAction = -1;" << endl << endl;
        cout << "switch ( thisState ) {" << endl;
    }
    
}

void FsmOraVisitor::visit(FsmNode& n) {
    
    if ( finalRun ) return;

    if ( getNew() ) {
        
        if ( not n.isInitial() ) {
            cout << endl << "default: break;" << endl << "}"
                 << endl << "break;" << endl;
        }
        
        cout << endl << "case " << n.getId() << ": "
        << "// State " << n.getName() << endl;
        cout << "switch ( rttIOPost->fsmEvent ) {" << endl;
        setNew(false);
    }
    else {
        
    }
    
}

void FsmOraVisitor::visit(FsmTransition& t) {
    
    if ( finalRun ) return;
    
    cout << "case " << t.getLabel()->getInput() << ": " << endl;
    cout << "// Transition -> " << t.getTarget()->getName() << endl;
    
    
    cout << "@rttAssert(rttIOPost->fsmAction == " << t.getLabel()->getOutput()
         << ");";
    
    cout << "thisState = " << t.getTarget()->getId() << ";" << endl;
    cout << "break;";
    
}


void FsmOraVisitor::visit(FsmLabel& /* x */) {

}
