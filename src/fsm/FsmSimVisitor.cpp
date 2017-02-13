//
//  FsmSimVisitor.cpp
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
#include "fsm/FsmSimVisitor.h"


using namespace std;


void FsmSimVisitor::visit(Fsm& f) {
    
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
        cout << "@func void sim_" << f.getName() << "() {" << endl << endl;
        cout << "int thisState = " << f.getInitStateIdx() << ";" << endl;
        cout << endl << "while ( @rttIsRunning ) {" << endl << endl;
        cout << "rttIOPost->fsmAction = -1;" << endl << endl;
        cout << "switch ( thisState ) {" << endl;
    }
    
}

void FsmSimVisitor::visit(FsmNode& n) {
    
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

void FsmSimVisitor::visit(FsmTransition& t) {
    
    if ( finalRun ) return;
    
    cout << "case " << t.getLabel()->getInput() << ": " << endl;
    cout << "// Transition -> " << t.getTarget()->getName() << endl;
    cout << "rttIOPost->fsmAction = " << t.getLabel()->getOutput() << ";" << endl;
    cout << "thisState = " << t.getTarget()->getId() << ";" << endl;
    cout << "break;";
    
}


void FsmSimVisitor::visit(FsmLabel& x) {

}
