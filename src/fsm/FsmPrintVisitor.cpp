//
//  FsmPrintVisitor.cpp
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
#include "fsm/FsmPrintVisitor.h"


using namespace std;


void FsmPrintVisitor::visit(Fsm& f) {
    cout << "FSM " << f.getName() << endl;
}

void FsmPrintVisitor::visit(FsmNode& n) {

    if ( getNew() ) cout << endl;
    cout << n.getName();
    setNew(false);
    
}

void FsmPrintVisitor::visit(FsmTransition& t) {
    cout << endl << "\t\tTransition " << t.getSource()->getName();
}


void FsmPrintVisitor::visit(FsmLabel& x) {
    cout << " --- "
    << x.getInput() << "/"
    << x.getOutput() << " ---> ";

}
