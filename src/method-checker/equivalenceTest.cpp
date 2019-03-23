#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"

using namespace std;

shared_ptr<FsmPresentationLayer> createPresentationLayer(const size_t maxInput, const size_t refMin, const size_t maxOutput)
{
    vector<string> in2string = vector<string>();
    for (size_t i = 0; i <= maxInput; i++)
    {
        in2string.push_back(to_string(i));
    }
    vector<string> out2string = vector<string>();
    for (size_t i = 0; i <= maxOutput; i++)
    {
        out2string.push_back(to_string(i));
    }
    vector<string>state2string = vector<string>();
    for (size_t i = 0; i < refMin; i++)
    {
        state2string.push_back(to_string(i));
    }
    shared_ptr<FsmPresentationLayer> pl{
          new FsmPresentationLayer(
            in2string,
            out2string,
            state2string)};

    return pl;
}

void runTestClass(shared_ptr<Dfsm> dfsm ,int numOutputFaults,int numTransitionFaults,
                  string faultClass,vector<IOTrace>& tMethodIOTraces,
                  vector<IOTrace>& wMethodIOTraces,vector<IOTrace>& wpMethodIOTraces,
                  vector<IOTrace>& hMethodIOTraces,vector<IOTrace>& dMethodIOTraces) {

    int count_equals = 0;
    int count_tMethod_pass = 0;
    int count_wMethod_pass = 0;
    int count_wpMethod_pass = 0;
    int count_hMethod_pass = 0;
    int count_dMethod_pass = 0;
    for(int i=0; i < 100;i++) {
        shared_ptr<Dfsm> mutant = dfsm->createMutant("Mutant",numOutputFaults,numTransitionFaults);

        bool result = dfsm->equivalenceCheck(*mutant);
        if(result) count_equals++;

        result = true;
        for(IOTrace io: tMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_tMethod_pass++;

        result = true;
        for(IOTrace io: wMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_wMethod_pass++;

        result = true;
        for(IOTrace io: wpMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_wpMethod_pass++;

        result = true;
        for(IOTrace io: hMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_hMethod_pass++;

        result = true;
        for(IOTrace io: dMethodIOTraces) {
          result &= mutant->pass(io);
        }
        if(result) count_dMethod_pass++;
    }
    cout << "Compare DFSM with Mutants of Class " << faultClass
       << "(numOutputFaults = " << numOutputFaults << ", "
       << "numTransitionFaults = " << numTransitionFaults << ")" << endl
       << "Number of equivalent DFSM : "  << count_equals << endl
       << "Number of T-Method pass : "  << count_tMethod_pass << endl
       << "Number of W-Method pass : "  << count_wMethod_pass << endl
       << "Number of Wp-Method pass : "  << count_wpMethod_pass << endl
       << "Number of D-Method pass : "  << count_dMethod_pass << endl
       << "Number of H-Method pass : "  << count_hMethod_pass << endl;
}

int main(int argc, char* argv[])
{
    while(true) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(6, 30, 6);

        shared_ptr<Dfsm> temp = make_shared<Dfsm>("Dfsm", 30, 6, 6, pl);
        shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(temp->minimise());

        IOListContainer e = dfsm->dMethodOnMinimisedDfsm(0);
        if(e.size() == 0) continue;
        vector<IOTrace> dMethodIOTraces;
        for (vector<int> v: *e.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            dMethodIOTraces.push_back(io);
        }
        cout << "dMethod size: " << dMethodIOTraces.size() << endl;

        IOListContainer a = dfsm->tMethod();
        vector<IOTrace> tMethodIOTraces;
        for (vector<int> v: *a.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            tMethodIOTraces.push_back(io);
        }
        cout << "TMethod size: " << tMethodIOTraces.size() << endl;

        IOListContainer b = dfsm->wMethod(0);
        vector<IOTrace> wMethodIOTraces;
        for (vector<int> v: *b.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            wMethodIOTraces.push_back(io);
        }
        cout << "WMethod size: " << wMethodIOTraces.size() << endl;

        IOListContainer c = dfsm->wpMethod(0);
        vector<IOTrace> wpMethodIOTraces;
        for (vector<int> v: *c.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            wpMethodIOTraces.push_back(io);
        }
        cout << "WpMethod size: " << wpMethodIOTraces.size() << endl;

        IOListContainer d = dfsm->hMethodOnMinimisedDfsm(0);
        vector<IOTrace> hMethodIOTraces;
        for (vector<int> v: *d.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            hMethodIOTraces.push_back(io);
        }
        cout << "hMethod size: " << hMethodIOTraces.size() << endl;

        runTestClass(dfsm, 1, 0, "Class 1", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 0, 1, "Class 2", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 2, 0, "Class 3", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 0, 2, "Class 4", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 1, 1, "Class 5", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 2, 1, "Class 6", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 1, 2, "Class 7", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 2, 2, "Class 8", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 3, 2, "Class 9", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        runTestClass(dfsm, 2, 3, "Class 10", tMethodIOTraces, wMethodIOTraces, wpMethodIOTraces, hMethodIOTraces,dMethodIOTraces);
        return 0;
    }
}
