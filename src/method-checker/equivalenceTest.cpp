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

shared_ptr<vector<IOTrace>> tMethodIOTraces = make_shared<vector<IOTrace>>();
shared_ptr<vector<IOTrace>> wMethodIOTraces = make_shared<vector<IOTrace>>();
shared_ptr<vector<IOTrace>> wpMethodIOTraces = make_shared<vector<IOTrace>>();
shared_ptr<vector<IOTrace>> hMethodIOTraces = make_shared<vector<IOTrace>>();
shared_ptr<vector<IOTrace>> dMethodIOTraces = make_shared<vector<IOTrace>>();
shared_ptr<vector<IOTrace>> adsMethodIOTraces = make_shared<vector<IOTrace>>();

void runTestClass(shared_ptr<Dfsm> dfsm ,int numOutputFaults,int numTransitionFaults,
                  string faultClass) {

    int count_equals = 0;
    int count_tMethod_pass = 0;
    int count_wMethod_pass = 0;
    int count_wpMethod_pass = 0;
    int count_hMethod_pass = 0;
    int count_dMethod_pass = 0;
    int count_adsMethod_pass = 0;

    int numberOfMutants = 100;

    for(int i=0; i < numberOfMutants;i++) {
        shared_ptr<Dfsm> mutant = dfsm->createMutant("Mutant",numOutputFaults,numTransitionFaults);

        bool result = dfsm->equivalenceCheck(*mutant);
        if(result) count_equals++;

        result = true;
        for(IOTrace io: *tMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_tMethod_pass++;

        result = true;
        for(IOTrace io: *wMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_wMethod_pass++;

        result = true;
        for(IOTrace io: *wpMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_wpMethod_pass++;

        result = true;
        for(IOTrace io: *hMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_hMethod_pass++;

        result = true;
        for(IOTrace io: *dMethodIOTraces) {
          result &= mutant->pass(io);
        }
        if(result) count_dMethod_pass++;

        result = true;
        for(IOTrace io: *adsMethodIOTraces) {
            result &= mutant->pass(io);
        }
        if(result) count_adsMethod_pass++;
    }
    cout << "Compare DFSM with Mutants of Class " << faultClass
       << "(numOutputFaults = " << numOutputFaults << ", "
       << "numTransitionFaults = " << numTransitionFaults << ")" << endl
       << "Number of equivalent DFSM : "  << count_equals << endl
       << "Number of T-Method pass : "  << count_tMethod_pass << endl
       << "Number of W-Method pass : "  << count_wMethod_pass << endl
       << "Number of Wp-Method pass : "  << count_wpMethod_pass << endl
       << "Number of D-Method pass : "  << count_dMethod_pass << endl
       << "Number of D-Method(ADS) pass : "  << count_adsMethod_pass << endl
       << "Number of H-Method pass : "  << count_hMethod_pass << endl;
}

int main(int argc, char* argv[])
{
    while(true) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(6, 30, 6);

        shared_ptr<Dfsm> temp = make_shared<Dfsm>("Dfsm", 30, 6, 6, pl);
        shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(temp->minimise());

        IOListContainer e = dfsm->dMethodOnMinimisedDfsm(0, false);
        if(e.size() == 0) continue;
        for (vector<int> v: *e.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            dMethodIOTraces->push_back(io);
        }
        cout << "dMethod size: " << e.getFlatSize() << endl;

        IOListContainer f = dfsm->dMethodOnMinimisedDfsm(0, true);
        for (vector<int> v: *e.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            adsMethodIOTraces->push_back(io);
        }
        cout << "dMethod(ADS) size: " << f.getFlatSize() << endl;

        IOListContainer a = dfsm->tMethod();
        for (vector<int> v: *a.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            tMethodIOTraces->push_back(io);
        }
        cout << "TMethod size: " << a.getFlatSize() << endl;

        IOListContainer b = dfsm->wMethod(0);
        for (vector<int> v: *b.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            wMethodIOTraces->push_back(io);
        }
        cout << "WMethod size: " << b.getFlatSize() << endl;

        IOListContainer c = dfsm->wpMethod(0);
        for (vector<int> v: *c.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            wpMethodIOTraces->push_back(io);
        }
        cout << "WpMethod size: " << c.getFlatSize() << endl;

        IOListContainer d = dfsm->hMethodOnMinimisedDfsm(0);
        for (vector<int> v: *d.getIOLists()) {
            InputTrace i(v, pl);
            IOTrace io = dfsm->applyDet(i);
            hMethodIOTraces->push_back(io);
        }
        cout << "hMethod size: " << d.getFlatSize() << endl;

        runTestClass(dfsm, 1, 0, "Class 1");
        runTestClass(dfsm, 0, 1, "Class 2");
        runTestClass(dfsm, 2, 0, "Class 3");
        runTestClass(dfsm, 0, 2, "Class 4");
        runTestClass(dfsm, 1, 1, "Class 5");
        runTestClass(dfsm, 2, 1, "Class 6");
        runTestClass(dfsm, 1, 2, "Class 7");
        runTestClass(dfsm, 2, 2, "Class 8");
        runTestClass(dfsm, 3, 2, "Class 9");
        runTestClass(dfsm, 2, 3, "Class 10");
        return 0;
    }
}
