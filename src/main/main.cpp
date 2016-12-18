#include <iostream>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <interface/FsmPresentationLayer.h>
#include <fsm/Dfsm.h>
#include <fsm/Fsm.h>
#include <fsm/IOTrace.h>
#include <trees/IOListContainer.h>
#include <trees/OutputTree.h>
#include <trees/TestSuite.h>


using namespace std;

void assert(string tc, bool verdict, string comment = "") {
    
    string sVerdict = (verdict) ? "PASS" : "FAIL";
    cout << sVerdict << ": " << tc
    << " : "
    << " assertion `" << comment << "'" <<  endl;
    
}

void test1() {
    
    cout << "TC-DFSM-0001 Show that Dfsm.applyDet() deals correctly with incomplete DFSMs "
    << endl;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    Dfsm d("TC-DFSM-0001.fsm",pl,"m1");
    d.toDot("TC-DFSM-0001");
 
    vector<int> inp;
    inp.push_back(1);
    inp.push_back(0);
    inp.push_back(0);
    inp.push_back(0);
    inp.push_back(1);


    InputTrace i(inp,pl);
    
    cout << "InputTrace = " << i << endl;

    
    IOTrace t = d.applyDet(i);
    
    cout << "IOTrace t = " << t << endl;
    
    vector<int> vIn = t.getInputTrace().get();
    vector<int> vOut = t.getOutputTrace().get();
    assert("TC-DFSM-0001",vIn.size() == 4
           and vOut.size() == 4
           and vOut[0] == 2
           and vOut[1] == 0
           and vOut[2] == 2
           and vOut[3] == 2,
           "For input trace 1.0.0.0.1, the output trace is 2.0.2.2");
    
    
    inp.insert(inp.begin(),9);
    InputTrace j(inp,pl);
    IOTrace u = d.applyDet(j);
    cout << "IOTrace u = " << u << endl;
    assert("TC-DFSM-0001",
           u.getOutputTrace().get().size() == 0 and
           u.getInputTrace().get().size() == 0,
           "For input trace 9, the output trace is empty.");

}


void test2() {
    
    cout << "TC-FSM-0001 Show that the copy constructor produces a deep copy of an FSM generated at random "
    << endl;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    shared_ptr<Fsm> f1 = Fsm::createRandomFsm("f1",3,5,10,pl);
    
    shared_ptr<Fsm> f2 = make_shared<Fsm>(*f1);
    
    f1->toDot("f1");
    f2->toDot("f1Copy");
    
    // Check using diff, that the dot-files of both FSMs
    // are identical
    assert("TC-FSM-0001", 0 == system("diff f1.dot f1Copy.dot"),
           "dot-files of original and copied FSM are identical");
    
    cout << "Show that original FSM and deep copy are equivalent, "
    << endl << "using the WpMethod";
    
    Fsm f1Obs = f1->transformToObservableFSM();
    Fsm f1Min = f1Obs.minimise();
    
    Fsm f2Obs = f2->transformToObservableFSM();
    Fsm f2Min = f2Obs.minimise();
    
    int m = (f2Min.getMaxNodes() > f1Min.getMaxNodes() ) ?
            (f2Min.getMaxNodes() - f1Min.getMaxNodes()) : 0;
    IOListContainer iolc = f1Min.wMethod(m);
    
    TestSuite t1 = f1Min.createTestSuite(iolc);
    TestSuite t2 = f2Min.createTestSuite(iolc);
    
    assert("TC-FSM-0001",
           t2.isEquivalentTo(t1),
           "Original FSM and its deep copy pass the same W-Method test suite");

    
    
}

void test3() {
    
    cout << "TC-FSM-0002 Show that createMutant() injects a fault into the original FSM" << endl;

    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    shared_ptr<Fsm> fsm = Fsm::createRandomFsm("F",3,3,5,pl);
    fsm->toDot("F");
    
    shared_ptr<Fsm> fsmMutant = fsm->createMutant("F_M",3,0);
    fsmMutant->toDot("FMutant");
    
    Fsm fsmMin = fsm->minimise();
    fsmMin.toDot("FM");
    
    Fsm fsmMutantMin = fsmMutant->minimise();
    
    int m = 0;
    if ( fsmMutantMin.getMaxNodes() > fsmMin.getMaxNodes() ) {
        m = (int)fsmMutantMin.getMaxNodes() - fsmMin.getMaxNodes();
    }
    
    IOListContainer iolc1 = fsmMin.wMethod(m);
    
    TestSuite t1 = fsmMin.createTestSuite(iolc1);
    TestSuite t2 = fsmMutantMin.createTestSuite(iolc1);
    
    assert("TC-FSM-0002", not t2.isEquivalentTo(t1),
           "Original FSM and mutant do not produce the same test suite results - tests are created by W-Method");
    
    IOListContainer iolc2 = fsmMin.wpMethod(m);
    
    TestSuite t1wp = fsmMin.createTestSuite(iolc2);
    TestSuite t2wp = fsmMutantMin.createTestSuite(iolc2);
    
    assert("TC-FSM-0002",
           not t2wp.isEquivalentTo(t1wp),
           "Original FSM and mutant do not produce the same test suite results - tests are created by W-Method");
    
    assert("TC-FSM-0002",
           t1wp.size() < t1.size(),
           "Wp-Method generates a smaller test suite than W-Method");

    
    

}

int main(int argc, char* argv[])
{
    
    test1();
    test2();
    test3();
    
    exit(0);
    
}



