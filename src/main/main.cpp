/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <stdlib.h>
#include <interface/FsmPresentationLayer.h>
#include <fsm/Dfsm.h>
#include <fsm/Fsm.h>
#include <fsm/FsmNode.h>
#include <fsm/IOTrace.h>
#include <fsm/FsmPrintVisitor.h>
#include <fsm/FsmSimVisitor.h>
#include <fsm/FsmOraVisitor.h>
#include <trees/IOListContainer.h>
#include <trees/OutputTree.h>
#include <trees/TestSuite.h>


using namespace std;

void assertInconclusive(string tc, string comment = "") {
    
    string sVerdict("INCONCLUSIVE");
    cout << sVerdict << ": " << tc << " : " << comment <<  endl;
    
}

void assert(string tc, bool verdict, string comment = "") {
    
    string sVerdict = (verdict) ? "PASS" : "FAIL";
    cout << sVerdict << ": " << tc
    << " : "
    << comment <<  endl;
    
}

void test1() {
    
    cout << "TC-DFSM-0001 Show that Dfsm.applyDet() deals correctly with incomplete DFSMs "
    << endl;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    Dfsm d("../../../resources/TC-DFSM-0001.fsm",pl,"m1");
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
    
    
    for ( size_t i = 0; i < 10; i++ ) {
        shared_ptr<FsmPresentationLayer> pl =
        make_shared<FsmPresentationLayer>();
        shared_ptr<Fsm> fsm = Fsm::createRandomFsm("F",5,5,8,pl,i);
        fsm->toDot("F");
        
        shared_ptr<Fsm> fsmMutant = fsm->createMutant("F_M",1,0);
        fsmMutant->toDot("FMutant");
        
        Fsm fsmMin = fsm->minimise();
        fsmMin.toDot("FM");
        
        Fsm fsmMutantMin = fsmMutant->minimise();
        
        unsigned int m = 0;
        if ( fsmMutantMin.getMaxNodes() > fsmMin.getMaxNodes() ) {
            m = fsmMutantMin.getMaxNodes() - fsmMin.getMaxNodes();
        }
        
        cout << "Call W-Method - additional states (m) = " << m << endl;
        
        IOListContainer iolc1 = fsmMin.wMethodOnMinimisedFsm(m);
        
        cout << "TS SIZE (W-Method): " << iolc1.size() << endl;
        
        if ( iolc1.size() > 100000) {
            cout << "Skip this test case, since size is too big" << endl;
            continue;
        }
        
        TestSuite t1 = fsmMin.createTestSuite(iolc1);
        TestSuite t2 = fsmMutantMin.createTestSuite(iolc1);
        
        assert("TC-FSM-0002", not t2.isEquivalentTo(t1),
               "Original FSM and mutant do not produce the same test suite results - tests are created by W-Method");
        
        IOListContainer iolc2 = fsmMin.wpMethod(m);
        
        cout << "TS SIZE (Wp-Method): " << iolc2.size() << endl;
        
        if ( iolc2.size() > iolc1.size() ) {
            
            ofstream outFile("fsmMin.fsm");
            fsmMin.dumpFsm(outFile);
            outFile.close();
             
            exit(1);
        }

        
        TestSuite t1wp = fsmMin.createTestSuite(iolc2);
        TestSuite t2wp = fsmMutantMin.createTestSuite(iolc2);
        
        assert("TC-FSM-0002",
               not t2wp.isEquivalentTo(t1wp),
               "Original FSM and mutant do not produce the same test suite results - tests are created by Wp-Method");
        
        assert("TC-FSM-0002",
               t1wp.size() <= t1.size(),
               "Wp-Method test suite size less or equal to W-Method size");
        
        if ( t1wp.size() > t1.size() ) {
            cout << "Test Suite Size (W-Method): " << t1.size()
            << endl << "Test Suite Size (Wp-Method): " << t1wp.size() << endl;
            cout << endl << "W-Method " << endl << iolc1 << endl;
            exit(1);
        }
        
        
    }
    
    
}


void test4() {
    
    cout << "TC-FSM-0004 Check correctness of state cover" << endl;
    
    const bool markAsVisited = true;
    bool havePassed = true;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
    for (size_t i = 0; i < 2000; i++) {
        
        // Create a random FSM
        std::shared_ptr<Fsm> f = Fsm::createRandomFsm("F",5,5,10,pl,i);
        std::shared_ptr<Tree> sc = f->getStateCover();
        
        if ( sc->size() != (size_t)f->getMaxNodes() + 1 ) {
            cout << "Size of state cover: " << sc->size()
            << " Number of states in FSM: " << f->getMaxNodes() + 1 << endl;
            assert("TC-FSM-0004",
                   sc->size() <= (size_t)f->getMaxNodes() + 1,
                   "Size of state cover must be less or equal than number of FSM states");
        }
        
        
        IOListContainer c = sc->getTestCases();
        std::shared_ptr<std::vector<std::vector<int>>> iols = c.getIOLists();
        
        for ( auto inLst : *iols ) {
            auto iTr = make_shared<InputTrace>(inLst,pl);
            f->apply(*iTr,true);
        }
        
        for ( std::shared_ptr<FsmNode> n : f->getNodes() ) {
            if ( not n->hasBeenVisited() ) {
                havePassed = false;
                assert("TC-FSM-0004",
                       n->hasBeenVisited(),
                       "State cover failed to visit node " + n->getName());
                
                f->toDot("FailedStateCoverFSM");
                
                filebuf fb;
                fb.open ("FailedStateCover.dot",std::ios::out);
                ostream os(&fb);
                sc->toDot(os);
                fb.close();
                
                int iCtr = 0;
                for ( auto inLst : *iols ) {
                    ostringstream oss;
                    oss << iCtr++;
                    auto iTr = make_shared<InputTrace>(inLst,pl);
                    filebuf fbot;
                    OutputTree ot = f->apply(*iTr,markAsVisited);
                    fbot.open ("FailedStateCover" + oss.str() + ".dot",
                               std::ios::out);
                    ostream osdot(&fbot);
                    sc->toDot(osdot);
                    fbot.close();
                }
                
                exit(1);
                
            }
        }
        
    }
    
    if ( havePassed ) {
        assert("TC-FSM-0004",
               true,
               "State cover reaches all states");
    }
    else {
        exit(0);
    }
    
    
}

void test5() {
    
    cout << "TC-FSM-0005 Check correctness of input " <<
    "equivalence classes" << endl;
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>();
    
    
    shared_ptr<Fsm> fsm =
    make_shared<Fsm>("../../../resources/TC-FSM-0005.fsm",pl,"F");
    fsm->toDot("TC-FSM-0005");
    
    vector< std::unordered_set<int> > v = fsm->getEquivalentInputs();
    
    for ( size_t s = 0; s < v.size(); s++ ) {
        cout << s << ": { ";
        bool isFirst = true;
        for ( auto x : v[s] ) {
            if ( isFirst ) {
                isFirst= false;
            }
            else   {
                cout << ", ";
            }
            cout << x;
        }
        cout << " }" << endl;
    }
    
    assert("TC-FSM-0005",
           v.size() == 3,
           "For TC-FSM-0005.fsm, there are 3 classes of equivalent inputs.");
    
    assert("TC-FSM-0005",
           v[0].size() == 1 and v[0].find(0) != v[0].end(),
           "Class 0 only contains input 0.");
    
    assert("TC-FSM-0005",
           v[1].size() == 1 and v[1].find(1) != v[1].end(),
           "Class 1 only contains input 1.");
    
    assert("TC-FSM-0005",
           v[2].size() == 2 and
           v[2].find(2) != v[2].end() and
           v[2].find(3) != v[2].end(),
           "Class 2 contains inputs 2 and 3.");
    
    
    // Check FSM without any equivalent inputs
    fsm = make_shared<Fsm>("../../../resources/fsmGillA7.fsm",pl,"F");
    fsm->toDot("fsmGillA7");
    v = fsm->getEquivalentInputs();
    
    assert("TC-FSM-0005",
           v.size() == 3,
           "For fsmGillA7, there are 3 input classes.");
    
    bool ok = true;
    for ( size_t s=0; s < v.size() and ok; s++ ) {
        if ( v[s].size() != 1 or
            v[s].find((int)s) == v[s].end() ) {
            ok =false;
        }
    }
    
    assert("TC-FSM-0005",
           ok,
           "For fsmGillA7, class x just contains input x.");
    
}

void test6() {
    
    cout << "TC-FSM-0006 Check correctness of FSM Print Visitor " << endl;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    Dfsm d("../../../resources/TC-DFSM-0001.fsm",pl,"m1");
    
    FsmPrintVisitor v;
    
    d.accept(v);
    
    cout << endl << endl;
    assertInconclusive("TC-FSM-0006",
                       "Output of print visitor has to be checked manually");
    
    
}

void test7() {
    
    //    cout << "TC-FSM-0007 Check correctness of FSM Simulation Visitor "
    //    << endl;
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
                                      "../../../resources/garageOut.txt",
                                      "../../../resources/garageState.txt");
    Dfsm d("../../../resources/garage.fsm",pl,"GC");
    d.toDot("GC");
    
    FsmSimVisitor v;
    
    d.accept(v);
    
    v.setFinalRun(true);
    d.accept(v);
    
    cout << endl << endl;
    //    assertInconclusive("TC-FSM-0007",
    //                       "Output of simulation visitor has to be checked manually");
}


void test8() {
    
    //    cout << "TC-FSM-0008 Check correctness of FSM Oracle Visitor "
    //    << endl;
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
                                      "../../../resources/garageOut.txt",
                                      "../../../resources/garageState.txt");
    Dfsm d("../../../resources/garage.fsm",pl,"GC");
    d.toDot("GC");
    
    FsmOraVisitor v;
    
    d.accept(v);
    
    v.setFinalRun(true);
    d.accept(v);
    
    cout << endl << endl;
    //    assertInconclusive("TC-FSM-0008",
    //                       "Output of oracle visitor has to be checked manually");
}


void fsbrts() {
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>("../../../resources/fsbrts.in",
                                      "../../../resources/fsbrts.out",
                                      "../../../resources/fsbrts.state");
    
    
    shared_ptr<Dfsm> fsmFsb =
    make_shared<Dfsm>("../../../resources/fsbrts.fsm",pl,"FSB");
    fsmFsb->toDot("FSBRTS");
    fsmFsb->toCsv("FSBRTS");
    
    cout << "FSBRTS is deterministic: " << fsmFsb->isDeterministic() << endl;
    
    Dfsm fsbMin = fsmFsb->minimise();
    fsbMin.toDot("FSBRTS_MIN");
    
    shared_ptr<Tree> scov = fsbMin.getStateCover();
    ofstream fscov("SCOV.dot");
    scov->toDot(fscov);
    fscov.close();
    
    
    shared_ptr<Tree> tcov = fsbMin.getTransitionCover();
    ofstream ftcov("TCOV.dot");
    tcov->toDot(ftcov);
    ftcov.close();
    
    
    IOListContainer w = fsbMin.getCharacterisationSet();
    cout << "W = " << endl << w << endl<<endl;
    
    
    
    IOListContainer ts =
    fsbMin.wMethodOnMinimisedDfsm(1);
    
    cout << "Test Suite Size = " << ts.size() << endl;
    
#if 1
    int n = 0;
    for ( auto lst : *ts.getIOLists()) {
        
        cout << " TC " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        OutputTree otree = fsmFsb->apply(*itr);
        
        cout << otree << endl << endl;
        
    }
#endif
    
}



void fsbrtssafe() {
    
    shared_ptr<FsmPresentationLayer> pl0 =
    make_shared<FsmPresentationLayer>("../../../resources/fsbrts.in",
                                      "../../../resources/fsbrts.out",
                                      "../../../resources/fsbrts.state");
    
    
    shared_ptr<Dfsm> fsmFsbOrig =
    make_shared<Dfsm>("../../../resources/fsbrts.fsm",pl0,"FSB");
    fsmFsbOrig->toDot("FSBRTS");
    
    cout << "FSBRTS is deterministic: " << fsmFsbOrig->isDeterministic() << endl;
    
    Dfsm fsbMinOrig = fsmFsbOrig->minimise();
    fsbMinOrig.toDot("FSBRTS_MIN");
    
    
    shared_ptr<Tree> tcovOrig = fsbMinOrig.getTransitionCover();
    ofstream ftcovOrig("TCOV.dot");
    tcovOrig->toDot(ftcovOrig);
    ftcovOrig.close();
    
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>("../../../resources/fsbrtssafe.in",
                                      "../../../resources/fsbrtssafe.out",
                                      "../../../resources/fsbrtssafe.state");
    
    
    shared_ptr<Dfsm> fsmFsb =
    make_shared<Dfsm>("../../../resources/fsbrtssafe.fsm",pl,"FSBSAFE");
    fsmFsb->toDot("FSBRTSSAFE");
    
    cout << "FSBRTSSAFE is deterministic: "
    << fsmFsb->isDeterministic() << endl;
    
    Dfsm fsbMin = fsmFsb->minimise();
    fsbMin.toDot("FSBRTSSAFE_MIN");
    
    shared_ptr<Tree> scov = fsbMin.getStateCover();
    ofstream fscov("SCOVSAFE.dot");
    scov->toDot(fscov);
    fscov.close();
    
    
    shared_ptr<Tree> tcov = fsbMin.getTransitionCover();
    ofstream ftcov("TCOVSAFE.dot");
    tcov->toDot(ftcov);
    ftcov.close();
    
    IOListContainer w = fsbMin.getCharacterisationSet();
    cout << "W_SAFE = " << endl << w << endl<<endl;
    
    IOListContainer ts =
    fsbMin.wMethodOnMinimisedDfsm(0);
    cout << "Test Suite Size = " << ts.size() << endl;
    
    int n = 0;
    for ( auto lst : *ts.getIOLists()) {
        
        cout << " TC " << ++n << endl;
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        OutputTree otree = fsbMin.apply(*itr);
        
        cout << otree << endl << endl;
        
        
    }
    
    
    cout << "NEW --------------------------" << endl;
    shared_ptr<Tree> iTree = tcovOrig;
    
    IOListContainer inputEnum = IOListContainer(6,
                                                1,
                                                1,
                                                pl0);
    iTree->add(inputEnum);
    iTree->add(w);
    
    cout << "iTree size = " << iTree->size() << endl;
    
    IOListContainer finalCont = iTree->getIOLists();
    
    vector<IOTrace> iotrLst;
    
    n = 0;
    for ( auto lst : *finalCont.getIOLists() ) {
        
        cout << "TCNEW " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        IOTrace iotr = fsbMinOrig.applyDet(*itr);
        cout << iotr << endl;
        iotrLst.push_back(iotr);
        
    }
    
    // Now run the test suite against an SUT FSM (mutant of the original
    // FSBRTS machine)
    shared_ptr<Dfsm> fsmFsbMutant =
    make_shared<Dfsm>("../../../resources/fsbrts_mutant.fsm",pl0,"FSB");
    fsmFsbMutant->toDot("FSBRTS_MUTANT");
    
    Dfsm fsmFsbMutantMin = fsmFsbMutant->minimise();
    cout << "Mutant state number = " << fsmFsbMutantMin.getNodes().size() <<endl;
    
    n = 0;
    for ( auto iot : iotrLst ) {
        
        ++n;
        if ( not fsmFsbMutant->pass(iot) ) {
            IOTrace passTrace = fsbMinOrig.applyDet(iot.getInputTrace());
            IOTrace failTrace = fsmFsbMutant->applyDet(iot.getInputTrace());

            
            cout << "FAIL TC " << n << endl << "Expected: " << passTrace
            << endl << "Observed: " << failTrace << endl<<endl;
            
            
        }
        else {
            cout << "PASS TC " << n << endl;
        }
        
    }
    
    
    
    
}

void check() {
    
    
    shared_ptr<FsmPresentationLayer> pl0 =
            make_shared<FsmPresentationLayer>();
    
    
    shared_ptr<Dfsm> fsmFsbOrig =
    make_shared<Dfsm>("../../../resources/fsbrts.fsm",pl0,"FSB");
    fsmFsbOrig->toDot("FSBRTS");
    
    cout << "FSBRTS is deterministic: " << fsmFsbOrig->isDeterministic() << endl;
    
    Dfsm fsbMinOrig = fsmFsbOrig->minimise();
    fsbMinOrig.toDot("FSBRTS_MIN");
    fsmFsbOrig->printTables();
    
    Fsm fsmIntersection = fsmFsbOrig->intersect(fsbMinOrig);
    fsmIntersection.toDot("INTER");
    
    cout << "intersection is completely specified: "
    << fsmIntersection.isCompletelyDefined() << endl;
    
}

void readCsv() {
    shared_ptr<Dfsm> fsmFsbOrig =
    make_shared<Dfsm>("plsautocnd.csv","FSB_AUTO_CND");
    
    fsmFsbOrig->toCsv("plsautocndOut");
    fsmFsbOrig->toDot("plsautocnd");

    shared_ptr<Dfsm> fsm2 =
    make_shared<Dfsm>("plsautocndOut.csv","FSB2");
    fsm2->toDot("plsautocndOut");
    
    Fsm inter = fsmFsbOrig->intersect(*fsm2);
    
    inter.toDot("inter");
    
    cout << endl << "Is fsmFsbOrig completely specified: "
    << fsmFsbOrig->isCompletelyDefined() << endl;
    
    cout << endl << "Is fsm2 completely specified: "
    << fsm2->isCompletelyDefined() << endl;

    
    cout << endl << "Is inter completely specified: "
    << inter.isCompletelyDefined() << endl;
    
}

void checkPlsAuto1() {
    
    shared_ptr<Dfsm> fsmFsbOrig =
    make_shared<Dfsm>("plsautocnd_I.csv","FSB_AUTO_CND");
    
    shared_ptr<Dfsm> fsm2 =
    make_shared<Dfsm>("plsautocnd_I_MUT.csv","FSB_AUTO_MUT");
    
    Fsm inter = fsmFsbOrig->intersect(*fsm2);
    
    cout << endl << "Is fsmFsbOrig completely specified: "
    << fsmFsbOrig->isCompletelyDefined() << endl;
    
    cout << endl << "Is fsm2 completely specified: "
    << fsm2->isCompletelyDefined() << endl;
    
    
    cout << endl << "Is inter completely specified: "
    << inter.isCompletelyDefined() << endl;
    
    Dfsm f1Min = fsmFsbOrig->minimise();
    Dfsm f2Min = fsm2->minimise();
    
    cout << "f1Min states: " << f1Min.getNodes().size() << endl
    << "f2Min states: " << f2Min.getNodes().size() << endl;
    
    unsigned int mMinusN =
    (f2Min.getNodes().size() > f1Min.getNodes().size()) ?
    (f2Min.getNodes().size() - f1Min.getNodes().size()) :
    0;

    IOListContainer ts = f1Min.wMethodOnMinimisedDfsm(mMinusN);
    
    std::shared_ptr<std::vector<std::vector<int>>> tcList = ts.getIOLists();
    
    int tcNo = 1;
    for ( auto tc : *tcList ) {
        
        // Get input trace from test case and calculate
        // associated IOTrace containing inputs plus
        // expected results
        cout << "TC " << tcNo++ << ": ";
        shared_ptr<InputTrace> inTrc = make_shared<InputTrace>(tc,
                                                               fsmFsbOrig->getPresentationLayer());
        IOTrace ioTrc = fsmFsbOrig->applyDet(*inTrc);
        if ( fsm2->pass(ioTrc) ) {
            cout << "PASS: " << ioTrc << endl;
        }
        else {
            IOTrace ioTrcMUT = fsm2->applyDet(*inTrc);
            cout << "FAIL" << endl
            << "Expected: " << ioTrc << endl
            << "Observed: " << ioTrcMUT << endl << endl;
        }
        
        
    }
    
    
    
}


void fsbrtsTestComplete() {
    
    shared_ptr<Dfsm> fsbrtsRef =
    make_shared<Dfsm>("FSBRTSX.csv","FSBRTS");
    Dfsm fsbrtsRefMin = fsbrtsRef->minimise();
    cout << "Reference model (minimised) state number = "
    << fsbrtsRefMin.getNodes().size() <<endl;
    
    shared_ptr<FsmPresentationLayer> pl = fsbrtsRef->getPresentationLayer();

    
    // Test mutant against reference model
    shared_ptr<Dfsm> fsbrts_m1 =
    make_shared<Dfsm>("FSBRTSX_M3.csv","FSBRTS_M3");
    Dfsm fsbrts_m1_min = fsbrts_m1->minimise();
    cout << "Mutant state number = " << fsbrts_m1_min.getNodes().size() <<endl;

    unsigned int mMinusN =
    (fsbrts_m1_min.getNodes().size() > fsbrtsRefMin.getNodes().size()) ?
    (fsbrts_m1_min.getNodes().size() - fsbrtsRefMin.getNodes().size()) : 0;
    
    IOListContainer tcInputs1
    = fsbrtsRefMin.wMethodOnMinimisedDfsm(mMinusN);
    cout << "Test Suite Size = " << tcInputs1.size() << endl;
    
    int n = 0;
    vector<IOTrace> iotrLst;
    for ( auto lst : *tcInputs1.getIOLists() ) {
        
        cout << "TCNEW " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        IOTrace iotr = fsbrtsRefMin.applyDet(*itr);
        cout << iotr << endl;
        iotrLst.push_back(iotr);
        
    }
    
    n = 0;
    for ( auto iot : iotrLst ) {
        
        ++n;
        if ( not fsbrts_m1->pass(iot) ) {
            IOTrace passTrace = fsbrtsRef->applyDet(iot.getInputTrace());
            IOTrace failTrace = fsbrts_m1->applyDet(iot.getInputTrace());
            
            cout << "FAIL TC " << n << endl << "Expected: " << passTrace
            << endl << "Observed: " << failTrace << endl<<endl;
            
        }
        else {
            cout << "PASS TC " << n << endl;
        }
        
    }
    
    // Test second mutant against reference model
    shared_ptr<Dfsm> fsbrts_m2 =
    make_shared<Dfsm>("FSBRTSX_M2.csv","FSBRTS_M2");
    Dfsm fsbrts_m2_min = fsbrts_m2->minimise();
    cout << "Mutant state number = " << fsbrts_m2_min.getNodes().size() <<endl;
    
    mMinusN =
    (fsbrts_m2_min.getNodes().size() > fsbrtsRefMin.getNodes().size()) ?
    (fsbrts_m2_min.getNodes().size() - fsbrtsRefMin.getNodes().size()) : 0;
    
    IOListContainer tcInputs2
    = fsbrtsRefMin.wMethodOnMinimisedDfsm(mMinusN);
    cout << "Test Suite Size = " << tcInputs2.size() << endl;
    
    n = 0;
    vector<IOTrace> iotrLst2;
    for ( auto lst : *tcInputs2.getIOLists() ) {
        
        cout << "TCNEW " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        IOTrace iotr = fsbrtsRefMin.applyDet(*itr);
        cout << iotr << endl;
        iotrLst2.push_back(iotr);
        
    }
    
    n = 0;
    for ( auto iot : iotrLst2 ) {
        
        ++n;
        if ( not fsbrts_m2->pass(iot) ) {
            IOTrace passTrace = fsbrtsRef->applyDet(iot.getInputTrace());
            IOTrace failTrace = fsbrts_m2->applyDet(iot.getInputTrace());
            
            cout << "FAIL TC " << n << endl << "Expected: " << passTrace
            << endl << "Observed: " << failTrace << endl<<endl;
            
        }
        else {
            cout << "PASS TC " << n << endl;
        }
        
    }
    
}


void fsbrtsTestCompleteSafe() {
    
    shared_ptr<Dfsm> fsbrtsRef =
    make_shared<Dfsm>("FSBRTSX.csv","FSBRTS");
    Dfsm fsbrtsRefMin = fsbrtsRef->minimise();
    cout << "Reference model (minimised) state number = "
    << fsbrtsRefMin.getNodes().size() <<endl;
    
    shared_ptr<FsmPresentationLayer> pl = fsbrtsRef->getPresentationLayer();
    
    // Get state cover of original model
    shared_ptr<Tree> scov = fsbrtsRefMin.getStateCover();
    
    // Get transition cover of original model
    shared_ptr<Tree> tcov = fsbrtsRefMin.getTransitionCover();

    // Get characterisation set of original model
    IOListContainer w = fsbrtsRefMin.getCharacterisationSet();
    cout << "W size = " << w.size() << endl;
    
    // Get characterisation set of reference model's safety abstraction
    shared_ptr<Dfsm> fsbrtsRefSafe =
    make_shared<Dfsm>("FSBRTSX_SAFE.csv","FSBRTS_SAFE");
    Dfsm fsbrtsRefSafeMin = fsbrtsRefSafe->minimise();
    cout << "Safety abstraction of reference model (minimised) state number = "
    << fsbrtsRefSafeMin.getNodes().size() <<endl;
    IOListContainer wSafe = fsbrtsRefSafeMin.getCharacterisationSet();
    cout << "WSafe size = " << wSafe.size() << endl;

    // Read mutant
    shared_ptr<Dfsm> fsbrts_m1 =
    make_shared<Dfsm>("FSBRTSX_M3.csv","FSBRTS_M3");
    Dfsm fsbrts_m1_min = fsbrts_m1->minimise();
    cout << "Mutant state number = " << fsbrts_m1_min.getNodes().size() <<endl;
    
    // Calculate W-Test Suite with safety reduction for first mutant
    int mMinusN =
    ( fsbrts_m1_min.getNodes().size() > fsbrtsRefMin.getNodes().size() )
    ? 0 : (fsbrts_m1_min.getNodes().size() - fsbrtsRefMin.getNodes().size());
    
    shared_ptr<Tree> tsTree = fsbrtsRefMin.getStateCover();
    tsTree->add(w);

    shared_ptr<Tree> vswsTree = fsbrtsRefMin.getTransitionCover();
    vswsTree->add(wSafe);
    tsTree->addToRoot(vswsTree->getIOLists());
    
    if ( mMinusN > 0 ) {
        shared_ptr<Tree> xTree = fsbrtsRefMin.getTransitionCover();
        IOListContainer inputEnum = IOListContainer(fsbrtsRef->getMaxInput(),
                                                    1,
                                                    mMinusN,
                                                    pl);
        xTree->add(inputEnum);
        xTree->add(wSafe);
        tsTree->addToRoot(xTree->getIOLists());
    }
    
    cout << "Total safety-reduced test suite size: " << tsTree->size() << endl;
    
    // Execute reduced test suite against mutant
    int n = 0;
    vector<IOTrace> iotrLst;
    IOListContainer tsCont = tsTree->getIOLists();
    cout << "x " << tsCont.size() << endl;
    auto iolCnt = tsCont.getIOLists();
    cout << "y " << iolCnt->size() << endl;

    for ( auto lst : *iolCnt ) {
        
        cout << "TCNEW " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        IOTrace iotr = fsbrtsRefMin.applyDet(*itr);
        cout << iotr << endl;
        iotrLst.push_back(iotr);
        
    }
    
    n = 0;
    for ( auto iot : iotrLst ) {
        
        ++n;
        if ( not fsbrts_m1->pass(iot) ) {
            IOTrace passTrace = fsbrtsRef->applyDet(iot.getInputTrace());
            IOTrace failTrace = fsbrts_m1->applyDet(iot.getInputTrace());
            
            cout << "FAIL TC " << n << endl << "Expected: " << passTrace
            << endl << "Observed: " << failTrace << endl<<endl;
            
        }
        else {
            cout << "PASS TC " << n << endl;
        }
        
    }

    
}


void fsbrtsTestCompleteSafe2() {
    
    shared_ptr<Dfsm> fsbrtsRef =
    make_shared<Dfsm>("FSBRTSX.csv","FSBRTS");
    Dfsm fsbrtsRefMin = fsbrtsRef->minimise();
    cout << "Reference model (minimised) state number = "
    << fsbrtsRefMin.getNodes().size() <<endl;
    
    shared_ptr<FsmPresentationLayer> pl = fsbrtsRef->getPresentationLayer();
    
    // Get state cover of original model
    shared_ptr<Tree> scov = fsbrtsRefMin.getStateCover();
    
    // Get transition cover of original model
    shared_ptr<Tree> tcov = fsbrtsRefMin.getTransitionCover();
    
    // Get characterisation set of original model
    IOListContainer w = fsbrtsRefMin.getCharacterisationSet();
    cout << "W size = " << w.size() << endl;
    
    // Get characterisation set of reference model's safety abstraction
    shared_ptr<Dfsm> fsbrtsRefSafe =
    make_shared<Dfsm>("FSBRTSX_SAFE.csv","FSBRTS_SAFE");
    Dfsm fsbrtsRefSafeMin = fsbrtsRefSafe->minimise();
    cout << "Safety abstraction of reference model (minimised) state number = "
    << fsbrtsRefSafeMin.getNodes().size() <<endl;
    IOListContainer wSafe = fsbrtsRefSafeMin.getCharacterisationSet();
    cout << "WSafe size = " << wSafe.size() << endl;
    
    // Read the first mutant
    shared_ptr<Dfsm> fsbrts_m2 =
    make_shared<Dfsm>("FSBRTSX_M2.csv","FSBRTS_M2");
    Dfsm fsbrts_m2_min = fsbrts_m2->minimise();
    cout << "Mutant state number = " << fsbrts_m2_min.getNodes().size() <<endl;
    
    // Calculate W-Test Suite with safety reduction for first mutant
    int mMinusN =
    ( fsbrts_m2_min.getNodes().size() > fsbrtsRefMin.getNodes().size() )
    ? 0 : (fsbrts_m2_min.getNodes().size() - fsbrtsRefMin.getNodes().size());
    
    shared_ptr<Tree> tsTree = fsbrtsRefMin.getStateCover();
    tsTree->add(w);
    
    shared_ptr<Tree> vswsTree = fsbrtsRefMin.getTransitionCover();
    vswsTree->add(wSafe);
    tsTree->addToRoot(vswsTree->getIOLists());
    
    if ( mMinusN > 0 ) {
        shared_ptr<Tree> xTree = fsbrtsRefMin.getTransitionCover();
        IOListContainer inputEnum = IOListContainer(fsbrtsRef->getMaxInput(),
                                                    1,
                                                    mMinusN,
                                                    pl);
        xTree->add(inputEnum);
        xTree->add(wSafe);
        tsTree->addToRoot(xTree->getIOLists());
    }
    
    cout << "Total safety-reduced test suite size: " << tsTree->size() << endl;
    
    // Execute reduced test suite against mutant
    int n = 0;
    vector<IOTrace> iotrLst;
    IOListContainer tsCont = tsTree->getIOLists();
    cout << "x " << tsCont.size() << endl;
    auto iolCnt = tsCont.getIOLists();
    cout << "y " << iolCnt->size() << endl;
    
    for ( auto lst : *iolCnt ) {
        
        cout << "TCNEW " << ++n << ": ";
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        
        IOTrace iotr = fsbrtsRefMin.applyDet(*itr);
        cout << iotr << endl;
        iotrLst.push_back(iotr);
        
    }
    
    n = 0;
    for ( auto iot : iotrLst ) {
        
        ++n;
        if ( not fsbrts_m2->pass(iot) ) {
            IOTrace passTrace = fsbrtsRef->applyDet(iot.getInputTrace());
            IOTrace failTrace = fsbrts_m2->applyDet(iot.getInputTrace());
            
            cout << "FAIL TC " << n << endl << "Expected: " << passTrace
            << endl << "Observed: " << failTrace << endl<<endl;
            
        }
        else {
            cout << "PASS TC " << n << endl;
        }
        
    }

}

void gdc_test1() {
    
    cout << "TC-GDC-0001 Check that the correct W-Method test suite "
    << endl << "is generated for the garage door controller example" << endl;

    
    shared_ptr<Dfsm> gdc =
    make_shared<Dfsm>("../../../resources/garage-door-controller.csv","GDC");
    
    shared_ptr<FsmPresentationLayer> pl = gdc->getPresentationLayer();
    
    gdc->toDot("GDC");
    gdc->toCsv("GDC");
    
    IOListContainer iolc =
        gdc->wMethod(2);
    
    shared_ptr< TestSuite > testSuite =
        make_shared< TestSuite >();
    for ( auto inVec : *iolc.getIOLists() ) {
        shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
        testSuite->push_back(gdc->apply(*itrc));
    }
    
    int tcNum = 0;
    for ( auto iotrc : *testSuite ) {
        cout << "TC-" << ++tcNum << ": " << iotrc;
    }
    
    testSuite->save("testsuite.txt");
    
    assert("TC-GDC-0001",
            0 == system("diff testsuite.txt ../../../resources/gdc-testsuite.txt"),
           "Expected GDC test suite and generated suite are identical");
    
    
}


int main()
{
    
    
    
#if 0
    int tc = atoi(argv[1]);
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    fsbrts();
    fsbrtssafe();
    check();
    readCsv();
    checkPlsAuto1();
    fsbrtsTestComplete();
    fsbrtsTestCompleteSafe();
    fsbrtsTestCompleteSafe2();
    
    gdc_test1();
    
    switch (tc) {
        case 1:
            checkPlsAuto1();
            break;
        case 2:
            fsbrtsTestComplete();
            break;
        case 3:
            fsbrtsTestCompleteSafe();
            break;
        case 4:
            fsbrtsTestCompleteSafe2();
            break;
        default:
            break;
    }
    
    

#endif
    
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    gdc_test1();

    exit(0);
    
}



