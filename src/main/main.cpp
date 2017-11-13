/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <random>
#include <stdlib.h>
#include <interface/FsmPresentationLayer.h>
#include <fsm/Dfsm.h>
#include <fsm/Fsm.h>
#include <fsm/FsmNode.h>
#include <fsm/IOTrace.h>
#include <fsm/IOTraceContainer.h>
#include <fsm/FsmPrintVisitor.h>
#include <fsm/FsmSimVisitor.h>
#include <fsm/FsmOraVisitor.h>
#include <trees/IOListContainer.h>
#include <trees/IOTreeContainer.h>
#include <trees/OutputTree.h>
#include <trees/TestSuite.h>
#include "json/json.h"
#include "logging/easylogging++.h"
#include "logging/Logging.h"

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

using namespace std;
using namespace Json;

void assertInconclusive(string tc, string comment = "") {
    
    string sVerdict("INCONCLUSIVE");
    CLOG(INFO, logging::globalLogger) << sVerdict << ": " << tc << " : " << comment <<  endl;
    
}

void assert(string tc, bool verdict, string comment = "") {
    
    string sVerdict = (verdict) ? "PASS" : "FAIL";
    CLOG(INFO, logging::globalLogger) << sVerdict << ": " << tc
    << " : "
    << comment <<  endl;
    
}

void assertOnFail(string tc, bool verdict, string comment = "") {

    string sVerdict = (verdict) ? "PASS: " + tc : "FAIL: " + tc + ": " + comment;
    if (verdict)
    {
        CLOG(INFO, logging::globalLogger) << sVerdict;
    }
    else
    {
        CLOG(ERROR, logging::globalLogger) << sVerdict;
    }


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
        shared_ptr<Fsm> fsm = Fsm::createRandomFsm("F",5,5,8,pl,false,i);
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
        std::shared_ptr<Fsm> f = Fsm::createRandomFsm("F",5,5,10,pl,false,i);
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

void test9() {
    
    cout << "TC-FSM-0009 Check correctness of method removeUnreachableNodes() "
         << endl;
    
    shared_ptr<Dfsm> d = nullptr;
    Reader jReader;
    Value root;
    stringstream document;
    ifstream inputFile("../../../resources/unreachable_gdc.fsm");
    document << inputFile.rdbuf();
    inputFile.close();
    
    if ( jReader.parse(document.str(),root) ) {
        d = make_shared<Dfsm>(root);
    }
    else {
        cerr << "Could not parse JSON model - exit." << endl;
        exit(1);
    }
    
    
    d->toDot("GU");
    
    size_t oldSize = d->size();
    
    vector<shared_ptr<FsmNode>> uNodes;
    if ( d->removeUnreachableNodes(uNodes) ) {
        
        d->toDot("G_all_reachable");
        
        for ( auto n : uNodes ) {
            cout << "Removed unreachable node: " << n->getName() << endl;
        }
        
        assert("TC-FSM-0009",
               uNodes.size() == 2 and (oldSize - d->size()) == 2,
               "All unreachable states have been removed");
    }
    else {
        assert("TC-FSM-0009",
               false,
               "Expected removeUnreachableNodes() to return FALSE");
    }
    
    
}


void test10() {
    
    cout << "TC-FSM-0010 Check correctness of Dfsm::minimise() "
    << endl;
    
    shared_ptr<Dfsm> d = nullptr;
    shared_ptr<FsmPresentationLayer> pl;
    Reader jReader;
    Value root;
    stringstream document;
    ifstream inputFile("../../../resources/unreachable_gdc.fsm");
    document << inputFile.rdbuf();
    inputFile.close();
    
    if ( jReader.parse(document.str(),root) ) {
        d = make_shared<Dfsm>(root);
        pl = d->getPresentationLayer();
    }
    else {
        cerr << "Could not parse JSON model - exit." << endl;
        exit(1);
    }
    
    
    Dfsm dMin = d->minimise();
    
    IOListContainer w = dMin.getCharacterisationSet();
    
    shared_ptr<std::vector<std::vector<int>>> inLst = w.getIOLists();
    
    bool allNodesDistinguished = true;
    for ( size_t n = 0; n < dMin.size(); n++ ) {
        
        shared_ptr<FsmNode> node1 = dMin.getNodes().at(n);
        
        for ( size_t m = n+1; m < dMin.size(); m++ ) {
            shared_ptr<FsmNode> node2 = dMin.getNodes().at(m);
            
            bool areDistinguished = false;
            
            for ( auto inputs : *inLst ) {
                
                shared_ptr<InputTrace> itr = make_shared<InputTrace>(inputs,pl);
                
                OutputTree o1 = node1->apply(*itr);
                OutputTree o2 = node2->apply(*itr);
                
                if ( o1 != o2 ) {
                    areDistinguished = true;
                    break;
                }
                
            }
            
            if ( not areDistinguished ) {
                
                assert("TC-FSM-0010",
                       false,
                       "All nodes of minimised DFSM must be distinguishable");
                cout << "Could not distinguish nodes "
                << node1->getName() << " and " << node2->getName() << endl;
                
                allNodesDistinguished = false;
            }
            
        }
        
    }
    
    if ( allNodesDistinguished ) {
        assert("TC-FSM-0010",
               true,
               "All nodes of minimised DFSM must be distinguishable");
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




vector<IOTrace> runAgainstRefModel(shared_ptr<Dfsm> refModel,
                                   IOListContainer& c) {
    
    shared_ptr<FsmPresentationLayer> pl = refModel->getPresentationLayer();
    
    auto iolCnt = c.getIOLists();
    
    // Register test cases in IO Traces
    vector<IOTrace> iotrLst;

    for ( auto lst : *iolCnt ) {
        
        shared_ptr<InputTrace> itr = make_shared<InputTrace>(lst,pl);
        IOTrace iotr = refModel->applyDet(*itr);
        iotrLst.push_back(iotr);
        
    }
    
    return iotrLst;
    
}

void runAgainstMutant(shared_ptr<Dfsm> mutant, vector<IOTrace>& expected) {
    
    for ( auto io : expected ) {
        
        InputTrace i = io.getInputTrace();
        
        if ( not mutant->pass(io) ) {
            cout << "FAIL: expected " << io << endl
            << "     : observed " << mutant->applyDet(i) << endl;
        }
        else {
            cout << "PASS: " << i << endl;
        }
        
    }
    
}

void wVersusT() {
    
    shared_ptr<Dfsm> refModel = make_shared<Dfsm>("FSBRTSX.csv","FSBRTS");

//    IOListContainer wTestSuite0 = refModel->wMethod(0);
//    IOListContainer wTestSuite1 = refModel->wMethod(1);
//    IOListContainer wTestSuite2 = refModel->wMethod(2);
//    IOListContainer wTestSuite3 = refModel->wMethod(3);
//    
  IOListContainer wpTestSuite0 = refModel->wpMethod(0);
//    IOListContainer wpTestSuite1 = refModel->wpMethod(1);
//    IOListContainer wpTestSuite2 = refModel->wpMethod(2);
//    IOListContainer wpTestSuite3 = refModel->wpMethod(3);
    
    //    IOListContainer tTestSuite = refModel->tMethod();
    
//    vector<IOTrace> expectedResultsW0 = runAgainstRefModel(refModel, wTestSuite0);
//    vector<IOTrace> expectedResultsW1 = runAgainstRefModel(refModel, wTestSuite1);
//    vector<IOTrace> expectedResultsW2 = runAgainstRefModel(refModel, wTestSuite2);
//    vector<IOTrace> expectedResultsW3 = runAgainstRefModel(refModel, wTestSuite3);
    vector<IOTrace> expectedResultsWp0 = runAgainstRefModel(refModel, wpTestSuite0);
//    vector<IOTrace> expectedResultsWp1 = runAgainstRefModel(refModel, wpTestSuite1);
//    vector<IOTrace> expectedResultsWp2  = runAgainstRefModel(refModel, wpTestSuite2);
//    vector<IOTrace> expectedResultsWp3 = runAgainstRefModel(refModel, wpTestSuite3);
//    vector<IOTrace> expectedResultsT = runAgainstRefModel(refModel, tTestSuite);


    
    for ( int i = 0; i < 10; i++ ) {
        
        cout << "Mutant No. " << (i+1) << ": " << endl;
        
        shared_ptr<Dfsm> mutant =
            make_shared<Dfsm>("FSBRTSX.csv","FSBRTS");
        mutant->createAtRandom();
        
//        runAgainstMutant(mutant,expectedResultsW0);
//        runAgainstMutant(mutant,expectedResultsW1);
//        runAgainstMutant(mutant,expectedResultsW2);
//        runAgainstMutant(mutant,expectedResultsW3);
        runAgainstMutant(mutant,expectedResultsWp0);
//        runAgainstMutant(mutant,expectedResultsWp1);
//        runAgainstMutant(mutant,expectedResultsWp2);
//        runAgainstMutant(mutant,expectedResultsWp3);
//        runAgainstMutant(mutant,expectedResultsT);

        
    }
    
    
    
    

    
}

//bool executeAdaptiveTest(const shared_ptr<Fsm>& spec, const shared_ptr<Fsm>& iut, bool& isReduction)
bool executeAdaptiveTest(
        const string prefix,
        const int numStates,
        const int numInput,
        const int numOutput,
        const size_t numOutFaults,
        const size_t numTransFaults,
        const unsigned int createRandomFsmSeed,
        const unsigned int createMutantSeed,
        const shared_ptr<FsmPresentationLayer>& pl,
        bool& isReduction)
{

    shared_ptr<FsmPresentationLayer> plCopy = make_shared<FsmPresentationLayer>(*pl);

    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating FSM.";
    shared_ptr<Fsm> spec = Fsm::createRandomFsm(prefix,
                                                numInput,
                                                numOutput,
                                                numStates,
                                                plCopy,
                                                true,
                                                createRandomFsmSeed);
    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating mutant.";

    shared_ptr<Fsm> iut;
    try {
        iut = spec->createMutant("mutant" + prefix,
                                 numOutFaults,
                                 numTransFaults,
                                 createMutantSeed);
    } catch (exception& e) {
        throw(e);
    }

    CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
    CLOG(INFO, logging::globalLogger) << "numInput: " << numInput + 1;
    CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutput + 1;
    CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
    CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
    CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed: " << createRandomFsmSeed;
    CLOG(INFO, logging::globalLogger) << "createMutantSeed: " << createMutantSeed;

    Fsm specMin = spec->minimise();
    Fsm iutMin = iut->minimise();

#ifdef ENABLE_DEBUG_MACRO
    const string dotPrefix = "../../../resources/adaptive-test/" + spec->getName() + "-";
    spec->toDot(dotPrefix + "spec");
    iut->toDot(dotPrefix + "iut");
    specMin.toDot(dotPrefix + "specMin");
    iutMin.toDot(dotPrefix + "iutMin");

    Fsm intersect = spec->intersect(*iut);
    intersect.toDot(dotPrefix + "intersect");

    Fsm product = spec->intersect(*iut);

    product.toDot(dotPrefix + "product");
#endif

    isReduction = !product.hasFailure();

    CLOG(INFO, logging::globalLogger) << "IUT is " + string((isReduction) ? "" : "NOT ") +
                                         "a reduction of the specification.";
    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Testing.";

    IOTraceContainer observedTraces;
    return isReduction == Fsm::adaptiveStateCounting(specMin, iutMin, static_cast<size_t>(iutMin.getMaxNodes()), observedTraces);
}

unsigned int getRandomSeed()
{
    std::random_device rd;
    return rd();
}

int getRandom(const int min, const int max, std::mt19937& gen)
{
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

int getRandom(const int max, std::mt19937& gen)
{
    return getRandom(0, max, gen);
}

int getRandom(std::mt19937& gen)
{
    std::uniform_int_distribution<int> dis;
    return dis(gen);
}

struct AdaptiveTestConfigDebug
{
    string prefix = "DEBUG";
    int numStates;
    int numInput;
    int numOutput;
    size_t numOutFaults;
    size_t numTransFaults;
    unsigned int createRandomFsmSeed;
    unsigned int createMutantSeed;
};

struct AdaptiveTestConfig
{

    // Required
    int numFsm = -1;
    int maxInput = -1;
    int maxOutput = -1;
    int maxStates = -1;
    int maxOutFaults = -1;
    int maxTransFaults = -1;

    // Optional
    int minStates = 2;
    int minInput = 2;
    int minOutput = 2;
    int minOutFaults = 1;
    int minTransFaults = 1;
    unsigned int seed = 0;
};

void adaptiveTest01(AdaptiveTestConfig& config)
{
    CLOG(INFO, logging::globalLogger) << "############## Adaptive Test 01 ##############";

    shared_ptr<FsmPresentationLayer> plTest =
    make_shared<FsmPresentationLayer>("../../../resources/adaptive-test-in.txt",
            + "../../../resources/adaptive-test-out.txt",
            + "../../../resources/adaptive-test-state.txt");

    if (config.numFsm < 0 ||
            config.maxInput < 0 ||
            config.maxOutput < 0 ||
            config.maxStates < 0 ||
            config.maxOutFaults < 0 ||
            config.maxTransFaults < 0)
    {
        CLOG(FATAL, logging::globalLogger) << "Please provide all required parameters.";
    }

    config.maxInput--;
    config.maxOutput--;
    config.maxStates--;
    config.minStates--;
    config.minInput--;
    config.minOutput--;

    if (config.seed == 0) {
        config.seed = getRandomSeed();
    }

    CLOG(INFO, logging::globalLogger) << "Seed: " << config.seed;

    std::mt19937 gen(config.seed);


    const int numberDigits = ((config.numFsm <= 1)? 1 : static_cast<int>(log10(config.numFsm)) + 1);


    const int diffInput = config.maxInput - config.minInput + 1;
    const int diffOutput = config.maxOutput - config.minOutput + 1;
    const int diffStates = config.maxStates - config.minStates + 1;

    if (diffInput <= 0 || diffOutput <= 0 || diffStates <= 0)
    {
        CLOG(FATAL, logging::globalLogger) << "Please check the test parameters.";
    }

    const int subLoopIterations = (config.numFsm < diffStates) ? 1 : 1 + ((config.numFsm - 1) / diffStates);

    CLOG(INFO, logging::globalLogger) << "numFsm: " << config.numFsm;
    CLOG(INFO, logging::globalLogger) << "maxInput: " << config.maxInput + 1;
    CLOG(INFO, logging::globalLogger) << "maxOutput: " << config.maxOutput + 1;
    CLOG(INFO, logging::globalLogger) << "maxStates: " << config.maxStates + 1;

    CLOG(INFO, logging::globalLogger) << "maxOutFaults: " << config.maxOutFaults;
    CLOG(INFO, logging::globalLogger) << "maxTransFaults: " << config.maxTransFaults;
    CLOG(INFO, logging::globalLogger) << "minOutFaults: " << config.minOutFaults;
    CLOG(INFO, logging::globalLogger) << "minTransFaults: " << config.minTransFaults;
    CLOG(INFO, logging::globalLogger) << "minInput: " << config.minInput + 1;
    CLOG(INFO, logging::globalLogger) << "minOutput: " << config.minOutput + 1;
    CLOG(INFO, logging::globalLogger) << "minStates: " << config.minStates + 1;
    CLOG(INFO, logging::globalLogger) << "seed: " << config.seed;

    CLOG(INFO, logging::globalLogger) << "diffStates: " << diffStates;
    CLOG(INFO, logging::globalLogger) << "subLoopIterations: " << subLoopIterations;

    TIMED_FUNC(timerObj);

    int executed = 0;
    int passed = 0;
    int i = 0;
    for (int numStates = config.minStates; numStates <= config.maxStates; ++numStates)
    {
        for (int j = 0; j < subLoopIterations; ++j)
        {

            shared_ptr<FsmPresentationLayer> plCopy = make_shared<FsmPresentationLayer>(*plTest);

            stringstream ss;
            ss << setw(numberDigits) << setfill('0') << i;
            string iteration = ss.str();
            logging::setLogfileSuffix(iteration);

            int numInput = getRandom(config.minInput, config.maxInput, gen);
            int numOutput = getRandom(config.minOutput, config.maxOutput, gen);
            size_t numOutFaults = static_cast<size_t>(getRandom(config.minOutFaults, config.maxOutFaults, gen));
            size_t numTransFaults = static_cast<size_t>(getRandom(config.minTransFaults, config.maxTransFaults, gen));

            const unsigned int createRandomFsmSeed = static_cast<unsigned int>(getRandom(gen));
            const unsigned int createMutantSeed = static_cast<unsigned int>(getRandom(gen));

            TIMED_SCOPE(timerBlkObj, "heavy-iter");
            CLOG(INFO, logging::globalLogger) << "------------------------------------------------------------------";
            CLOG(INFO, logging::globalLogger) << "i: " << iteration;

            bool isReduction;
            bool result;

            long durationMS;
            long durationMin;

            bool couldCreateMutant = false;
            while (!couldCreateMutant)
            {
                try {
                    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
                    result = executeAdaptiveTest(
                                iteration,
                                numStates,
                                numInput,
                                numOutput,
                                numOutFaults,
                                numTransFaults,
                                createRandomFsmSeed,
                                createMutantSeed,
                                plTest,
                                isReduction);
                    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                    durationMS = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    durationMin = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
                    couldCreateMutant = true;
                } catch (exception& e) {
                    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create mutant. Decreasing faults.";
                    int random = getRandom(0, 1, gen);
                    if (random == 0 && numOutFaults > 0)
                    {
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing output faults.";
                        --numOutFaults;
                        continue;
                    }
                    else if (numTransFaults > 0)
                    {
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing transition faults.";
                        --numTransFaults;
                        continue;
                    }
                    else
                    {
                        CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
                        CLOG(INFO, logging::globalLogger) << "numInput: " << numInput + 1;
                        CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutput + 1;
                        CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
                        CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
                        CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed: " << createRandomFsmSeed;
                        CLOG(INFO, logging::globalLogger) << "createMutantSeed: " << createMutantSeed;
                        CLOG(WARNING, logging::globalLogger) << "Could not create mutant. Skipping.";
                        break;
                    }
                }
            }
            if (!couldCreateMutant)
            {
                ++i;
                continue;
            }

            if (result)
            {
                ++passed;
            }

            assertOnFail("TC-AT-01-" + iteration, result);
            CLOG(INFO, logging::globalLogger) << "Calculation took " << durationMS << " ms (" << durationMin << " minutes).";

            ++executed;
            ++i;
        }
    }

    CLOG(INFO, logging::globalLogger) << "";
    CLOG(INFO, logging::globalLogger) << "#################### SUMMARY ####################";
    CLOG(INFO, logging::globalLogger) << "# Total tests  : " << executed;
    CLOG(INFO, logging::globalLogger) << "# Passed       : " << passed;
    CLOG(INFO, logging::globalLogger) << "# Failed       : " << executed - passed;
    CLOG(INFO, logging::globalLogger) << "# Not executed : " << i - executed;
    CLOG(INFO, logging::globalLogger) << "#################################################";
}

std::string getcwd() {
    std::string result(1024,'\0');
    while( getcwd(&result[0], result.size()) == 0) {
        if( errno != ERANGE ) {
          throw std::runtime_error(strerror(errno));
        }
        result.resize(result.size()*2);
    }
    result.resize(result.find('\0'));
    return result;
}



int main(int argc, char* argv[])
{
    START_EASYLOGGINGPP(argc, argv);
    logging::initLogging();

    CLOG(INFO, logging::globalLogger) << "############## Starting Application ##############";

    AdaptiveTestConfig config;
    config.numFsm = 1000;

    config.minInput = 1;
    config.maxInput = 10;

    config.minOutput = 1;
    config.maxOutput = 10;

    config.minStates = 1;
    config.maxStates = 10;

    config.minTransFaults = 0;
    config.maxTransFaults = 5;

    config.minOutFaults = 0;
    config.maxOutFaults = 5;

    config.seed = 4240152224;


    bool debug = true;
    if (debug)
    {
        CLOG(INFO, logging::globalLogger) << "############## Debugging ##############";

        AdaptiveTestConfigDebug debugConfig;
        debugConfig.numStates = 2;
        debugConfig.numInput = 5;
        debugConfig.numOutput = 5;
        debugConfig.numOutFaults = 1;
        debugConfig.numTransFaults = 0;
        debugConfig.createRandomFsmSeed = 430243142;
        debugConfig.createMutantSeed = 1425131424;

        shared_ptr<FsmPresentationLayer> plTest =
        make_shared<FsmPresentationLayer>("../../../resources/adaptive-test-in.txt",
                + "../../../resources/adaptive-test-out.txt",
                + "../../../resources/adaptive-test-state.txt");
        bool isReduction;
        bool result = executeAdaptiveTest(
                    debugConfig.prefix,
                    debugConfig.numStates - 1,
                    debugConfig.numInput - 1,
                    debugConfig.numOutput - 1,
                    debugConfig.numOutFaults,
                    debugConfig.numTransFaults,
                    debugConfig.createRandomFsmSeed,
                    debugConfig.createMutantSeed,
                    plTest,
                    isReduction);
        assertOnFail("TC-AT-DEBUG", result);
    }
    else
    {
        adaptiveTest01(config);
    }

	cout << endl << endl;
}



