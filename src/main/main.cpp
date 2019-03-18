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
#include "json/json.h"

#include "sets/HittingSet.h"
#include "sets/HsTreeNode.h"
#include <algorithm>
#include <cmath>
#include "fsm/PkTableRow.h"
#include "fsm/PkTable.h"
#include "fsm/DFSMTable.h"
#include "fsm/DFSMTableRow.h"
#include "fsm/OFSMTableRow.h"
#include "fsm/OFSMTable.h"
#include "fsm/FsmTransition.h"
#include "fsm/FsmLabel.h"

#include <tuple>

#include "main/Tests.h"


using namespace std;
using namespace Json;

void assertInconclusive(string tc, string comment = "") {
    
    string sVerdict("INCONCLUSIVE");
    cout << sVerdict << ": " << tc << " : " << comment <<  endl;
    
}

void fsmlib_assert(string tc, bool verdict, string comment = "") {
    
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
    fsmlib_assert("TC-DFSM-0001",vIn.size() == 4
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
    fsmlib_assert("TC-DFSM-0001",
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
    fsmlib_assert("TC-FSM-0001", 0 == system("diff f1.dot f1Copy.dot"),
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
    
    fsmlib_assert("TC-FSM-0001",
           t2.isEquivalentTo(t1),
           "Original FSM and its deep copy pass the same W-Method test suite");
    
    
    
}

void test3() {
    
    cout << "TC-FSM-0002 Show that createMutant() injects a fault into the original FSM" << endl;
    
    
    for ( size_t i = 0; i < 4; i++ ) {
        shared_ptr<FsmPresentationLayer> pl =
        make_shared<FsmPresentationLayer>();
        shared_ptr<Fsm> fsm = Fsm::createRandomFsm("F",5,5,8,pl,(unsigned)i);
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
        
        if ( iolc1.size() > 1000) {
            cout << "Skip this test case, since size is too big" << endl;
            continue;
        }
        
        TestSuite t1 = fsmMin.createTestSuite(iolc1);
        TestSuite t2 = fsmMutantMin.createTestSuite(iolc1);
        
        fsmlib_assert("TC-FSM-0002", not t2.isEquivalentTo(t1),
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
        
        fsmlib_assert("TC-FSM-0002",
               not t2wp.isEquivalentTo(t1wp),
               "Original FSM and mutant do not produce the same test suite results - tests are created by Wp-Method");
        
        fsmlib_assert("TC-FSM-0002",
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
        std::shared_ptr<Fsm> f = Fsm::createRandomFsm("F",5,5,10,pl,(unsigned)i);
        std::shared_ptr<Tree> sc = f->getStateCover();
        
        if ( sc->size() != (size_t)f->getMaxNodes() + 1 ) {
            cout << "Size of state cover: " << sc->size()
            << " Number of states in FSM: " << f->getMaxNodes() + 1 << endl;
            fsmlib_assert("TC-FSM-0004",
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
                fsmlib_assert("TC-FSM-0004",
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
        fsmlib_assert("TC-FSM-0004",
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
    
    fsmlib_assert("TC-FSM-0005",
           v.size() == 3,
           "For TC-FSM-0005.fsm, there are 3 classes of equivalent inputs.");
    
    fsmlib_assert("TC-FSM-0005",
           v[0].size() == 1 and v[0].find(0) != v[0].end(),
           "Class 0 only contains input 0.");
    
    fsmlib_assert("TC-FSM-0005",
           v[1].size() == 1 and v[1].find(1) != v[1].end(),
           "Class 1 only contains input 1.");
    
    fsmlib_assert("TC-FSM-0005",
           v[2].size() == 2 and
           v[2].find(2) != v[2].end() and
           v[2].find(3) != v[2].end(),
           "Class 2 contains inputs 2 and 3.");
    
    
    // Check FSM without any equivalent inputs
    fsm = make_shared<Fsm>("../../../resources/fsmGillA7.fsm",pl,"F");
    fsm->toDot("fsmGillA7");
    v = fsm->getEquivalentInputs();
    
    fsmlib_assert("TC-FSM-0005",
           v.size() == 3,
           "For fsmGillA7, there are 3 input classes.");
    
    bool ok = true;
    for ( size_t s=0; s < v.size() and ok; s++ ) {
        if ( v[s].size() != 1 or
            v[s].find((int)s) == v[s].end() ) {
            ok =false;
        }
    }
    
    fsmlib_assert("TC-FSM-0005",
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
        
        fsmlib_assert("TC-FSM-0009",
               uNodes.size() == 2 and (oldSize - d->size()) == 2,
               "All unreachable states have been removed");
    }
    else {
        fsmlib_assert("TC-FSM-0009",
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
                
                fsmlib_assert("TC-FSM-0010",
                       false,
                       "All nodes of minimised DFSM must be distinguishable");
                cout << "Could not distinguish nodes "
                << node1->getName() << " and " << node2->getName() << endl;
                
                allNodesDistinguished = false;
            }
            
        }
        
    }
    
    if ( allNodesDistinguished ) {
        fsmlib_assert("TC-FSM-0010",
               true,
               "All nodes of minimised DFSM must be distinguishable");
    }
    
}


void test10b() {
    
    cout << "TC-FSM-1010 Check correctness of Dfsm::minimise() with DFSM huang201711"
    << endl;
    
    shared_ptr<FsmPresentationLayer> pl =
        make_shared<FsmPresentationLayer>("../../../resources/huang201711in.txt",
                                          "../../../resources/huang201711out.txt",
                                          "../../../resources/huang201711state.txt");
    
    
    shared_ptr<Dfsm> d = make_shared<Dfsm>("../../../resources/huang201711.fsm",
                                           pl,
                                           "F");
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
                
                fsmlib_assert("TC-FSM-1010",
                       false,
                       "All nodes of minimised DFSM must be distinguishable");
                cout << "Could not distinguish nodes "
                << node1->getName() << " and " << node2->getName() << endl;
                
                allNodesDistinguished = false;
            }
            
        }
        
    }
    
    if ( allNodesDistinguished ) {
        fsmlib_assert("TC-FSM-1010",
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
    
    Dfsm gdcMin = gdc->minimise();
    
    gdcMin.toDot("GDC_MIN");
    
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
    
    fsmlib_assert("TC-GDC-0001",
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

void test11() {
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
                                                                            "../../../resources/garageOut.txt",
                                                                            "../../../resources/garageState.txt");
    
    shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/garage.fsm",pl,"GDC");
    
    
    gdc->toDot("GDC");
    
    Fsm gdcMin = gdc->minimise();
    
    gdcMin.toDot("GDC_MIN");
    
}

void test12() {
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
                                                                            "../../../resources/garageOut.txt",
                                                                            "../../../resources/garageState.txt");
    
    shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/garage.fsm",pl,"GDC");
    
    
    gdc->toDot("GDC");
    
    Dfsm gdcMin = gdc->minimise();
    
    gdcMin.toDot("GDC_MIN");
    
}

void test13() {
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
    shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/garage.fsm",pl,"GDC");
    
    
    gdc->toDot("GDC");
    
    Dfsm gdcMin = gdc->minimise();
    
    gdcMin.toDot("GDC_MIN");
    
}


void test14() {
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
    shared_ptr<Fsm> fsm = make_shared<Fsm>("../../../resources/NN.fsm",pl,"NN");
    
    fsm->toDot("NN");
    
    Fsm fsmMin = fsm->minimiseObservableFSM();
    
    fsmMin.toDot("NN_MIN");
    
}


void test15() {
    
    cout << "TC-DFSM-0015 Show that Fsm::transformToObservableFSM() produces an "
    << "equivalent observable FSM"
    << endl;
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
    shared_ptr<Fsm> nonObs = make_shared<Fsm>("../../../resources/nonObservable.fsm",pl,"NON_OBS");
    
    
    nonObs->toDot("NON_OBS");
    
    Fsm obs = nonObs->transformToObservableFSM();
    
    obs.toDot("OBS");
    
    fsmlib_assert("TC-DFSM-0015",
           obs.isObservable(),
           "Transformed FSM is observable");
    
    // Show that nonObs and obs have the same language.
    // We use brute force test that checks all traces of length n*m
    int n = (int)nonObs->size();
    int m = (int)obs.size();
    int theLen = n+m-1;
    
    IOListContainer allTrc = IOListContainer(nonObs->getMaxInput(),
                                             1,
                                             theLen,
                                             pl);
    
    shared_ptr<vector<vector<int>>> allTrcLst = allTrc.getIOLists();
    
    for ( auto trc : *allTrcLst ) {
        
        // Run the test case against both FSMs and compare
        // the (nondeterministic) result
        shared_ptr<InputTrace> iTr =
        make_shared<InputTrace>(trc,pl);
        OutputTree o1 = nonObs->apply(*iTr);
        OutputTree o2 = obs.apply(*iTr);
        
        if ( o1 != o2 ) {
            
            fsmlib_assert("TC-DFSM-0015",
                   o1 == o2,
                   "Transformed FSM has same language as original FSM");
            
            cout << "o1 = " << o1 << endl;
            cout << "o2 = " << o2 << endl;
            return;
            
        }
       
    }
    
}

void faux() {
    
    
    shared_ptr<FsmPresentationLayer> pl =
    make_shared<FsmPresentationLayer>("../../../resources/gillIn.txt",
                                      "../../../resources/gillOut.txt",
                                      "../../../resources/gillState.txt");
    
    shared_ptr<Dfsm> d = make_shared<Dfsm>("../../../resources/gill.fsm",
                                           pl,
                                           "G0");
    
    d->toDot("G0");
    
    d->toCsv("G0");
    
    Dfsm dMin = d->minimise();
    
    dMin.toDot("G0_MIN");
    
    
    
}

void test16() {
    
    shared_ptr<Dfsm> exp1 = nullptr;
    Reader jReader;
    Value root;
    stringstream document;
    ifstream inputFile("../../../resources/exp1.fsm");
    document << inputFile.rdbuf();
    inputFile.close();
    
    if ( jReader.parse(document.str(),root) ) {
        exp1 = make_shared<Dfsm>(root);
    }
    else {
        cerr << "Could not parse JSON model - exit." << endl;
        exit(1);
    }
    
    exp1->toDot("exp1");
    
    shared_ptr<Dfsm> exp2 = nullptr;
    Reader jReader2;
    Value root2;
    stringstream document2;
    ifstream inputFile2("../../../resources/exp2.fsm");
    document2 << inputFile2.rdbuf();
    inputFile2.close();
    
    if ( jReader2.parse(document2.str(),root) ) {
        exp2 = make_shared<Dfsm>(root);
    }
    else {
        cerr << "Could not parse JSON model - exit." << endl;
        exit(1);
    }
    
    exp2->toDot("exp2");
    
    Fsm prod = exp1->intersect(*exp2);
    
    cout << endl << "NEW PL STATES" << endl ;
    prod.getPresentationLayer()->dumpState(cout);
    
    
    prod.toDot("PRODexp1exp2");
    
    
    
    
}

// ===============================================================================================================================

///*
//	Calculates the set of transitions labels of outgoing transitions from given nodes set.
//*/
//unordered_set<FsmLabel> calcLblSet(std::unordered_set < shared_ptr<FsmNode>> &nodes) {
//	unordered_set<FsmLabel> lblSet;
//	for (auto n : nodes) {
//		for (auto tr : n->getTransitions()) {
//			lblSet.insert(*tr->getLabel());
//		}
//	}
//	return lblSet;
//}
//
///*
//	Checks if processed contains front.
//*/
////bool containsPair(std::vector<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> &processed,
////	std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>> &front) {
////	for (auto p : processed) {
////		if (p.first == front.first && p.second == front.second) {
////			return true;
////		}
////	}
////	return false;
////}
//
///*
//	Checks if wl contains pair.
//*/
////bool containsPair(std::deque<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> &wl,
////	std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>> &pair) {
////	for (auto p : wl) {
////		if (p.first == pair.first && p.second == pair.second) {
////			return true;
////		}
////	}
////	return false;
////}
//
//template<typename T>
//bool containsPair(T &lst, std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>> &pair)
//{
//	typename T::const_iterator it;
//	for (it = lst.begin(); it != lst.end(); ++it)
//	{                                      
//		if (it->first == pair.first && it->second == pair.second) {
//			return true;
//		}
//	}
//	return false;
//}
//
//
//
///*
//	Calculates and returns the set of target nodes reached from states contained in given node list with transitions labeled with given lbl.
//*/
//unordered_set<shared_ptr<FsmNode>> calcTargetNodes(unordered_set<shared_ptr<FsmNode>> &nodes, FsmLabel &lbl) {
//	unordered_set<shared_ptr<FsmNode>> tgtNds;
//	for (auto n : nodes) {
//		for (auto tr : n->getTransitions()) {
//			if (*tr->getLabel() == lbl) tgtNds.insert(tr->getTarget());
//		}
//	}
//	return tgtNds;
//}
//
///**
//	This function returns true iff L(q) = L(u). Otherwise the function returns false.
//*/
//bool ioEquivalenceCheck(const std::shared_ptr<FsmNode> q, const std::shared_ptr<FsmNode> u) {
//	// init wl with (q,u)
//	std::deque<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> wl;
//	std::unordered_set<std::shared_ptr<FsmNode>> l = std::unordered_set<std::shared_ptr<FsmNode>>{q};
//	std::unordered_set<std::shared_ptr<FsmNode>> r = std::unordered_set<std::shared_ptr<FsmNode>>{u};	
//	//wl.push_back(std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>(first, second));
//	wl.push_back({ l,r });
//	std::vector<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> processed;
//
//	while (not wl.empty()) {
//		auto front = wl.front();
//		wl.pop_front();
//		//calculate lblSet_l
//		unordered_set<FsmLabel> lblSet_l = calcLblSet(front.first);
//		//calculate lblSet_r
//		unordered_set<FsmLabel> lblSet_r = calcLblSet(front.second);
//		if (lblSet_l != lblSet_r) return false;
//
//		// insert front to processed if processed does not contain front already
//		if (not containsPair(processed, front)) {
//			processed.push_back(front);
//		}
//
//		for (auto lbl : lblSet_l) {
//			// calc tgtNds_l
//			unordered_set<shared_ptr<FsmNode>> tgtNds_l = calcTargetNodes(front.first, lbl);
//			// calc tgtNds_r
//			unordered_set<shared_ptr<FsmNode>> tgtNds_r = calcTargetNodes(front.second, lbl);
//			pair<unordered_set<shared_ptr<FsmNode>>, unordered_set<shared_ptr<FsmNode>>> pair{tgtNds_l, tgtNds_r};
//			if (not containsPair(wl, pair) and not containsPair(processed, pair)) {
//				wl.push_back(pair);
//			}
//		}
//	}
//	return true;
//}
//
//
//
///*
//	Checks if given fsm contains a pair of states with the same language. Returns true iff fsm contains states q, q' with
//	q != q' and L(q) = L(q').
//*/
//bool hasEquivalentStates(const Fsm &fsm) {
//	for (int c1 = 0; c1 < fsm.getNodes().size(); c1++) {
//		for (int c2 = c1 + 1; c2 < fsm.getNodes().size(); c2++) {
//			if (ioEquivalenceCheck(fsm.getNodes().at(c1), fsm.getNodes().at(c2))) {
//				return true;
//			}
//		}
//	}
//	return false;
//}
//
///*
//	Check if fsm1 and fsm2 have the same structure (same labeled transitions between nodes with the same indices).
//	This method can be used to check if some Fsms structure was changed by some method call.
//	Its faster than checking for isomorphism, because of restrictiveness.
//	In other words this function tests if fsm1 and fsm2 have identical nodes lists.
//*/
//bool checkForEqualStructure(const Fsm &fsm1, const Fsm &fsm2) {
//	// fsm1 and fsm2 need to be the same size
//	if (fsm1.getNodes().size() != fsm2.getNodes().size()) return false;
//	for (int i = 0; i < fsm1.getNodes().size(); ++i) {
//		// each node should have the same number of transitions
//		if (fsm1.getNodes().at(i)->getTransitions().size() != fsm2.getNodes().at(i)->getTransitions().size()) return false;
//		for (int j = 0; j < fsm1.getNodes().at(i)->getTransitions().size(); ++j) {
//			auto fsm1Tr = fsm1.getNodes().at(i)->getTransitions().at(j);
//			auto fsm2Tr = fsm2.getNodes().at(i)->getTransitions().at(j);
//			// compare fsm1Tr and fsm2Tr
//			if (fsm1Tr->getSource()->getId() != fsm2Tr->getSource()->getId()
//				|| fsm1Tr->getTarget()->getId() != fsm2Tr->getTarget()->getId()
//				|| not (*fsm1Tr->getLabel() == *fsm2Tr->getLabel())) {
//				return false;
//			}
//		}
//	}
//	return true;
//}
//
///*
//	Returns a set of all the states of the given Fsm that are reachable from the initial state of that Fsm.
//*/
//unordered_set<shared_ptr<FsmNode>> getReachableStates(const Fsm &fsm) {
//	unordered_set<shared_ptr<FsmNode>> reached{fsm.getInitialState()};
//	deque<shared_ptr<FsmNode>> wl{ fsm.getInitialState() };
//	while (not wl.empty()) {
//		shared_ptr<FsmNode> q = wl.front();
//		wl.pop_front();
//		for (auto tr : q->getTransitions()) {
//			// add reached target to wl if this target wasn't reached before
//			if (reached.insert(tr->getTarget()).second) {
//				wl.push_back(tr->getTarget());
//			}
//		}
//	}
//	return reached;
//}
//
///*
//	Returns true iff nodes[i].id == i for all 0 <= i < fsm.getNodes().size()
//*/
//bool checkNodeIds(const Fsm &fsm) {
//	for (size_t i = 0; i < fsm.getNodes().size(); ++i) {
//		if (fsm.getNodes().at(i)->getId() != i) return false;
//	}
//	return true;
//}
//
///*
//	Returns true iff fsm.getNodes() contains given node pointer.
//*/
//bool contains(const Fsm &fsm, const shared_ptr<FsmNode> node) {
//	for (auto n : fsm.getNodes()) {
//		if (n == node) return true;
//	}
//	return false;
//}
//
///*
//	Returns true iff fsm.getNodes() contains any of the given node pointers in nodes.
//*/
//bool contains(const Fsm &fsm, const unordered_set<shared_ptr<FsmNode>> nodes) {
//	for (auto n : nodes) {
//		if (contains(fsm, n)) return true;
//	}
//	return false;
//}
//
//
///*
//	Checks the transitions and return false iff any transitions hurts the invariant of Fsm.
//*/
//bool checkAllTransitions(const Fsm &fsm) { 
//	for (auto n : fsm.getNodes()) {
//		for (auto tr : n->getTransitions()) {
//			if (tr == nullptr || tr->getLabel() == nullptr || tr->getLabel()->getInput() > fsm.getMaxInput()
//				|| tr->getLabel()->getOutput() > fsm.getMaxOutput()  || tr->getLabel()->getInput() < 0 || tr->getLabel()->getOutput() < 0
//				|| tr->getSource() != n
//				|| not contains(fsm, tr->getTarget())) {
//				return false;
//			}
//		}
//	}
//	return true;
//}
//
///*
//	This function checks the Fsm class invariant for the given Fsm object.
//*/
//bool checkFsmClassInvariant(const Fsm &fsm) {
//	if (fsm.getMaxInput() < 0) return false;
//	if (fsm.getMaxOutput() < 0) return false;	
//	if (fsm.getNodes().size() < 1) return false;	
//	if (not checkNodeIds(fsm)) return false;
//	if (contains(fsm, nullptr)) return false;		
//	if (not checkAllTransitions(fsm)) return false;	
//	if (fsm.getMaxState() != fsm.getNodes().size() - 1) return false;	
//	if (not(0 <= fsm.getInitStateIdx() and fsm.getInitStateIdx() <= fsm.getMaxState())) return false;	
//	return true;
//}
//
//bool checkDfsmClassInvariant(Dfsm &dfsm) {
//	return checkFsmClassInvariant(dfsm) and dfsm.isDeterministic();
//}
//
///*
//	Checks if the given fsm is initial connected.
//*/
//bool isInitialConnected(const Fsm &fsm) {
//	auto reachable = getReachableStates(fsm);
//	auto nodes = fsm.getNodes();
//	unordered_set < shared_ptr<FsmNode> >nodeSet(nodes.cbegin(), nodes.cend());
//	return reachable == nodeSet;
//}
//
///*
//	Checks if unreachableNodesAfter contains all elements from unreachableNodesBefore and unreachable but no other element.
//*/
//bool checkUnreachableNodesList(const vector<shared_ptr<FsmNode>> &unreachableNodesBefore, const vector<shared_ptr<FsmNode>> &unreachableNodesAfter,
//	unordered_set<shared_ptr<FsmNode>> &unreachable) {
//	// check the size
//	if (unreachableNodesAfter.size() != unreachableNodesBefore.size() + unreachable.size()) return false;
//
//	// check if each node in unreachableNodesBefore is in unreachableNodesAfter	
//	for (auto n : unreachableNodesBefore) {
//		bool found = false;
//		for (auto n2 : unreachableNodesAfter) {
//			if (n == n2) {
//				found = true;
//				break;
//			}
//		}
//		if (not found) return false;
//	}
//
//	// check if each node in unreachable is in unreachableNodesAfter
//	for (auto n : unreachable) {
//		bool found = false;
//		for (auto n2 : unreachableNodesAfter) {
//			if (n == n2) {
//				found = true;
//				break;
//			}
//		}
//		if (not found) return false;
//	}
//
//	// check if each node of unreachableNodesAfter is in unreachableNodesBefore or unreachable
//	for (auto n : unreachableNodesAfter) {
//		bool found = false;
//		for (auto n2 : unreachableNodesBefore) {
//			if (n == n2) {
//				found = true;
//				break;
//			}
//		}
//		if (not found) {
//			for (auto n2 : unreachable) {
//				if (n == n2) {
//					found = true;
//					break;
//				}
//			}
//		}
//		if (not found) return false;
//	}
//	return true;
//}
//
///*
//	This function creates and returns a randomly created Dfsm object. It is needed because the Dfsm Constructor that create randomized Dfsms does not set
//	the maxState member correctly, so no Dfsm Object created would fullfill the invariant.
//*/
////Dfsm createRandomDfsm(const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer> presentationLayer) {
////	Dfsm dfsm(fsmName, maxNodes, maxInput, maxOutput, presentationLayer);
////	dfsm.setMaxState(dfsm.getNodes().size() - 1);
////	return dfsm;
////}
//
////TODO Other random creation wrapper (createRandomMinimisedFsm) if needed (some may be needed because constructor or some transformation method
//// does not set maxState correctly, which matters if Fsm::createMutant() will be used on the created fsm)
//
//
////void testFsmClassInvariant() {
////	for (int i = 0; i < 30; i++) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////		//auto fsm2 = fsm->minimise();
////		fsmlib_assert("TC", checkFsmClassInvariant(*fsm), "Random FSM fullfills invariant");
////		//fsmlib_assert("TC", checkFsmClassInvariant(fsm2), "Minimised Random FSM fullfills invariant");
////	}
////
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> q1 = make_shared<FsmNode>(1, pl);
////		shared_ptr<FsmNode> q2 = make_shared<FsmNode>(2, pl);
////		shared_ptr<FsmNode> q3 = make_shared<FsmNode>(3, pl);
////
////		shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr0);
////		shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(q0, q2, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr1);
////		shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(q2, q1, make_shared<FsmLabel>(2, 2, pl));
////		q2->addTransition(tr2);
////		Fsm fsm("M", 2, 2, { q0,q1,q2,q3 }, pl);
////		fsmlib_assert("TC", checkFsmClassInvariant(fsm), "FSM fullfills invariant");
////	}
////
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		Fsm fsm("M", 0, 0, { q0 }, pl);
////		fsmlib_assert("TC", checkFsmClassInvariant(fsm), "FSM fullfills invariant");
////	}
////
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> q1 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(1, 0, pl));
////		q0->addTransition(tr0);
////		Fsm fsm("M", 0, 0, { q0, q1 }, pl);
////		fsmlib_assert("TC", not checkFsmClassInvariant(fsm), "FSM does not fullfill invariant if some input is greater than maxInput");
////	}
////
////
////}
//
////void testCheckDfsmClassInvariant() {
////	for (int i = 0; i < 10; ++i) {
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		Dfsm dfsm = createRandomDfsm("M", 10, 3, 4, pl);
////		fsmlib_assert("TC", checkDfsmClassInvariant(dfsm), "Random DFSM fullfills invariant.");
////	}
////}
//
//
////void testIOEquivalenceCheck() {
////	{
////	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////	shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////	shared_ptr<FsmNode> q1 = make_shared<FsmNode>(1, pl);
////
////	shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(0, 1, pl));
////	q0->addTransition(tr0);
////	shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(q1, q1, make_shared<FsmLabel>(1, 1, pl));
////	q1->addTransition(tr1);
////
////	shared_ptr<FsmNode> u0 = make_shared<FsmNode>(0, pl);
////	shared_ptr<FsmNode> u1 = make_shared<FsmNode>(1, pl);
////	shared_ptr<FsmNode> u2 = make_shared<FsmNode>(2, pl);
////
////	shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(u0, u1, make_shared<FsmLabel>(0, 1, pl));
////	u0->addTransition(tr2);
////	shared_ptr<FsmTransition> tr3 = make_shared<FsmTransition>(u1, u1, make_shared<FsmLabel>(1, 1, pl));
////	u1->addTransition(tr3);
////	shared_ptr<FsmTransition> tr4 = make_shared<FsmTransition>(u0, u2, make_shared<FsmLabel>(1, 1, pl));
////	u0->addTransition(tr4);
////
////	cout << ioEquivalenceCheck(q0, u0) << endl;
////	}
////
////	for (int i = 0; i < 30; i++) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////		cout << "fsm.size: " << fsm->size() << endl;
////		auto fsm2 = fsm->minimise();
////		cout << "fsm2.size: " << fsm2.size() << endl;
////		cout << ioEquivalenceCheck(fsm->getInitialState(), fsm2.getInitialState()) << endl;
////	}
////
////	cout << "-------------------------------" << endl;
////
////	for (int i = 0; i < 5; i++) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////		cout << "fsm.size: " << fsm->size() << endl;
////		auto fsm2 = fsm->minimise();
////		cout << "fsm2.size: " << fsm2.size() << endl;
////		if (hasEquivalentStates(fsm2)) {
////			cout << "FAULT" << endl;
////		}
////	}
////
////	cout << "-------------------------------" << endl;
////
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> q1 = make_shared<FsmNode>(1, pl);
////
////		shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(0, 1, pl));
////		q0->addTransition(tr0);
////		shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(q1, q1, make_shared<FsmLabel>(1, 1, pl));
////		q1->addTransition(tr1);
////		shared_ptr<FsmTransition> tr5 = make_shared<FsmTransition>(q0, q0, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr5);
////
////		shared_ptr<FsmNode> u0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> u1 = make_shared<FsmNode>(1, pl);
////		shared_ptr<FsmNode> u2 = make_shared<FsmNode>(2, pl);
////
////		shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(u0, u1, make_shared<FsmLabel>(0, 1, pl));
////		u0->addTransition(tr2);
////		shared_ptr<FsmTransition> tr3 = make_shared<FsmTransition>(u1, u1, make_shared<FsmLabel>(1, 1, pl));
////		u1->addTransition(tr3);
////		shared_ptr<FsmTransition> tr4 = make_shared<FsmTransition>(u0, u2, make_shared<FsmLabel>(1, 1, pl));
////		u0->addTransition(tr4);
////
////		cout << ioEquivalenceCheck(q0, u0) << endl;
////	}
////
////	cout << "-------------------------------" << endl;
////
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> q1 = make_shared<FsmNode>(1, pl);
////		shared_ptr<FsmNode> q2 = make_shared<FsmNode>(2, pl);
////		shared_ptr<FsmNode> q3 = make_shared<FsmNode>(3, pl);
////
////		shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr0);
////		shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(q0, q2, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr1);
////		shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(q2, q1, make_shared<FsmLabel>(2, 2, pl));
////		q2->addTransition(tr2);
////		shared_ptr<FsmTransition> tr3 = make_shared<FsmTransition>(q1, q3, make_shared<FsmLabel>(0, 1, pl));
////		q1->addTransition(tr3);
////		shared_ptr<FsmTransition> tr4 = make_shared<FsmTransition>(q3, q0, make_shared<FsmLabel>(0, 0, pl));
////		q3->addTransition(tr4);
////
////		shared_ptr<FsmNode> u0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> u1 = make_shared<FsmNode>(1, pl);
////		shared_ptr<FsmNode> u2 = make_shared<FsmNode>(2, pl);
////		shared_ptr<FsmNode> u3 = make_shared<FsmNode>(3, pl);
////
////		shared_ptr<FsmTransition> tr5 = make_shared<FsmTransition>(u0, u1, make_shared<FsmLabel>(1, 1, pl));
////		u0->addTransition(tr5);
////		shared_ptr<FsmTransition> tr6 = make_shared<FsmTransition>(u1, u2, make_shared<FsmLabel>(2, 2, pl));
////		u1->addTransition(tr6);
////		shared_ptr<FsmTransition> tr7 = make_shared<FsmTransition>(u1, u3, make_shared<FsmLabel>(0, 1, pl));
////		u1->addTransition(tr7);
////		shared_ptr<FsmTransition> tr8 = make_shared<FsmTransition>(u2, u3, make_shared<FsmLabel>(0, 1, pl));
////		u2->addTransition(tr8);
////		shared_ptr<FsmTransition> tr9 = make_shared<FsmTransition>(u3, u0, make_shared<FsmLabel>(0, 0, pl));
////		u3->addTransition(tr9);
////		shared_ptr<FsmTransition> tr10 = make_shared<FsmTransition>(u2, u1, make_shared<FsmLabel>(0, 1, pl));
////		u2->addTransition(tr10);
////		
////
////		cout << ioEquivalenceCheck(q0, u0) << endl;
////	}
////
////}
//
///*
//	This function is used to test the checkForEqualStructure function
//*/
////void testCheckForEqualStructure() {
////	cout << "testCheckForEqualStructure" << endl;
////
////	cout << "positive cases:" << endl;
////	for (int i = 0; i < 10; ++i) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////
////		cout << checkForEqualStructure(*fsm, *fsm) << endl;
////
////		Fsm copy = Fsm(*fsm);
////		cout << checkForEqualStructure(*fsm, copy) << endl;
////		Fsm ofsm = fsm->transformToObservableFSM();
////
////		cout << checkForEqualStructure(*fsm, copy) << endl;
////
////		Fsm copy2 = Fsm(ofsm);
////		Fsm minOfsm = ofsm.minimiseObservableFSM();
////		cout << checkForEqualStructure(copy2, ofsm) << endl;
////
////		cout << "-----------------------------------" << endl;
////	}
////
////	cout << "negative cases:" << endl;
////
////	for (int i = 0; i < 10; ++i) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////
////		auto mutant = fsm->createMutant("mutant", 1, 1);
////		
////		cout << checkForEqualStructure(*fsm, *mutant) << endl;
////
////		cout << "-----------------------------------" << endl;
////	}
////
////}
//
///*
//	This function is used to test the getReachableStates function
//*/
////void testGetReachableStates() {
////	for (int i = 0; i < 10; ++i) {
////		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
////
////		unordered_set<shared_ptr<FsmNode>> reachable = getReachableStates(*fsm);
////
////		vector<shared_ptr<FsmNode>> nodes = fsm->getNodes();
////		unordered_set<shared_ptr<FsmNode>> nodeSet(nodes.begin(), nodes.end());
////
////		fsmlib_assert("TC", reachable == nodeSet, "getReachableStates returns set containing each reachable state");
////	}
////	
////	{
////		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////		shared_ptr<FsmNode> q0 = make_shared<FsmNode>(0, pl);
////		shared_ptr<FsmNode> q1 = make_shared<FsmNode>(1, pl);
////		shared_ptr<FsmNode> q2 = make_shared<FsmNode>(2, pl);
////		shared_ptr<FsmNode> q3 = make_shared<FsmNode>(3, pl);
////
////		shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(q0, q1, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr0);
////		shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(q0, q2, make_shared<FsmLabel>(1, 1, pl));
////		q0->addTransition(tr1);
////		shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(q2, q1, make_shared<FsmLabel>(2, 2, pl));
////		q2->addTransition(tr2);
////		//shared_ptr<FsmTransition> tr3 = make_shared<FsmTransition>(q1, q3, make_shared<FsmLabel>(0, 1, pl));
////		//q1->addTransition(tr3);
////		/*shared_ptr<FsmTransition> tr4 = make_shared<FsmTransition>(q3, q0, make_shared<FsmLabel>(0, 0, pl));
////		q3->addTransition(tr4);*/
////
////		Fsm fsm("M",2,2,{q0,q1,q2,q3},pl);
////		unordered_set<shared_ptr<FsmNode>> reachable = getReachableStates(fsm);
////
////		vector<shared_ptr<FsmNode>> nodes = fsm.getNodes();
////		unordered_set<shared_ptr<FsmNode>> nodeSet(nodes.begin(), nodes.end());
////
////		fsmlib_assert("TC", reachable != nodeSet, "getReachableStates returns set containing only reachable state");
////	}
////}
//
//// ====================================================================================================
//// Prüfverfahren "Spracherhaltende FSM-Transformationen"
//
///**
// * Test function: Fsm::removeUnreachableNodes()
// */
//void testRemoveUnreachableNodes(Fsm &m1) {
//	// get copy of m1
//	Fsm copyOfM1 = Fsm(m1);
//	vector<shared_ptr<FsmNode>> unreachableNodes;
//
//	// determine set of unreachable nodes in m1
//	auto reachable = getReachableStates(m1);
//	unordered_set<shared_ptr<FsmNode>> unreachable;
//	for (auto n : m1.getNodes()) {
//		if (reachable.count(n) == 0) unreachable.insert(n);
//	}
//
//	//vector<shared_ptr<FsmNode>> copyOfUnreachableNodes(unreachableNodes.begin(), unreachableNodes.end());
//
//	// use algorithm to transform m1
//	bool b = m1.removeUnreachableNodes(unreachableNodes);
//
//	// first check invariant of m1
//	bool invariantViolation = not checkFsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m1	
//	fsmlib_assert("TC", isInitialConnected(m1), "Result of removeUnreachableNodes() is initial connected");
//	fsmlib_assert("TC", not contains(m1, unreachable), "Resulting FSM of removeUnreachableNodes() contains none of the nodes that were unreachable before.");
//
//	// check if L(m1) = L(copyOfM1)
//	fsmlib_assert("TC", ioEquivalenceCheck(m1.getInitialState(), copyOfM1.getInitialState()), "removeUnreachableNodes() does not change language of the FSM");
//
//	unordered_set<shared_ptr<FsmNode>> unreachableNodesSet{unreachableNodes.cbegin(), unreachableNodes.cend()};
//	// check b and unreachableNodes
//	fsmlib_assert("TC", (b and (not unreachable.empty())) || (not b and unreachable.empty()), "removeUnreachableNodes() returns true iff FSM contains some unreachable node");
//	fsmlib_assert("TC", (unreachableNodes.size() == unreachable.size()) and (unreachableNodesSet == unreachable), "unreachableNodes contains each unreachable node that was removed");
//	/*fsmlib_assert("TC", checkUnreachableNodesList(copyOfUnreachableNodes, unreachableNodes, unreachable), "unreachableNodes contains all unreachable nodes that were removed and all nodes from before");*/
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkFsmClassInvariant(m1), "FSM still fullfills class invariants after transformation");
//}
//
/////**
//// * Test function: Fsm::removeUnreachableNodes()
//// */
////void testRemoveUnreachableNodes(Fsm &m1, vector<shared_ptr<FsmNode>> &unreachableNodes) {
////	// determine set of unreachable nodes in m1
////	auto reachable = getReachableStates(m1);
////	unordered_set<shared_ptr<FsmNode>> unreachable;
////	for (auto n : m1.getNodes()) {
////		if (reachable.count(n) == 0) unreachable.insert(n);
////	}
////
////	// get copy of m1 and unreachableNodes
////	Fsm copyOfM1 = Fsm(m1);
////	vector<shared_ptr<FsmNode>> copyOfUnreachableNodes(unreachableNodes.begin(), unreachableNodes.end());
////
////	// use algorithm to transform m1
////	bool b = m1.removeUnreachableNodes(unreachableNodes);
////
////	// first check invariant of m1
////	bool invariantViolation = not checkFsmClassInvariant(m1);
////	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
////	// stop test execution at this point if invariant of m does not hold anymore
////	if (invariantViolation) return;
////
////	// check properties of m1
////	fsmlib_assert("TC", not contains(m1,unreachable), "Resulting FSM of removeUnreachableNodes() contains none of the nodes that were unreachable before.");
////	fsmlib_assert("TC", isInitialConnected(m1), "Result of removeUnreachableNodes() is initial connected");
////
////	// check if L(m1) = L(copyOfM1)
////	fsmlib_assert("TC", ioEquivalenceCheck(m1.getInitialState(), copyOfM1.getInitialState()), "removeUnreachableNodes() does not change language of the FSM");
////
////	// check b and unreachableNodes
////	fsmlib_assert("TC", (b and (not unreachable.empty())) || (not b and unreachable.empty()), "removeUnreachableNodes() returns true iff FSM contains some unreachable node");
////	fsmlib_assert("TC", checkUnreachableNodesList(copyOfUnreachableNodes, unreachableNodes, unreachable), "unreachableNodes contains all unreachable nodes that were removed and all nodes from before");
////
////	//// check unexpected side effects
////	//fsmlib_assert("TC", checkFsmClassInvariant(m1), "FSM still fullfills class invariants after transformation");
////}
//
///**
// * Test function: Fsm::removeUnreachableNodes()
// */
//void testRemoveUnreachableNodes(Dfsm &m1, vector<shared_ptr<FsmNode>> &unreachableNodes) {
//	// determine set of unreachable nodes in m1
//	auto reachable = getReachableStates(m1);
//	unordered_set<shared_ptr<FsmNode>> unreachable;
//	for (auto n : m1.getNodes()) {
//		if (reachable.count(n) == 0) unreachable.insert(n);
//	}
//
//	// get copy of m1 and unreachableNodes
//	Dfsm copyOfM1 = m1;	
//	vector<shared_ptr<FsmNode>> copyOfUnreachableNodes(unreachableNodes.begin(), unreachableNodes.end());
//
//	// use algorithm to transform m1
//	bool b = m1.removeUnreachableNodes(unreachableNodes);
//
//	// first check invariant of m1
//	bool invariantViolation = not checkDfsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m1
//	fsmlib_assert("TC", not contains(m1, unreachable), "Resulting FSM of removeUnreachableNodes() contains none of the nodes that were unreachable before.");
//	fsmlib_assert("TC", isInitialConnected(m1), "Result of removeUnreachableNodes() is initial connected");
//
//	// check if L(m1) = L(copyOfM1)
//	fsmlib_assert("TC", ioEquivalenceCheck(m1.getInitialState(), copyOfM1.getInitialState()), "removeUnreachableNodes() does not change language of the FSM");
//
//	// check b and unreachableNodes
//	fsmlib_assert("TC", (b and (not unreachable.empty())) || (not b and unreachable.empty()), "removeUnreachableNodes() returns true iff FSM contains some unreachable node");
//	fsmlib_assert("TC", checkUnreachableNodesList(copyOfUnreachableNodes, unreachableNodes, unreachable), "unreachableNodes contains all unreachable nodes that were removed and all nodes from before");
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkDfsmClassInvariant(m1), "DFSM still fullfills class invariants after transformation");
//}
//
///**
// * Test function: Fsm::transformToObservableFSM()
// */
//void testTransformToObservableFSM(Fsm &m1) {
//	// get copy of m1
//	Fsm copyOfM1 = Fsm(m1);
//
//	// use algorithm to transform m1
//	Fsm m2 = m1.transformToObservableFSM();
//
//	// first check invariant of m2
//	bool invariantViolation = not checkFsmClassInvariant(m2);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M2 after transformation");
//	// stop test execution at this point if invariant of m2 does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m2
//	fsmlib_assert("TC", m2.isObservable(), "M2 is observable after transformToObservable()");
//
//	// check if L(m1) = L(m2)
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), "transformToObservable() does not change the language");
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkFsmClassInvariant(m1), "M1 still fullfills class invariants after transformation");
//	//fsmlib_assert("TC", checkFsmClassInvariant(m2), "M2 still fullfills class invariants after transformation");
//
//	// check invariant of m1
//	invariantViolation = not checkFsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m1 does not hold anymore
//	if (invariantViolation) return;
//
//	fsmlib_assert("TC", checkForEqualStructure(m1, copyOfM1), "M1 was not changed by algorithm");
//}
//
///**
// * Test function: Dfsm::minimise()
// */
//void testMinimise_Dfsm(Dfsm &m1) {
//	// get copy of m1
//	Dfsm copyOfM1 = m1;//Dfsm(m1);
//
//	// use algorithm to transform m1
//	Dfsm m2 = m1.minimise();
//
//	// first check invariant of m2
//	bool invariantViolation = not checkDfsmClassInvariant(m2);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M2 after transformation");
//	// stop test execution at this point if invariant of m2 does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m2
//	fsmlib_assert("TC", m2.isDeterministic(), "M2 is deterministic after minimise()");
//	fsmlib_assert("TC", isInitialConnected(m2), "M2 is initial connected after minimise()");
//	fsmlib_assert("TC", not hasEquivalentStates(m2), "M2 has no equivalent states after minimise()");
//
//	// check if L(m1) = L(m2)
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), "minimise() does not change the language");
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkDfsmClassInvariant(m1), "M1 still fullfills class invariants after transformation");
//	//fsmlib_assert("TC", checkDfsmClassInvariant(m2), "M2 still fullfills class invariants after transformation");
//
//	// check invariant of m1
//	invariantViolation = not checkDfsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m1 does not hold anymore
//	if (invariantViolation) return;
//	fsmlib_assert("TC", isInitialConnected(m1), "M1 is initial connected after minimise()");
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m1.getInitialState()), "Language of M1 was not changed by algorithm");
//}
//
///**
// * Test function: Fsm::minimiseObservableFSM()
// */
//void testMinimiseObservableFSM(Fsm &m1) {
//	// get copy of m1
//	Fsm copyOfM1 = Fsm(m1);
//
//	// use algorithm to transform m1
//	Fsm m2 = m1.minimiseObservableFSM();
//
//	// first check invariant of m2
//	bool invariantViolation = not checkFsmClassInvariant(m2);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M2 after transformation");
//	// stop test execution at this point if invariant of m2 does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m2
//	fsmlib_assert("TC", m2.isObservable(), "M2 is observable after minimiseObservable()");
//	fsmlib_assert("TC", not hasEquivalentStates(m2), "M2 has no equivalent states after minimise()");
//
//	// check if L(m1) = L(m2)
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), "minimiseObservableFSM() does not change the language");
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkFsmClassInvariant(m1), "M1 still fullfills class invariants after transformation");
//	//fsmlib_assert("TC", checkFsmClassInvariant(m2), "M2 still fullfills class invariants after transformation");
//
//	// check invariant of m1
//	invariantViolation = not checkFsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m1 does not hold anymore
//	if (invariantViolation) return;
//	fsmlib_assert("TC", checkForEqualStructure(m1, copyOfM1), "M1 was not changed by algorithm");
//}
//
///**
// * Test function: Fsm::minimise()
// */
//void testMinimise_Fsm(Fsm &m1) {
//	// get copy of m1
//	Fsm copyOfM1 = Fsm(m1);
//
//	// use algorithm to transform m1
//	Fsm m2 = m1.minimise();
//	
//	// first check invariant of m2
//	bool invariantViolation = not checkFsmClassInvariant(m2);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M2 after transformation");
//	// stop test execution at this point if invariant of m2 does not hold anymore
//	if (invariantViolation) return;
//
//	// check properties of m2	
//	fsmlib_assert("TC", m2.isObservable(), "M2 is observable after minimise()");
//	fsmlib_assert("TC", not hasEquivalentStates(m2), "M2 has no equivalent states after minimise()");
//	fsmlib_assert("TC", isInitialConnected(m2), "M2 is initial connected after minimise()");
//
//	// check if L(m1) = L(m2)
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), "minimise() does not change the language");
//
//	//// check unexpected side effects
//	//fsmlib_assert("TC", checkFsmClassInvariant(m1), "M1 still fullfills class invariants after transformation");
//	//fsmlib_assert("TC", checkFsmClassInvariant(m2), "M2 still fullfills class invariants after transformation");
//
//	// check invariant of m1
//	invariantViolation = not checkFsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolation, "class invariant holds for M1 after transformation");
//	// stop test execution at this point if invariant of m1 does not hold anymore
//	if (invariantViolation) return;
//	
//	// check unexpected sideeffects
//	fsmlib_assert("TC", isInitialConnected(m1), "M1 is initial connected after minimise()");
//	fsmlib_assert("TC", ioEquivalenceCheck(copyOfM1.getInitialState(), m1.getInitialState()), "Language of M1 was not changed by algorithm");
//}
//
//template <typename T>
//struct FsmTransformationTestCase {
//	string id;
//	shared_ptr<T> m;
//};
//
//template <typename T>
//void parseFsmTransformationTSFile(const string &testSuitePath, vector<FsmTransformationTestCase<T>> &testSuite) {
//	string fname = testSuitePath;
//	//vector<FsmTransformationTestCase> testSuite;
//	ifstream inputFile(fname);
//	if (inputFile.is_open())
//	{
//		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//		string line;
//		while (getline(inputFile, line))
//		{
//			vector<string> lineContent;
//			stringstream ss(line);
//			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
//			FsmTransformationTestCase<T> tc;
//			tc.id = lineContent.at(0);
//			string m1Path = lineContent.at(1);
//			if (m1Path == "../../../resources/TestSuites/FSM-Transformations/FSM002.fsm") {
//				shared_ptr<FsmNode> n = make_shared<FsmNode>(0, pl);
//				//vector<shared_ptr<FsmNode>> lst{ n };
//				T m("M", 0, 0, { n }, pl);
//				tc.m = make_shared<T>(m);
//				testSuite.push_back(tc);
//			}
//			else {
//				shared_ptr<T> m = make_shared<T>(m1Path, pl, "M");
//				tc.m = m;
//				testSuite.push_back(tc);
//			}
//		}
//		inputFile.close();
//
//		//return make_shared<vector<FsmTransformationTestCase>>(testSuite);
//	}
//	else
//	{
//		cout << "Unable to open input file" << endl;
//		exit(EXIT_FAILURE);
//	}
//}
//
///**
// * Test Suite: Fsm::removeUnreachableNodes()
// */
//void removeUnreachableNodes_TS() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	//random tests
//	for (int i = 0; i < 100; ++i) {
//		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, pl);		
//		testRemoveUnreachableNodes(*fsm);
//	}
//
//	// practical examples
//	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	testRemoveUnreachableNodes(*csm);
//	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	testRemoveUnreachableNodes(*fsb);
//	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "FSB");
//	testRemoveUnreachableNodes(*gdc);
//
//	// partition test cases
//	vector<FsmTransformationTestCase<Fsm>> testSuite;
//	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_removeUnreachableNodes.testsuite", testSuite);
//	for (auto tc : testSuite) {
//		cout << "Start TC " << tc.id << endl;
//		testRemoveUnreachableNodes(*tc.m);
//	}
//}
//
///**
// * Test Suite: Fsm::transformToObservableFSM()
// */
//void transformToObservableFSM_TS() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	//random tests
//	for (int i = 0; i < 100; ++i) {
//		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, pl);
//		testTransformToObservableFSM(*fsm);
//	}
//
//	// practical examples
//	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	testTransformToObservableFSM(*csm);
//	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	testTransformToObservableFSM(*fsb);
//	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "FSB");
//	testTransformToObservableFSM(*gdc);
//
//	// partition test cases
//	vector<FsmTransformationTestCase<Fsm>> testSuite;
//	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_transformToObservable.testsuite", testSuite);
//	for (auto tc : testSuite) {
//		cout << "Start TC " << tc.id << endl;
//		testTransformToObservableFSM(*tc.m);
//	}
//}
//
///**
// * Test Suite: Dfsm::minimise()
// */
//void minimise_Dfsm_TS() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	//random tests
//	for (int i = 0; i < 100; ++i) {
//		//auto dfsm = createRandomDfsm("M", 10, 4, 4, pl);
//		//testMinimise_Dfsm(dfsm);
//	}
//
//	// practical examples
//	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	testMinimise_Dfsm(*csm);
//	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	testMinimise_Dfsm(*fsb);
//	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "FSB");
//	testMinimise_Dfsm(*gdc);
//
//	// partition test cases
//	vector<FsmTransformationTestCase<Dfsm>> testSuite;
//	parseFsmTransformationTSFile<Dfsm>("../../../resources/TestSuites/FSM-Transformations/Dfsm_minimise.testsuite", testSuite);
//	for (auto tc : testSuite) {
//		cout << "Start TC " << tc.id << endl;
//		testMinimise_Dfsm(*tc.m);
//	}
//}
//
///**
// * Test Suite: Fsm::minimiseObservableFSM()
// */
//void minimiseObservableFSM_TS() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	//random tests
//	for (int i = 0; i < 100; ++i) {
//		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, pl);
//		Fsm ofsm = fsm->transformToObservableFSM();
//		testMinimiseObservableFSM(ofsm);
//	}
//
//	// practical examples
//	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	testMinimiseObservableFSM(*csm);
//	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	testMinimiseObservableFSM(*fsb);
//	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "FSB");
//	testMinimiseObservableFSM(*gdc);
//
//	// partition test cases
//	vector<FsmTransformationTestCase<Fsm>> testSuite;
//	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_minimiseObservable.testsuite", testSuite);
//	for (auto tc : testSuite) {
//		cout << "Start TC " << tc.id << endl;
//		testMinimiseObservableFSM(*tc.m);
//	}
//}
//
///**
// * Test Suite: Fsm::minimise()
// */
//void minimise_Fsm_TS() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	//random tests
//	for (int i = 0; i < 100; ++i) {
//		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 10, pl);
//		testMinimise_Fsm(*fsm);
//	}
//
//	// practical examples
//	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	testMinimise_Fsm(*csm);
//	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	testMinimise_Fsm(*fsb);
//	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "FSB");
//	testMinimise_Fsm(*gdc);
//
//	//		// these make the program crash
//    //"../../../resources/TestSuites/FSM029.fsm"
//	//"../../../resources/TestSuites/FSM053.fsm") {
//
//	// partition test cases
//	vector<FsmTransformationTestCase<Fsm>> testSuite;
//	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_minimise.testsuite", testSuite);
//	for (auto tc : testSuite) {
//		cout << "Start TC " << tc.id << endl;
//		testMinimise_Fsm(*tc.m);
//	}
//}
//// ====================================================================================================
//
//
//// ====================================================================================================
//// Prüfverfahren "Konstruktion des Produkts"
//
//typedef std::unordered_set<std::shared_ptr<FsmNode>> reachedStates_t;
//typedef std::tuple<reachedStates_t, reachedStates_t, reachedStates_t> reachedStatesTuple_t;
//
//template<typename T>
//bool containsTuple(T &lst, reachedStatesTuple_t &tuple)
//{
//	typename T::const_iterator it;
//	for (it = lst.begin(); it != lst.end(); ++it)
//	{
//		if (std::get<0>(*it) == std::get<0>(tuple) && std::get<1>(*it) == std::get<1>(tuple) && std::get<2>(*it) == std::get<2>(tuple)) {
//			return true;
//		}
//	}
//	return false;
//}
//
///**
//	This function returns true iff L(intersection) = L(m1) ∩ L(m2). Otherwise the function returns false.
//*/
//bool languageIntersectionCheck(const Fsm &m1, const Fsm &m2, const Fsm &intersection) {
//	// init wl with ({m1.q_0},{m2.q_0},{intersection.q_0})	
//	deque<reachedStatesTuple_t> wl;
//	reachedStates_t a = reachedStates_t{ m1.getInitialState() };
//	reachedStates_t b = reachedStates_t{ m2.getInitialState() };
//	reachedStates_t c = reachedStates_t{ intersection.getInitialState() };
//
//	wl.push_back({ a,b,c });
//	std::vector<reachedStatesTuple_t> processed;
//
//	while (not wl.empty()) {
//		auto front = wl.front();
//		wl.pop_front();
//		//calculate lblSet_a
//		unordered_set<FsmLabel> lblSet_a = calcLblSet(std::get<0>(front));
//		//calculate lblSet_b
//		unordered_set<FsmLabel> lblSet_b = calcLblSet(std::get<1>(front));
//		//calculate lblSet_c
//		unordered_set<FsmLabel> lblSet_c = calcLblSet(std::get<2>(front));
//
//		// calculate intersection of lblSet_a and lblSet_b
//		unordered_set<FsmLabel> lblSet_i;
//		for (auto &lbl : lblSet_a) {
//			if (lblSet_b.count(lbl) > 0) lblSet_i.insert(lbl);
//		}
//		if (lblSet_c != lblSet_i) return false;
//
//		// insert front to processed if processed does not contain front already
//		if (not containsTuple(processed, front)) {
//			processed.push_back(front);
//		}
//
//		for (auto lbl : lblSet_c) {
//			// calc tgtNds_a
//			reachedStates_t tgtNds_a = calcTargetNodes(std::get<0>(front), lbl);
//			// calc tgtNds_b
//			reachedStates_t tgtNds_b = calcTargetNodes(std::get<1>(front), lbl);
//			// calc tgtNds_c
//			reachedStates_t tgtNds_c = calcTargetNodes(std::get<2>(front), lbl);
//			reachedStatesTuple_t tuple{ tgtNds_a, tgtNds_b, tgtNds_c };
//			if (not containsTuple(wl, tuple) and not containsTuple(processed, tuple)) {
//				wl.push_back(tuple);
//			}
//		}
//	}
//	return true;
//}
//
///**
// * Test function for Fsm::intersect(const Fsm & f).
// */
//void testIntersection(Fsm &m1, const Fsm &m2) {
//	// get copy of m1 and m2
//    const Fsm copyOfM1 = Fsm(m1);
//
//    // use Algorithm to calculate result
//	const Fsm intersection = m1.intersect(m2);
//
//	// first check invariant for m1 and intersection   (we don't need to check invariant for m2 because it's const)
//	bool invariantViolationOfM1 = not checkFsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolationOfM1, "Fsm class invariant still holds for M1 after calculation.");
//	bool invariantViolationOfIntersection = not checkFsmClassInvariant(intersection);
//	fsmlib_assert("TC", not invariantViolationOfIntersection, "Fsm class invariant holds for intersection after calculation.");
//	// stop test execution at this point if invariant of m or intersection does not hold anymore
//	if (invariantViolationOfM1 || invariantViolationOfIntersection) return;
//
//	// check language intersection
//	fsmlib_assert("TC", languageIntersectionCheck(m1,m2,intersection), "Language of the result is intersection of L(M1) and L(M2)");
//
//	// check for forbidden side effects
//	fsmlib_assert("TC", checkForEqualStructure(m1, copyOfM1), "M1 was not changed by algorithm");
//}
//
///**
// * Test function for Fsm::intersect(const Fsm & f). (Dfsm Context)
// */
//void testIntersection(Dfsm &m1, const Fsm &m2) {
//	// get copy of m1 and m2
//	const Dfsm copyOfM1 = Dfsm(m1);
//
//	// use Algorithm to calculate result
//	const Fsm intersection = m1.intersect(m2);
//
//	// first check invariant for m1 and intersection   (we don't need to check invariant for m2 because it's const)
//	bool invariantViolationOfM1 = not checkDfsmClassInvariant(m1);
//	fsmlib_assert("TC", not invariantViolationOfM1, "Dfsm class invariant still holds for M1 after calculation.");
//	bool invariantViolationOfIntersection = not checkFsmClassInvariant(intersection);
//	fsmlib_assert("TC", not invariantViolationOfIntersection, "Fsm class invariant holds for intersection after calculation.");
//	// stop test execution at this point if invariant of m or intersection does not hold anymore
//	if (invariantViolationOfM1 || invariantViolationOfIntersection) return;
//
//	// check language intersection
//	fsmlib_assert("TC", languageIntersectionCheck(m1, m2, intersection), "Language of the result is intersection of L(M1) and L(M2)");
//
//	// check for forbidden side effects
//	fsmlib_assert("TC", checkForEqualStructure(m1, copyOfM1), "M1 was not changed by algorithm");
//}
//
//struct IntersectTestCase {
//	string id;
//	string m1Path;
//	string m2Path;
//};
//
//shared_ptr<vector<IntersectTestCase>> parseIntersectTSFile(const string &testSuitePath) {
//	string fname = testSuitePath;
//	vector<IntersectTestCase> testSuite;
//	ifstream inputFile(fname);
//	if (inputFile.is_open())
//	{
//		string line;
//
//		while (getline(inputFile, line))
//		{
//			vector<string> lineContent;
//			stringstream ss(line);
//
//			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
//
//			IntersectTestCase tc;
//			tc.id = lineContent.at(0);
//			tc.m1Path = lineContent.at(1);
//			tc.m2Path = lineContent.at(2);			
//
//			testSuite.push_back(tc);
//		}
//		inputFile.close();
//
//		return make_shared<vector<IntersectTestCase>>(testSuite);
//	}
//	else
//	{
//		cout << "Unable to open input file" << endl;
//		exit(EXIT_FAILURE);
//	}
//}
//
///*
// *	Random Test Suite for test of Fsm::getCharacterisationSet().
//*/
//void intersection_TS_Random() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	const int seed = 8950;
//	srand(seed);
//	// random tests
//	for (int i = 0; i < 100; ++i) {
//		cout << "i:" << i << endl;
//		auto m1 = Fsm::createRandomFsmRepeatable("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
//		const auto m2 = m1->createMutantRepeatable("M2", 2, 2);
//		testIntersection(*m1, *m2);
//	}
//
//	for (int i = 0; i < 100; ++i) {
//		cout << "i:" << i << endl;
//		auto m1 = Fsm::createRandomFsmRepeatable("M1", 4, 4, 10, make_shared<FsmPresentationLayer>());
//		const auto m2 = Fsm::createRandomFsmRepeatable("M2", 4, 4, 10, make_shared<FsmPresentationLayer>());
//		testIntersection(*m1, *m2);
//	}
//
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	auto m1 = Dfsm("M", 15, 4, 4, pl);
//	//	cout << m1.getMaxState() << endl;
//	//	//auto m1 = createRandomDfsm("M1", 15, 4, 4, pl);
//	//	const auto m2 = m1.createMutant("M2", 2, 2);
//	//	testIntersection(m1, *m2);
//	//}
//
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	auto m1 = Dfsm("M1", 15, 4, 4, pl, true);
//	//	//auto m1 = createRandomDfsm("M1", 15, 4, 4, pl);
//	//	const auto m2 = Dfsm("M1", 15, 4, 4, pl, true);
//	//	testIntersection(m1, m2);
//	//}
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//	for(int i = 0; i < 100; ++i){
//		testIntersection(*csm, *csm->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1));
//	}
//	cout << "FSB Tests" << endl;
//	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//	for (int i = 0; i < 100; ++i) {
//		testIntersection(*fsb, *fsb->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1));
//	}
//	cout << "GDC Tests" << endl;
//	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//	for (int i = 0; i < 100; ++i) {
//		testIntersection(*gdc, *gdc->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1));
//	}
//
//	// partition test cases
//	auto testSuite = parseIntersectTSFile("../../../resources/TestSuites/Intersection/Fsm_intersect.testsuite");
//	for (auto tc : *testSuite) {
//		cout << "Start Test Case : " << tc.id << endl;
//		shared_ptr<Fsm> m1 = make_shared<Fsm>(tc.m1Path, pl, "M1");
//		shared_ptr<Fsm> m2 = make_shared<Fsm>(tc.m2Path, pl, "M2");
//		testIntersection(*m1, *m2);
//	}
//}
//
//
//// ====================================================================================================
//
//// ====================================================================================================
//// Prüfverfahren "Berechnung von Distinguishing Traces"
//
////set<vector<int>> calcCompleteOutputTraces(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc) {	
////	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
////	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
////
////	for (int x : inputTrc) {
////		for (auto reachedNode : wl) {
////			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {
////				if (transition->getLabel()->getInput() != x) continue; 
////				vector<int> outputTrc = get<1>(reachedNode);
////				outputTrc.push_back(transition->getLabel()->getOutput());
////				wl_next.insert({ transition->getTarget(), outputTrc });
////			}
////		}
////		wl = wl_next;
////		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
////	}
////
////	set<vector<int>> outputTrcs;
////	for (auto reachedNode : wl) {
////		outputTrcs.insert(std::get<1>(reachedNode));
////	}
////
////	return outputTrcs;
////}
//
///*
//	Second Version of Algorithm. Applies inputTrc to startNode and produces set of outputtraces. If some FsmNode is reached with a prefix
//	of inputTrc in which the next input of inputTrc is undefined, the corresponding output trace will be expanded by an 'NULL' output
//	(not contained in the output alphabet) and algorithm stays in this FsmNode. Then the next input is applied.
//*/
//set<vector<int>> calcCompleteOutputTraces2(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc, const int maxOutput) {
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
//	const int nullOutput = maxOutput + 1; // or nullOutput = -1
//
//	for (int x : inputTrc) {
//		for (auto reachedNode : wl) {
//			bool defined = false;
//			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {				
//				if (transition->getLabel()->getInput() != x) continue;
//				defined = true;
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(transition->getLabel()->getOutput());
//				wl_next.insert({ transition->getTarget(), outputTrc });
//			}
//			if (not defined) {
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(nullOutput);
//				wl_next.insert({ std::get<0>(reachedNode), outputTrc });
//			}
//		}
//		wl = wl_next;
//		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
//	}
//
//	set<vector<int>> outputTrcs;
//	for (auto reachedNode : wl) {
//		//for (auto i : std::get<1>(reachedNode)) cout << i << ",";
//		//cout << "\n";
//		outputTrcs.insert(std::get<1>(reachedNode));
//	}
//	return outputTrcs;
//}
//
///**
// *Third Version of Algorithm. Applies inputTrc to startNode and produces set of outputtraces.If some FsmNode is reached with a prefix
// *of inputTrc in which the next input of inputTrc is undefined, the corresponding output trace will be expanded by an 'NULL' output
// *(not contained in the output alphabet) and algorithm stays in this FsmNode.Then the next input is applied.
// */
//set<vector<int>> calcCompleteOutputTraces3(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc, const int nullOutput, bool b) {
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
//	cout << "start" << endl;
//
//	for (int x : inputTrc) {
//		for (auto reachedNode : wl) {
//			bool defined = false;
//			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {
//				if (transition->getLabel()->getInput() != x) continue;
//				defined = true;
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(transition->getLabel()->getOutput());
//				wl_next.insert({ transition->getTarget(), outputTrc });
//			}
//			if (not defined) {
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(nullOutput);
//				wl_next.insert({ std::get<0>(reachedNode), outputTrc });
//			}
//		}
//		wl = wl_next;
//		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
//	}
//
//	set<vector<int>> outputTrcs;
//	for (auto reachedNode : wl) {
//		//for (auto i : std::get<1>(reachedNode)) cout << i << ",";
//		//cout << "\n";
//		outputTrcs.insert(std::get<1>(reachedNode));
//	}
//	cout << "end" << endl;
//	return outputTrcs;
//}
//
///*
//	Checks if given inTrc is a Distinguishing Trace for q1 and q2. 
//	Returns true iff inTrc produces some outTrc of the same length
//	which is only contained in the language of one of these FsmNodes.
//*/
//bool isDistTrc(const shared_ptr<FsmNode> q1, const shared_ptr<FsmNode> q2, const vector<int> &inTrc, const int maxOutput) {
//	//return calcCompleteOutputTraces(q1, inTrc) != calcCompleteOutputTraces(q2, inTrc);
//	return calcCompleteOutputTraces2(q1, inTrc, maxOutput) != calcCompleteOutputTraces2(q2, inTrc, maxOutput);
//}
//
///*
//	Returns true iff w contains a Distinguishing Trace for q1 and q2.
//*/
//bool containsDistTrcForPair(const shared_ptr<FsmNode> q1, const shared_ptr<FsmNode> q2, const IOListContainer &w, const int maxOutput) {
//	for (auto inTrc : *w.getIOLists()) {
//		/*if (isDistTrc(q1, q2, inTrc)) return true;*/
//		if (isDistTrc(q1, q2, inTrc, maxOutput)) return true;
//	}
//	return false;
//}
//
///*
//	m has to be minimal and observable.
//	Returns true iff w is a Characterisation Set of m.
//*/
//bool isCharaterisationSet(const Fsm &m, const IOListContainer w) {
//	for (int q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
//		for (int q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
//			if (not containsDistTrcForPair(m.getNodes().at(q1Idx), m.getNodes().at(q2Idx), w, m.getMaxOutput())) return false;
//		}
//	}
//	return true;
//}
//
///*
//	Checks if trc is a non empty prefix of some element contained in w.
//*/
//bool isPrefixOfElement(const vector<int> &trc, const std::shared_ptr<const Tree> w) {
//	//w->addToRoot() may be faster
//	if (trc.size() == 0) return false;
//
//	for (auto elem : *w->getIOLists().getIOLists()) {
//		if (elem.size() < trc.size()) continue;
//
//		for (int i = 0; i < trc.size(); ++i) {
//			if (trc.at(i) != elem.at(i)) continue;
//		}
//		return true;
//	}
//	return false;
//}
//
///*
//	m has to be minimal and observable. qi ist expected to be a state of m.
//	wi is expected to contain traces over the input alphabet of m. w is expected to be a characterisation set of m.
//
//	Returns true iff wi is a State Identification Set for qi in m with Characterisation Set w.
//*/
//bool isStateIdentificationSet(const Fsm &m, const shared_ptr<FsmNode> qi, const std::shared_ptr<Tree> wi, const std::shared_ptr<Tree> w) {
//	for (auto trc : *wi->getIOLists().getIOLists()) {
//		if (not isPrefixOfElement(trc, w)) return false;
//	}
//
//	for (auto q : m.getNodes()) {
//		if (q == qi) continue;
//
//		if (not containsDistTrcForPair(q, qi, wi->getIOLists(), m.getMaxOutput())) return false;
//	}
//
//	return true;
//}
//
///*
//	m has to be minimal and observable. qi ist expected to be a state of m.
//	wi is expected to be a State Identification Set of qi in m. w is expected to be a characterisation set of m.
//
//	Returns true iff there is no subset of wi that is a State Identification Set of qi in m.
//*/
//bool isMinimalStateIdentificationSet(const Fsm &m, const shared_ptr<FsmNode> qi, const std::shared_ptr<Tree> wi, const std::shared_ptr<Tree> w) {
//	auto pl = make_shared<FsmPresentationLayer>();
//	for (int i = 0; i < wi->getIOLists().getIOLists()->size(); ++i) {
//		IOListContainer iolc(pl);
//		for (int j = 0; j < wi->getIOLists().getIOLists()->size(); ++j) {
//			if (i == j) continue;
//			iolc.add({ wi->getIOLists().getIOLists()->at(j), pl });
//		}
//		shared_ptr<Tree> alternativeWi = make_shared<Tree>(make_shared<TreeNode>(), pl);
//		alternativeWi->addToRoot(iolc);
//		if (isStateIdentificationSet(m, qi, alternativeWi, w)) return false;
//	}
//	return true;
//}
//
//
//
///**
// * Test function for Dfsm::getCharacterisationSet().
// * Parameter m is expected to be a minimal and complete Dfsm.
// */
//void testGetCharacterisationSet_Dfsm(Dfsm &m) {
//	// get copy of m
//	const Dfsm copyOfM = Dfsm(m);
//
//	// use Algorithm to calculate result
//	const auto w = m.getCharacterisationSet();
//
//	// first check invariant of m
//	bool invariantViolation = not checkDfsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// check definition of 'Characterisation Set' for w
//	fsmlib_assert("TC", isCharaterisationSet(m, w), "Result is a Characterisation Set for M.");
//
//	fsmlib_assert("TC", *m.characterisationSet->getIOLists().getIOLists() == *w.getIOLists(), "Result is stored in attribute.");
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//}
//
///**
// * Test function for Fsm::getCharacterisationSet().
// * Parameter m is expected to be a minimal and observable Fsm.
// */
//void testGetCharacterisationSet_Fsm(Fsm &m) {
//	// get copy of m
//	const Fsm copyOfM = Fsm(m);
//
//	// use Algorithm to calculate result
//	const auto w = m.getCharacterisationSet();
//
//	// first check invariant of m
//	bool invariantViolation = not checkFsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// check definition of 'Characterisation Set' for w
//	fsmlib_assert("TC", isCharaterisationSet(m, w), "Result is a Characterisation Set for M.");
//
//	fsmlib_assert("TC", *m.characterisationSet->getIOLists().getIOLists() == *w.getIOLists(), "Result is stored in attribute.");
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//}
//
///**
// * Test function for FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
// *                                           const vector<shared_ptr<PkTable>>& pktblLst,
// *                                           const int maxInput)
// * Parameter m is expected to be a minimal and complete Dfsm.
// */
//void testCalcDistinguishingTrace1(Dfsm &m) {
//	// get copy of m
//	const Dfsm copyOfM = Dfsm(m);
//
//	// calculate the needed parameters from m
//	m.calcPkTables();
//	const auto tables = m.pktblLst;
//
//	// test each pair of different nodes
//	for (int q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
//		for (int q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
//			if (q1Idx == q2Idx) continue;
//			const auto q1 = m.getNodes().at(q1Idx);
//			const auto q2 = m.getNodes().at(q2Idx);
//
//			// use Algorithm to calculate result
//			InputTrace inTrc = q1->calcDistinguishingTrace(q2, tables, m.getMaxInput());
//
//			// first check invariant of m
//			bool invariantViolation = not checkDfsmClassInvariant(m);
//			fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//			// stop test execution at this point if invariant of m does not hold anymore
//			if (invariantViolation) return;
//
//			// check definition of 'Distinguishing Trace' for inTrc
//			fsmlib_assert("TC", isDistTrc(q1, q2, inTrc.get(), m.getMaxOutput()), "Calculated Trace is a Distinguishing Trace for q1 and q2.");
//
//			// check if structure of m has changed
//			fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//		}
//	}
//}
//
///**
// * Test function for FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
// *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
// *                                           const int maxInput,
// *                                           const int maxOutput)
// * Parameter m is expected to be a minimal and observable Fsm.
// */
//void testCalcDistinguishingTrace2(Fsm &m) {
//	// get copy of m
//	const Fsm copyOfM = Fsm(m);
//
//	// calculate the needed parameters from m
//	m.calcOFSMTables();
//	const auto tables = m.ofsmTableLst;
//
//	// test each pair of different nodes
//	for (int q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
//		for (int q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
//			if (q1Idx == q2Idx) continue;
//			const auto q1 = m.getNodes().at(q1Idx);
//			const auto q2 = m.getNodes().at(q2Idx);
//
//			// use Algorithm to calculate result
//			InputTrace inTrc = q1->calcDistinguishingTrace(q2, tables, m.getMaxInput(), m.getMaxOutput());
//
//			// first check invariant of m
//			bool invariantViolation = not checkFsmClassInvariant(m);
//			fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
//			// stop test execution at this point if invariant of m does not hold anymore
//			if (invariantViolation) return;
//
//			// check definition of 'Distinguishing Trace' for inTrc
//			fsmlib_assert("TC", isDistTrc(q1, q2, inTrc.get(), m.getMaxOutput()), "Calculated Trace is a Distinguishing Trace for q1 and q2.");
//
//			// check if structure of m has changed
//			fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//		}
//	}
//}
//
///**
// * Test function for Fsm::calcStateIdentificationSets().
// * m is expected to be a minimal and observable Fsm.
// */
//void testCalcStateIdentificationSets(Fsm &m) {
//	// get copy of m
//	const Fsm copyOfM = Fsm(m);
//
//	// calculate the needed parameters from m
//	m.getCharacterisationSet();
//	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();
//	
//
//	// use Algorithm to calculate result
//	m.calcStateIdentificationSets();
//	const auto stateIdSets = m.stateIdentificationSets;
//
//	// first check invariant of m
//	bool invariantViolation = not checkFsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// Check Definition of minimal State Identification Set for each element in stateIdSets
//	fsmlib_assert("TC", stateIdSets.size() == m.getNodes().size(), "Number of calculated State Identification Sets matches the number of states of M.");
//	for (int i = 0; i < stateIdSets.size(); ++i) {
//		fsmlib_assert("TC", isStateIdentificationSet(m,m.getNodes().at(i),stateIdSets.at(i),m.characterisationSet), "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");
//		fsmlib_assert("TC", isMinimalStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), "M.stateIdentificationSets[i] is a minimal State Identification Set for M.nodes[i].");
//	}
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//
//	// check if m.characterisationSet has changed
//	fsmlib_assert("TC", *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), "characterisation set of M has not changed");
//}
//
///**
// * Test function for Fsm::calcStateIdentificationSets(). Test in context of Dfsm.
// * m is expected to be a minimal Dfsm.
// */
//void testCalcStateIdentificationSets(Dfsm &m) {
//	// get copy of m
//	const Dfsm copyOfM = Dfsm(m);
//
//	// calculate the needed parameters from m
//	m.Fsm::getCharacterisationSet();
//	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();
//
//
//	// use Algorithm to calculate result
//	m.calcStateIdentificationSets();
//	const auto stateIdSets = m.stateIdentificationSets;
//
//	// first check invariant of m
//	bool invariantViolation = not checkDfsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// Check Definition of minimal State Identification Set for each element in stateIdSets
//	fsmlib_assert("TC", stateIdSets.size() == m.getNodes().size(), "Number of calculated State Identification Sets matches the number of states of M.");
//	for (int i = 0; i < stateIdSets.size(); ++i) {
//		fsmlib_assert("TC", isStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");
//		fsmlib_assert("TC", isMinimalStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), "M.stateIdentificationSets[i] is a minimal State Identification Set for M.nodes[i].");
//	}
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//
//	// check if m.characterisationSet has changed
//	fsmlib_assert("TC", *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), "characterisation set of M has not changed");
//}
//
///**
// * Test function for Fsm::calcStateIdentificationSetsFast().
// * m is expected to be a minimal and observable Fsm.
// */
//void testCalcStateIdentificationSetsFast(Fsm &m) {
//	// get copy of m
//	const Fsm copyOfM = Fsm(m);
//
//	// calculate the needed parameters from m
//	m.getCharacterisationSet();
//	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();
//
//
//	// use Algorithm to calculate result
//	m.calcStateIdentificationSetsFast();
//	const auto stateIdSets = m.stateIdentificationSets;
//
//	// first check invariant of m
//	bool invariantViolation = not checkFsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// Check Definition of State Identification Set for each element in stateIdSets
//	fsmlib_assert("TC", stateIdSets.size() == m.getNodes().size(), "Number of calculated State Identification Sets matches the number of states of M.");
//	for (int i = 0; i < stateIdSets.size(); ++i) {
//		fsmlib_assert("TC", isStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");
//	}
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//
//	// check if m.characterisationSet has changed
//	fsmlib_assert("TC", *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), "characterisation set of M has not changed");
//}
//
///**
// * Test function for Fsm::calcStateIdentificationSetsFast().
// * m is expected to be a minimal and complete Dfsm.
// */
//void testCalcStateIdentificationSetsFast(Dfsm &m) {
//	// get copy of m
//	const Dfsm copyOfM = Dfsm(m);
//
//	// calculate the needed parameters from m
//	m.Fsm::getCharacterisationSet();
//	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();
//
//
//	// use Algorithm to calculate result
//	m.calcStateIdentificationSetsFast();
//	const auto stateIdSets = m.stateIdentificationSets;
//
//	// first check invariant of m
//	bool invariantViolation = not checkDfsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	// Check Definition of State Identification Set for each element in stateIdSets
//	fsmlib_assert("TC", stateIdSets.size() == m.getNodes().size(), "Number of calculated State Identification Sets matches the number of states of M.");
//	for (int i = 0; i < stateIdSets.size(); ++i) {
//		fsmlib_assert("TC", isStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");
//	}
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//
//	// check if m.characterisationSet has changed
//	fsmlib_assert("TC", *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), "characterisation set of M has not changed");
//}
//
//struct DistinguishingTraceTestCase {
//	string id;
//	string fsmPath;
//	string csp;
//	string sisp;
//};
//
//shared_ptr<vector<DistinguishingTraceTestCase>> parseDistinguishingTraceTSFile(const string &testSuitePath) {
//	string fname = testSuitePath;
//	vector<DistinguishingTraceTestCase> testSuite;
//	ifstream inputFile(fname);
//	if (inputFile.is_open())
//	{
//		string line;
//
//		while (getline(inputFile, line))
//		{
//			vector<string> lineContent;
//			stringstream ss(line);
//
//			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
//
//			DistinguishingTraceTestCase tc;
//			tc.id = lineContent.at(0);
//			tc.fsmPath = lineContent.at(1);
//			tc.csp = lineContent.at(2);
//			tc.sisp = lineContent.at(3);
//
//			testSuite.push_back(tc);
//		}
//		inputFile.close();
//
//		return make_shared<vector<DistinguishingTraceTestCase>>(testSuite);
//	}
//	else
//	{
//		cout << "Unable to open input file" << endl;
//		exit(EXIT_FAILURE);
//	}
//}
//
///*
// *	Random Test Suite for test of Fsm::getCharacterisationSet().
//*/
//void getCharacterisationSet_Dfsm_TS_Random() {
//	const int seed = 261250;
//	srand(seed);
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	auto m = Dfsm("M", 15, 4, 4, pl);
//	//	auto minM = m.minimise();
//	//	cout << "minFsm size: " << minM.size() << endl;
//	//	testGetCharacterisationSet_Dfsm(minM);
//	//}
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testGetCharacterisationSet_Dfsm(csm);
//
//	cout << "FSB Tests" << endl;
//	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testGetCharacterisationSet_Dfsm(fsb);
//
//	cout << "GDC Tests" << endl;
//	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testGetCharacterisationSet_Dfsm(gdc);
//
//
//	// partition test cases
///*	auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Dfsm_getCharacterisationSet.testsuite");
//	for (auto tc : *testSuite) {
//		cout << "Start Test Case : " << tc.id << endl;
//		shared_ptr<Dfsm> m = make_shared<Dfsm>(tc.fsmPath, pl, "M");
//		if (tc.csp == "true") {
//			cout << "csp = true" << endl;
//			m->getCharacterisationSet();
//		}
//		if (tc.sisp == "true") {
//			cout << "sisp = true" << endl;
//			m->calcStateIdentificationSetsFast();
//		}
//		testGetCharacterisationSet_Dfsm(*m);
//	}*/     
//}
//
///*
// *	Random Test Suite for test of Fsm::getCharacterisationSet().
//*/
//void getCharacterisationSet_Fsm_TS_Random() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	const int seed = 64162;
//	srand(seed);
//
//	// random tests
//	/*for (int i = 0; i < 100; ++i) {
//		cout << "i:" << i << endl;
//		auto fsm = Fsm::createRandomFsm("M1", 4, 4, 5, pl);
//		auto minFsm = fsm->minimise();
//		cout << "minFsm size: " << minFsm.size() << endl;
//		testGetCharacterisationSet_Fsm(minFsm);
//	}*/
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testGetCharacterisationSet_Fsm(csm);
//
//	cout << "FSB Tests" << endl;
//	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testGetCharacterisationSet_Fsm(fsb);
//
//	cout << "GDC Tests" << endl;
//	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testGetCharacterisationSet_Fsm(gdc);
//
//	// partition test cases
//	//auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_getCharacterisationSet.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> m = make_shared<Fsm>(tc.fsmPath, pl, "M");
//	//	if (tc.csp == "true") {
//	//		cout << "csp = true" << endl;
//	//		m->getCharacterisationSet();
//	//	}
//	//	if (tc.sisp == "true") {
//	//		cout << "sisp = true" << endl;
//	//		m->calcStateIdentificationSetsFast();
//	//	}
//	//	testGetCharacterisationSet_Fsm(*m);
//	//}
//}
//
///*
// *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
// *                                           const vector<shared_ptr<PkTable>>& pktblLst,
// *                                           const int maxInput)
//*/
//void calcDistinguishingTrace1_TS_Random() {
//	const int seed = 329132;
//	srand(seed);
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	auto m = Dfsm("M", 15, 4, 2, pl);
//	//	auto minM = m.minimise();
//	//	cout << "minFsm size: " << minM.size() << endl;
//	//	testCalcDistinguishingTrace1(minM);
//	//}
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testCalcDistinguishingTrace1(csm);
//
//	cout << "FSB Tests" << endl;
//	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testCalcDistinguishingTrace1(fsb);
//
//	cout << "GDC Tests" << endl;
//	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testCalcDistinguishingTrace1(gdc);
//
//	// partition test cases
//	//auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/FsmNode_calcDistinguishingTrace_pk.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> m = make_shared<Dfsm>(tc.fsmPath, pl, "M");
//	//	if (tc.csp == "true") {
//	//		cout << "csp = true" << endl;
//	//		m->getCharacterisationSet();
//	//	}
//	//	if (tc.sisp == "true") {
//	//		cout << "sisp = true" << endl;
//	//		m->calcStateIdentificationSetsFast();
//	//	}
//	//	testCalcDistinguishingTrace1(*m);
//	//}
//}
//
///*
// *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
// *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
// *                                           const int maxInput,
// *                                           const int maxOutput)
//*/
//void calcDistinguishingTrace2_TS_Random() {
//	const int seed = 36185;
//	srand(seed);
//	
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	auto fsm = Fsm::createRandomFsm("M1", 4, 4, 5, pl);
//	//	auto minFsm = fsm->minimise();
//	//	cout << "minFsm size: " << minFsm.size() << endl;
//	//	testCalcDistinguishingTrace2(minFsm);
//	//}
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testCalcDistinguishingTrace2(csm);
//
//	cout << "FSB Tests" << endl;
//	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testCalcDistinguishingTrace2(fsb);
//
//	cout << "GDC Tests" << endl;
//	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testCalcDistinguishingTrace2(gdc);
//
//	// partition tests
//	//auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/FsmNode_calcDistinguishingTrace_ofsm.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> m = make_shared<Fsm>(tc.fsmPath, pl, "M");
//	//	if (tc.csp == "true") {
//	//		cout << "csp = true" << endl;
//	//		m->getCharacterisationSet();
//	//	}
//	//	if (tc.sisp == "true") {
//	//		cout << "sisp = true" << endl;
//	//		m->calcStateIdentificationSetsFast();
//	//	}
//	//	testCalcDistinguishingTrace2(*m);
//	//}
//}
//
///*
// *	Random Test Suite for test of Fsm::calcStateIdentificationSets().
//*/
//void calcStateIdentificationSets_TS_Random() {
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	
//	//const int seed = 3447;
//	//srand(seed);
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 10 + 1;
//	//	int mI = rand() % 7;
//	//	int mO = mI;
//	//	auto fsm = Fsm::createRandomFsmRepeatable("M1", mI, mO, size, make_shared<FsmPresentationLayer>());
//	//	auto minFsm = fsm->minimise();
//	//	cout << "minFsm size: " << minFsm.size() << endl;
//	//	testCalcStateIdentificationSets(minFsm);
//	//}
//
//	//shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = rand() % 6;
//	//	auto m = Dfsm("M", size, mI, mO, pl, true);
//	//	auto minM = m.minimise();
//	//	cout << "minFsm size: " << minM.size() << endl;
//	//	cout << minM << endl;
//	//	testCalcStateIdentificationSets(minM);
//	//}
//
//	// practical examples
//	cout << "CSM Tests" << endl;
//	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testCalcStateIdentificationSets(csm);
//
//	cout << "FSB Tests" << endl;
//	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testCalcStateIdentificationSets(fsb);
//
//	cout << "GDC Tests" << endl;
//	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testCalcStateIdentificationSets(gdc);
//
//	// partition tests
//	//auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_calcStateIdentificationSets.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> m = make_shared<Fsm>(tc.fsmPath, pl, "M");
//	//	if (tc.csp == "true") {
//	//		cout << "csp = true" << endl;
//	//		m->getCharacterisationSet();
//	//	}
//	//	if (tc.sisp == "true") {
//	//		cout << "sisp = true" << endl;
//	//		m->calcStateIdentificationSetsFast();
//	//	}
//	//	testCalcStateIdentificationSets(*m);
//	//}
//}
//
//shared_ptr<Fsm> createPartialMutant(shared_ptr<Fsm> m) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
//	}
//
//	// Now add transitions that correspond exactly to the transitions in
//	// this FSM
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m->getNodes()[n];
//		int ignoreInput = rand() % (m->getMaxInput() + 1);
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			if (tr->getLabel()->getInput() == ignoreInput) continue;
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//	}
//	return make_shared<Fsm>(m->getName(), m->getMaxInput(), m->getMaxOutput(), lst, m->getPresentationLayer());
//}
//
///**
// * Transform m to a complete Fsm by adding self loops in states for undefined inputs producing some nullouput not contained in the
// * regular output alphabet.
// */
//shared_ptr<Fsm> transformToComplete(const shared_ptr<const Fsm> m, const size_t nullOutput) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
//	}
//
//	bool partial = false;
//
//	// Now add transitions that correspond exactly to the transitions in
//	// this FSM
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m->getNodes()[n];
//		set<int> definedInputs;
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			definedInputs.insert(tr->getLabel()->getInput());
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//		// add self loops with nullOutputs for undefined inputs
//		for (int input = 0; input <= m->getMaxInput(); ++input) {
//			if (definedInputs.count(input) == 0) {
//				partial = true;
//				shared_ptr<FsmTransition> newTr =
//					make_shared<FsmTransition>(theNewFsmNodeSrc, lst[n], make_shared<FsmLabel>(input, nullOutput, m->getPresentationLayer()));
//				theNewFsmNodeSrc->addTransition(newTr);
//			}
//		}
//	}
//	if (partial) {
//		auto completeM = make_shared<Fsm>(m->getName(), m->getMaxInput(), nullOutput, lst, m->getPresentationLayer());
//		completeM->initStateIdx = m->initStateIdx;
//		return completeM;
//	}
//	else {
//		auto completeM = make_shared<Fsm>(m->getName(), m->getMaxInput(), m->getMaxOutput(), lst, m->getPresentationLayer());
//		completeM->initStateIdx = m->initStateIdx;
//		return completeM;
//	}
//
//}
//
///**
// * Transform m to a complete Fsm by adding self loops in states for undefined inputs producing some nullouput not contained in the
// * regular output alphabet. This function specifies the nulloutput as m.maxOutput + 1.
// */
//shared_ptr<Fsm> transformToComplete(shared_ptr<Fsm> m) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
//	}
//
//	int nullOutput = m->getMaxOutput() + 1;
//
//	// Now add transitions that correspond exactly to the transitions in
//	// this FSM
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m->getNodes()[n];	
//		set<int> definedInputs;
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			definedInputs.insert(tr->getLabel()->getInput());
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//		// add self loops with nullOutputs for undefined inputs
//		for (int input = 0; input <= m->getMaxInput(); ++input) {
//			if (definedInputs.count(input) == 0) {
//				shared_ptr<FsmTransition> newTr =
//					make_shared<FsmTransition>(theNewFsmNodeSrc, lst[n], make_shared<FsmLabel>(input, nullOutput, m->getPresentationLayer()));
//				theNewFsmNodeSrc->addTransition(newTr);
//			}
//		}
//	}
//	return make_shared<Fsm>(m->getName(), m->getMaxInput(), m->getMaxOutput() + 1, lst, m->getPresentationLayer());
//}
//
//
//
///*
// *	Random Test Suite for test of Fsm::calcStateIdentificationSetsFast().
//*/
//void calcStateIdentificationSetsFast_TS_Random() {
//	//const int seed = 1376;
//	//srand(seed);
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 10 + 1;
//	//	int mI = rand() % 7;
//	//	int mO = mI;
//	//	auto fsm = Fsm::createRandomFsmRepeatable("M1", mI, mO, size, make_shared<FsmPresentationLayer>());
//	//	auto tmp = createPartialMutant(fsm);
//	//	auto minFsm = tmp->minimise();
//	//	cout << "minFsm size: " << minFsm.size() << endl;
//	//	cout << "minFsm cs? " << minFsm.isCompletelyDefined() << endl;
//	//	testCalcStateIdentificationSetsFast(minFsm);
//	//}
//
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 10 + 1;
//	//	int mI = rand() % 7;		
//	//	int mO = mI;
//	//	auto fsm = Fsm::createRandomFsmRepeatable("M1", mI, mO, size, make_shared<FsmPresentationLayer>());
//	//	auto minFsm = fsm->minimise();
//	//	cout << "minFsm size: " << minFsm.size() << endl;
//	//	testCalcStateIdentificationSetsFast(minFsm);
//	//}
//
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = rand() % 6;
//	//	auto m = Dfsm("M", size, mI, mO, pl, true);
//	//	auto minM = m.minimise();
//	//	cout << "minFsm size: " << minM.size() << endl;
//	//	cout << minM << endl;
//	//	testCalcStateIdentificationSetsFast(minM);
//	//}
//
//    // practical examples
//	cout << "CSM Tests" << endl;
//	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//	testCalcStateIdentificationSetsFast(csm);
//
//	cout << "FSB Tests" << endl;
//	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//	testCalcStateIdentificationSetsFast(fsb);
//
//	cout << "GDC Tests" << endl;
//	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//	testCalcStateIdentificationSetsFast(gdc);
//
//
//	// partition tests
//	//auto testSuite = parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_calcStateIdentificationSetsFast.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> m = make_shared<Fsm>(tc.fsmPath, pl, "M");
//	//	if (tc.csp == "true") {
//	//		cout << "csp = true" << endl;
//	//		m->getCharacterisationSet();
//	//	}
//	//	if (tc.sisp == "true") {
//	//		cout << "sisp = true" << endl;
//	//		m->calcStateIdentificationSetsFast();
//	//	}
//	//	testCalcStateIdentificationSetsFast(*m);
//	//}
//}
//
//// ====================================================================================================
//
///**
// * Calculate and return the maxOutput of the Fsms m and mutants.
// */
//int getMaxOutput(const Fsm & m, const vector<shared_ptr<const Fsm>>& mutants) {
//	int maxOutput = m.getMaxOutput();
//	for (shared_ptr<const Fsm > mut : mutants) {
//		maxOutput < mut->getMaxOutput() ? maxOutput = mut->getMaxOutput() : maxOutput;
//	}
//	return maxOutput;
//}
//
//class TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates) = 0;
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates) = 0;
//};
//
//class WMethodGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates) 
//	{
//		cout << "WMethodGenerator fsm variant" << endl;
//		return make_shared<IOListContainer>(m.wMethod(numAddStates)); 
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "WMethodGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.wMethod(numAddStates));
//	}
//};
//
//class WMethodOnMinimisedFsmGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "WMethodOnMinimisedFsmGenerator fsm variant" << endl;
//		return make_shared<IOListContainer>(m.wMethodOnMinimisedFsm(numAddStates));
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "WMethodOnMinimisedFsmGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.wMethodOnMinimisedFsm(numAddStates));
//	}
//};
//
//class WMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "WMethodOnMinimisedDfsmGenerator fsm variant" << endl;
//		return nullptr;
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "WMethodOnMinimisedDfsmGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.wMethodOnMinimisedDfsm(numAddStates));
//	}
//};
//
//class WpMethodGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "WpMethodGenerator fsm variant" << endl;
//		return make_shared<IOListContainer>(m.wpMethod(numAddStates));
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "WpMethodGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.wpMethod(numAddStates));
//	}
//};
//
//class WpMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "WpMethodOnMinimisedDfsmGenerator fsm variant" << endl;
//		return nullptr;
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "WpMethodOnMinimisedDfsmGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.wpMethodOnMinimisedDfsm(numAddStates));
//	}
//};
//
//class HsiMethodGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "HsiMethodGenerator fsm variant" << endl;
//		return make_shared<IOListContainer>(m.hsiMethod(numAddStates));
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "HsiMethodGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.hsiMethod(numAddStates));
//	}
//};
//
//class HMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
//public:
//	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
//	{
//		cout << "HMethodOnMinimisedDfsmGenerator fsm variant" << endl;
//		return nullptr;
//	}
//	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
//	{
//		cout << "HMethodOnMinimisedDfsmGenerator dfsm variant" << endl;
//		return make_shared<IOListContainer>(m.hMethodOnMinimisedDfsm(numAddStates));
//	}
//};
//
///**
// * Test function for Fsm::wMethod(const unsigned int numAddStates)
// * Each element of mutants is expected to be complete and normalized and to have the same input alphabet as m.
// * nullOutput is expected to be the output that was used as nulloutputs in the completion of the mutants.
// */
//void testTestTheory(Fsm & m, const vector<shared_ptr<const Fsm>>& mutants, const size_t nullOutput, const shared_ptr<TestSuiteGenerator> tsGen) {
//	const Fsm copyOfM = Fsm(m);
//	// calculate numAddStates 
//	vector<shared_ptr<const Fsm>> filteredMutants;
//	int numAddStates = 0;
//	auto completeM = transformToComplete(make_shared<Fsm>(m), nullOutput);
//	auto minComplM = completeM->minimise();
//	if (minComplM.size() > 30) {
//		cout << "FSM too big. Stop Test Case." << endl;
//		return;
//	}
//	const int maxAddStates = 2; //3;	
//
//	// filter out Fsms with too many states
///*	for (const auto mutant : mutants) {
//		int sizeDiff = mutant->size() - minComplM.size();	
//		cout << "maxAddStates: " <<  maxAddStates << endl;
//		cout << "sizeDiff: " << sizeDiff << endl;
//		if (sizeDiff > maxAddStates) continue;		
//		filteredMutants.push_back(mutant);
//		if (sizeDiff > numAddStates) {
//			numAddStates = sizeDiff;
//		}
//	}	*/	
//
//	for (const auto mutant : mutants) {
//		if (mutant->size() > minComplM.size()) {
//			int sizeDiff = mutant->size() - minComplM.size();
//			cout << "sizeDiff: " << sizeDiff << endl;
//			if (sizeDiff > maxAddStates) continue;
//			filteredMutants.push_back(mutant);
//			if (sizeDiff > numAddStates) {
//				numAddStates = sizeDiff;
//			}
//		}
//		else {
//			filteredMutants.push_back(mutant);
//		}
//	}
//
//	// Calculate complete test suite
//	const auto ts = tsGen->generateTestSuite(m, numAddStates);
//
//	//int nullOutput = getMaxOutput(m, mutants) + 1;
//
//	// first check invariant of m
//	bool invariantViolation = not checkFsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//	
//	cout << "ts.size: " << ts->size() << endl;
//	cout << "ts[0].size: " << ts->getIOLists()->at(0).size() << endl;
//	int maxLength = 0;
//	for (auto tc : *ts->getIOLists()) { if (tc.size() > maxLength) maxLength = tc.size(); }
//	cout << "maxLength: " << maxLength << endl;
//	cout << "m.size: " << m.size() << endl;
//	cout << "completeM.size: " << completeM->size() << endl;
//	cout << "minComplM.size: " << minComplM.size() << endl;
//	cout << "filteredMutants.size: " << filteredMutants.size() << endl;
//	cout << "numAddStates: " << numAddStates << endl;
//
//	// stop test case is testcases are too long or test suite contains too many test cases
//	if (maxLength > 12 or ts->size() > 2000) {
//		cout << "Test Suite too big or test cases too long. Stop test case." << endl;
//		return;
//	}
//
//	// save outputs of completeM in hash table
//	//unordered_map < int, set<vector<int>>> outputsOfCompleteM;
//	//for (int i = 0; i < ts->getIOLists()->size(); ++i) {
//	//	outputsOfCompleteM[i] = calcCompleteOutputTraces2(completeM->getInitialState(), ts->getIOLists()->at(i), nullOutput);
//	//}
//
//	// Check completeness of test suite with help of the mutants
//	for (const auto mutant : filteredMutants) {
//		bool diff = false;		
//		for (const auto &tc : *ts->getIOLists()) {   			
//			if (calcCompleteOutputTraces2(completeM->getInitialState(), tc, nullOutput)
//				!= calcCompleteOutputTraces2(mutant->getInitialState(), tc, nullOutput)) {
//				diff = true;
//				break;
//			}
//	/*		if (outputsOfCompleteM.at(i) != calcCompleteOutputTraces2(mutant->getInitialState(), ts->getIOLists()->at(i), nullOutput)) {
//				diff = true;
//				break;
//			}*/
//		}
//		fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()) != diff, "M and mutant are i/o-equivalent iff mutant passed test suite.");
//		//if (not diff) {
//		//	cout << "calcs equivalence" << endl;
//		//	fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()), "M and mutant are i/o-equivalent if mutant passed test suite.");
//		//}
//	}
//
//	// check if language of m has changed
//	/*fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");*/
//	fsmlib_assert("TC", ioEquivalenceCheck(m.getInitialState(), copyOfM.getInitialState()), "Language of M has not changed");
//}
//
///**
// * Test function for Fsm::wMethod(const unsigned int numAddStates)
// * Each element of mutants is expected to be complete and normalized and to have the same input alphabet as m.
// * nullOutput is expected to be the output that was used as nulloutputs in the completion of the mutants.
// */
//void testTestTheory(Dfsm & m, const vector<shared_ptr<const Fsm>>& mutants, const size_t nullOutput, const shared_ptr<TestSuiteGenerator> tsGen) {
//	const Dfsm copyOfM = Dfsm(m);
//	// calculate numAddStates 
//	vector<shared_ptr<const Fsm>> filteredMutants;
//	int numAddStates = 0;
//	auto completeM = transformToComplete(make_shared<Dfsm>(m), nullOutput);
//	auto minComplM = completeM->minimise();
//	if (minComplM.size() > 50) {
//		cout << "FSM too big. Stop Test Case." << endl;
//		return;
//	}
//	const int maxAddStates = 3;
//	//for (const auto mutant : mutants) {
//	//	int sizeDiff = mutant->size() - minComplM.size();
//	//	if (sizeDiff > maxAddStates) continue;
//	//	filteredMutants.push_back(mutant);
//	//	if (sizeDiff > numAddStates) {
//	//		numAddStates = sizeDiff;
//	//	}
//	//}
//	for (const auto mutant : mutants) {
//		if (mutant->size() > minComplM.size()) {
//			int sizeDiff = mutant->size() - minComplM.size();
//			if (sizeDiff > maxAddStates) continue;
//			filteredMutants.push_back(mutant);
//			if (sizeDiff > numAddStates) {
//				numAddStates = sizeDiff;
//			}
//		}
//		else {
//			filteredMutants.push_back(mutant);
//		}
//	}
//	const auto ts = tsGen->generateTestSuite(m, numAddStates);
//
//	// first check invariant of m
//	bool invariantViolation = not checkDfsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	cout << "ts.size: " << ts->size() << endl;
//	cout << "ts[0].size: " << ts->getIOLists()->at(0).size() << endl;
//	int maxLength = 0;
//	for (auto tc : *ts->getIOLists()) { if (tc.size() > maxLength) maxLength = tc.size(); }
//	//cout << "maxLength: " << maxLength << endl;
//	//cout << "m.size: " << m.size() << endl;
//	//cout << "completeM.size: " << completeM->size() << endl;
//	//cout << "minComplM.size: " << minComplM.size() << endl;
//	//cout << "filteredMutants.size: " << filteredMutants.size() << endl;
//
//	for (const auto mutant : filteredMutants) {
//		bool diff = false;
//		for (const auto &tc : *ts->getIOLists()) {
//			if (calcCompleteOutputTraces2(completeM->getInitialState(), tc, nullOutput)
//				!= calcCompleteOutputTraces2(mutant->getInitialState(), tc, nullOutput)) {
//				diff = true;
//				break;
//			}
//		}
//		//fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()) != diff, "M and mutant are i/o-equivalent iff mutant passed test suite.");
//		if (not diff) {
//			fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()), "M and mutant are i/o-equivalent if mutant passed test suite.");
//		}
//	}
//
//	// check if structure of m has changed
//	//fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//	fsmlib_assert("TC", ioEquivalenceCheck(m.getInitialState(), copyOfM.getInitialState()), "Language of M has not changed");
//}
//
//
///**
// * Test function for Fsm::wMethod(const unsigned int numAddStates)
// * Each element of mutants is expected to be complete and normalized and to have the same input alphabet as m.
// * nullOutput is expected to be the output that was used as nulloutputs in the completion of the mutants.
// */
////void testWMethod(Fsm & m, const vector<shared_ptr<const Fsm>>& mutants, const size_t nullOutput) {
////	const Fsm copyOfM = Fsm(m);
////	// calculate numAddStates 
////	vector<shared_ptr<const Fsm>> filteredMutants;
////	size_t numAddStates = 0;
////	auto completeM = transformToComplete(make_shared<Fsm>(m),nullOutput);
////	auto minComplM = completeM->minimise();
////	if (minComplM.size() > 50) {
////		cout << "FSM too big. Stop Test Case." << endl;
////		return;
////	}
////	const size_t maxAddStates = 3;
////	for (const auto mutant : mutants) {
////		int sizeDiff = mutant->size() - minComplM.size();
////		if (sizeDiff > maxAddStates) continue;
////		filteredMutants.push_back(mutant);
////		if (sizeDiff > numAddStates) {
////			numAddStates = sizeDiff;
////		}
////	}
////
////	cout << "vor" << endl;
////	cout << "m.size: " << m.size() << endl;
////	cout << "minComplM.size: " << minComplM.size() << endl;
////	cout << "numAddStates: " << numAddStates << endl;
////	cout << "mI" << m.getMaxInput() << endl;
////	const auto ts = m.wMethod(numAddStates);
////	cout << "nach" << endl;
////
////	//int nullOutput = getMaxOutput(m, mutants) + 1;
////
////	// first check invariant of m
////	bool invariantViolation = not checkFsmClassInvariant(m);
////	fsmlib_assert("TC", not invariantViolation, "Fsm class invariant still holds for M after calculation.");
////	// stop test execution at this point if invariant of m does not hold anymore
////	if (invariantViolation) return;
////
////	for (const auto mutant : filteredMutants) {
////		bool diff = false;
////		for (const auto &tc : *ts.getIOLists()) {
////			if (calcCompleteOutputTraces2(completeM->getInitialState(), tc, nullOutput)
////				!= calcCompleteOutputTraces2(mutant->getInitialState(), tc, nullOutput)) {
////				diff = true;
////				break;
////			}
////		}
////		fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()) != diff, "M and mutant are i/o-equivalent iff mutant passed test suite.");
////		//if (not diff) {
////		//	fsmlib_assert("TC", ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState()), "M and mutant are i/o-equivalent if mutant passed test suite.");
////		//}
////	}
////
////	// check if structure of m has changed
////	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
////}
//
//
////void wMethod_TS_Random() {
////	const int seed = time(NULL);
////	srand(seed);
////
////	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
////	for (int i = 0; i < 200; ++i) {
////		cout << "i:" << i << endl;
////		int size = rand() % 6 + 1; // = 6; 
////		int mI = rand() % 6;
////		int mO = (rand() % 6) + 1; 
////		//auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl);
////		auto m = createPartialMutant( Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl));
////		
////		cout << m << endl;
////
////		vector<shared_ptr<const Fsm>> mutants;
////		for (int j = 0; j < 20; ++j) {
////			size_t numOutFaults = (rand() % 2);
////			size_t numTrFaults = (rand() % 2);
////			if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
////			cout << "numOutFaults: " << numOutFaults << endl;
////			cout << "numTrFaults: " << numTrFaults << endl;
////			cout << "m.mI: " << m->getMaxInput() << endl;
////			cout << "m.mO: " << m->getMaxOutput() << endl;
////			cout << "m.size: " << m->size() << endl;
////			cout << "pre" << endl;
////			auto tmp = m->createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults);
////			cout << "nach mutate" << endl;
////			/*auto minMut = m->createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();*/
////			auto minMut = tmp->minimise();
////			cout << "post" << endl;
////			//mutants.push_back(make_shared<Fsm>(minMut));
////			mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), m->getMaxOutput()));
////		}
////
////		testWMethod(*m, mutants,m->getMaxOutput() +1);
////	}
////}
//
//struct TestTheoryTestCase {
//	string id;
//	string rmPath;
//	vector<string> mutantPaths;
//};
//
//shared_ptr<vector<TestTheoryTestCase>> parseTestTheoryTSFile(const string &testSuitePath) {
//	string fname = testSuitePath;
//	vector<TestTheoryTestCase> testSuite;
//	ifstream inputFile(fname);
//	if (inputFile.is_open())
//	{
//		string line;
//
//		while (getline(inputFile, line))
//		{
//			vector<string> lineContent;
//			stringstream ss(line);
//
//			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
//
//			TestTheoryTestCase tc;
//			tc.id = lineContent.at(0);
//			tc.rmPath = lineContent.at(1);
//			//tc.mutantPaths 
//
//
//			for (int i = 2; i < lineContent.size(); ++i) {
//				tc.mutantPaths.push_back(lineContent.at(i));
//				//cout << lineContent.at(i) << endl;
//				//shared_ptr<Fsm> fsm = make_shared<Fsm>(lineContent.at(i), make_shared<FsmPresentationLayer>(), "M");
//				//fsmlib_assert("TC", checkFsmClassInvariant(*fsm), "inv check");
//				//cout << fsm->size() << endl;
//			}
//			testSuite.push_back(tc);
//		}
//		inputFile.close();
//
//		return make_shared<vector<TestTheoryTestCase>>(testSuite);
//	}
//	else
//	{
//		cout << "Unable to open input file" << endl;
//		exit(EXIT_FAILURE);
//	}
//}
//
//// Test Fsm::wMethod(...)
//void wMethod_TS_Random2() {
//	const int seed = 12792;
//	srand(seed);
//	shared_ptr<WMethodGenerator> tsGenerator = make_shared<WMethodGenerator>();
//
//	// random tests
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	for (int i = 0; i < 1000; ++i) {
//		cout << "i:" << i << endl;
//		int size = (rand() % 6) + 1; // = 6; 
//		int mI = rand() % 4; // (rand() % 5) + 1;
//		int mO = (rand() % 6) + 1;
//		auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl);
//		//auto m = createPartialMutant(Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl));
//		const size_t nullOutput = m->getMaxOutput() + 1;
//
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int j = 0; j < 20; ++j) {
//			size_t numOutFaults = (rand() % 2);
//			size_t numTrFaults = (rand() % 2);
//			if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//			auto minMut = m->createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//			//mutants.push_back(make_shared<Fsm>(minMut));
//			mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput)); // m->getMaxOutput()
//		}
//		testTestTheory(*m, mutants, nullOutput, tsGenerator);
//	}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> ref = make_shared<Fsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Dfsm::wMethod(...)
//void wMethod_Dfsm_TS_Random() {
//	const int seed = 412725;
//	srand(seed);
//
//	shared_ptr<WMethodGenerator> tsGenerator = make_shared<WMethodGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1; 
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//	//	cout << "ic: " << isInitialConnected(m) << endl;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Fsm::wMethodOnMinimisedFsm(...)
//void wMethodOnMinimisedFsm_Fsm_TS_Random() {
//	const int seed = 27161;
//	srand(seed);
//	
//	shared_ptr<WMethodOnMinimisedFsmGenerator> tsGenerator = make_shared<WMethodOnMinimisedFsmGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 50; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 6 + 1; // = 6; 
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	//auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl)->minimise();
//	//	auto m = createPartialMutant(Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl))->minimise();
//	//	//Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//
//	//	cout << m << endl;
//
//	//	// filter fsm that have too many states
//	//	if (m.size() > 30) {
//	//		cout << "FSM too big. Stop Test Case." << endl;
//	//		continue;
//	//	}
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		//mutants.push_back(make_shared<Fsm>(minMut));
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput)); // m->getMaxOutput()
//	//		cout << "added" << endl;
//	//	}
//
//	//	//testWMethod(*m, mutants, m->getMaxOutput() + 1);
//	//	cout << "invoked" << endl;
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	//// hier ne schleife mit dfsm testfällen
//	//for (int i = 0; i < 10; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//	//	m = m.minimise();
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wMethodOnMinimisedFsm.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> ref = make_shared<Fsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// test Dfsm::wMethodOnMinimisedDfsm(...)
//void wMethodOnMinimisedDfsm_Dfsm_TS_Random() {
//	const int seed = 625;
//	srand(seed);
//
//	shared_ptr<WMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<WMethodOnMinimisedDfsmGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 2;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//	//	m = m.minimise();
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	/*auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wMethodOnMinimisedDfsm.testsuite");
//	for (auto tc : *testSuite) {
//		cout << "Start Test Case : " << tc.id << endl;
//		shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//		vector<shared_ptr<const Fsm>> partialMutants;
//		for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//			shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//			partialMutants.push_back(mutant);
//		}
//		size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//		vector<shared_ptr<const Fsm>> completeMutants;
//		for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//		testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	}*/
//}
//
//// test Fsm::wpMethod(...)
//void wpMethod_Fsm_TS_Random() {
//	const int seed = 5217;
//	srand(seed);
//	shared_ptr<WpMethodGenerator> tsGenerator = make_shared<WpMethodGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 50; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 6 + 1; // = 6; 
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl)->minimise();
//	//	//auto m = createPartialMutant(Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl))->minimise();
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//	//	cout << m << endl;
//
//	//	// filter fsm that have too many states
//	//	if (m.size() > 30) {
//	//		cout << "FSM too big. Stop Test Case." << endl;
//	//		continue;
//	//	}
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		cout << "numOutFaults: " << numOutFaults << endl;
//	//		cout << "numTrFaults: " << numTrFaults << endl;
//	//		cout << "m.mI: " << m.getMaxInput() << endl;
//	//		cout << "m.mO: " << m.getMaxOutput() << endl;
//	//		cout << "m.size: " << m.size() << endl;
//	//		cout << "pre" << endl;
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		cout << "post" << endl;
//	//		//mutants.push_back(make_shared<Fsm>(minMut));
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wpMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> ref = make_shared<Fsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Dfsm::wpMethod(...)
//void wpMethod_Dfsm_TS_Random() {
//	const int seed = 712871;
//	srand(seed);
//	shared_ptr<WpMethodGenerator> tsGenerator = make_shared<WpMethodGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wpMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Dfsm::wpMethodOnMinimisedDfsm(...)
//void wpMethodOnMinimisedDfsm_Dfsm_TS_Random() {
//	const int seed = 625166;
//	srand(seed);
//	shared_ptr<WpMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<WpMethodOnMinimisedDfsmGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//	//	m = m.minimise();
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wpMethodOnMinimisedDfsm.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// test Fsm::hsiMethod(...)
//void hsiMethod_Fsm_TS_Random() {
//	const int seed = 12;
//	srand(seed);
//	shared_ptr<HsiMethodGenerator> tsGenerator = make_shared<HsiMethodGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 50; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 6 + 1; // = 6; 
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	//auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl);
//	//	auto m = createPartialMutant(Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl))->minimise();
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	// filter fsm that have too many states
//	//	if (m.size() > 30) {
//	//		cout << "FSM too big. Stop Test Case." << endl;
//	//		continue;
//	//	}
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		//mutants.push_back(make_shared<Fsm>(minMut));
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_hsiMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Fsm> ref = make_shared<Fsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Dfsm::hsiMethod(...)
//void hsiMethod_Dfsm_TS_Random() {
//	const int seed = 907;
//	srand(seed);
//	shared_ptr<HsiMethodGenerator> tsGenerator = make_shared<HsiMethodGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_hsiMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
//// Test Dfsm::hMethodOnMinimisedDfsm(...)
//void hMethodOnMinimisedDfsm_Dfsm_TS_Random() {
//	const int seed = 3234;
//	srand(seed);
//	shared_ptr<HMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<HMethodOnMinimisedDfsmGenerator>();
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 20 + 2;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1;
//	//	Dfsm m("M", size, mI, mO, pl, true); m.setMaxState(m.getNodes().size() - 1);
//	//	m = m.minimise();
//
//	//	if (m.size() < 2) {
//	//		cout << "M has less than 2 states. Stop test case." << endl;
//	//	}
//
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 2);
//	//		size_t numTrFaults = (rand() % 2);
//	//		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
//	//		auto minMut = m.createMutantRepeatable("Mutant_" + j, numOutFaults, numTrFaults)->minimise();
//	//		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
//	//	}
//	//	testTestTheory(m, mutants, nullOutput, tsGenerator);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// CSM is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// FSB is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			// GDC is complete => created mutants are complete => minimised mutants are complete
//			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, rand() % 6 + 1)->minimise()));
//		}
//		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator);
//	}
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_hMethodOnMinimisedDfsm.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> partialMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		partialMutants.push_back(mutant);
//	//	}
//	//	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;
//
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
//	//	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator);
//	//}
//}
//
///**
// * Test function for Dfsm::tMethod().
// * m is expected to be complete. Each element of mutants is expected to differ from m only by zero or more output faults.
// */
//void testTMethod(Dfsm & m, const vector<shared_ptr<const Fsm>>& mutants) {
//	const Dfsm copyOfM = Dfsm(m);
//
//	const auto ts = m.tMethod();
//
//	//int nullOutput = getMaxOutput(m, mutants) + 1;
//
//	// first check invariant of m
//	bool invariantViolation = not checkDfsmClassInvariant(m);
//	fsmlib_assert("TC", not invariantViolation, "Dfsm class invariant still holds for M after calculation.");
//	// stop test execution at this point if invariant of m does not hold anymore
//	if (invariantViolation) return;
//
//	for (const auto mutant : mutants) {
//		bool diff = false;
//		for (const auto &tc : *ts.getIOLists()) {
//			if (calcCompleteOutputTraces2(m.getInitialState(), tc, m.getMaxOutput())
//				!= calcCompleteOutputTraces2(mutant->getInitialState(), tc, mutant->getMaxOutput())) {
//				diff = true;
//				break;
//			}
//		}
//		fsmlib_assert("TC", ioEquivalenceCheck(m.getInitialState(), mutant->getInitialState()) != diff, "M and mutant are i/o-equivalent iff mutant passed test suite.");
//		//if (not diff) {
//		//	fsmlib_assert("TC", ioEquivalenceCheck(m.getInitialState(), mutant->getInitialState()), "M and mutant are i/o-equivalent if mutant passed test suite.");
//		//}
//	}
//
//	// check if structure of m has changed
//	fsmlib_assert("TC", checkForEqualStructure(m, copyOfM), "M was not changed by algorithm");
//}
//
///*
// *	Random Test Suite for test of Dfsm::tMethod().
//*/
//void tMethod_TS_Random() {
//	const int seed = 876298;
//	srand(seed);
//
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//
//	// random tests
//	//for (int i = 0; i < 1000; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 15 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1; // cases with maxOutput = 0 are trivial
//	//	auto m = Dfsm("M", size, mI, mO, pl, true);
//	//	m.setMaxState(m.getNodes().size() - 1);
//	//	cout << m << endl;
//	//	vector<shared_ptr<const Fsm>> mutants;
//	//	for (int j = 0; j < 20; ++j) {
//	//		size_t numOutFaults = (rand() % 4) + 1;
//	//		cout << "numOutFaults: " << numOutFaults << endl;
//	//		mutants.push_back(m.createMutantRepeatable("Mutant_" + j, numOutFaults, 0));
//	//	}
//	//	cout << "here" << endl;
//	//	testTMethod(m, mutants);
//	//}
//
//	// practical examples
//	{
//		cout << "CSM Tests" << endl;
//		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			mutants.push_back(csm->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, 0));
//		}
//		testTMethod(*csm, mutants);
//	}
//	{
//		cout << "FSB Tests" << endl;
//		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			mutants.push_back(fsb->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, 0));			
//		}
//		testTMethod(*fsb, mutants);
//	}
//	{
//		cout << "GDC Tests" << endl;
//		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
//		vector<shared_ptr<const Fsm>> mutants;
//		for (int i = 0; i < 100; ++i) {
//			mutants.push_back(gdc->createMutantRepeatable("Mutant_" + i, rand() % 6 + 1, 0));			
//		}
//		testTMethod(*gdc, mutants);
//	}
//
//
//	// partition tests
//	//auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_tMethod.testsuite");
//	//for (auto tc : *testSuite) {
//	//	cout << "Start Test Case : " << tc.id << endl;
//	//	shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
//	//	vector<shared_ptr<const Fsm>> completeMutants;
//	//	for (int i = 0; i < tc.mutantPaths.size(); ++i) {
//	//		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + i);
//	//		completeMutants.push_back(mutant);
//	//	}
//	//	testTMethod(*ref, completeMutants);
//	//}
//}
//
//// ====================================================================================================
//
//
////void createMutants() {
////	// set accordingly
////	const string path = "../../../resources/TestSuites/TestTheories/RM31.fsm";
////	const string writePathPref = "../../../resources/TestSuites/TestTheories/Mut32_";
////	const int numAddStates = 1;
////	const int sizeOfMinRef = 4;
////
////	shared_ptr<Fsm> fsm = make_shared<Fsm>(path, make_shared<FsmPresentationLayer>(), "M");
////    cout << "initIdx: " << fsm->getInitStateIdx() << endl;
////    cout << "mI: " << fsm->getMaxInput() << endl;
////	cout << "mO:" << fsm->getMaxOutput() << endl;
////	cout << "size: " << fsm->getNodes().size() << endl;
////	cout << "comp. spec: " << fsm->isCompletelyDefined() << endl;
////	cout << "det: " << fsm->isDeterministic() << endl;
////	cout << "obs: " << fsm->isObservable() << endl;
////	cout << checkFsmClassInvariant(*fsm) << endl;
////
////
////	// create Mutants
////	srand(time(NULL));
////	cout << rand()<< endl;
////	int numMutants = 0;
////	vector<shared_ptr<Fsm>> mutants;
////	while (numMutants < 5) {
////		size_t numOutFaults = rand() % 3;
////		size_t numTrFaults = rand() % 3;
////		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;
////		cout << "numOutFaults: " << numOutFaults << endl;
////		cout << "numTrFaults: " << numTrFaults << endl;
////		auto mutant = fsm->createMutant("Mut", numOutFaults, numTrFaults)->minimise();
////		if (mutant.size() == sizeOfMinRef + numAddStates) {//(mutant.size() <= sizeOfMinRef + numAddStates) {
////			cout << "created mutant with size " << mutant.size() << endl;
////			++numMutants;
////			mutants.push_back(make_shared<Fsm>(mutant));
////		}
////		else {
////			cout << "too big" << endl;
////		}
////	}
////
////	//vector<shared_ptr<Fsm>> remaining;
////	//for (auto m : mutants) {
////	//	bool found = false;
////	//	for (auto mu : remaining) {
////	//		if (not checkForEqualStructure(*m, *mu)) {
////	//			remaining.push_back(m)
////	//		}
////	//	}
////	//}
////
////	// write mutants 
////	for (int i = 0; i < mutants.size(); ++i) {
////		std::stringstream writePath;
////		writePath << writePathPref << i << ".fsm";
////		ofstream outFile(writePath.str());
////		mutants.at(i)->dumpFsm(outFile);
////		outFile.close();
////	}
////}
//
////void createMutantsForTMethod() {
////	// set accordingly
////	const string path = "../../../resources/TestSuites/TestTheories/RM25.fsm";
////	const string writePathPref = "../../../resources/TestSuites/TestTheories/Mut25_t_";
////	const int numAddStates = 0;
////	const int sizeOfMinRef = 3;
////
////	shared_ptr<Fsm> fsm = make_shared<Fsm>(path, make_shared<FsmPresentationLayer>(), "M");
////    cout << "initIdx: " << fsm->getInitStateIdx() << endl;
////    cout << "mI: " << fsm->getMaxInput() << endl;
////	cout << "mO:" << fsm->getMaxOutput() << endl;
////	cout << "size: " << fsm->getNodes().size() << endl;
////	cout << "comp. spec: " << fsm->isCompletelyDefined() << endl;
////	cout << "det: " << fsm->isDeterministic() << endl;
////	cout << "obs: " << fsm->isObservable() << endl;
////	cout << checkFsmClassInvariant(*fsm) << endl;
////
////
////	// create Mutants
////	srand(time(NULL));
////	cout << rand()<< endl;
////	int numMutants = 0;
////	vector<shared_ptr<Fsm>> mutants;
////	while (numMutants < 5) {
////		size_t numOutFaults = 1;
////		size_t numTrFaults = 0;
////		
////		cout << "numOutFaults: " << numOutFaults << endl;
////		cout << "numTrFaults: " << numTrFaults << endl;
////		auto mutant = fsm->createMutant("Mut", numOutFaults, numTrFaults);
////		
////			
////			++numMutants;
////			mutants.push_back(mutant);
////
////	}
////
////	//vector<shared_ptr<Fsm>> remaining;
////	//for (auto m : mutants) {
////	//	bool found = false;
////	//	for (auto mu : remaining) {
////	//		if (not checkForEqualStructure(*m, *mu)) {
////	//			remaining.push_back(m)
////	//		}
////	//	}
////	//}
////
////	// write mutants 
////	for (int i = 0; i < mutants.size(); ++i) {
////		std::stringstream writePath;
////		writePath << writePathPref << i << ".fsm";
////		ofstream outFile(writePath.str());
////		mutants.at(i)->dumpFsm(outFile);
////		outFile.close();
////	}
////}
//
//
//// Selects two nodes of m randomly and changes the transitions of one of those nodes in a way that makes both equivalent.
//// It is expected that srand() was called before.
//shared_ptr<Fsm> makeStatesEquivalent(const Fsm &m) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m.getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m.getName(), m.getPresentationLayer()));
//	}
//
//	unsigned int first = rand() % m.getNodes().size();
//	unsigned int second = rand() % m.getNodes().size();
//	if (first == second) second = (second + 1) % m.getNodes().size();
//	cout << "first: " << first << endl;
//	cout << "second: " << second << endl;
//	
//
//	// Now add transitions that correspond exactly to the transitions in m,
//	for (int n = 0; n <= m.getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m.getNodes()[n];
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//	}
//	// replaced transitions of lst[second] by transitions with the same labels and targets as the transitions of 
//    // lst[first]. This results in two equivalent nodes iff first != second.
//	if (first != second) {
//		lst[second]->getTransitions().clear();
//		for (const auto tr : lst[first]->getTransitions()) {
//			lst[second]->getTransitions().push_back(make_shared<FsmTransition>(lst[second], tr->getTarget(), make_shared<FsmLabel>(*(tr->getLabel()))));
//		}
//	}
//	unsigned int mI = 0;
//	unsigned int mO = 0;
//	for (const auto n : lst) {
//		for (const auto tr : n->getTransitions()) {
//			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
//			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
//		}
//	}
//	cout << "mI: " << mI << endl;
//	cout << "mO: " << mO << endl;
//	return make_shared<Fsm>(m.getName(), mI, mO, lst, m.getPresentationLayer());
//}
//
//// Selects some node of m at random and deletes all incoming transitions of this node in the new Fsm.
//// The returned Fsm contains unreachable nodes if the selected node is not the initial node.
//// It is expected that srand() was called before.
//shared_ptr<Fsm> makeStatesUnreachable(const Fsm &m) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m.getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m.getName(), m.getPresentationLayer()));
//	}
//
//	unsigned int nodeIdx = rand() % m.getNodes().size();
//	cout << "n: " << nodeIdx << endl;
//
//	// Now add transitions that correspond exactly to the transitions in m, but ignore transitions
//	// to node at nodeIdx
//	for (int n = 0; n <= m.getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m.getNodes()[n];
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			int tgtId = tr->getTarget()->getId();
//			if (tgtId == nodeIdx) continue;
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//	}
//	unsigned int mI = 0;
//	unsigned int mO = 0;
//	for (const auto n : lst) {
//		for (const auto tr : n->getTransitions()) {
//			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
//			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
//		}
//	}
//	cout << "mI: " << mI << endl;
//	cout << "mO: " << mO << endl;
//	return make_shared<Fsm>(m.getName(), mI, mO, lst, m.getPresentationLayer());
//}
//
//// Selects some input in each state of m and removes each transition for this input.
//// It is expected that srand() was called before.
//shared_ptr<Fsm> makeStatesPartial(shared_ptr<Fsm> m) {
//	vector<shared_ptr<FsmNode> > lst;
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
//	}
//
//	// Now add transitions that correspond exactly to the transitions in
//	// m, but ignore transitions with some randomly selected input
//	for (int n = 0; n <= m->getMaxState(); n++) {
//		auto theNewFsmNodeSrc = lst[n];
//		auto theOldFsmNodeSrc = m->getNodes()[n];
//		int ignoreInput = rand() % (m->getMaxInput() + 1);
//		cout << "ignoreInput: " << ignoreInput << endl;
//		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
//			if (tr->getLabel()->getInput() == ignoreInput) continue;
//			int tgtId = tr->getTarget()->getId();
//			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
//			shared_ptr<FsmTransition> newTr =
//				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
//			theNewFsmNodeSrc->addTransition(newTr);
//		}
//	}
//	unsigned int mI = 0;
//	unsigned int mO = 0;
//	for (const auto n : lst) {
//		for (const auto tr : n->getTransitions()) {
//			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
//			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
//		}
//	}
//	cout << "mI: " << mI << endl;
//	cout << "mO: " << mO << endl;
//	return make_shared<Fsm>(m->getName(), mI, mO, lst, m->getPresentationLayer());	
//}
//
//void randomFSMTestData() {
//	const int seed = time(NULL);//1234;
//	srand(seed);
//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
//	//size_t obsC = 0, icC = 0, reducedC = 0, minimalC = 0, normC = 0, complC = 0, obsAndMiZero = 0, size1C = 0;
//	//for (int i = 0; i < 1000; ++i) {
//	//	int size = rand() % 10 + 1;
//	//	int mI = rand() % 6;
//	//	int mO = (rand() % 6) + 1; // cases with maxOutput = 0 are trivial
//	//	auto m = Dfsm("M", size, mI, mO, pl, true);
//	//	cout << "mI: " << m.getMaxInput() << endl;
//	//	cout << "mO: " << m.getMaxOutput() << endl;
//	//	cout << "size: " << m.getNodes().size() << endl;
//	//	if (m.isObservable()) ++obsC;
//	//	if (isInitialConnected(m)) ++icC;
//	//	if (not hasEquivalentStates(m)) ++reducedC;
//	//	if (isInitialConnected(m) && not hasEquivalentStates(m)) ++minimalC;
//	//	if (m.isObservable() && isInitialConnected(m) && not hasEquivalentStates(m)) ++normC;
//	//	if (m.isCompletelyDefined()) ++complC;
//	//	if (m.isObservable() and m.getMaxInput() == 0) ++obsAndMiZero;
//	//	if (m.getNodes().size() == 1) ++size1C;
//
//	//	cout << "=====================================================" << endl;
//	//}
//	//cout << "obsC: " << obsC << endl;
//	//cout << "icC: " << icC << endl;
//	//cout << "reducedC: " << reducedC << endl;
//	//cout << "minimalC: " << minimalC << endl;
//	//cout << "normC: " << normC << endl;
//	//cout << "complC: " << complC << endl;
//	//cout << "obsAndMiZero: " << obsAndMiZero << endl;
//	//cout << "size1C: " << size1C << endl;
//
//
//	//vector<shared_ptr<Fsm>> fsmList;
//	//for (int i = 0; i < 100; ++i) {
//	//	cout << "i:" << i << endl;
//	//	int size = rand() % 6 + 1; // = 6; 
//	//	int mI = rand() % 5 + 1;
//	//	int mO = (rand() % 6) + 1;
//	//	//auto m = Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl);
//	//	auto m = createPartialMutant(Fsm::createRandomFsmRepeatable("M", mI, mO, size, pl))->minimise();
//	//	for (auto fsm : fsmList) {
//	//		if (ioEquivalenceCheck(fsm->getInitialState(), m.getInitialState())) {
//	//			cout << "equal language found" << endl;
//	//		}
//	//	}
//	//fsmList.push_back(make_shared<Fsm>(m));
//	//	const size_t nullOutput = m.getMaxOutput() + 1;
//	//	cout << m << endl;
//	//	cout << "mI: " << m.getMaxInput() << endl;
//	//	cout << "mO: " << m.getMaxOutput() << endl;
//	//	cout << "size: " << m.getNodes().size() << endl;
//	//	cout << "=====================================================" << endl;
//	//}
//
//	size_t obsC = 0, icC = 0, reducedC = 0, minimalC = 0, normC = 0, complC = 0, obsAndMiZero = 0, size1C = 0;
//	for (int i = 0; i < 10; ++i) {
//		cout << "i:" << i << endl;
//		int size = (rand() % 6) + 1; // = 6; 
//		int mI = rand() % 4 + 1; //rand() % 4;
//		int mO = (rand() % 6) + 1;
//		auto m = Fsm::createRandomFsmRepeatable("M",mI, mO, size, pl);
//		cout << "before_mI: " << m->getMaxInput() << endl;
//		cout << "before_mO: " << m->getMaxOutput() << endl;
//		cout << "before_size: " << m->getNodes().size() << endl;
//
//		//m = makeStatesPartial(makeStatesEquivalent(*m));
//		m = makeStatesEquivalent(*makeStatesPartial(m));
//		cout << "Fsm Inv: " << checkFsmClassInvariant(*m) << endl;
//		cout << *m << endl;
//
//		cout << "mI: " << m->getMaxInput() << endl;
//		cout << "mO: " << m->getMaxOutput() << endl;
//		cout << "size: " << m->getNodes().size() << endl;
//		if (m->isObservable()) ++obsC;
//		if (isInitialConnected(*m)) ++icC;
//		if (not hasEquivalentStates(*m)) ++reducedC;
//		if (isInitialConnected(*m) && not hasEquivalentStates(*m)) ++minimalC;
//		if (m->isObservable() && isInitialConnected(*m) && not hasEquivalentStates(*m)) ++normC;
//		if (m->isCompletelyDefined()) ++complC;
//		if (m->isObservable() and m->getMaxInput() == 0) ++obsAndMiZero;
//		if (m->getNodes().size() == 1) ++size1C;
//
//		cout << "=====================================================" << endl;
//	}
//	cout << "obsC: " << obsC << endl;
//	cout << "icC: " << icC << endl;
//	cout << "reducedC: " << reducedC << endl;
//	cout << "minimalC: " << minimalC << endl;
//	cout << "normC: " << normC << endl;
//	cout << "complC: " << complC << endl;
//	cout << "obsAndMiZero: " << obsAndMiZero << endl;
//	cout << "size1C: " << size1C << endl;
//}




int main(int argc, char** argv)
{


	std::cout << "test start" << std::endl;
	cout << argc << endl;
	cout << argv[0] << endl;
	//if (argc == 2) {
	//	string input(argv[1]);

	vector<TestResult> results;

	// FSM Transformation Tests:
	//if(argc == 2 and input == "1") //
	
	//removeUnreachableNodes_TS().printResults();
	//transformToObservableFSM_TS().printResults();
	minimise_Dfsm_TS().printResults();
	//minimiseObservableFSM_TS().printResults();
	//minimise_Fsm_TS().printResults();

	// Intersection Test:
	//intersect_TS().printResults();

	//// Calculation of Distinguishing Traces Test:
	//getCharacterisationSet_Dfsm_TS().printResults();
	//getCharacterisationSet_Fsm_TS().printResults();
	//calcDistinguishingTrace_PkTables_TS().printResults();
	//calcDistinguishingTrace_OFSMTables_TS().printResults();
	//calcStateIdentificationSets_TS().printResults();
	//calcStateIdentificationSetsFast_TS().printResults();

	//// Complete Test Theories Test:
	//tMethod_TS().printResults();
	//wMethod_Fsm_TS().printResults();
	//wMethod_Dfsm_TS().printResults();
	//wMethodOnMinimisedFsm_TS().printResults();
	//wMethodOnMinimisedDfsm_TS().printResults();
	//wpMethod_Fsm_TS().printResults();
	//wpMethod_Dfsm_TS().printResults();
	//wpMethodOnMinimisedDfsm_TS().printResults();
	//hsiMethod_Fsm_TS().printResults();
	//hsiMethod_Dfsm_TS().printResults();
	//hMethodOnMinimisedDfsm_TS().printResults();

	for (auto &result : results) result.printResults();



	//testIOEquivalenceCheck();
	//testCheckForEqualStructure();
	//testGetReachableStates();
	//testFsmClassInvariant();
	//testCheckDfsmClassInvariant();
    
#if 0
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    test9();
    test10();
    test10b();
    test11();
    test13();
    test14();
    test15();

    faux();

    
    gdc_test1();
    
    

    wVersusT();

    if ( argc < 6 ) {
        cerr << endl <<
        "Missing file names - exit." << endl;
        exit(1);
    }
    
    
    
    string fsmName(argv[1]);
    string fsmFile(argv[2]);
    
    /*
    string inputFile(argv[3]);
    string outputFile(argv[4]);
    string stateFile(argv[5]);
    */
    
    /* Create the presentation layer */
    //shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>(inputFile,outputFile,stateFile);
    
    shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
    
    /* Create an Fsm instance, using the transition relation file,
     * the presentation layer, and the FSM name
     */
    shared_ptr<Fsm> fsm = make_shared<Fsm>(fsmFile,pl,fsmName);
    
    /* Produce a GraphViz (.dot) representation of the created FSM */
    fsm->toDot(fsmName);
    
    /* Transform the FSM into an equivalent observable one */
    Fsm fsmObs = fsm->transformToObservableFSM();
    
    /* Output the observable FSM to a GraphViz file (.dot-file) */
    fsmObs.toDot(fsmObs.getName());
    
#endif
	//testTreeNodeAddConstInt1();
	//testTreeNodeAddConstInt2();
	//testTreeNodeAddConstInt3();
	//testTreeNodeEqualOperator1();
	//testTreeNodeEqualOperator2();
	//testTreeNodeCalcLeaves();
	//testTreeNodeClone();
	//testTreeNodeGetPath();
	//testTreeNodeSuperTreeOf1();
	//testTreeNodeSuperTreeOf2();
	//testTreeNodeTraverse();
	//testTreeNodeDeleteNode();
	//testTreeNodeDeleteSingleNode();
	//testAddToThisNode();
	//testTreeNodeAddIOListContainer();
	//testTreeNodeTentativeAddToThisNode();
	//testTreeNodeAfter();
	//testTreeNodeAddToThisNodeIOListContainer();

	//testTreeRemove();
	//testTreeToDot();
	//testTreeGetPrefixRelationTree();
	//testTreeTentativeAddToRoot();

	//testFsmPresentationLayerFileConstructor();
	//testFsmPresentationLayerDumpIn();
	//testFsmPresentationLayerComparePositive();
	//testFsmPresentationLayerCompareNegative();
	
	//testTraceEquals1Positive();
	//testTraceEquals1Negative();
	//testTraceEquals2Positive();
	//testTraceEquals2Negative();
	//testTraceOutputOperator();

	//testInputTraceOutputOperator();
	
	//testOutputTraceOutputOperator();

	//testOutputTreeContainsNegative();
	//testOutputTreeContainsPositive();
	//testOutputTreeToDot();
	//testOutputTreeGetOutputTraces();
	//testOutputTreeOutputOperator();

	//testTestSuiteIsEquivalentToPositive();
	//testTestSuiteIsEquivalentToNegative();
	//testTestSuiteIsReductionOfNegative();
	//testTestSuiteIsReductionOfPositive();

	//testIOListContainerConstructor();

	//testTraceSegmentGetCopy();

	//testSegmentedTraceEqualOperatorPositive();
	//testSegmentedTraceEqualOperatorNegative();

	//testHittingSetConstructor();
	//testHittingSetCalcMinCardHittingSet();

	//testHsTreeNodeIsHittingSetPositive();
	//testHsTreeNodeIsHittingSetNegative();
	//testHsTreeNodeExpandNode();
	//testHsTreeNodeToDot();

	//testInt2IntMapConstructor();

	//testPkTableRowIsEquivalentPositive();
	//testPkTableRowIsEquivalentNegative();

	//testPkTableMaxClassId();
	//testPkTableGetPkPlusOneTable();
	//testPkTableGetMembers();
	//testPkTableToFsm();

	//testDFSMTableGetP1Table();

	//testOFSMTableRowConstructor();
	//testOFSMTableIoEqualsPositive();
	//testOFSMTableIoEqualsNegative();
	//testOFSMTableClassEqualsPositive();
	//testOFSMTableClassEqualsNegative();

	//testOFSMTableConstructor();
	//testOFSMTableMaxClassId();
	//testOFSMTableNext();
	//testOFSMTableToFsm();

	//testFsmLabelOperatorLessThan();

	//testFsmNodeAddTransition();
	//testFsmNodeApply();
	//testFsmNodeAfter1();
	//testFsmNodeAfter2();
	//testFsmNodeGetDFSMTableRow();
	//testFsmNodeDistinguishedPositive();
	//testFsmNodeDistinguishedNegative();
	//testFsmNodeCalcDistinguishingTrace1();
	//testFsmNodeCalcDistinguishingTrace2();
	//testFsmNodeIsObservable();
	//testFsmNodeIsDeterministic();

    //testFsmAccept();
	//testFsmDeepCopyConstructor();
	//testFsmConstructor1();
	//testFsmConstructor2();
	//testFsmCreateRandomFsm();
	//testFsmCreateMutant();
	//testFsmDumpFsm();


	/*testMinimise();
	testWMethod();*/
	//testCharacterisationSet();
	//testGetDistTraces();
	//testHMethod();
	//testWpMethodWithDfsm(); 
	//testIntersectionCharacteristics();
    //test1();
    //test2();
    //test3();
    //test4();
    //test5();
    //test6();
    //test7();
    //test8();
    //test9();
    //test10();
    //test10b();
    //test11();
    //test13();
    //test14();
    //test15();
    exit(0);
    
}



