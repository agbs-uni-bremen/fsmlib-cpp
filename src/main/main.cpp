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

bool checkDistinguishingCond(Dfsm &minimized) {
	const std::vector<shared_ptr<FsmNode>> nodes = minimized.getNodes();
	for (size_t i = 0; i < nodes.size() - 1; ++i) {
		for (size_t j = i + 1; j < nodes.size(); ++j) {
			if (!minimized.distinguishable(*nodes[i], *nodes[j])) {
				return false;
			}
		}
	}
	return true;
}


void testMinimise() {
	cout << "TC-DFSM-0017 Show that Dfsm::minimise() produces an "
		<< "equivalent minimal FSM"
		<< endl;

	auto pl = make_shared<FsmPresentationLayer>();
	auto dfsm = make_shared<Dfsm>("DFSM", 50, 5, 5, pl);
	Dfsm minimized = dfsm->minimise();
	std::vector<shared_ptr<FsmNode>> unreachableNodes;

	// check for unreachable nodes
	fsmlib_assert("TC-DFSM-0017",
		not minimized.removeUnreachableNodes(unreachableNodes),
		"Minimized Dfsm doesn't contain unreachable nodes");

	// check if states are distinguishable
	fsmlib_assert("TC-DFSM-0017",
		checkDistinguishingCond(minimized),
		"Each node pair of the minimized Dfsm is distinguishable");

	// check language equality
	fsmlib_assert("TC-DFSM-0017",
		minimized.intersect(*dfsm).isCompletelyDefined(),
		"Language of minimized Dfsm equals language of unminimized Dfsm");
}

void testWMethod() {
	cout << "TC-DFSM-0018 Show that Dfsm implModel only passes W-Method Testsuite "
		<< "if intersection is completely defined"
		<< endl;

	auto pl = make_shared<FsmPresentationLayer>();
	auto refModel = make_shared<Dfsm>("refModel", 50, 5, 5, pl);
	auto implModel = make_shared<Dfsm>("implModel", 50, 5, 5, pl);
	IOListContainer iolc = refModel->wMethod(0);

	// check language equality with W-Method Testsuite
	bool equal = true;
	for (auto trc : *(iolc.getIOLists())) {
		shared_ptr<InputTrace> iTr =
			make_shared<InputTrace>(trc, pl);
		if (not implModel->pass(refModel->applyDet(*iTr))) {
			equal = false;
			break;
		}
	}

	fsmlib_assert("TC-DFSM-0018",
		refModel->intersect(*implModel).isCompletelyDefined() == equal,
		"implModel passes W-Method Testsuite if and only if intersection is completely defined");
}

// Checks if tr1 is a prefix of tr2.
bool isPrefix(const std::vector<int> &tr1, const std::vector<int> &tr2) {
	if (tr1.size() > tr2.size()) {
		return false;
	}
	for (size_t i = 0; i < tr1.size(); ++i) {
		if (tr1[i] != tr2[i]) {
			return false;
		}
	}
	return true;
}

//Checks if ot1 is part of ot2. This is only the case if the InputTrace of ot1 is a prefix
//of the InputTrace of ot2 and if every OutputTrace of ot1 is a prefix of an OutputTrace
//of ot2.
bool containsOuputTree(OutputTree &ot1, OutputTree &ot2) {
	InputTrace it1 = ot1.getInputTrace();
	InputTrace it2 = ot2.getInputTrace();
	if (not isPrefix(it1.get(), it2.get())) {
		return false;
	}

	for (OutputTrace &outTr1 : ot1.getOutputTraces()) {
		bool prefix(false);
		for (OutputTrace &outTr2 : ot2.getOutputTraces()) {
			if (isPrefix(outTr1.get(), outTr2.get())) {
				prefix = true;
				break;
			}
		}
		if (not prefix) {
			return false;
		}
	}
	return true;
}

//void testIntersect() {
//	cout << "TC-DFSM-0019 Show that Fsm::intersect() produces FSM which accepts intersection "
//		<< "of the languages from the original FSMs"
//		<< endl;
//
//	auto pl = make_shared<FsmPresentationLayer>();
//	auto m1 = make_shared<Dfsm>("refModel", 3, 3, 3, pl)->minimise();
//	auto m2 = m1.createMutant("mutant", 2, 2)->minimise();
//	Fsm intersection = m1.intersect(m2).minimise();
//
//	std::cout << "m1.isDeterministic(): " << m1.isDeterministic() << std::endl;
//	std::cout << "m1.isCompletelyDefined(): " << m1.isCompletelyDefined() << std::endl;
//
//	std::cout << "m2.isDeterministic(): " << m2.isDeterministic() << std::endl;
//	std::cout << "m2.isCompletelyDefined(): " << m2.isCompletelyDefined() << std::endl;
//
//	std::cout << "intersection.isDeterministic(): " << intersection.isDeterministic() << std::endl;
//	std::cout << "intersection.isCompletelyDefined(): " << intersection.isCompletelyDefined() << std::endl;
//
//	std::cout << "m1.size(): " << m1.size();
//	std::cout << "intersection.size(): " << intersection.size();
//
//	// show that every trace of length n+m-1 in language of intersection is in language of m1 and m2 
//	IOListContainer iolc = IOListContainer(m1.getMaxInput(),
//			1,
//			intersection.size() + m1.size() - 1,
//			pl);
//
//	std::cout << "iolc created" << std::endl;
//
//	int c = 0;
//	for (auto trc : *(iolc.getIOLists())) {
//		if (c++ >= 10) break;
//		shared_ptr<InputTrace> iTr =
//			make_shared<InputTrace>(trc, pl);
//		OutputTree ot1 = intersection.apply(*iTr);
//		OutputTree ot2 = m1.apply(*iTr);
//		OutputTree ot3 = m2.apply(*iTr);		
//
//		std::cout << "ot2.contains(ot1): " << ot2.contains(ot1) << std::endl;
//		std::cout << "ot3.contains(ot1): " << ot3.contains(ot1) << std::endl;
//
//		std::cout << "containsOuputTree(ot1, ot2)" << containsOuputTree(ot1, ot2) << std::endl;
//		std::cout << "containsOuputTree(ot1, ot3)" << containsOuputTree(ot1, ot3) << std::endl;
//
//		std::cout << "intersection tr:\t" << ot1 << std::endl;
//		std::cout << std::endl;
//		std::cout << "ot2:\t\t\t" << ot2 << std::endl;
//		std::cout << std::endl;
//		std::cout << "ot3:\t\t\t" << ot3 << std::endl;
//	}
//}

void testIntersectionCharacteristics() {
	auto pl = make_shared<FsmPresentationLayer>();
	auto m1 = make_shared<Dfsm>("m1", 10, 3, 3, pl)->minimise();
	auto m2 = m1.createMutant("m2", 2, 2);

	fsmlib_assert("TC-DFSM-0019b",
		m1.intersect(*m2).isDeterministic(),
		"m1 or m2 deterministic => product automata deterministic");

	auto m3 = Fsm::createRandomFsm("m3", 3, 3, 3, pl);
	auto m4 = Fsm::createRandomFsm("m4", 3, 3, 3, pl);
	Fsm intersection = m3->intersect(*m4);
	if (not intersection.isDeterministic()) {
		fsmlib_assert("TC-DFSM-0019b",
			(not m3->isDeterministic()) and (not m4->isDeterministic()),
			"product automata of m3 and m4 nondeterministic => m3 and m4 nondeterministic");
	}
	if (intersection.isCompletelyDefined()) {
		fsmlib_assert("TC-DFSM-0019b",
			m3->isCompletelyDefined() and m4->isCompletelyDefined(),
			"product automata of m3 and m4 completely specified => m3 and m4 completely specified");
	}

	// TODO: m1 or m2 incomplete specified => product automata incomplete specified
}

bool equalSetOfOutputTrees(std::vector<OutputTree> &otv1, std::vector<OutputTree> &otv2) {
	if (otv1.size() != otv2.size()) {
		return false;
	}
	for (size_t i = 0; i < otv1.size(); ++i) {
		if (otv1[i] != otv2[i]) {
			return false;
		}
	}
	return true;
}

void testCharacterisationSet() {
	cout << "TC-FSM-0019 Show that calculated characterisation set "
		<< "distinguishes each pair of FSM states"
		<< endl;
	auto pl = make_shared<FsmPresentationLayer>();
	auto m1 = Fsm::createRandomFsm("M1", 3, 3, 10, pl)->minimise();
	IOListContainer iolc = m1.getCharacterisationSet();

	// calculate output trees for every node
	std::vector<std::vector<OutputTree>> outputTrees;
	for (const auto &node : m1.getNodes()) {		
		std::vector<OutputTree> traces;
		for (const auto trc : *(iolc.getIOLists())) {
			shared_ptr<InputTrace> iTr =
				make_shared<InputTrace>(trc, pl);
			traces.push_back(node->apply(*iTr));
		}
		outputTrees.push_back(traces);
	}

	// check if vector contains equal sets of output trees
	for (size_t i = 0; i < m1.getNodes().size() - 1; ++i) {
		for (size_t j = i + 1; j < m1.getNodes().size(); ++j) {
			if (equalSetOfOutputTrees(outputTrees[i], outputTrees[j])) {
				std::cout << "============= FAIL ==============" << std::endl;
			}
		}
	}
	std::cout << "============= PASS ============" << std::endl;
}

bool checkDistTracesForEachNodePair(Dfsm &m) {
	m.calculateDistMatrix();
	for (size_t i = 0; i < m.size() - 1; ++i) {
		for (size_t j = i + 1; j < m.size(); ++j) {
			auto ni = m.getNodes().at(i);
			auto nj = m.getNodes().at(j);
			auto distTraces = m.getDistTraces(*ni, *nj);
			for (auto trc : distTraces) {
				shared_ptr<InputTrace> iTr =
					make_shared<InputTrace>(*trc, m.getPresentationLayer());
				OutputTree oti = ni->apply(*iTr);
				OutputTree otj = nj->apply(*iTr);
				if (oti == otj) {
					return false;
				}
			}
		}
	}
	return true;
}

void testGetDistTraces() {
	cout << "TC-DFSM-0020 Show that calculated distinguishing traces "
		<< "in fact distinguish states"
		<< endl;

	auto pl = make_shared<FsmPresentationLayer>();
	auto m = make_shared<Dfsm>("M", 50, 5, 5, pl);

	fsmlib_assert("TC-DFSM-0020",
		checkDistTracesForEachNodePair(*m),
		"Each calculated distinguishing trace produces unequal set of output traces");

}

void testHMethod() {
	cout << "TC-DFSM-0021 Show that Dfsm implModel only passes H-Method Testsuite "
		<< "if intersection is completely defined"
		<< endl;

	auto pl = make_shared<FsmPresentationLayer>();
	auto refModel = make_shared<Dfsm>("refModel", 50, 5, 5, pl)->minimise();
	Fsm implModel = refModel.createMutant("mutant", 2, 2)->minimise();

	// refModel and implModel have to be compl. specified, deterministic and minimal
	// implModel should have at most the same size as refModel
	IOListContainer iolc = refModel.hMethodOnMinimisedDfsm(0);
	TestSuite ts1 = refModel.createTestSuite(iolc);
	TestSuite ts2 = implModel.createTestSuite(iolc);

	fsmlib_assert("TC-DFSM-0021",
		refModel.intersect(implModel).isCompletelyDefined() == ts1.isEquivalentTo(ts2),
		"implModel passes H-Method Testsuite if and only if intersection is completely defined");
}

void testWpMethodWithDfsm() {
	cout << "TC-DFSM-0022 Show that Dfsm implModel only passes Wp-Method Testsuite "
		<< "if intersection is completely defined"
		<< endl;

	auto pl = make_shared<FsmPresentationLayer>();
	auto refModel = make_shared<Dfsm>("refModel", 50, 5, 5, pl)->minimise();
	Fsm implModel = refModel.createMutant("mutant", 1, 1)->minimise();

	// refModel required to be minimised and observable
	// implModel required to have at most 0 additional states compared to refModel (both are prime machines)
	IOListContainer iolc = refModel.wpMethodOnMinimisedDfsm(0);
	TestSuite ts1 = refModel.createTestSuite(iolc);
	TestSuite ts2 = implModel.createTestSuite(iolc);

	// refModel and implModel required to be deterministic and completely specified
	fsmlib_assert("TC-DFSM-0022",
		refModel.intersect(implModel).isCompletelyDefined() == ts1.isEquivalentTo(ts2),
		"implModel passes Wp-Method Testsuite if and only if intersection is completely defined");
}


//===================================== TreeNode Tests ===================================================

// tests TreeNode::add(const int x). Checks if new TreeEdge is created for given input.
void testTreeNodeAddConstInt1(){
	int io = 1;
	shared_ptr<TreeNode> n1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> ref = n1->add(io);
	fsmlib_assert("TC-TreeNode-NNNN",
		static_cast<shared_ptr<TreeNode>>(ref->getParent()) == n1,
		"parent of new node is old node");

	bool containedInChildren = false;
	for (shared_ptr<TreeEdge> e : *(n1->getChildren())) {
		if (e->getIO() == io && e->getTarget() == ref) {
			containedInChildren = true;
		}
	}
	fsmlib_assert("TC-TreeNode-NNNN",
		containedInChildren,
		"after call to TreeNode::add(x) there has to be a child labeled with x");
}

// tests TreeNode::add(const int x). Checks if no new TreeEdge is created if TreeEdge with matching io label
// already exists.
void testTreeNodeAddConstInt2() {
	int io = 1;
	shared_ptr<TreeNode> n1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = n1->add(io);
	int oldNumChilds = n1->getChildren()->size();
	shared_ptr<TreeNode> child2 = n1->add(io);
	int newNumChilds = n1->getChildren()->size();
	fsmlib_assert("TC-TreeNode-NNNN",
		child2 == child1,
		"TreeNode::add(x) returns reference to target node of existing TreeEdge with matching io");
	fsmlib_assert("TC-TreeNode-NNNN",
		oldNumChilds == newNumChilds,
		"TreeNode::add(x) doesn't add new TreeEdge if TreeEdge with matching io already exists");
}

// tests TreeNode::add(const int x). Checks if new TreeEdge is created for given input. TreeNode already has children, but none with matching
// io label.
void testTreeNodeAddConstInt3() {
	shared_ptr<TreeNode> n1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = n1->add(1);
	shared_ptr<TreeNode> child2 = n1->add(2);
	fsmlib_assert("TC-TreeNode-NNNN",
		child1 != child2,
		"calling TreeNode::add(x) and TreeNode::add(y) with x != y returns two different nodes");


	fsmlib_assert("TC-TreeNode-NNNN",
		n1->getChildren()->size() == 2,
		"number of TreeEdges contained in children attribute matches number of actually added values");
}

// tests addToThisNode(const std::vector<int> &lst) ( add(std::vector<int>::const_iterator lstIte, const std::vector<int>::const_iterator end) respectivly)
void testAddToThisNode() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	//shared_ptr<TreeNode> copy = root->clone();
	std::vector<int> inputs = {};
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		root->isLeaf(),
		"addToThisNode() doesn't change tree if vector is empty");

	// root is leaf. input vector contains only one element
	inputs = { 1 };
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		root->getChildren()->size() == 1
		&& root->getChildren()->at(0)->getIO() == 1
		&& root->getChildren()->at(0)->getTarget()->isLeaf(),
		"addToThisNode({1}) called on root (leaf) adds only one edge labeled by correct input. Target node is leaf");

	// root is no leaf (one edge labeled with 1; child is leaf). input vector contains two elements and shares no prefix with 
	// path starting at root.
	inputs = { 2,3 };
	shared_ptr<TreeEdge> e1 = make_shared<TreeEdge>(1, make_shared<TreeNode>());
	shared_ptr<TreeEdge> e2 = make_shared<TreeEdge>(2, make_shared<TreeNode>());
	shared_ptr<TreeEdge> e3 = make_shared<TreeEdge>(3, make_shared<TreeNode>());
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		root->getChildren()->size() == 2
		&& root->hasEdge(e1) != nullptr
		&& root->hasEdge(e2) != nullptr
		&& root->hasEdge(e2)->getTarget()->getChildren()->size() == 1
		&& root->hasEdge(e2)->getTarget()->hasEdge(e3) != nullptr
		&& root->hasEdge(e2)->getTarget()->hasEdge(e3)->getTarget()->isLeaf(),
		"addToThisNode({2,3}) called on root (non leaf) adds only new edge to root labeled by 2. Target node gets edge with label 3 targeting leaf node");

	// root has one edge labeled with 1. child is leaf. vector contains two elements and shares prefix with path starting at root
	root = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, make_shared<TreeNode>()));
	inputs = { 1,2 };
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		root->getChildren()->size() == 1
		&& root->getChildren()->at(0)->getIO() == 1
		&& root->getChildren()->at(0)->getTarget()->getChildren()->size() == 1
		&& root->getChildren()->at(0)->getTarget()->getChildren()->at(0)->getIO() == 2
		&& root->getChildren()->at(0)->getTarget()->getChildren()->at(0)->getTarget()->isLeaf(),
		"addToThisNode() called with vector sharing prefix with path starting at root doesn't add new edges for this prefix. "
		"result contains path represented by vector.");

	// root has two edges (1-edge, 2-edge). target of 1-edge has one edge (2-edge). vector equals path contained in this tree. 
	root = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	shared_ptr<TreeNode> grandChild = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	root->add(make_shared<TreeEdge>(2, child2));
	child1->add(make_shared<TreeEdge>(2, grandChild));
	shared_ptr<TreeNode> copy = root->clone();
	inputs = { 1,2 };
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *copy),
		"addToThisNode() doesn't change tree if vector equals path contained in tree emanating from root");

	// root is leaf. vector contains two elements
	root = make_shared<TreeNode>();
	inputs = { 1,2 };
	root->addToThisNode(inputs);
	fsmlib_assert("TC-TreeNode-NNNN",
		root->getChildren()->size() == 1
		&& root->getChildren()->at(0)->getIO() == inputs.at(0)
		&& root->getChildren()->at(0)->getTarget()->getChildren()->size() == 1
		&& root->getChildren()->at(0)->getTarget()->getChildren()->at(0)->getIO() == inputs.at(1)
		&& root->getChildren()->at(0)->getTarget()->getChildren()->at(0)->getTarget()->isLeaf(),
		"addToThisNode() called on leaf adds path represented by vector.");
}

// helper prints vector<int>
void printVector(vector<int> & vec) {
	std::cout << "[";
	for (int i : vec) {
		std::cout << i << " ";
	}
	std::cout << "]\n";
}

// helper to print vector<vector<int>>
void printVectors(shared_ptr<vector<vector<int>>> vectors) {
	for (vector<int> &vec : *vectors) {
		printVector(vec);
	}
	std::cout << "\n";
}

// checks if every path contained in cloneIoll extended by every path in ioLstPtr is contained in rootIoll.
bool containsExpectedPaths(shared_ptr<vector<vector<int>>> cloneIoll, shared_ptr<vector<vector<int>>> rootIoll, shared_ptr<std::vector<std::vector<int>>> iolLstPtr) {
	printVectors(cloneIoll);
	printVectors(rootIoll);
	printVectors(iolLstPtr);
	for (vector<int> &clonePath : *cloneIoll) {
		for (vector<int> &extendingPath : *iolLstPtr) {
			vector<int> expectedPath(clonePath);
			std::cout << " expected vorher:";
			printVector(expectedPath);
			expectedPath.insert(expectedPath.cend(), extendingPath.cbegin(), extendingPath.cend());
			std::cout << " expected nachher:";
			printVector(expectedPath);
			if (find(rootIoll->cbegin(), rootIoll->cend(), expectedPath) == rootIoll->cend()) {
				return false;
			}
		}
	}
	return true;
}

// checks if every path of leaves is a result of a concatenation of a path from cloneIoll and a path from iolLstPtr.
bool containsNoUnexpectedPath(shared_ptr<vector<vector<int>>> cloneIoll, std::vector<std::shared_ptr<TreeNode>>& leaves, shared_ptr<std::vector<std::vector<int>>> iolLstPtr) {
	std::cout << "============================================================================" << std::endl;
	printVectors(cloneIoll);
	printVectors(iolLstPtr);

	// first build all possible concatenations and save them
	vector<vector<int>> concatenations;
	for (vector<int> &clonePath : *cloneIoll) {
		for (vector<int> &extendingPath : *iolLstPtr) {
			std::vector<int> concatenation(clonePath);
			concatenation.insert(concatenation.cend(), extendingPath.cbegin(), extendingPath.cend());
			concatenations.push_back(concatenation);
		}
	}

	std::cout << "concatenations: " << std::endl;
	printVectors(make_shared<vector<vector<int>>>(concatenations));

	// now check if every path from leaves is contained in concatenations
	for (shared_ptr<TreeNode> leave : leaves) {
		std::vector<int> path = leave->getPath();
		std::cout << "path: ";
		printVector(path);
		if (find(concatenations.cbegin(), concatenations.cend(), path) == concatenations.cend()) {
			return false;
		}
	}

	std::cout << "============================================================================" << std::endl;
	return true;
}

// tests TreeNode::add(const IOListContainer & tcl)
void testTreeNodeAddIOListContainer() {
	// root is leaf. iolc contains only one trace which is empty.
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	shared_ptr<TreeNode> clone = root->clone();
	std::vector<int> ioTrace1 = {};
	std::vector<std::vector<int>> ioLst = { ioTrace1 };
	shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();	
	root->add(IOListContainer(iolLstPtr, presentationLayer));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& (*root == *clone),
		"add(IOListContainer &iolc) called on leaf doesn't change tree if iolc only contains empty traces");

	// root is leaf. iolc contains only one trace which consists only of one input.
	root = make_shared<TreeNode>();
	clone = root->clone();
	ioTrace1 = { 1 };
	ioLst = { ioTrace1 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	shared_ptr<vector<vector<int>>> cloneIoll = make_shared<vector<vector<int>>>();
	std::vector<int> cloneV;
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	shared_ptr<vector<vector<int>>> rootIoll = make_shared<vector<vector<int>>>();
	std::vector<int> rootV;
	root->traverse(rootV, rootIoll);
	vector<shared_ptr<TreeNode>> leaves;
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr), //containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr)
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");

	// root is leaf. iolc contains two traces. Each trace contains only one element (differing)
	root = make_shared<TreeNode>();
	clone = root->clone();
	ioTrace1 = { 1 };
	std::vector<int> ioTrace2;
	ioTrace2 = { 2 };
	ioLst = { ioTrace1, ioTrace2 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	cloneIoll = make_shared<vector<vector<int>>>();
	cloneV.clear();
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	rootIoll = make_shared<vector<vector<int>>>();
	rootV.clear();
	root->traverse(rootV, rootIoll);
	leaves.clear();
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");


	// root is leaf. iolc contains two paths. first path consists of two elements. second path consists of one element.
	root = make_shared<TreeNode>();
	clone = root->clone();
	ioTrace1 = { 1,2 };
	ioTrace2 = { 2 };
	ioLst = { ioTrace1, ioTrace2 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	cloneIoll = make_shared<vector<vector<int>>>();
	cloneV.clear();
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	rootIoll = make_shared<vector<vector<int>>>();
	rootV.clear();
	root->traverse(rootV, rootIoll);
	leaves.clear();
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");

	//// root is leaf. iolc contains two paths. first path consists of two elements. second path consists of one element.
	//root = make_shared<TreeNode>();
	//clone = root->clone();
	//ioTrace1 = { 1,2 };
	//ioTrace2 = { 1 };
	//ioLst = { ioTrace1, ioTrace2 };
	//iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	//presentationLayer = make_shared<FsmPresentationLayer>();
	//
	//cloneIoll = make_shared<vector<vector<int>>>();
	//cloneV.clear();
	//clone->traverse(cloneV, cloneIoll);
	//
	//root->add(IOListContainer(iolLstPtr, presentationLayer));
	//rootIoll = make_shared<vector<vector<int>>>();
	//rootV.clear();
	//root->traverse(rootV, rootIoll);
	//leaves.clear();
	//root->calcLeaves(leaves);
    //
	//assert("TC-TreeNode-NNNN",
	//	root->superTreeOf(clone)
	//	&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
	//	"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");

	// root has two childs (leaves). iolc contains one paths, which doesn't share a prefix with a path already contained in tree.
	root = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	root->add(make_shared<TreeEdge>(2, child2));
	clone = root->clone();
	ioTrace1 = { 3 };
	ioLst = { ioTrace1 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	cloneIoll = make_shared<vector<vector<int>>>();
	cloneV.clear();
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	rootIoll = make_shared<vector<vector<int>>>();
	rootV.clear();
	root->traverse(rootV, rootIoll);
	leaves.clear();
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");

	// root has two childs (leaves). iolc contains one path, which is equals to a path already contained in tree.
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	root->add(make_shared<TreeEdge>(2, child2));
	clone = root->clone();
	ioTrace1 = { 1 };
	ioLst = { ioTrace1 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	cloneIoll = make_shared<vector<vector<int>>>();
	cloneV.clear();
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	rootIoll = make_shared<vector<vector<int>>>();
	rootV.clear();
	root->traverse(rootV, rootIoll);
	leaves.clear();
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");

	// root has two childs (leaves). iolc contains one path, which has a prefix already contained in tree.
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	root->add(make_shared<TreeEdge>(2, child2));
	clone = root->clone();
	ioTrace1 = { 1, 2 };
	ioLst = { ioTrace1 };
	iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
	presentationLayer = make_shared<FsmPresentationLayer>();

	cloneIoll = make_shared<vector<vector<int>>>();
	cloneV.clear();
	clone->traverse(cloneV, cloneIoll);

	root->add(IOListContainer(iolLstPtr, presentationLayer));
	rootIoll = make_shared<vector<vector<int>>>();
	rootV.clear();
	root->traverse(rootV, rootIoll);
	leaves.clear();
	root->calcLeaves(leaves);

	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(clone)
		&& containsExpectedPaths(cloneIoll, rootIoll, iolLstPtr),
		"add(IOListContainer &iolc) result is super tree of old tree and result contains the expected paths");
	fsmlib_assert("TC-TreeNode-NNNN",
		containsNoUnexpectedPath(cloneIoll, leaves, iolLstPtr),
		"add(IOListContainer &iolc) result contains only expected paths");
}

// tests operator==(TreeNode const & treeNode1, TreeNode const & treeNode2)  (positive case)
void testTreeNodeEqualOperator1() {
	shared_ptr<TreeNode> n1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n2 = make_shared<TreeNode>();
	fsmlib_assert("TC-TreeNode-NNNN",
		*n1 == *n2,
		"operator== returns true if both nodes are equal");

	shared_ptr<TreeNode> n11 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n21 = make_shared<TreeNode>();
	n1->add(make_shared<TreeEdge>(1, n11));
	n2->add(make_shared<TreeEdge>(1, n21));
	fsmlib_assert("TC-TreeNode-NNNN",
		*n1 == *n2,
		"operator== returns true if both nodes are equal");

	shared_ptr<TreeNode> n12 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n22 = make_shared<TreeNode>();
	n1->add(make_shared<TreeEdge>(2, n12));
	n2->add(make_shared<TreeEdge>(2, n22));
	fsmlib_assert("TC-TreeNode-NNNN",
		*n1 == *n2,
		"operator== returns true if both nodes are equal");

	shared_ptr<TreeNode> n111 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n112 = make_shared<TreeNode>();
	n11->add(make_shared<TreeEdge>(1, n111));
	n11->add(make_shared<TreeEdge>(2, n112));
	shared_ptr<TreeNode> n211 = make_shared<TreeNode>();	
	shared_ptr<TreeNode> n212 = make_shared<TreeNode>();
	n21->add(make_shared<TreeEdge>(1, n211));
	n21->add(make_shared<TreeEdge>(2, n212));

	fsmlib_assert("TC-TreeNode-NNNN",
		*n1 == *n2,
		"operator== returns true if both nodes are equal");
}

// tests operator==(TreeNode const & treeNode1, TreeNode const & treeNode2)  (negative case)
void testTreeNodeEqualOperator2() {
	shared_ptr<TreeNode> n1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n2 = make_shared<TreeNode>();
	n1->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2),
		"operator== returns false if only one of the TreeNode instances is marked as deleted");

	n1 = make_shared<TreeNode>();
	n2 = make_shared<TreeNode>();
	n2->add(make_shared<TreeEdge>(1, make_shared<TreeNode>()));
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2),
		"operator== returns false if the compared TreeNode instances have different number of children");

	n1->add(make_shared<TreeEdge>(2, make_shared<TreeNode>()));
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2) && n1->getChildren()->size() == n2->getChildren()->size(),
		"operator== returns false if both TreeNode instances have same number of children but edges are labeled differently");

	n1 = make_shared<TreeNode>();
	n2 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n11 = make_shared<TreeNode>();
	shared_ptr<TreeNode> n21 = make_shared<TreeNode>();
	n1->add(make_shared<TreeEdge>(1, n11));
	n2->add(make_shared<TreeEdge>(1, n21));
	n11->add(make_shared<TreeEdge>(1, make_shared<TreeNode>()));
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2) && n11->getChildren()->size() != n21->getChildren()->size(),
		"operator== returns false if two corresponding childs of both TreeNode instances differ in the number of children");

	n21->add(make_shared<TreeEdge>(2, make_shared<TreeNode>()));
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2) && n11->getChildren()->size() == n21->getChildren()->size(),
		"operator== returns false if two corresponding childs of both TreeNode instances differ in the labeling of their children");

	n11->add(make_shared<TreeEdge>(2, make_shared<TreeNode>()));
	n21->add(make_shared<TreeEdge>(1, make_shared<TreeNode>()));
	n11->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(*n1 == *n2),
		"operator== returns false if two corresponding childs differ in beeing marked as deleted");

	n21->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*n1 == *n2),
		"operator== returns true if both instances are equal");
}

// tests TreeNode::calcLeaves(vector<shared_ptr<TreeNode const>>& leaves) const
void testTreeNodeCalcLeaves() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	std::vector<shared_ptr<TreeNode>> leaves;
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 1 && leaves[0] == root,
		"calcLeaves() called on leave adds this leave");

	leaves = std::vector<shared_ptr<TreeNode>>();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 1 && leaves[0] == child1,
		"calcLeaves() called on parent with leave-child adds this leave-child");

	leaves = std::vector<shared_ptr<TreeNode>>();
	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 2 && std::find(leaves.cbegin(), leaves.cend(), child1) != leaves.cend() 
		&& std::find(leaves.cbegin(), leaves.cend(), child2) != leaves.cend(),
		"calcLeaves() called on parent with leave-childs adds all leave-childs");

	leaves = std::vector<shared_ptr<TreeNode>>();
	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(1, grandChild1));
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 2 && std::find(leaves.cbegin(), leaves.cend(), grandChild1) != leaves.cend()
		&& std::find(leaves.cbegin(), leaves.cend(), child2) != leaves.cend(),
		"calcLeaves() called on parent with leave-child and leave-grandchild adds leave-child and leave-grandchild");

	leaves = std::vector<shared_ptr<TreeNode>>();
	shared_ptr<TreeNode> grandChild2 = make_shared<TreeNode>();
	child2->add(make_shared<TreeEdge>(1, grandChild2));
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 2 && std::find(leaves.cbegin(), leaves.cend(), grandChild1) != leaves.cend()
		&& std::find(leaves.cbegin(), leaves.cend(), grandChild2) != leaves.cend(),
		"calcLeaves() called on root with two leave-grandchilds adds both leave-grandchilds");

	leaves = std::vector<shared_ptr<TreeNode>>();
	shared_ptr<TreeNode> grandChild3 = make_shared<TreeNode>();
	child2->add(make_shared<TreeEdge>(3, grandChild3));
	root->calcLeaves(leaves);
	fsmlib_assert("TC-TreeNode-NNNN",
		leaves.size() == 3 && std::find(leaves.cbegin(), leaves.cend(), grandChild1) != leaves.cend()
		&& std::find(leaves.cbegin(), leaves.cend(), grandChild2) != leaves.cend()
		&& std::find(leaves.cbegin(), leaves.cend(), grandChild3) != leaves.cend(),
		"calcLeaves() called on root with three leave-grandchilds adds these leave-grandchilds");
}

// extracts all TreeNodes reachable from given node.
void extractAllTreeNodes(shared_ptr<TreeNode> node, std::vector<shared_ptr<TreeNode>> & nodes) {
	nodes.push_back(node);
	for (const auto &edge : *(node->getChildren())) {
		extractAllTreeNodes(edge->getTarget(), nodes);
	}
}

// extracts all TreeEdges reachable from given node.
void extractAllTreeEdges(shared_ptr<TreeNode> node, std::vector<shared_ptr<TreeEdge>> & edges) {
	for (const auto &edge : *(node->getChildren())) {
		edges.push_back(edge);
		extractAllTreeEdges(edge->getTarget(), edges);
	}
}

// only called with TreeNodes which are known to be equal (*original == *copy).
// checks if no TreeNode and no TreeEdge is contained in both trees.
bool isDeepCopyOfEqualNode(shared_ptr<TreeNode> original, shared_ptr<TreeNode> copy) {
	std::vector<shared_ptr<TreeNode>> originalNodes;
	extractAllTreeNodes(original, originalNodes);
	std::vector<shared_ptr<TreeNode>> copyNodes;
	std::cout << "originalNodes.size(): " << originalNodes.size() << std::endl;
	
	extractAllTreeNodes(copy, copyNodes);
	std::cout << "copyNodes.size(): " << copyNodes.size() << std::endl;
	for (shared_ptr<TreeNode> node : originalNodes) {
		for (shared_ptr<TreeNode> cnode : copyNodes) {
			if (node == cnode) {
				return false;
			}
		}
	}

	//std::vector<shared_ptr<TreeEdge>> originalEdges;
	//extractAllTreeEdges(original, originalEdges);
	//std::vector<shared_ptr<TreeEdge>> copyEdges;
	//extractAllTreeEdges(copy, copyEdges);
	//for (shared_ptr<TreeEdge> edge : originalEdges) {
	//	for (shared_ptr<TreeEdge> cedge : copyEdges) {
	//		if (edge == cedge) {
	//			return false;
	//		}
	//	}
	//}
}

// tests TreeNode::clone() const
void testTreeNodeClone() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	shared_ptr<TreeNode> clone = root->clone();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *clone) && isDeepCopyOfEqualNode(root, clone),    //(root != clone),
		"clone equals original and clone is deep copy");

	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	clone = root->clone();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *clone) && isDeepCopyOfEqualNode(root, clone),
		"clone equals original and clone is deep copy");

	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	clone = root->clone();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *clone) && isDeepCopyOfEqualNode(root, clone),
		"clone equals original and clone is deep copy");

	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(1, grandChild1));
	clone = root->clone();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *clone) && isDeepCopyOfEqualNode(root, clone),
		"clone equals original and clone is deep copy");

	shared_ptr<TreeNode> grandChild2 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(2, grandChild2));
	clone = root->clone();
	fsmlib_assert("TC-TreeNode-NNNN",
		(*root == *clone) && isDeepCopyOfEqualNode(root, clone),
		"clone equals original and clone is deep copy");
}

// tests TreeNode::getPath() const
void testTreeNodeGetPath() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	std::vector<int> expected = {};
	fsmlib_assert("TC-TreeNode-NNNN",
		expected == root->getPath(),
		"getPath invoked on root returns empty list");

	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	expected = {1};
	fsmlib_assert("TC-TreeNode-NNNN",
		expected == child1->getPath(),
		"getPath invoked on child returns list containing only the input needed to reach it");

	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	expected = { 2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		expected == child2->getPath(),
		"getPath invoked on child returns list containing only the input needed to reach it");

	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> grandChild2 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(1, grandChild1));
	child1->add(make_shared<TreeEdge>(2, grandChild2));
	expected = { 1, 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		expected == grandChild1->getPath(),
		"getPath invoked on grandchild returns list containing only the two inputs needed to reach it");

	expected = { 1, 2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		expected == grandChild2->getPath(),
		"getPath invoked on grandchild returns list containing only the two inputs needed to reach it");
}

// tests TreeNode::superTreeOf(const shared_ptr<TreeNode> otherNode) const.
// negative case
void testTreeNodeSuperTreeOf1() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	shared_ptr<TreeNode> rootOther = make_shared<TreeNode>();
	shared_ptr<TreeNode> childOther1 = make_shared<TreeNode>();
	rootOther->add(make_shared<TreeEdge>(1, childOther1));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has more children than root");

	shared_ptr<TreeNode> childOther2 = make_shared<TreeNode>();
	rootOther->add(make_shared<TreeEdge>(2, childOther2));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has more children than root");

	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has more children than root");

	shared_ptr<TreeNode> grandChildOther1 = make_shared<TreeNode>();
	childOther2->add(make_shared<TreeEdge>(1, grandChildOther1));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has more children than root");

	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if one corresponding child of root and rootOther has different number of childs");

	//-----------------------------------------------------------------------------------------------

	root = make_shared<TreeNode>();
	rootOther = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	childOther1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	rootOther->add(make_shared<TreeEdge>(2, childOther1));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has TreeEdge with label not existent in root");

	child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	childOther2 = make_shared<TreeNode>();
	rootOther->add(make_shared<TreeEdge>(3, childOther2));
	fsmlib_assert("TC-TreeNode-NNNN",
		!root->superTreeOf(rootOther),
		"superTreeOf() returns false if rootOther has TreeEdge with label not existent in root");
}

// tests TreeNode::superTreeOf(const shared_ptr<TreeNode> otherNode) const.
// positive case
void testTreeNodeSuperTreeOf2() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	shared_ptr<TreeNode> rootOther = make_shared<TreeNode>();
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther),
		"superTreeOf() returns true if root and rootOther are equal");

	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther),
		"superTreeOf() returns true if root contains rootOther");

	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	shared_ptr<TreeNode> childOther1 = make_shared<TreeNode>();
	rootOther->add(make_shared<TreeEdge>(1, childOther1));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther),
		"superTreeOf() returns true if root contains rootOther");

	shared_ptr<TreeNode> childOther2 = make_shared<TreeNode>();
	rootOther->add(make_shared<TreeEdge>(2, childOther2));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther) && *root == *rootOther,
		"superTreeOf() returns true if root and rootOther are equal");

	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(1, grandChild1));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther),
		"superTreeOf() returns true if root contains rootOther");

	shared_ptr<TreeNode> grandChild2 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(2, grandChild2));
	shared_ptr<TreeNode> grandChildOther1 = make_shared<TreeNode>();
	childOther1->add(make_shared<TreeEdge>(2, grandChildOther1));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->superTreeOf(rootOther),
		"superTreeOf() returns true if root contains rootOther");
}

// tests TreeNode::traverse(vector<int>& v,shared_ptr<vector<vector<int>>> ioll)
void testTreeNodeTraverse() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	std::vector<int> v;
	std::shared_ptr<std::vector<std::vector<int>>> ioll = make_shared<std::vector<std::vector<int>>>();
	root->traverse(v, ioll);
	fsmlib_assert("TC-TreeNode-NNNN",
		v.size() == 0 && ioll->size() == 1 && ioll->at(0) == v,
		"traverse() called on leave doesn't add inputs to current input vector v but adds v to ioll");

	v.clear();
	ioll->clear();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, child1));
	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(2, child2));
	root->traverse(v, ioll);
	/*const std::vector<int> e1 = {};
	const std::vector<int> e2 = { 1 };
	const std::vector<int> e3 = { 2 };*/
	fsmlib_assert("TC-TreeNode-NNNN",
		ioll->size() == 3 && std::find(ioll->cbegin(), ioll->cend(), root->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), child1->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), child2->getPath()) != ioll->cend(),
		"traverse() called on node n adds all int paths from n to leaves of the tree");

	v.clear();
	ioll->clear();
	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(3, grandChild1));
	shared_ptr<TreeNode> grandChild2 = make_shared<TreeNode>();
	child1->add(make_shared<TreeEdge>(4, grandChild2));
	root->traverse(v, ioll);
	//const std::vector<int> e4 = { 1,3 };
	//const std::vector<int> e5 = { 1,4 };
	fsmlib_assert("TC-TreeNode-NNNN",
		ioll->size() == 5 && std::find(ioll->cbegin(), ioll->cend(), root->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), child1->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), child2->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), grandChild1->getPath()) != ioll->cend()
		&& std::find(ioll->cbegin(), ioll->cend(), grandChild2->getPath()) != ioll->cend(),
		"traverse() called on node n adds all int paths from n to leaves of the tree");
}

// tests TreeNode::deleteNode()
void testTreeNodeDeleteNode() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	root->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		root->isDeleted(),
		"deleteNode() called on root marks root as deleted");

	// parent (root) with one child (leaf). child gets deleted
	root = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	shared_ptr<TreeEdge> rootToChild1 = make_shared<TreeEdge>(1, child1);
	root->add(rootToChild1);
	child1->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) == root->getChildren()->cend(),
		"deleteNode() called on child marks child as deleted, removes child from children list of parent and doesn't mark parent as deleted");

	// parent (root) with one child (leaf). parent gets deleted
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	root->add(rootToChild1);
	root->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		root->isDeleted()
		&& !(child1->isDeleted())
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend(),
		"deleteNode() called on parent of undeleted child marks parent as deleted but doesn't change children");

	// parent (root) with two children (leaves) child1 and child2. child2 gets deleted
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> child2 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	shared_ptr<TreeEdge> rootToChild2 = make_shared<TreeEdge>(2, child2);
	root->add(rootToChild1);
	root->add(rootToChild2);
	child2->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& !(child1->isDeleted())
		&& child2->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild2) == root->getChildren()->cend(),
		"deleteNode() called on child2 (leaf) marks child2 as deleted, removes it from child list of its parent and doesn't change the other child");

	// root has child1 as only child. child1 has grandChild1 as only child, which is a leaf. child1 gets deleted.
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	shared_ptr<TreeEdge> child1ToGrandChild1 = make_shared<TreeEdge>(1, grandChild1);
	root->add(rootToChild1);
	child1->add(child1ToGrandChild1);
	child1->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& !(grandChild1->isDeleted())
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend()
		&& std::find(child1->getChildren()->cbegin(), child1->getChildren()->cend(), child1ToGrandChild1) != child1->getChildren()->cend(),
		"deleteNode() called on child (non leaf) marks it as deleted but doesn't change parent or grandchilds of parent");
	
	// root has child1 as only child. child1 has grandChild1 as only child, which is a leaf. child1 is already deleted.
	// In the next step grandChild1 gets deleted
	grandChild1->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& grandChild1->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) == root->getChildren()->cend()
		&& std::find(child1->getChildren()->cbegin(), child1->getChildren()->cend(), child1ToGrandChild1) == child1->getChildren()->cend(),
		"deleteNode() called on child (leaf) with already deleted parent marks this child as deleted, removes it from "
		"list of childs from its parent and removes parent from child list of parents parent");

	// root has child1 as only child. child1 has grandChild1 as only child, which is a leaf. grandChild1 gets deleted.
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	grandChild1 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	child1ToGrandChild1 = make_shared<TreeEdge>(1, grandChild1);
	root->add(rootToChild1);
	child1->add(child1ToGrandChild1);
	grandChild1->deleteNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& !(child1->isDeleted())
		&& grandChild1->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend()
		&& std::find(child1->getChildren()->cbegin(), child1->getChildren()->cend(), child1ToGrandChild1) == child1->getChildren()->cend(),
		"deleteNode() called on child (leaf) with already non deleted parent marks this child as deleted, removes it from "
		"list of childs from its parent and doesn't delete parent");
}

// tests TreeNode::deleteSingleNode()
void testTreeNodeDeleteSingleNode() {
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	root->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		root->isDeleted(),
		"deleteSingleNode() called on root marks root as deleted");

	// parent (root) with one child (leaf). child gets deleted
	root = make_shared<TreeNode>();
	shared_ptr<TreeNode> child1 = make_shared<TreeNode>();
	shared_ptr<TreeEdge> rootToChild1 = make_shared<TreeEdge>(1, child1);
	root->add(rootToChild1);
	child1->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) == root->getChildren()->cend(),
		"deleteSingleNode() called on child marks child as deleted, removes child from children list of parent and doesn't mark parent as deleted");

	// parent (root) with one child (leaf). parent gets deleted
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	root->add(rootToChild1);
	root->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		root->isDeleted()
		&& !(child1->isDeleted())
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend(),
		"deleteSingleNode() called on parent of undeleted child marks parent as deleted but doesn't change children");

	// root has child1 as only child. child1 has grandChild1 as only child, which is a leaf. child1 gets deleted.
	root = make_shared<TreeNode>();
	child1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> grandChild1 = make_shared<TreeNode>();
	rootToChild1 = make_shared<TreeEdge>(1, child1);
	shared_ptr<TreeEdge> child1ToGrandChild1 = make_shared<TreeEdge>(1, grandChild1);
	root->add(rootToChild1);
	child1->add(child1ToGrandChild1);
	child1->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& !(grandChild1->isDeleted())
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend()
		&& std::find(child1->getChildren()->cbegin(), child1->getChildren()->cend(), child1ToGrandChild1) != child1->getChildren()->cend(),
		"deleteSingleNode() called on child (non leaf) marks it as deleted but doesn't change parent or grandchilds of parent");

	// root has child1 as only child. child1 has grandChild1 as only child, which is a leaf. child1 is already deleted.
	// In the next step grandChild1 gets deleted
	grandChild1->deleteSingleNode();
	fsmlib_assert("TC-TreeNode-NNNN",
		!(root->isDeleted())
		&& child1->isDeleted()
		&& grandChild1->isDeleted()
		&& std::find(root->getChildren()->cbegin(), root->getChildren()->cend(), rootToChild1) != root->getChildren()->cend()
		&& std::find(child1->getChildren()->cbegin(), child1->getChildren()->cend(), child1ToGrandChild1) == child1->getChildren()->cend(),
		"deleteNode() called on child (leaf) with already deleted parent marks this child as deleted, removes it from "
		"list of childs of its parent but doesn't remove parent from child list of parents parent");
}

// tests TreeNode::tentativeAddToThisNode(std::vector<int>::const_iterator start, std::vector<int>::const_iterator stop, std::shared_ptr<TreeNode>& n)
void testTreeNodeTentativeAddToThisNode() {
	// case 1: path is already contained in tree -> returns 0

	// root is leaf and inpPath is empty
	shared_ptr<TreeNode> ref = make_shared<TreeNode>();
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	vector<int> inpPath = {};
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 0
		&& ref == root,
		"tentativeAddToThisNode() returns 0 if path is already contained in tree");

	// root has two children (leaves). inpPath contains one element (already contained as a label)
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	inpPath = { 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 0
		&& ref == c1,
		"tentativeAddToThisNode() returns 0 if path is already contained in tree");

	// root has two childs (c1 and c2). c1 is a leaf. c2 has two childs (gc1 and gc2). Both are leaves. inpPath contains one element (already contained as a label)
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	inpPath = { 2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 0
		&& ref == c2,
		"tentativeAddToThisNode() returns 0 if path is already contained in tree");

	// root has two childs (c1 and c2). c1 is a leaf. c2 has two childs (gc1 and gc2). Both are leaves. inpPath (2 elements) is already contained in tree.
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	gc1 = make_shared<TreeNode>();
	gc2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	inpPath = { 2, 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 0
		&& ref == gc1,
		"tentativeAddToThisNode() returns 0 if path is already contained in tree");

	// case 2: prefix of path reaches leaf of tree (no new branch needed) -> returns 1

	// root is leaf. path contains one element
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	inpPath = { 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 1
		&& ref == root,
		"tentativeAddToThisNode() returns 1 if a prefix of the path reaches a leaf "
		"and adding the path would require to extend the tree");

	// root has two childs (c1 and c2). c1 and c2 are leaves. inpPath contains two elements.
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	//std::cout << "root leaf:" << root->isLeaf() << std::endl;
	//std::cout << "root children size: " << root->getChildren()->size() << std::endl;
	//std::cout << "c1 children size: " << c1->getChildren()->size() << std::endl;
	//std::cout << "c2 children size: " << c2->getChildren()->size() << std::endl;
	inpPath = { 1,2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 1
		&& ref == c1,
		"tentativeAddToThisNode() returns 1 if a prefix of the path reaches a leaf "
		"and adding the path would require to extend the tree");

	// root has two childs (c1 and c2). c1 is a leaf. c2 has two childs (gc1 and gc2). Both are leaves. inpPath (3 elements) 
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	gc1 = make_shared<TreeNode>();
	gc2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	inpPath = { 2, 1, 2 };

	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 1
		&& ref == gc1,
		"tentativeAddToThisNode() returns 1 if a prefix of the path reaches a leaf "
		"and adding the path would require to extend the tree");

	// case 3: path not fully contained and no prefix of path reaches leaf of tree (new branch needed) -> returns 2

	// root is has child c1. c1 is a leaf. path contains one element
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	inpPath = { 2 };

	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 2
		&& ref == root,
		"tentativeAddToThisNode() returns 2 if no prefix of the path reaches a leaf "
		"and adding the path would require to create a new branch");

	// root is has childs c1 and c2. c1 is a leaf. c2 has two childs gc1 and gc2. gc1 and gc2 are leaves.
	// inpPath contains 2 elements.
	ref = make_shared<TreeNode>();
	root = make_shared<TreeNode>();
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	gc1 = make_shared<TreeNode>();
	gc2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	inpPath = { 2, 3 };

	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 2
		&& ref == c2,
		"tentativeAddToThisNode() returns 2 if no prefix of the path reaches a leaf "
		"and adding the path would require to create a new branch");

	// root is has childs c1 and c2. c1 is a leaf. c2 has two childs gc1 and gc2. gc1 and gc2 are leaves.
	// inpPath contains 3 elements.
	inpPath = { 2, 3, 2};
	fsmlib_assert("TC-TreeNode-NNNN",
		root->tentativeAddToThisNode(inpPath.cbegin(), inpPath.cend(), ref) == 2
		&& ref == c2,
		"tentativeAddToThisNode() returns 2 if no prefix of the path reaches a leaf "
		"and adding the path would require to create a new branch");
}

// tests TreeNode::after(std::vector<int>::const_iterator lstIte, const std::vector<int>::const_iterator end)
void testTreeNodeAfter() {
	// case 1: result is no nullptr

	// root is a leaf. path is empty.
	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	vector<int> path = {};
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == root,
		"after() called with empty path returns root");

	// root has two children (c1 and c2). path is empty
	shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == root,
		"after() called with empty path returns root");

	// root has two children (c1 and c2). c1 is a leaf. c2 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// path contains one element and is a prefix of a path which is contained in the tree.
	shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
	shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	path = { 2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == c2,
		"after() called with non empty path which is a prefix of a contained path returns the node, which can be reached with this prefix");

	// root has two children (c1 and c2). c1 is a leaf. c2 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// path contains two elements and equals a path which is contained in the tree.
	path = { 2, 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == gc1,
		"after() called with non empty path which equals a contained path returns the node, which can be reached with this path");


	// case 2: result is nullptr

	// root is leaf. path is not empty
	root = make_shared<TreeNode>();
	path = { 1 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == nullptr,
		"after() called with path that can't be completely matched against tree returns nullptr.");

	// root has two children (c1 and c2). path contains one element which doesn't match any edge label.
	c1 = make_shared<TreeNode>();
	c2 = make_shared<TreeNode>();
	root->add(make_shared<TreeEdge>(1, c1));
	root->add(make_shared<TreeEdge>(2, c2));
	path = { 3 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == nullptr,
		"after() called with path that can't be completely matched against tree returns nullptr.");

	// root has two children (c1 and c2). c1 and c2 are leaves. path contains two elements (path is longer than any contained path). 
	path = { 1,2 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == nullptr,
		"after() called with path that can't be completely matched against tree returns nullptr.");

	// root has two children (c1 and c2). c1 is a leaf. c2 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// path contains two elements.
	gc1 = make_shared<TreeNode>();
	gc2 = make_shared<TreeNode>();
	c2->add(make_shared<TreeEdge>(1, gc1));
	c2->add(make_shared<TreeEdge>(2, gc2));
	path = { 2,3 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == nullptr,
		"after() called with path that can't be completely matched against tree returns nullptr.");

	// root has two children (c1 and c2). c1 is a leaf. c2 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// path contains three elements (path is longer than any contained path).
	path = { 2, 1, 3 };
	fsmlib_assert("TC-TreeNode-NNNN",
		root->after(path.cbegin(), path.cend()) == nullptr,
		"after() called with path that can't be completely matched against tree returns nullptr.");
}

// checks if every path (from root to leaf) in newNode is a path in oldNode or in iolc.
// This is a condition the result of TreeNode::addToThisNode(const IOListContainer & tcl) has to fullfill
bool resultContainsOnlyExpectedPaths(shared_ptr<TreeNode> newNode, shared_ptr<TreeNode> oldNode, IOListContainer & iolc) {
	std::cout << "==========================================================" << std::endl;
	std::cout << "resultContainsOnlyExpectedPaths:" << std::endl;
	std::vector<shared_ptr<TreeNode>> newLeaves;
	newNode->calcLeaves(newLeaves);

	std::vector<shared_ptr<TreeNode>> oldLeaves;
	oldNode->calcLeaves(oldLeaves);

	// store all paths from oldNode (from root to a leaf) and from iolc in vector paths.
	std::vector<std::vector<int>> paths;
	for (shared_ptr<TreeNode> oldLeaf : oldLeaves) {
		paths.push_back(oldLeaf->getPath());
	}
	paths.insert(paths.cend(), iolc.getIOLists()->cbegin(), iolc.getIOLists()->cend());

	printVectors(make_shared<std::vector<std::vector<int>>>(paths));

	// check if every path in newNode is contained in paths
	for (shared_ptr<TreeNode> newLeaf : newLeaves) {
		std::vector<int> p = newLeaf->getPath();
		printVector(p);
		if (find(paths.cbegin(), paths.cend(), newLeaf->getPath()) == paths.cend()) {
			return false;
		}
	}
	std::cout << "==========================================================" << std::endl;
	return true;
}

// checks if newNode contains each path contained in iolc.
// This is a condition the result of TreeNode::addToThisNode(const IOListContainer & tcl) has to fullfill
bool resultContainsEachAddedPath(shared_ptr<TreeNode> newNode, IOListContainer & iolc) {	
	std::cout << "==========================================================" << std::endl;
	std::cout << "resultContainsEachAddedPath:" << std::endl;
	std::vector<shared_ptr<TreeNode>> reachable;
	extractAllTreeNodes(newNode, reachable);

	// store each path from newNode in vector paths
	std::vector<std::vector<int>> paths;
	for (shared_ptr<TreeNode> node : reachable) {
		paths.push_back(node->getPath());
	}
	printVectors(make_shared<std::vector<std::vector<int>>>(paths));

	// check if each path from iolc is contained in paths
	for (std::vector<int> path : *(iolc.getIOLists())) {
		if (find(paths.cbegin(), paths.cend(), path) == paths.cend()) {
			return false;
		}
	}
	std::cout << "==========================================================" << std::endl;
	return true;
}

// checks if each TreeNode in the tree emanating from newNode is observable. (There are no two TreeEdges emanating from the same TreeNode, which
// have the same label)
bool treeIsObservable(shared_ptr<TreeNode> newNode) {	
	std::vector<shared_ptr<TreeNode>> reachable;
	extractAllTreeNodes(newNode, reachable);

	// check the children of each node
	for (shared_ptr<TreeNode> node : reachable) {
		std::unordered_set<int> nodeLabels;
		for (shared_ptr<TreeEdge> edge : *(node->getChildren())) {
			// if label of edge couldn't be inserted, it is already contained, so there is another edge with same label
			if (!nodeLabels.insert(edge->getIO()).second) {
				return false;
			}
		}
	}
	return true;
}

// tests TreeNode::addToThisNode(const IOListContainer & tcl)
void testTreeNodeAddToThisNodeIOListContainer() {
	// root is a leaf. IOListContainer is empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> old = root->clone();
		std::vector<std::vector<int>> ioLst = { };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root is a leaf. IOListContainer contains only empty vectors.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = {};
		std::vector<std::vector<int>> ioLst = {inputs1};
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root is a leaf. IOListContainer contains two paths ({1} and {2}).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = { 1 };
		std::vector<int> inputs2 = { 2 };
		std::vector<std::vector<int>> ioLst = { inputs1, inputs2 };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root has one child (c1). c1 is a leaf. IOListContainer contains two paths ({1} and {2}). TreeEdge matches one contained path.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = { 1 };
		std::vector<int> inputs2 = { 2 };
		std::vector<std::vector<int>> ioLst = { inputs1, inputs2 };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root has two childs (c1 and c2). c2 is a leaf. c1 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// IOListContainer contains two paths ({1, 3} and {1, 3}). (both are equal and match a prefix already contained in tree)
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		root->add(make_shared<TreeEdge>(2, c2));
		c1->add(make_shared<TreeEdge>(1, gc1));
		c1->add(make_shared<TreeEdge>(2, gc2));
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = { 1, 3 };
		std::vector<int> inputs2 = { 1, 3 };
		std::vector<std::vector<int>> ioLst = { inputs1, inputs2 };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root has two childs (c1 and c2). c2 is a leaf. c1 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// IOListContainer contains two paths ({1, 3} and {1, 3, 1}). (first is prefix of second, both have prefix in tree)
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		root->add(make_shared<TreeEdge>(2, c2));
		c1->add(make_shared<TreeEdge>(1, gc1));
		c1->add(make_shared<TreeEdge>(2, gc2));
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = { 1, 3 };
		std::vector<int> inputs2 = { 1, 3, 1 };
		std::vector<std::vector<int>> ioLst = { inputs1, inputs2 };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}

	// root has two childs (c1 and c2). c2 is a leaf. c1 has two children (gc1 and gc2). gc1 and gc2 are leaves.
	// IOListContainer contains three paths ({1, 1} and {1, 2} and {2}). (each path is already contained in tree)
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		root->add(make_shared<TreeEdge>(2, c2));
		c1->add(make_shared<TreeEdge>(1, gc1));
		c1->add(make_shared<TreeEdge>(2, gc2));
		shared_ptr<TreeNode> old = root->clone();
		std::vector<int> inputs1 = { 1, 1 };
		std::vector<int> inputs2 = { 1, 2 };
		std::vector<int> inputs3 = { 2 };
		std::vector<std::vector<int>> ioLst = { inputs1, inputs2, inputs3 };
		shared_ptr<std::vector<std::vector<int>>> iolLstPtr = make_shared < std::vector<std::vector<int>>>(ioLst);
		shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
		IOListContainer iolc1(iolLstPtr, presentationLayer);
		root->addToThisNode(iolc1);
		fsmlib_assert("TC-TreeNode-NNNN",
			root->superTreeOf(old),
			"result of addToThisNode(const IOListContainer & tcl) is a super tree of the original TreeNode, so every path from original is still contained in result");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsOnlyExpectedPaths(root, old, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains only expected paths (paths from original tree or added paths)");
		fsmlib_assert("TC-TreeNode-NNNN",
			resultContainsEachAddedPath(root, iolc1),
			"result of addToThisNode(const IOListContainer & tcl) contains each added path");
		fsmlib_assert("TC-TreeNode-NNNN",
			treeIsObservable(root),
			"result of addToThisNode(const IOListContainer & tcl) contains no redundant prefixes");
	}
}

//===================================== Tree Tests ===================================================

// tests Tree::remove(const std::shared_ptr<Tree> otherTree)
void testTreeRemove() {
	
	// thisTree is a leaf. otherTree is a leaf.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());

		thisTree.remove(make_shared<Tree>(otherTree));

		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted(),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isLeaf(),
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists");
	}

	// thisTree is a leaf. otherTree has root otherRoot. otherRoot has one child (otherC1). otherC1 is a leaf.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());

		thisTree.remove(make_shared<Tree>(otherTree));

		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted(),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isLeaf(),
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists");
	}

	// thisTree has root thisRoot. thisRoot has one child (thisC1). thisC1 is a leaf. otherTree is a leaf.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());

		thisTree.remove(make_shared<Tree>(otherTree));
		
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted()
			&& !(thisC1->isDeleted()),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted. Non corresponding nodes aren't deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			!(thisRoot->isLeaf())
			&& thisC1->isLeaf(),
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists. Each non "
			"deleted leaf isn't removed.");		
	}
	
	// thisTree is super tree of otherTree. thisTree contains two edges emanating from thisRoot. otherTree contains one edge.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());
		
		thisTree.remove(make_shared<Tree>(otherTree));
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted()
			&& thisC1->isDeleted()
			&& !(thisC2->isDeleted()),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted. Non corresponding nodes aren't deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->getChildren()->size() == 1
			&& thisRoot->getChildren()->at(0)->getTarget() == thisC2,
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists. Each non "
			"deleted leaf isn't removed.");
	}

	// otherTree is super tree of thisTree. thisTree contains two edges emanating from thisRoot. otherTree contains three edges emanating from otherRoot.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC3 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		otherRoot->add(make_shared<TreeEdge>(2, otherC2));
		otherRoot->add(make_shared<TreeEdge>(3, otherC3));
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());

		thisTree.remove(make_shared<Tree>(otherTree));

		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted()
			&& thisC1->isDeleted()
			&& thisC2->isDeleted(),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted. Non corresponding nodes aren't deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isLeaf(),
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists. Each non "
			"deleted leaf isn't removed.");
	}	

	// thisTree is super tree of otherTree. height of both trees is 2.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC3 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		thisC1->add(make_shared<TreeEdge>(1, thisGC1));
		thisC1->add(make_shared<TreeEdge>(2, thisGC2)); 
		thisC2->add(make_shared<TreeEdge>(1, thisGC3));
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC2 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		otherRoot->add(make_shared<TreeEdge>(2, otherC2));
		otherC1->add(make_shared<TreeEdge>(2, otherGC1));
		otherC2->add(make_shared<TreeEdge>(1, otherGC2));
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());

		thisTree.remove(make_shared<Tree>(otherTree));

		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->isDeleted()
			&& thisC1->isDeleted()
			&& thisC2->isDeleted()
			&& !(thisGC1->isDeleted())
			&& thisGC2->isDeleted()
			&& thisGC3->isDeleted(),
			"remove(const std::shared_ptr<Tree> otherTree): All corresponding nodes (source and target of matching edges) "
			"of thisTree and otherTree are marked as deleted. Non corresponding nodes aren't deleted.");
		fsmlib_assert("TC-Tree-NNNN",
			thisRoot->getChildren()->size() == 1
			&& thisRoot->getChildren()->at(0)->getTarget() == thisC1
			&& thisC1->getChildren()->size() == 1
			&& thisC1->getChildren()->at(0)->getTarget() == thisGC1,
			"remove(const std::shared_ptr<Tree> otherTree): Each deleted leaf is removed from children lists. Each non "
			"deleted leaf isn't removed.");
	}
}

// gets string representation of a Tree::toDot() result and counts the declarations of edges.
int countEdgesInToDotResult(std::string content) {
	int counter = 0;
	int pos = content.find(" -> ", 0);
	while (pos != string::npos) {
		++counter;
		pos = content.find(" -> ", pos + 1);
	}
	return counter;
}

// tests Tree::toDot(ostream & out)
void testTreeToDot() {
	// tree contains only the root (leaf). result contains no edge.
	{
		std::ostringstream stream;
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		Tree tree(root, make_shared<FsmPresentationLayer>());
		tree.toDot(stream);
		std::string content = stream.str();
		fsmlib_assert("TC-Tree-NNNN",
			content.find(" -> ") == string::npos,
			"result of toDot(ostream & out) contains only expected edges");
	}

	
	{
		// root of the tree has one child (c1). c1 is a leaf. result contains 1 edge.
		std::ostringstream stream;
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		Tree tree(root, make_shared<FsmPresentationLayer>());
		tree.toDot(stream);
		std::string content = stream.str();
		fsmlib_assert("TC-Tree-NNNN",
			content.find("0 -> 1[label = \"1\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-Tree-NNNN",
			countEdgesInToDotResult(content) == 1,
			"result of toDot(ostream & out) contains only expected edges");


		// root has two children (c1 and c2). Both are leaves. result contains two edges.
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(2, c2));
		stream.str("");
		stream.clear();
		tree.toDot(stream);
		content = stream.str();

		fsmlib_assert("TC-Tree-NNNN",
			content.find("0 -> 1[label = \"1\" ];") != string::npos
			&& content.find("0 -> 2[label = \"2\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-Tree-NNNN",
			countEdgesInToDotResult(content) == 2,
			"result of toDot(ostream & out) contains only expected edges");

		// root has two children (c1 and c2). c1 has two children (gc1 and gc2). gc1, gc2 and c2 are leaves. result contains 4 edges.
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		c1->add(make_shared<TreeEdge>(1, gc1));
		c1->add(make_shared<TreeEdge>(2, gc2));
		stream.str("");
		stream.clear();
		tree.toDot(stream);
		content = stream.str();

		fsmlib_assert("TC-Tree-NNNN",
			content.find("0 -> 1[label = \"1\" ];") != string::npos
			&& content.find("1 -> 2[label = \"1\" ];") != string::npos
			&& content.find("1 -> 3[label = \"2\" ];") != string::npos
			&& content.find("0 -> 4[label = \"2\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-Tree-NNNN",
			countEdgesInToDotResult(content) == 4,
			"result of toDot(ostream & out) contains only expected edges");
	}

}

// tests Tree::getPrefixRelationTree(const shared_ptr<Tree> & b)
void testTreeGetPrefixRelationTree() {
	// thisTree contains no test case (root without children). otherTree contains two test cases.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		Tree thisTree(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC2 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		otherRoot->add(make_shared<TreeEdge>(2, otherC2));
		Tree otherTree(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> otherTreePtr = make_shared<Tree>(otherTree);
		shared_ptr<Tree> result = thisTree.getPrefixRelationTree(otherTreePtr);
		fsmlib_assert("TC-Tree-NNNN",
			result == otherTreePtr,
			"getPrefixRelationTree(const shared_ptr<Tree> & b) returns pointer to b if thisTree contains no test case (root without childs) "
			"and b contains at least one test case.");		
	}

	// thisTree contains two test cases. otherTree contains no test case (root without children).
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);
		fsmlib_assert("TC-Tree-NNNN",
			result == thisTree,
			"getPrefixRelationTree(const shared_ptr<Tree> & b) returns pointer to thisTree if b contains no test case (root without childs) "
			"and thisTree contains at least one test case.");
	}

	// thisTree and otherTree both contain no test case (root without children).
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);
		fsmlib_assert("TC-Tree-NNNN",
			result != thisTree
			&& result != otherTree
			&& result->size() == 1,
			"getPrefixRelationTree(const shared_ptr<Tree> & b) returns pointer to a new empty Tree if b and thisTree are empty.");
	}

	// thisTree contains two test cases. otherTree contains one testcase, which equals one of the test cases of thisTree.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);
		vector<int> expected = { 1 };
		fsmlib_assert("TC-Tree-NNNN",
			result->getTestCases().getIOLists()->size() == 1
			&& result->getTestCases().getIOLists()->at(0) == expected,
			"result of getPrefixRelationTree(const shared_ptr<Tree> & b) contains each expected test case and no unexpected test case.");
	}

	// thisTree contains two test cases. otherTree contains only one testcase. Trees don't share any prefixes.
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(3, otherC1));
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);
		fsmlib_assert("TC-Tree-NNNN",
			result->size() == 1,
			"result of getPrefixRelationTree(const shared_ptr<Tree> & b) contains each expected test case and no unexpected test case.");
	}

	// thisTree contains two test cases (thisTC1 and thisTC2). otherTree contains three testcase (otherTC1, otherTC2 and otherTC3). 
	// otherTC1 is a prefix of thisTC1. thisTC2 is a prefix of otherTC2 and otherTC3.
	// (Each test case is either a prefix of a test case contained in the other tree or has a prefix in the other tree)
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC1 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		thisC1->add(make_shared<TreeEdge>(1, thisGC1));
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC2 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		otherRoot->add(make_shared<TreeEdge>(2, otherC2));
		otherC2->add(make_shared<TreeEdge>(2, otherGC1));
		otherC2->add(make_shared<TreeEdge>(1, otherGC2));
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);

		vector<int> thisTC1 = { 1,1 };
		vector<int> otherTC2 = { 2,2 };
		vector<int> otherTC3 = { 2,1 };

		shared_ptr<vector<vector<int>>> testCases = result->getTestCases().getIOLists();

		fsmlib_assert("TC-Tree-NNNN",
			testCases->size() == 3
			&& find(testCases->cbegin(), testCases->cend(), thisTC1) != testCases->cend()
			&& find(testCases->cbegin(), testCases->cend(), otherTC2) != testCases->cend()
			&& find(testCases->cbegin(), testCases->cend(), otherTC3) != testCases->cend(),
			"result of getPrefixRelationTree(const shared_ptr<Tree> & b) contains each expected test case and no unexpected test case.");
	}

	// thisTree contains two test cases (thisTC1 and thisTC2). otherTree contains three testcase (otherTC1, otherTC2 and otherTC3). 
	// otherTC1 is a prefix of thisTC1. thisTC2, otherTC2 and otherTC3 aren't prefixes of each other (but share prefixes).
	{
		shared_ptr<TreeNode> thisRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> thisGC2 = make_shared<TreeNode>();
		thisRoot->add(make_shared<TreeEdge>(1, thisC1));
		thisRoot->add(make_shared<TreeEdge>(2, thisC2));
		thisC1->add(make_shared<TreeEdge>(1, thisGC1));
		thisC2->add(make_shared<TreeEdge>(3, thisGC2));
		shared_ptr<Tree> thisTree = make_shared<Tree>(thisRoot, make_shared<FsmPresentationLayer>());

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherC2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> otherGC2 = make_shared<TreeNode>();
		otherRoot->add(make_shared<TreeEdge>(1, otherC1));
		otherRoot->add(make_shared<TreeEdge>(2, otherC2));
		otherC2->add(make_shared<TreeEdge>(2, otherGC1));
		otherC2->add(make_shared<TreeEdge>(1, otherGC2));
		shared_ptr<Tree> otherTree = make_shared<Tree>(otherRoot, make_shared<FsmPresentationLayer>());
		shared_ptr<Tree> result = thisTree->getPrefixRelationTree(otherTree);

		vector<int> thisTC1 = { 1,1 };

		shared_ptr<vector<vector<int>>> testCases = result->getTestCases().getIOLists();

		fsmlib_assert("TC-Tree-NNNN",
			testCases->size() == 1
			&& find(testCases->cbegin(), testCases->cend(), thisTC1) != testCases->cend(),
			"result of getPrefixRelationTree(const shared_ptr<Tree> & b) contains each expected test case and no unexpected test case.");
	}
}

// tests Tree::tentativeAddToRoot(SegmentedTrace& alpha)
void testTreeTentativeAddToRoot() {
	// case 1: alpha is already contained in Tree.

	// root of Tree is a leaf. alpha = <>
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		vector<int> seg1vec = { };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = {seg1};
		SegmentedTrace alpha(segments);
	
		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");
	}

	// Tree contains one path (2 edges). alpha matches this path (with one, two or three segments)
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		c1->add(make_shared<TreeEdge>(2, gc1));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 2 segments
		vector<int> seg1vec = { 1 };				
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 2 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1, seg2 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");

		// alpha = 1 segment
		seg1vec = { 1, 2 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		segments = { seg1 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");

		// alpha = 3 segments
		seg1vec = { 1 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		seg2vec = { 2 };
		seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		vector<int> seg3vec = {  };
		shared_ptr<TraceSegment> seg3 = make_shared<TraceSegment>(make_shared<vector<int>>(seg3vec));
		segments = { seg1, seg2, seg3 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");
	}

	// Bigger Tree with longer paths. alpha is already contained.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> ggc1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		root->add(make_shared<TreeEdge>(2, c2));
		c1->add(make_shared<TreeEdge>(2, gc1));
		c2->add(make_shared<TreeEdge>(3, gc2));
		gc1->add(make_shared<TreeEdge>(3, ggc1));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 2 segments. alpha is prefix of a path in Tree
		vector<int> seg1vec = { 1 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 2 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1, seg2 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");

		/// alpha = 2 segments. alpha is prefix of a path in Tree
		seg1vec = { 1, 2 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		seg2vec = { 3 };
		seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		segments = { seg1, seg2 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 0,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 0 if alpha is already contained in Tree.");
	}

	// case 2: prefix of path alpha reaches leaf of tree (no new branch needed) -> returns 1

	// root of the Tree is a leaf. alpha is not empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 2 segments.
		vector<int> seg1vec = { 1 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 2 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1, seg2 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 1,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 1 if prefix of alpha reaches leaf of Tree.");
	}

	// root of the Tree has two childs (c1 and c2). c1 and c2 are leaves. Prefix of alpha reaches c1. 
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		root->add(make_shared<TreeEdge>(2, c2));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 1 segment (<1,2>).
		vector<int> seg1vec = { 1, 2 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 1,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 1 if prefix of alpha reaches leaf of Tree.");

		// alpha = 2 segments (<1>.<2>).
		seg1vec = { 1 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 2 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		segments = { seg1, seg2 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 1,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 1 if prefix of alpha reaches leaf of Tree.");
	}

	// root of the Tree has two childs (c1 and c2). c1 has one child (gc1).  gc1 is a leaf. c2 has a child (gc2). gc2 is a leaf.
	// Prefix of alpha reaches gc1. 
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> c2 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc2 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		c1->add(make_shared<TreeEdge>(2, gc1));
		root->add(make_shared<TreeEdge>(2, c2));
		c2->add(make_shared<TreeEdge>(3, gc2));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 2 segments (<1,2>.<3>).
		vector<int> seg1vec = { 1, 2 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 3 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1, seg2 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 1,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 1 if prefix of alpha reaches leaf of Tree.");

		// alpha = 2 segments (<1>.<2,3>).
		seg1vec = { 1 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		seg2vec = { 2, 3 };
		seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		segments = { seg1, seg2 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 1,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 1 if prefix of alpha reaches leaf of Tree.");
	}

	// case 3: path not fully contained and no prefix of path reaches leaf of tree (new branch needed) -> returns 2

	// root of the Tree has one childs (c1). c1 is a leaf. First element of alpha has no corresponding edge in the Tree.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 1 segment (<2>).
		vector<int> seg1vec = { 2 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 2,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 2 if prefix of alpha reaches leaf of Tree.");

		// alpha = 2 segments (<2>.<1>)
		vector<int> seg2vec = { 1 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		segments = { seg1, seg2 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 2,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 2 if prefix of alpha reaches leaf of Tree.");

	}

	// root of the Tree has one child (c1). c1 has one child (gc1). gc1 is a leaf.
	// Prefix of alpha matches a prefix of a path contained in the Tree, but alpha contains additional elements (new Branch needed).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		shared_ptr<TreeNode> c1 = make_shared<TreeNode>();
		shared_ptr<TreeNode> gc1 = make_shared<TreeNode>();
		root->add(make_shared<TreeEdge>(1, c1));
		c1->add(make_shared<TreeEdge>(2, gc1));
		shared_ptr<Tree> tree = make_shared<Tree>(root, make_shared<FsmPresentationLayer>());

		// alpha = 1 segment (<1, 1>).
		vector<int> seg1vec = { 1, 1 };
		shared_ptr<TraceSegment> seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		std::deque< std::shared_ptr<TraceSegment> > segments = { seg1 };
		SegmentedTrace alpha(segments);

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 2,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 2 if prefix of alpha reaches leaf of Tree.");

		// alpha = 2 segments (<1>.<1>)
		seg1vec = { 1 };
		seg1 = make_shared<TraceSegment>(make_shared<vector<int>>(seg1vec));
		vector<int> seg2vec = { 1 };
		shared_ptr<TraceSegment> seg2 = make_shared<TraceSegment>(make_shared<vector<int>>(seg2vec));
		segments = { seg1, seg2 };
		alpha = segments;

		fsmlib_assert("TC-Tree-NNNN",
			tree->tentativeAddToRoot(alpha) == 2,
			"Tree::tentativeAddToRoot(SegmentedTrace& alpha) return 2 if prefix of alpha reaches leaf of Tree.");

	}
}

//===================================== FsmPresentationLayer Tests ===================================================

// tests FsmPresentationLayer(const std::string & inputs, const std::string & outputs, const std::string & states);
void testFsmPresentationLayerFileConstructor() {
	// Correct file names. No file is empty.
	shared_ptr<FsmPresentationLayer> pl =
		make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
			"../../../resources/garageOut.txt",
			"../../../resources/garageState.txt");

	vector<string> inputs{ "e1", "e2", "e3", "e4" };
	vector<string> outputs{ "a0", "a1", "a2", "a3", "a4"};
	vector<string> states{ "Door Up", "Door Down", "Door stopped going down", "Door stopped going up", "Door closing", "Door opening" };
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl->getIn2String() == inputs
		&& pl->getOut2String() == outputs
		&& pl->getState2String() == states,
		"FsmPresentationLayer(const std::string & inputs, const std::string & outputs, const std::string & states) correctly initializes "
		"contents of in2String, out2String and state2String");

	// input file doesn't exist
	pl =
		make_shared<FsmPresentationLayer>("nonExistingFileName.txt",
			"../../../resources/garageOut.txt",
			"../../../resources/garageState.txt");
	inputs = {  };
	outputs = { "a0", "a1", "a2", "a3", "a4" };
	states = { "Door Up", "Door Down", "Door stopped going down", "Door stopped going up", "Door closing", "Door opening" };
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl->getIn2String() == inputs
		&& pl->getOut2String() == outputs
		&& pl->getState2String() == states,
		"FsmPresentationLayer(const std::string & inputs, const std::string & outputs, const std::string & states) correctly initializes "
		"contents of in2String, out2String and state2String");

	// Using the same file name for each parameter.
	pl =
		make_shared<FsmPresentationLayer>("../../../resources/garageIn.txt",
			"../../../resources/garageIn.txt",
			"../../../resources/garageIn.txt");
	inputs = { "e1", "e2", "e3", "e4" };
	outputs = { "e1", "e2", "e3", "e4" };
	states = { "e1", "e2", "e3", "e4" };
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl->getIn2String() == inputs
		&& pl->getOut2String() == outputs
		&& pl->getState2String() == states,
		"FsmPresentationLayer(const std::string & inputs, const std::string & outputs, const std::string & states) correctly initializes "
		"contents of in2String, out2String and state2String");

	// One file is completely empty
	pl =
		make_shared<FsmPresentationLayer>("../../../resources/emptyIn.txt",
			"../../../resources/garageOut.txt",
			"../../../resources/garageState.txt");
	inputs = {  };
	outputs = { "a0", "a1", "a2", "a3", "a4" };
	states = { "Door Up", "Door Down", "Door stopped going down", "Door stopped going up", "Door closing", "Door opening" };
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl->getIn2String() == inputs
		&& pl->getOut2String() == outputs
		&& pl->getState2String() == states,
		"FsmPresentationLayer(const std::string & inputs, const std::string & outputs, const std::string & states) correctly initializes "
		"contents of in2String, out2String and state2String");
}

// tests FsmPresentationLayer::dumpIn(std::ostream & out)
void testFsmPresentationLayerDumpIn() {
	// in2String is empty
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	ostringstream stream;
	pl->dumpIn(stream);
	string s = stream.str();
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		s.empty(),
		"FsmPresentationLayer::dumpIn(std::ostream & out) writes each element of in2String to out.");

	stream.str("");
	stream.clear();

	// in2String contains only one element
	vector<string> in2String{ "e1" };
	vector<string> out2String{ "o1", "o2" };
	vector<string> state2String{ "s1" };
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	pl->dumpIn(stream);
	s = stream.str();
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		s == "e1",
		"FsmPresentationLayer::dumpIn(std::ostream & out) writes each element of in2String to out.");

	stream.str("");
	stream.clear();

	// in2String contains two elements
	in2String = { "e1", "e2" };
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	pl->dumpIn(stream);
	s = stream.str();
	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		s == "e1\ne2",
		"FsmPresentationLayer::dumpIn(std::ostream & out) writes each element of in2String to out.");
}

// tests FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer> otherPresentationLayer)
// Positive Case
void testFsmPresentationLayerComparePositive() {
	// in2String1 and in2String2 are empty. out2String1 and out2String2 are empty. state2String1 == state2String
	vector<string> in2String1{};
	vector<string> out2String1{};
	vector<string> state2String1{};
	shared_ptr<FsmPresentationLayer> pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);

	vector<string> in2String2{};
	vector<string> out2String2{};
	vector<string> state2String2{};
	shared_ptr<FsmPresentationLayer> pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl1->compare(pl2)
		&& pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"true if both in2String and out2String lists are equal.");

	// in2String and out2String lists are equals but not empty. state2String lists differ in size.
	in2String1 = { "e1" };
	in2String2 = { "e1" };
	out2String1 = { "o1", "o2" };
	out2String2 = { "o1", "o2" };
	state2String1 = {};
	state2String2 = { "s1" };

	pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);
	pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		pl1->compare(pl2)
		&& pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"true if both in2String and out2String lists are equal.");

}

// tests FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer> otherPresentationLayer)
// Positive Case
void testFsmPresentationLayerCompareNegative() {
	// in2String lists differ in size. out2String lists are equal.
	vector<string> in2String1{};
	vector<string> out2String1{"o1"};
	vector<string> state2String1{};
	shared_ptr<FsmPresentationLayer> pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);

	vector<string> in2String2{ "e1" };
	vector<string> out2String2{"o1"};
	vector<string> state2String2{};
	shared_ptr<FsmPresentationLayer> pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		!pl1->compare(pl2)
		&& !pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"false if both in2String lists differ in size");

	// in2String lists are equals. out2String lists differ in size.
	in2String1 = { "e1" };
	in2String2 = { "e1" };
	out2String1 = { "o1" };
	out2String2 = {  };
	state2String1 = {};
	state2String2 = {};

	pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);
	pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		!pl1->compare(pl2)
		&& !pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"false if both out2String lists differ in size.");

	// in2String lists are equal. out2String lists have the same size but contain at least one different element.
	in2String1 = { "e1", "e2" };
	in2String2 = { "e1", "e2" };
	out2String1 = { "o1", "o2" };
	out2String2 = { "o1", "o3" };
	state2String1 = {};
	state2String2 = {};

	pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);
	pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		!pl1->compare(pl2)
		&& !pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"false if both out2String lists contain different elements.");

	// out2String lists are equal. in2String lists have the same size but contain at least one different element.
	in2String1 = { "e1", "e2" };
	in2String2 = { "e2", "e1" };
	out2String1 = { "o1", "o2" };
	out2String2 = { "o1", "o2" };
	state2String1 = {};
	state2String2 = {};

	pl1 = make_shared<FsmPresentationLayer>(in2String1, out2String1, state2String1);
	pl2 = make_shared<FsmPresentationLayer>(in2String2, out2String2, state2String2);

	fsmlib_assert("TC-FsmPresentationLayer-NNNN",
		!pl1->compare(pl2)
		&& !pl2->compare(pl1),
		"FsmPresentationLayer::compare(std::shared_ptr<FsmPresentationLayer>) returns "
		"false if both out2String lists contain different elements.");
}

//===================================== Trace Tests ===================================================

// tests operator==(Trace const & trace1, Trace const & trace2)
// Positive case.
void testTraceEquals1Positive() {
	// tr1 and tr2 are empty.
	vector<int> v1 = {};
	Trace tr1{ v1, make_shared<FsmPresentationLayer>() };

	vector<int> v2 = {};
	Trace tr2{ v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");

	// tr1 and tr2 both contain one element.
	v1 = {1};
	v2 = {1};
	tr1 = { v1, make_shared<FsmPresentationLayer>() };
	tr2 = { v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");

	// tr1 and tr2 both contain two elements.
	v1 = { 1, 2 };
	v2 = { 1, 2 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };
	tr2 = { v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");
}

// tests operator==(Trace const & trace1, Trace const & trace2)
// Negative case.
void testTraceEquals1Negative() {
	// tr1 is empty. tr2 isn't empty.
	vector<int> v1 = {};
	Trace tr1{ v1, make_shared<FsmPresentationLayer>() };

	vector<int> v2 = { 1 };
	Trace tr2{ v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr2 is empty. tr1 isn't empty.
	v1 = { 1 };
	v2 = {  };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };
	tr2 = { v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr1 and tr2 have the same size but contain different elements.
	v1 = { 1 };
	v2 = { 2 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };
	tr2 = { v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr1 and tr2 have the same size but contain different elements.
	v1 = { 1, 2 };
	v2 = { 1, 3 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };
	tr2 = { v2, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");
}

// tests operator==(Trace const & trace1, std::vector<int> const & trace2)
// Positive case.
void testTraceEquals2Positive() {
	// tr1 and tr2 are empty.
	vector<int> v1 = {};
	Trace tr1{ v1, make_shared<FsmPresentationLayer>() };

	vector<int> tr2 = {};

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");

	// tr1 and tr2 both contain one element.
	v1 = { 1 };
	tr2 = { 1 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");

	// tr1 and tr2 both contain two elements.
	v1 = { 1, 2 };
	tr2 = { 1, 2 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		tr1 == tr2,
		"tr1 == tr2 if the underlying vectors are equal.");
}

// tests operator==(Trace const & trace1, std::vector<int> const & trace2)
// Negative case.
void testTraceEquals2Negative() {
	// tr1 is empty. tr2 isn't empty.
	vector<int> v1 = {};
	Trace tr1{ v1, make_shared<FsmPresentationLayer>() };

	vector<int> tr2 = { 1 };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr2 is empty. tr1 isn't empty.
	v1 = { 1 };
	tr2 = {};
	tr1 = { v1, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr1 and tr2 have the same size but contain different elements.
	v1 = { 1 };
	tr2 = { 2 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");

	// tr1 and tr2 have the same size but contain different elements.
	v1 = { 1, 2 };
	tr2 = { 1, 3 };
	tr1 = { v1, make_shared<FsmPresentationLayer>() };

	fsmlib_assert("TC-Trace-NNNN",
		not (tr1 == tr2),
		"tr1 == tr2 is false if the underlying vectors are unequal.");
}

// tests operator<<(std::ostream & out, const Trace & trace)
void testTraceOutputOperator() {
	// trace is empty.
	vector<int> v1 = {};
	Trace tr1{ v1, make_shared<FsmPresentationLayer>() };
	ostringstream out;
	out << tr1;
	string result = out.str();

	fsmlib_assert("TC-Trace-NNNN",
		result == "",
		"operator<<(std::ostream & out, const Trace & trace) writes every element of trace to out in the right order.");

	out.str("");
	out.clear();

	// trace contains one element.
	v1 = { 1 };
	tr1 = { v1,  make_shared<FsmPresentationLayer>() };
	out << tr1;
	result = out.str();

	fsmlib_assert("TC-Trace-NNNN",
		result == "1",
		"operator<<(std::ostream & out, const Trace & trace) writes every element of trace to out in the right order.");

	out.str("");
	out.clear();

	// trace contains two elements.
	v1 = { 1,2 };
	tr1 = { v1,  make_shared<FsmPresentationLayer>() };
	out << tr1;
	result = out.str();

	fsmlib_assert("TC-Trace-NNNN",
		result == "1.2",
		"operator<<(std::ostream & out, const Trace & trace) writes every element of trace to out in the right order.");
}

//===================================== InputTrace Tests ===================================================

// tests operator<<(std::ostream & out, const InputTrace & trace)
void testInputTraceOutputOperator() {
	// in2String of pl is empty. inTrace contains one element.
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	vector<int> inVec{0};
	InputTrace inTrace{ inVec, pl };
	ostringstream out;
	out << inTrace;
	string result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "0",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// in2String of pl is empty. inTrace contains two elements.
	inVec = { 0, 3 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "0.3",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// in2String of pl contains one element. inTrace contains one element.
	vector<string> in2String{ "e0" };
	vector<string> out2String{ };
	vector<string> state2String{ };
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	inVec = { 0 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "e0",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// in2String contains one element. inTrace contains one element >= in2String.size()
	inVec = { 1 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "1",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// in2String contains one element. inTrace contains two elements. One elements is smaller than in2String.size() and one is equal to in2String.size().
	inVec = { 0,1 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "e0.1",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// in2String contains two elements. inTrace contains three elements. All elements are smaller than in2String.size().
	in2String = { "e0", "e1" };
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	inVec = { 0,1,0 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "e0.e1.e0",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();


	// in2String contains two elements. inTrace contains three elements. Two elements are smaller than in2String.size(). One element is equal to this size.
	inVec = { 1,2,1 };
	inTrace = { inVec, pl };
	out << inTrace;
	result = out.str();
	fsmlib_assert("TC-InputTrace-NNNN",
		result == "e1.2.e1",
		"operator<<(std::ostream & out, const InputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();
}

//===================================== OutputTrace Tests ===================================================

// tests operator<<(std::ostream & out, const OutputTrace & trace)
void testOutputTraceOutputOperator() {
	// out2String of pl is empty. outTrace contains one element.
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	vector<int> outVec{ 0 };
	OutputTrace outTrace{ outVec, pl };
	ostringstream out;
	out << outTrace;
	string result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "0",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// out2String of pl is empty. outTrace contains two elements.
	outVec = { 0, 3 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "0.3",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// out2String of pl contains one element. outTrace contains one element.
	vector<string> in2String{};
	vector<string> out2String{ "e0" };
	vector<string> state2String{};
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	outVec = { 0 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "e0",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// out2String contains one element. outTrace contains one element >= out2String.size()
	outVec = { 1 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "1",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// out2String contains one element. outTrace contains two elements. One elements is smaller than out2String.size() and one is equal to out2String.size().
	outVec = { 0,1 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "e0.1",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();

	// out2String contains two elements. outTrace contains three elements. All elements are smaller than out2String.size().
	out2String = { "e0", "e1" };
	pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
	outVec = { 0,1,0 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "e0.e1.e0",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();


	// out2String contains two elements. outTrace contains three elements. Two elements are smaller than out2String.size(). One element is equal to this size.
	outVec = { 1,2,1 };
	outTrace = { outVec, pl };
	out << outTrace;
	result = out.str();
	fsmlib_assert("TC-OutputTrace-NNNN",
		result == "e1.2.e1",
		"operator<<(std::ostream & out, const OutputTrace & trace) writes every element of trace to out in the right order and with the right name");

	out.str("");
	out.clear();
}

//===================================== OutputTree Tests ===================================================

// tests OutputTree::contains(OutputTree& ot)
// Negative Case.
void testOutputTreeContainsNegative() {
	// inputTrace of this Tree is empty. otherInputTrace of otherTree isn't empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{1};
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if the InputTraces differ.");
	}

	// inputTrace of this Tree isn't empty. otherInputTrace of otherTree is empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if the InputTraces differ.");
	}

	// inputTraces both aren't empty. inputTraces are unequal. 
	// Both Trees contain the same output trace.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 1 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{1,2};
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 1 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if the InputTraces differ.");
	}

	// inputTraces are equal (both empty).
	// tree contains empty trace. otherTree contains one output trace which is not empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{  };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{  };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 1 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}

	// inputTraces are equal (both empty).
	// tree contains only output trace [2]. otherTree contains only output trace [1].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 2 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{};
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 1 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}

	// inputTraces are equal.
	// tree contains only output trace [1,2]. otherTree contains only output trace [1] (prefix)).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{1, 2};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 1, 2 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1, 2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 1 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}

	// inputTraces are equal.
	// tree contains output traces [1,2] and [2,3]. otherTree contains output traces [1,2] and [2,2] (they share one Trace)).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1, 2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 1, 2 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 2, 3 };
		tree.addToRoot(outVec2);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1, 2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec1{ 1, 2 };
		otherTree.addToRoot(otherOutVec1);
		vector<int> otherOutVec2{ 2, 2 };
		otherTree.addToRoot(otherOutVec2);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}

	// inputTraces are equal.
	// tree contains output trace [1,2]. otherTree contains output trace [1,2,1].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1, 2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 1, 2 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1, 2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 1, 2, 1 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}

	// inputTraces are equal.
	// tree contains output trace [0]. otherTree contains only empty trace [].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 0 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);

		fsmlib_assert("TC-OutputTree-NNNN",
			not tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns false if otherTree contains an output trace which is not contained in tree.");
	}
}

// tests OutputTree::contains(OutputTree& ot)
// Positive Case.
void testOutputTreeContainsPositive() {
	// all Traces are empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{};
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);

		fsmlib_assert("TC-OutputTree-NNNN",
			tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns true if the InputTraces are the same and each output trace contained in otherTree is "
			"also contained in tree.");
	}

	// inputTraces are equal.
	// tree contains only output trace [0]. otherTree contains only output trace [0].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 0 };
		tree.addToRoot(outVec);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 0 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns true if the InputTraces are the same and each output trace contained in otherTree is "
			"also contained in tree.");
	}

	// inputTraces are equal.
	// tree contains output traces [0,0] and [0,2]. otherTree contains only output trace [0,2].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1,2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 0,0 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 0,2 };
		tree.addToRoot(outVec2);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1,2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 0,2 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns true if the InputTraces are the same and each output trace contained in otherTree is "
			"also contained in tree.");
	}

	// inputTraces are equal.
	// tree contains output traces [0,0] and [0,2]. otherTree contains only output trace [0,0].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1,2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 0,0 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 0,2 };
		tree.addToRoot(outVec2);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1,2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec{ 0,0 };
		otherTree.addToRoot(otherOutVec);

		fsmlib_assert("TC-OutputTree-NNNN",
			tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns true if the InputTraces are the same and each output trace contained in otherTree is "
			"also contained in tree.");
	}

	// inputTraces are equal.
	// Both Trees are equal and not empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 1,2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 0,0 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 0,2 };
		tree.addToRoot(outVec2);

		shared_ptr<TreeNode> otherRoot = make_shared<TreeNode>();
		vector<int> otherInVec{ 1,2 };
		InputTrace otherInputTrace(otherInVec, pl);
		OutputTree otherTree(otherRoot, otherInputTrace, pl);
		vector<int> otherOutVec1{ 0,0 };
		otherTree.addToRoot(otherOutVec1);
		vector<int> otherOutVec2{ 0,2 };
		otherTree.addToRoot(otherOutVec2);

		fsmlib_assert("TC-OutputTree-NNNN",
			tree.contains(otherTree),
			"OutputTree::contains(OutputTree& ot) returns true if the InputTraces are the same and each output trace contained in otherTree is "
			"also contained in tree.");
	}
}

// tests OutputTree::toDot(std::ostream& out) 
void testOutputTreeToDot() {
	// root is a leaf. InputTrace is empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		std::ostringstream stream;
		tree.toDot(stream);
		string content = stream.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			content.find(" -> ") == string::npos,
			"result of toDot(ostream & out) contains only expected edges");
	}

	// OutputTree contains only one path [1]. InputTrace is [0].
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{0};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec{ 1 };
		tree.addToRoot(outVec); 

		std::ostringstream stream;
		tree.toDot(stream);
		string content = stream.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			content.find("0 -> 1[label = \"0/1\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-OutputTree-NNNN",
			countEdgesInToDotResult(content) == 1,
			"result of toDot(ostream & out) contains only expected edges");
	}

	// OutputTree contains two paths [1] and [2]. InputTrace is [0] (empty). Length of InputTrace is equal to the length of each Trace contained in 
	// OutputTree.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 1 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec2);

		std::ostringstream stream;
		tree.toDot(stream);
		string content = stream.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			content.find("0 -> 1[label = \"0/1\" ];") != string::npos
			&& content.find("0 -> 2[label = \"0/2\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-OutputTree-NNNN",
			countEdgesInToDotResult(content) == 2,
			"result of toDot(ostream & out) contains only expected edges");
	}

	// OutputTree contains two paths [1] and [2]. InputTrace is [0,1]. Length of InputTrace is greater than the length of at least one Trace
	// contained in OutputTree.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0,1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 1 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec2);

		std::ostringstream stream;
		tree.toDot(stream);
		string content = stream.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			content.find("0 -> 1[label = \"0/1\" ];") != string::npos
			&& content.find("0 -> 2[label = \"0/2\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-OutputTree-NNNN",
			countEdgesInToDotResult(content) == 2,
			"result of toDot(ostream & out) contains only expected edges");
	}

	// OutputTree contains paths [1,1], [1,2] and [2]. InputTrace is [0,1]. Length of InputTrace is greater than the length of at least one Trace
	// contained in OutputTree.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0,1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 1, 1 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 1, 2 };
		tree.addToRoot(outVec2);
		vector<int> outVec3{ 2 };
		tree.addToRoot(outVec3);

		std::ostringstream stream;
		tree.toDot(stream);
		string content = stream.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			content.find("0 -> 1[label = \"0/1\" ];") != string::npos
			&& content.find("1 -> 2[label = \"1/1\" ];") != string::npos
			&& content.find("1 -> 3[label = \"1/2\" ];") != string::npos
			&& content.find("0 -> 4[label = \"0/2\" ];") != string::npos,
			"result of toDot(ostream & out) contains each expected edge");

		fsmlib_assert("TC-OutputTree-NNNN",
			countEdgesInToDotResult(content) == 4,
			"result of toDot(ostream & out) contains only expected edges");
	}
}

// tests OutputTree::getOutputTraces()
void testOutputTreeGetOutputTraces() {
	// root is a leaf. InputTrace is empty.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		vector<int> outVec{};
		OutputTrace outputTrace(outVec, pl);

		vector<OutputTrace> result = tree.getOutputTraces();

		fsmlib_assert("TC-OutputTree-NNNN",
			result.size() == 1
			&& result.at(0) == outputTrace,
			"OutputTree::getOutputTraces() called on an OutputTree which consists only of a root, returns only an empty OutputTrace.");
	}

	// OutputTree contains two Traces ([1] and [2]).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{0};
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		vector<int> outVec1{ 1 };
		tree.addToRoot(outVec1);
		OutputTrace outputTrace1(outVec1, pl);
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec2);
		OutputTrace outputTrace2(outVec2, pl);
		
		vector<OutputTrace> result = tree.getOutputTraces();

		fsmlib_assert("TC-OutputTree-NNNN",
			find(result.cbegin(), result.cend(), outputTrace1) != result.cend()
			&& find(result.cbegin(), result.cend(), outputTrace2) != result.cend(),
			"result of OutputTree::getOutputTraces() contains each expected OutputTrace");
		fsmlib_assert("TC-OutputTree-NNNN",
			result.size() == 2,
			"result of OutputTree::getOutputTraces() contains only expected OutputTraces");
	}

	// OutputTree contains three Traces ([1,1], [1,2] and [2]).
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0, 1 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);		

		vector<int> outVec1{ 1,1 };
		tree.addToRoot(outVec1);
		OutputTrace outputTrace1(outVec1, pl);
		vector<int> outVec2{ 1,2 };
		tree.addToRoot(outVec2);
		OutputTrace outputTrace2(outVec2, pl);
		vector<int> outVec3{ 2 };
		tree.addToRoot(outVec3);
		OutputTrace outputTrace3(outVec3, pl);

		vector<OutputTrace> result = tree.getOutputTraces();

		fsmlib_assert("TC-OutputTree-NNNN",
			find(result.cbegin(), result.cend(), outputTrace1) != result.cend()
			&& find(result.cbegin(), result.cend(), outputTrace2) != result.cend()
			&& find(result.cbegin(), result.cend(), outputTrace3) != result.cend(),
			"result of OutputTree::getOutputTraces() contains each expected OutputTrace");
		fsmlib_assert("TC-OutputTree-NNNN",
			result.size() == 3,
			"result of OutputTree::getOutputTraces() contains only expected OutputTraces");
	}
}

// tests operator<<(std::ostream& out, OutputTree& ot)
void testOutputTreeOutputOperator() {
	// root of OutputTree is a leaf (tree contains only the empty trace). FsmPresentationLayer contains no names.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{  };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);

		ostringstream out;
		out << tree;
		string result = out.str();
		fsmlib_assert("TC-OutputTree-NNNN",
			result == "\n",   
			"operator<<(std::ostream& out, OutputTree& ot) writes each trace contained in OutputTree with the right names to out.");
	}

	// OutputTree contains Traces [0] and [1]. InputTrace = [0].  FsmPresentationLayer contains names for each Input/Output.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		vector<string> in2String{ "e0", "e1" };
		vector<string> out2String{ "o0", "o1" };
		vector<string> state2String;
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
		vector<int> inVec{ 0 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 0 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 1 };
		tree.addToRoot(outVec2);

		ostringstream out;
		out << tree;
		string result = out.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			result.find("(e0/o0)") != string::npos
			&& result.find("(e0/o1)") != string::npos,
			"operator<<(std::ostream& out, OutputTree& ot) writes each trace contained in OutputTree with the right names to out.");
	}

	// OutputTree contains Traces [0,2] and [1]. InputTrace = [0,2].  FsmPresentationLayer contains names for some but not for each Input/Output.
	{
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		vector<string> in2String{ "e0", "e1" };
		vector<string> out2String{ "o0", "o1" };
		vector<string> state2String;
		std::shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>(in2String, out2String, state2String);
		vector<int> inVec{ 0, 2 };
		InputTrace inputTrace(inVec, pl);
		OutputTree tree(root, inputTrace, pl);
		vector<int> outVec1{ 0, 2 };
		tree.addToRoot(outVec1);
		vector<int> outVec2{ 1 };
		tree.addToRoot(outVec2);

		ostringstream out;
		out << tree;
		string result = out.str();

		fsmlib_assert("TC-OutputTree-NNNN",
			result.find("(e0/o0).(2/2)") != string::npos
			&& result.find("(e0/o1)") != string::npos,
			"operator<<(std::ostream& out, OutputTree& ot) writes each trace contained in OutputTree with the right names to out.");
	}
}

//===================================== TestSuite Tests ===================================================

// tests TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput)
// Positive Case
// Additional tests for TestSuite::isReductionOf(TestSuite& theOtherTs,bool writeOutput): if ts1 equals ts2 then ts1 is a reduction
// of ts2 and vice versa
void testTestSuiteIsEquivalentToPositive() {
	// empty TestSuites
	{
		TestSuite ts1;
		TestSuite ts2;
		fsmlib_assert("TC-TestSuite-NNNN",
			ts1.isEquivalentTo(ts2)
			&& ts2.isEquivalentTo(ts1),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns true if both TestSuites are empty.");

		fsmlib_assert("TC-TestSuite-NNNN",
			ts1.isReductionOf(ts2)
			&& ts2.isReductionOf(ts1),
			"If TestSuites are equivalent they are also reductions of each other");
	}

	// both TestSuites contain one test case (OutputTree).
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec1);
		tree.addToRoot(outVec2);
		TestSuite ts;
		ts.push_back(tree);

		shared_ptr<FsmPresentationLayer> o_pl = make_shared<FsmPresentationLayer>();
		vector<int> o_inVec{ 0 };
		InputTrace o_inTrc{ o_inVec, o_pl };
		shared_ptr<TreeNode> o_root = make_shared<TreeNode>();
		OutputTree o_tree{ o_root, o_inTrc, o_pl };
		vector<int> o_outVec1{ 1 };
		vector<int> o_outVec2{ 2 };
		o_tree.addToRoot(o_outVec1);
		o_tree.addToRoot(o_outVec2);
		TestSuite o_ts;
		o_ts.push_back(o_tree);

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isEquivalentTo(o_ts)
			&& o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns true if both "
			"TestSuites contain the same OutputTrees in the same order.");

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isReductionOf(o_ts)
			&& o_ts.isReductionOf(ts),
			"If TestSuites are equivalent they are also reductions of each other");
	}

	// both TestSuites contain two test cases (OutputTrees).
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0,1 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1,1 };
		vector<int> outVec2Tree1{ 1,2 };
		vector<int> outVec3Tree1{ 2,3 };
		tree1.addToRoot(outVec1Tree1);
		tree1.addToRoot(outVec2Tree1);
		tree1.addToRoot(outVec3Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 1,1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 3,3 };
		vector<int> outVec2Tree2{ 4,4 };
		tree2.addToRoot(outVec1Tree2);
		tree2.addToRoot(outVec2Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		TestSuite o_ts;
		o_ts.push_back(tree1);
		o_ts.push_back(tree2);

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isEquivalentTo(o_ts)
			&& o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns true if both "
			"TestSuites contain the same OutputTrees in the same order.");

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isReductionOf(o_ts)
			&& o_ts.isReductionOf(ts),
			"If TestSuites are equivalent they are also reductions of each other");
	}
}

// tests TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput)
// Negative Case
// Additional tests for TestSuite::isReductionOf(TestSuite& theOtherTs,bool writeOutput): if not(ts1isEquivalentTo(ts2) then ts1 and ts2 can't be
// reductions of each other
void testTestSuiteIsEquivalentToNegative() {
	// TestSuites contain different number of OutputTrees.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec1);
		tree.addToRoot(outVec2);
		TestSuite ts;
		ts.push_back(tree);

		TestSuite o_ts;

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isEquivalentTo(o_ts)
			&& not o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns false if both TestSuites have different sizes.");

		fsmlib_assert("TC-TestSuite-NNNN",
			not (ts.isReductionOf(o_ts) && o_ts.isReductionOf(ts)),
			"If TestSuites aren't equal they can't be reductions of each other.");
	}

	// both TestSuites have the same number of OutputTrees (1), but they aren't equal. 
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec1);
		tree.addToRoot(outVec2);
		TestSuite ts;
		ts.push_back(tree);

		shared_ptr<FsmPresentationLayer> o_pl = make_shared<FsmPresentationLayer>();
		vector<int> o_inVec{ 0 };
		InputTrace o_inTrc{ o_inVec, o_pl };
		shared_ptr<TreeNode> o_root = make_shared<TreeNode>();
		OutputTree o_tree{ o_root, o_inTrc, o_pl };
		vector<int> o_outVec1{ 1 };
		o_tree.addToRoot(o_outVec1);
		TestSuite o_ts;
		o_ts.push_back(o_tree);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isEquivalentTo(o_ts)
			&& not o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns false if one "
			"TestSuite contains an OutputTree that is not contained in the other TestSuite.");

		fsmlib_assert("TC-TestSuite-NNNN",
			not (ts.isReductionOf(o_ts) && o_ts.isReductionOf(ts)),
			"If TestSuites aren't equal they can't be reductions of each other.");
	}

	// both TestSuites contain two test cases (OutputTrees). The TestSuites differ in the second OutputTree (in the InputTraces).
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1 };
		vector<int> outVec2Tree1{ 2 };
		tree1.addToRoot(outVec1Tree1);
		tree1.addToRoot(outVec2Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 0,1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 1,3 };
		vector<int> outVec2Tree2{ 2,4 };
		tree2.addToRoot(outVec1Tree2);
		tree2.addToRoot(outVec2Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		vector<int> o_inVecTree2{ 0,0 };
		InputTrace o_inTrcTree2{ o_inVecTree2, pl };
		shared_ptr<TreeNode> o_rootTree2 = make_shared<TreeNode>();
		OutputTree o_tree2{ o_rootTree2, o_inTrcTree2, pl };
		vector<int> o_outVec1Tree2{ outVec1Tree2 };
		vector<int> o_outVec2Tree2{ outVec2Tree2 };
		o_tree2.addToRoot(o_outVec1Tree2);
		o_tree2.addToRoot(o_outVec2Tree2);

		TestSuite o_ts;
		o_ts.push_back(tree1);
		o_ts.push_back(o_tree2);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isEquivalentTo(o_ts)
			&& not o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns false if one "
			"TestSuite contains an OutputTree that is not contained in the other TestSuite.");

		fsmlib_assert("TC-TestSuite-NNNN",
			not (ts.isReductionOf(o_ts) && o_ts.isReductionOf(ts)),
			"If TestSuites aren't equal they can't be reductions of each other.");
	}

	// both TestSuites contain the same two test cases (OutputTrees) in different order.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1 };
		vector<int> outVec2Tree1{ 2 };
		tree1.addToRoot(outVec1Tree1);
		tree1.addToRoot(outVec2Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 0,1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 1,3 };
		vector<int> outVec2Tree2{ 2,4 };
		tree2.addToRoot(outVec1Tree2);
		tree2.addToRoot(outVec2Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		TestSuite o_ts;
		o_ts.push_back(tree2);
		o_ts.push_back(tree1);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isEquivalentTo(o_ts)
			&& not o_ts.isEquivalentTo(ts),
			"TestSuite::isEquivalentTo(TestSuite& theOtherTs,bool writeOutput) returns false if one "
			"TestSuite contains an OutputTree that is not contained in the other TestSuite.");

		fsmlib_assert("TC-TestSuite-NNNN",
			not (ts.isReductionOf(o_ts) && o_ts.isReductionOf(ts)),
			"If TestSuites aren't equal they can't be reductions of each other.");
	}
}

// tests TestSuite::isReductionOf(TestSuite& theOtherTs, bool writeOutput)
// Negative case.
void testTestSuiteIsReductionOfNegative() {
	// TestSuites contain different number of OutputTrees.
	{
		TestSuite ts;

		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec1);
		tree.addToRoot(outVec2);

		TestSuite o_ts;
		o_ts.push_back(tree);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isReductionOf(o_ts),
			"TestSuite::isReductionOf(TestSuite& theOtherTs, bool writeOutput) returns false if both TestSuites have different sizes.");
	}

	// Both TestSuites contain only one OutputTree. The OutputTree of ts contains a path that is not contained in the OutputTree of o_ts.
	{		
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		vector<int> outVec2{ 2 };
		tree.addToRoot(outVec1);
		tree.addToRoot(outVec2);

		TestSuite ts;
		ts.push_back(tree);

		vector<int> o_inVec{ 0 };
		InputTrace o_inTrc{ o_inVec, pl };
		shared_ptr<TreeNode> o_root = make_shared<TreeNode>();
		OutputTree o_tree{ o_root, o_inTrc, pl };
		vector<int> o_outVec1{ 1 };
		o_tree.addToRoot(o_outVec1);

		TestSuite o_ts;
		o_ts.push_back(o_tree);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isReductionOf(o_ts),
			"ts1.isReductionOf(ts2) returns false if ts1[i] contains a path that is not contained in ts2[i], for any 0 <= i < size.");
	}

	// Both TestSuites contain two OutputTrees. Two corresponding OutputTrees have different InputTraces.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1 };
		vector<int> outVec2Tree1{ 2 };
		tree1.addToRoot(outVec1Tree1);
		tree1.addToRoot(outVec2Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 1,1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 3,4 };
		tree2.addToRoot(outVec1Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		vector<int> o_inVecTree2{ 1,2 };
		InputTrace o_inTrcTree2{ o_inVecTree2, pl };
		shared_ptr<TreeNode> o_rootTree2 = make_shared<TreeNode>();
		OutputTree o_tree2{ o_rootTree2, o_inTrcTree2, pl };
		vector<int> o_outVec1Tree2{ outVec1Tree2 };
		o_tree2.addToRoot(o_outVec1Tree2);

		TestSuite o_ts;
		o_ts.push_back(tree1);
		o_ts.push_back(o_tree2);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isReductionOf(o_ts),
			"ts1.isReductionOf(ts2) returns false if ts1[i] and ts2[i] differ in the InputTrace, for any 0 <= i < size.");
	}

	// Both TestSuites contain two OutputTrees. Each OutputTree of ts is a contained in some OutputTree of o_ts, but ts has the wrong order.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1 };
		vector<int> outVec2Tree1{ 2 };
		tree1.addToRoot(outVec1Tree1);
		tree1.addToRoot(outVec2Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 3 };
		tree2.addToRoot(outVec1Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		vector<int> o_inVecTree2{ 1 };
		InputTrace o_inTrcTree2{ o_inVecTree2, pl };
		shared_ptr<TreeNode> o_rootTree2 = make_shared<TreeNode>();
		OutputTree o_tree2{ o_rootTree2, o_inTrcTree2, pl };
		vector<int> o_outVec1Tree2{ 3 };
		vector<int> o_outVec2Tree2{ 4 };
		o_tree2.addToRoot(o_outVec1Tree2);
		o_tree2.addToRoot(o_outVec2Tree2);

		TestSuite o_ts;
		o_ts.push_back(o_tree2);
		o_ts.push_back(tree1);

		fsmlib_assert("TC-TestSuite-NNNN",
			not ts.isReductionOf(o_ts),
			"ts1.isReductionOf(ts2) returns false if ts1[i] contains a path that is not contained in ts2[i], for any 0 <= i < size.");
	}
}

// tests TestSuite::isReductionOf(TestSuite& theOtherTs, bool writeOutput)
// Positive case.
void testTestSuiteIsReductionOfPositive() {
	// both TestSuites are empty.
	{
		TestSuite ts;
		TestSuite o_ts;
		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isReductionOf(o_ts)
			&& o_ts.isReductionOf(ts),
			"ts1.isReductionOf(ts2) returns true if ts1 and ts2 are equivalent.");
	}

	// both TestSuites contain one OutputTree. Each path of the OutputTree of ts is also a path in the OutputTree of o_ts.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> inVec{ 0 };
		InputTrace inTrc{ inVec, pl };
		shared_ptr<TreeNode> root = make_shared<TreeNode>();
		OutputTree tree{ root, inTrc, pl };
		vector<int> outVec1{ 1 };
		tree.addToRoot(outVec1);

		TestSuite ts;
		ts.push_back(tree);

		vector<int> o_inVec{ 0 };
		InputTrace o_inTrc{ o_inVec, pl };
		shared_ptr<TreeNode> o_root = make_shared<TreeNode>();
		OutputTree o_tree{ o_root, o_inTrc, pl };
		vector<int> o_outVec1{ 1 };
		vector<int> o_outVec2{ 2 };
		o_tree.addToRoot(o_outVec1);
		o_tree.addToRoot(o_outVec2);

		TestSuite o_ts;
		o_ts.push_back(o_tree);

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isReductionOf(o_ts),
			"ts1.isReductionOf(ts2) returns true if each path of ts1[i] is contained in ts2[i], for each 0 <= i < size.");
	}

	// Both TestSuites contain two OutputTrees. Each OutputTree of ts is contained in the corresponding OutputTree of o_ts (o_ts[i].contains(ts[i])).
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// constructing first tree
		vector<int> inVecTree1{ 0 };
		InputTrace inTrcTree1{ inVecTree1, pl };
		shared_ptr<TreeNode> rootTree1 = make_shared<TreeNode>();
		OutputTree tree1{ rootTree1, inTrcTree1, pl };
		vector<int> outVec1Tree1{ 1 };
		tree1.addToRoot(outVec1Tree1);

		// constructing second tree
		vector<int> inVecTree2{ 1,1 };
		InputTrace inTrcTree2{ inVecTree2, pl };
		shared_ptr<TreeNode> rootTree2 = make_shared<TreeNode>();
		OutputTree tree2{ rootTree2, inTrcTree2, pl };
		vector<int> outVec1Tree2{ 2,3 };
		tree2.addToRoot(outVec1Tree2);

		TestSuite ts;
		ts.push_back(tree1);
		ts.push_back(tree2);

		vector<int> o_inVecTree2{ 1, 1 };
		InputTrace o_inTrcTree2{ o_inVecTree2, pl };
		shared_ptr<TreeNode> o_rootTree2 = make_shared<TreeNode>();
		OutputTree o_tree2{ o_rootTree2, o_inTrcTree2, pl };
		vector<int> o_outVec1Tree2{ 2,2 };
		vector<int> o_outVec2Tree2{ 2,3 };
		o_tree2.addToRoot(o_outVec1Tree2);
		o_tree2.addToRoot(o_outVec2Tree2);

		TestSuite o_ts;
		o_ts.push_back(tree1);
		o_ts.push_back(o_tree2);

		fsmlib_assert("TC-TestSuite-NNNN",
			ts.isReductionOf(o_ts),
			"ts1.isReductionOf(ts2) returns true if each path of ts1[i] is contained in ts2[i], for each 0 <= i < size.");
	}
}

//===================================== IOListContainer Tests ===================================================

// This function can be used to check if the given element is a valid one (size between minLength and maxLength,
// each int contained in element is a value between 0 and maxInput). 
bool checkIOListContainerElement(const int maxInput, const int minLength, const int maxLength, const vector<int> &element) {
	if (element.size() < minLength || element.size() > maxLength) {
		return false;
	}
	for (const auto &i : element) {
		if (i < 0 || i > maxInput) {
			return false;
		}
	}
	return true;
}

// This function can be used to calculate the expected size of a constructed IOListContainer (using the 
// IOListContainer(const int maxInput, const int minLength, const int maxLength, const std::shared_ptr<FsmPresentationLayer>)
// constructor). (maxInput + 1)^minLength + (maxInput + 1)^(minLength+1) + ... + (maxInput + 1)^maxLength
int calculateExpectedIOListContainerSize(int maxInput, int minLength, int maxLength) {
	int expectedSize = 0;
	for (int i = minLength; i <= maxLength; ++i) {
		expectedSize += pow(maxInput + 1, i);
	}
	return expectedSize;
}

// tests IOListContainer::IOListContainer(const int maxInput, const int minLength, 
//             const int maxLenght, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
void testIOListContainerConstructor() {
	// maxInput = minLength = maxLength = 0.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		IOListContainer iolc(0, 0, 0, pl);
		shared_ptr<vector<vector<int>>> ioLists = iolc.getIOLists();

		fsmlib_assert("TC-IOListContainer-NNNN",
			ioLists->size() == 1
			&& ioLists->at(0).size() == 0,
			"IOListContainer::IOListContainer(const int maxInput, const int minLength, "
			"const int maxLenght, const std::shared_ptr<FsmPresentationLayer>) generates only the "
			"empty list if minLength = maxLength = 0");
	}

	// maxInput = minLength = 0, maxLength = 1
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxInput = 0;
		const int minLength = 0;
		const int maxLength = 1;
		IOListContainer iolc(maxInput, minLength, maxLength, pl);
		shared_ptr<vector<vector<int>>> ioLists = iolc.getIOLists();		

		int expectedSize = calculateExpectedIOListContainerSize(maxInput, minLength, maxLength);

		fsmlib_assert("TC-IOListContainer-NNNN",
			ioLists->size() == expectedSize,
			"Size of the constructed IOListContainer matches the expected size.");

		set<vector<int>> traces(ioLists->cbegin(), ioLists->cend());
		fsmlib_assert("TC-IOListContainer-NNNN",
			traces.size() == expectedSize,
			"Constructed IOListContainer only contains unique elements.");

		for (const vector<int> &trace : traces) {
			fsmlib_assert("TC-IOListContainer-NNNN",
				checkIOListContainerElement(maxInput, minLength, maxLength, trace),
				"Only valid elements are contained in the constructed IOListContainer.");
		}
	}

	//  maxInput = 1, minLength = 0, maxLength = 1
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxInput = 1;
		const int minLength = 0;
		const int maxLength = 1;
		IOListContainer iolc(maxInput, minLength, maxLength, pl);
		shared_ptr<vector<vector<int>>> ioLists = iolc.getIOLists();

		int expectedSize = calculateExpectedIOListContainerSize(maxInput, minLength, maxLength);

		fsmlib_assert("TC-IOListContainer-NNNN",
			ioLists->size() == expectedSize,
			"Size of the constructed IOListContainer matches the expected size.");

		set<vector<int>> traces(ioLists->cbegin(), ioLists->cend());
		fsmlib_assert("TC-IOListContainer-NNNN",
			traces.size() == expectedSize,
			"Constructed IOListContainer only contains unique elements.");

		for (const vector<int> &trace : traces) {
			fsmlib_assert("TC-IOListContainer-NNNN",
				checkIOListContainerElement(maxInput, minLength, maxLength, trace),
				"Only valid elements are contained in the constructed IOListContainer.");
		}
	}

	//  maxInput = 2, minLength = 1, maxLength = 1
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxInput = 2;
		const int minLength = 1;
		const int maxLength = 1;
		IOListContainer iolc(maxInput, minLength, maxLength, pl);
		shared_ptr<vector<vector<int>>> ioLists = iolc.getIOLists();

		int expectedSize = calculateExpectedIOListContainerSize(maxInput, minLength, maxLength);

		fsmlib_assert("TC-IOListContainer-NNNN",
			ioLists->size() == expectedSize,
			"Size of the constructed IOListContainer matches the expected size.");

		set<vector<int>> traces(ioLists->cbegin(), ioLists->cend());
		fsmlib_assert("TC-IOListContainer-NNNN",
			traces.size() == expectedSize,
			"Constructed IOListContainer only contains unique elements.");

		for (const vector<int> &trace : traces) {
			fsmlib_assert("TC-IOListContainer-NNNN",
				checkIOListContainerElement(maxInput, minLength, maxLength, trace),
				"Only valid elements are contained in the constructed IOListContainer.");
		}
	}

	//  maxInput = 2, minLength = 1, maxLength = 2
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxInput = 2;
		const int minLength = 1;
		const int maxLength = 2;
		IOListContainer iolc(maxInput, minLength, maxLength, pl);
		shared_ptr<vector<vector<int>>> ioLists = iolc.getIOLists();

		int expectedSize = calculateExpectedIOListContainerSize(maxInput, minLength, maxLength);

		fsmlib_assert("TC-IOListContainer-NNNN",
			ioLists->size() == expectedSize,
			"Size of the constructed IOListContainer matches the expected size.");

		set<vector<int>> traces(ioLists->cbegin(), ioLists->cend());
		fsmlib_assert("TC-IOListContainer-NNNN",
			traces.size() == expectedSize,
			"Constructed IOListContainer only contains unique elements.");

		for (const vector<int> &trace : traces) {
			fsmlib_assert("TC-IOListContainer-NNNN",
				checkIOListContainerElement(maxInput, minLength, maxLength, trace),
				"Only valid elements are contained in the constructed IOListContainer.");
		}
	}
}

//===================================== TraceSegment Tests ===================================================

// tests TraceSegment::getCopy()
void testTraceSegmentGetCopy() {
	// segment is empty. prefix = string::npos
	{
		vector<int> seg{};
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			seg == result,
			"TraceSegment::getCopy() returns a copy of the whole segment stored in the "
			"TraceSegment object if prefix equals string::npos.");
	}

	// segment.size() = 1 . prefix = string::npos
	{
		vector<int> seg{ 1 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			seg == result,
			"TraceSegment::getCopy() returns a copy of the whole segment stored in the "
			"TraceSegment object if prefix equals string::npos.");
	}

	// segment.size() = 0 . prefix = 0 (prefix >= segment.size() case)
	{
		vector<int> seg{ };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 0 };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			seg == result,
			"TraceSegment::getCopy() returns a copy of the whole segment stored in the "
			"TraceSegment object if prefix equals segment.size().");
	}

	// segment.size() = 1 . prefix = 1 (prefix >= segment.size() case)
	{
		vector<int> seg{ 1 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 1 };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			seg == result,
			"TraceSegment::getCopy() returns a copy of the whole segment stored in the "
			"TraceSegment object if prefix equals segment.size().");
	}

	// segment.size() = 2 . prefix = 2 (prefix >= segment.size() case)
	{
		vector<int> seg{ 1, 2 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 2 };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			seg == result,
			"TraceSegment::getCopy() returns a copy of the whole segment stored in the "
			"TraceSegment object if prefix equals segment.size().");
	}

	// segment.size() != 0 . prefix = 0 
	{
		vector<int> seg{ 1 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 0 };
		vector<int> result = trcSeg.getCopy();
		fsmlib_assert("TC-TraceSegment-NNNN",
			result.empty(),
			"TraceSegment::getCopy() returns empty vector if prefix is set to 0. ");
	}

	// segment.size() = 2 . prefix = 1 ( 0 < prefix < segment.size() case) 
	{
		vector<int> seg{ 0, 1 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 1 };
		vector<int> result = trcSeg.getCopy();
		vector<int> expected{ 0 };
		fsmlib_assert("TC-TraceSegment-NNNN",
			expected == result,
			"TraceSegment::getCopy() returns a copy of the prefix of the stored segment if "
			"0 < prefix < segment.size()");
	}

	// segment.size() = 3 . prefix = 1 ( 0 < prefix < segment.size() case) 
	{
		vector<int> seg{ 0, 1, 2 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 1 };
		vector<int> result = trcSeg.getCopy();
		vector<int> expected{ 0 };
		fsmlib_assert("TC-TraceSegment-NNNN",
			expected == result,
			"TraceSegment::getCopy() returns a copy of the prefix of the stored segment if "
			"0 < prefix < segment.size()");
	}

	// segment.size() = 3 . prefix = 2 ( 0 < prefix < segment.size() case) 
	{
		vector<int> seg{ 0, 1, 2 };
		std::shared_ptr< std::vector<int> > segment = make_shared<vector<int>>(seg);
		TraceSegment trcSeg{ segment, 2 };
		vector<int> result = trcSeg.getCopy();
		vector<int> expected{ 0, 1 };
		fsmlib_assert("TC-TraceSegment-NNNN",
			expected == result,
			"TraceSegment::getCopy() returns a copy of the prefix of the stored segment if "
			"0 < prefix < segment.size()");
	}
}

//===================================== SegmentedTrace Tests ===================================================

// tests operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2)
// positive case
void testSegmentedTraceEqualOperatorPositive() {
	// both SegmentedTraces are empty.
	{
		deque< std::shared_ptr<TraceSegment> > segments;
		SegmentedTrace segTrc1(segments);
		SegmentedTrace segTrc2(segments);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are empty.");
	}

	// segTrc1 = <<1>> , segTrc2 = <<> <1>>
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1 = make_shared<TraceSegment>(v1);
		segments1.push_back(trcSeg1);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1> <>>, segTrc2 = <<1>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1>>, segTrc2 = <<1>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1,2>,<3>>, segTrc2 = <<1>,<2>,<3>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 2});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		shared_ptr<vector<int>> v2_3 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg2_3 = make_shared<TraceSegment>(v2_3);
		segments2.push_back(trcSeg2_3);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1,2>,<3>,<45>>, segTrc2 = <<1>,<2,3,4>,<5>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 2});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		shared_ptr<vector<int>> v1_3 = make_shared<vector<int>>(vector<int>{4,5});
		shared_ptr<TraceSegment> trcSeg1_3 = make_shared<TraceSegment>(v1_3);
		segments1.push_back(trcSeg1_3);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2,3,4});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		shared_ptr<vector<int>> v2_3 = make_shared<vector<int>>(vector<int>{5});
		shared_ptr<TraceSegment> trcSeg2_3 = make_shared<TraceSegment>(v2_3);
		segments2.push_back(trcSeg2_3);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1,3>> (prefixTrcSeg1 = 1), segTrc2 = <<1>> (prefix = string::npos)
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 3});
		int prefixTrcSeg1 = 1;
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1, prefixTrcSeg1);
		segments1.push_back(trcSeg1_1);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");
	}

	// segTrc1 = <<1,2,3>,<3>> (prefixTrcSeg1_1 = 2, prefixTrcSeg1_2 = 1), segTrc2 = <<1,3>, <2,4>, <4> , <3>> 
	// (prefixTrcSeg2_1 = 1, prefixTrcSeg2_2 = 1, prefixTrcSeg2_3 = 0, prefixTrcSeg2_4 = 1)
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 2, 3});
		int prefixTrcSeg1_1 = 2;
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1, prefixTrcSeg1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3});
		int prefixTrcSeg1_2 = 1;
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2, prefixTrcSeg1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1,3});
		int prefixTrcSeg2_1 = 1;
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1, prefixTrcSeg2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2, 4});
		int prefixTrcSeg2_2 = 1;
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2, prefixTrcSeg2_2);
		segments2.push_back(trcSeg2_2);
		shared_ptr<vector<int>> v2_3 = make_shared<vector<int>>(vector<int>{4});
		int prefixTrcSeg2_3 = 0;
		shared_ptr<TraceSegment> trcSeg2_3 = make_shared<TraceSegment>(v2_3, prefixTrcSeg2_3);
		segments2.push_back(trcSeg2_3);
		shared_ptr<vector<int>> v2_4 = make_shared<vector<int>>(vector<int>{3});
		int prefixTrcSeg2_4 = 1;
		shared_ptr<TraceSegment> trcSeg2_4 = make_shared<TraceSegment>(v2_4, prefixTrcSeg2_4);
		segments2.push_back(trcSeg2_4);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			segTrc1 == segTrc2,
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns true if both traces are equal.");


		std::cout << "copy equals..." << std::endl;
		std::cout << (segTrc1.getCopy() == segTrc2.getCopy()) << std::endl;
	}

}

// tests operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2)
// negative case
void testSegmentedTraceEqualOperatorNegative() {
	// segTrc1 = <<>> , segTrc2 = <<1>>
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1 = make_shared<vector<int>>(vector<int>{});
		shared_ptr<TraceSegment> trcSeg1 = make_shared<TraceSegment>(v1);
		segments1.push_back(trcSeg1);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1>, <2,3>>, segTrc2 = <<1>,<2>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{2,3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1>, <2,3>>, segTrc2 = <<1>,<2,4>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{2, 3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2,4});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1>, <2>>, segTrc2 = <<2>,<1>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{2});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{2});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1>> prefixTrcSeg1_1 = 0, segTrc2 = <<1>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});
		int prefixTrcSeg1_1 = 0;
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1, prefixTrcSeg1_1);
		segments1.push_back(trcSeg1_1);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1>,<2>> prefixTrcSeg1_2 = 0, segTrc2 = <<1>,<2>>
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1});		
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{2});
		int prefixTrcSeg1_2 = 0;
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2, prefixTrcSeg1_2);
		segments1.push_back(trcSeg1_2);		
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1});
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{2});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1,2>,<3>>, segTrc2 = <<1,2>,<3>> prefixTrcSeg2_1 = 1
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1,2});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1,2});
		int prefixTrcSeg2_1 = 1;
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1, prefixTrcSeg2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1,2>,<3>>, segTrc2 = <<1,2>,<3>> prefixTrcSeg2_1 = 0
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 2});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2);
		segments1.push_back(trcSeg1_2);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1, 2});
		int prefixTrcSeg2_1 = 0;
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1, prefixTrcSeg2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{3});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}

	// segTrc1 = <<1,2>,<3,4,5>,<6>> prefixTrcSeg1_2 = 2, segTrc2 = <<1,2>,<3,4>,<5,6>> 
	{
		deque< std::shared_ptr<TraceSegment> > segments1;
		shared_ptr<vector<int>> v1_1 = make_shared<vector<int>>(vector<int>{1, 2});
		shared_ptr<TraceSegment> trcSeg1_1 = make_shared<TraceSegment>(v1_1);
		segments1.push_back(trcSeg1_1);
		shared_ptr<vector<int>> v1_2 = make_shared<vector<int>>(vector<int>{3,4,5});
		int prefixTrcSeg1_2 = 2;
		shared_ptr<TraceSegment> trcSeg1_2 = make_shared<TraceSegment>(v1_2, prefixTrcSeg1_2);
		segments1.push_back(trcSeg1_2);
		shared_ptr<vector<int>> v1_3 = make_shared<vector<int>>(vector<int>{6});
		shared_ptr<TraceSegment> trcSeg1_3 = make_shared<TraceSegment>(v1_3);
		segments1.push_back(trcSeg1_3);
		SegmentedTrace segTrc1(segments1);

		deque< std::shared_ptr<TraceSegment> > segments2;
		shared_ptr<vector<int>> v2_1 = make_shared<vector<int>>(vector<int>{1, 2});		
		shared_ptr<TraceSegment> trcSeg2_1 = make_shared<TraceSegment>(v2_1);
		segments2.push_back(trcSeg2_1);
		shared_ptr<vector<int>> v2_2 = make_shared<vector<int>>(vector<int>{3,4});
		shared_ptr<TraceSegment> trcSeg2_2 = make_shared<TraceSegment>(v2_2);
		segments2.push_back(trcSeg2_2);
		shared_ptr<vector<int>> v2_3 = make_shared<vector<int>>(vector<int>{5,6});
		shared_ptr<TraceSegment> trcSeg2_3 = make_shared<TraceSegment>(v2_3);
		segments2.push_back(trcSeg2_3);
		SegmentedTrace segTrc2(segments2);

		fsmlib_assert("TC-SegmentedTrace-NNNN",
			not(segTrc1 == segTrc2),
			"operator==(SegmentedTrace const & trace1, SegmentedTrace const & trace2) returns false if both traces are unequal.");
	}
}

//===================================== HittingSet Tests ===================================================

// tests HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s)
void testHittingSetConstructor() {
	// empty vector (no set contained)
	{
		vector<unordered_set<int>> s;
		HittingSet hs(s);
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest.empty(),
			"HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s) called with "
			"empty s initializes HsTreeNode::hSmallest as empty set");
	}

	// vector contains only the empty set.
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{});
		HittingSet hs(s);
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest.empty(),
			"HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s) initializes "
			"HsTreeNode::hSmallest as empty set, if s contains only empty sets.");
	}

	// vector contains one non-empty set.
	{
		vector<unordered_set<int>> s;
		unordered_set<int> set1{ 1, 2 };
		s.push_back(set1);
		HittingSet hs(s);
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest.size() == 2
			&& HsTreeNode::hSmallest.count(1) != 0
			&& HsTreeNode::hSmallest.count(2) != 0,
			"HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s) initializes  "
			"HsTreeNode::hSmallest as a set containing all elements contained in sets of s.");
	}

	// vector contains two non-empty sets. The intersection of both sets isn't empty.
	{
		vector<unordered_set<int>> s;
		unordered_set<int> set1{ 1, 2 };
		unordered_set<int> set2{ 2, 3 };
		s.push_back(set1);
		s.push_back(set2);
		HittingSet hs(s);
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest.size() == 3
			&& HsTreeNode::hSmallest.count(1) != 0
			&& HsTreeNode::hSmallest.count(2) != 0
			&& HsTreeNode::hSmallest.count(3) != 0,
			"HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s) initializes  "
			"HsTreeNode::hSmallest as a set containing all elements contained in sets of s.");
	}
}

// tests HittingSet::calcMinCardHittingSet()
void testHittingSetCalcMinCardHittingSet() {
	// empty vector (no set contained)
	{
		vector<unordered_set<int>> s;
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 0;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		// this returns true (ok because s is empty, so one could argue that every possible set is a hitting set )
		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains only the empty set.
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 0;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		// this returns false (ok because no element of minHittingSet is contained in the only set of s in this special case)
		//HsTreeNode node(minHittingSet, s);
		//fsmlib_assert("TC-HittingSet-NNNN",
		//	node.isHittingSet(),
		//	"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains {1,2} => expectedSize = 1
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{1,2});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 1;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains {1,2} and {3} => expectedSize = 2
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{1, 2});
		s.push_back(unordered_set<int>{3});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 2;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains {1,2} and {2,3} => expectedSize = 1
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{1, 2});
		s.push_back(unordered_set<int>{2, 3});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 1;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains {1,2}, {2,3} and {3,1} => expectedSize = 2
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{1, 2});
		s.push_back(unordered_set<int>{2, 3});
		s.push_back(unordered_set<int>{3, 1});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 2;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}

	// vector contains {1,2,3}, {2,3} and {3,1} => expectedSize = 1
	{
		vector<unordered_set<int>> s;
		s.push_back(unordered_set<int>{1, 2, 3});
		s.push_back(unordered_set<int>{2, 3});
		s.push_back(unordered_set<int>{3, 1});
		HittingSet hs(s);
		unordered_set<int> minHittingSet = hs.calcMinCardHittingSet();
		fsmlib_assert("TC-HittingSet-NNNN",
			HsTreeNode::hSmallest == minHittingSet,
			"The result of HittingSet::calcMinCardHittingSet() is assigned to HsTreeNode::hSmallest.");

		int expectedSize = 1;
		fsmlib_assert("TC-HittingSet-NNNN",
			minHittingSet.size() == expectedSize,
			"Result of HittingSet::calcMinCardHittingSet() is minimal.");

		HsTreeNode node(minHittingSet, s);
		fsmlib_assert("TC-HittingSet-NNNN",
			node.isHittingSet(),
			"Result of HittingSet::calcMinCardHittingSet() is a hitting set.");
	}
}

//===================================== HsTreeNode Tests ===================================================

// tests HsTreeNode::isHittingSet() 
// Positive case
void testHsTreeNodeIsHittingSetPositive() {
	// s = <> , x = {}
	{
		unordered_set<int> x{};
		vector<unordered_set<int>> s{};
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns true if x is a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}> , x = {2}
	{
		unordered_set<int> x{2};
		vector<unordered_set<int>> s{ unordered_set<int>{1,2} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns true if x is a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}, {1,3}> , x = {1}
	{
		unordered_set<int> x{ 1 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns true if x is a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}, {1,3}> , x = {2,3,4}
	{
		unordered_set<int> x{ 2,3,4 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns true if x is a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}
}

// tests HsTreeNode::isHittingSet() 
// Negative case
void testHsTreeNodeIsHittingSetNegative() {
	// s = <{}> , x = {}
	{
		unordered_set<int> x{};
		vector<unordered_set<int>> s{ unordered_set<int>{} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1}> , x = {}
	{
		unordered_set<int> x{};
		vector<unordered_set<int>> s{ unordered_set<int>{1} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}> , x = {3}
	{
		unordered_set<int> x{3};
		vector<unordered_set<int>> s{ unordered_set<int>{1,2} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}, {1,3}> , x = {2}
	{
		unordered_set<int> x{ 2 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}, {1,3}> , x = {3,4}
	{
		unordered_set<int> x{ 3,4 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}

	// s = <{1,2}, {1,3}, {3,4}> , x = {2,4}
	{
		unordered_set<int> x{ 2,4 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3}, unordered_set<int>{3,4} };
		HsTreeNode n(x, s);
		int oldMaxNumNode = HsTreeNode::maxNodeNum;

		fsmlib_assert("TC-HsTreeNode-NNNN",
			not n.isHittingSet(),
			"HsTreeNode::isHittingSet() returns false if x is not a hitting set of s.");

		fsmlib_assert("TC-HsTreeNode-NNNN",
			oldMaxNumNode == HsTreeNode::maxNodeNum,
			"HsTreeNode::isHittingSet() doesn't change HsTreeNode::maxNodeNum.");
	}
}

// tests HsTreeNode::expandNode()
void testHsTreeNodeExpandNode() {
	// s = <{3},{2}>, x = {1,2} -> x isn't a Hitting Set of s.
	{
		unordered_set<int> x = { 1,2 };
		vector<unordered_set<int>> s{ unordered_set<int>{3}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode n(x, s);
		n.expandNode();

		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest == x,
			"HsTreeNode::expandNode() doesn't change HsTreeNode::hSmallest if x isn't a hitting set of s.");

	}

	// s = <>, x = {1,2} -> s is empty so no elements are needed in x. (A created hitting set can't be empty, so the expected size is 1.)
	{
		unordered_set<int> x = { 1,2 };
		vector<unordered_set<int>> s{ };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();

		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest.size() == 1,
			"HsTreeNode::expandNode() calculates all possible hitting sets of s, consisting of elements of x, "
			"and stores the minimal one in HsTreeNode::hSmallest.");
	}

	// s = <{1},{2}>, x = {1,2} -> No child is a hitting set of s. (hitting set  = {1,2})
	{
		unordered_set<int> x = { 1,2 };
		vector<unordered_set<int>> s{ unordered_set<int>{1}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();

		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest == x,
			"HsTreeNode::expandNode() calculates all possible hitting sets of s, consisting of elements of x, "
			"and stores the minimal one in HsTreeNode::hSmallest.");
	}

	// s = <{1,2},{2,3}>, x = {1,2,3} -> Each direct child of n is a hitting set. Minimal hitting set contains 1 element.
	{
		unordered_set<int> x = { 1,2,3 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{2,3} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();

		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest.size() == 1,
			"HsTreeNode::expandNode() calculates all possible hitting sets of s, consisting of elements of x, "
			"and stores the minimal one in HsTreeNode::hSmallest.");
	}

	// s = <{1,2},{1,3}>, x = {1} -> Each direct child of n is empty (no child is hitting set). Minimal hitting set contains 1 element.
	{
		unordered_set<int> x = { 1 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{1,3} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();

		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest.size() == 1,
			"HsTreeNode::expandNode() calculates all possible hitting sets of s, consisting of elements of x, "
			"and stores the minimal one in HsTreeNode::hSmallest.");
	}

	// s = <{4},{1,4},{2,4}>, x = {1,2,4} -> Not every child of n is a hitting set (but some). 
	// Minimal hitting set contains 1 element.
	{
		unordered_set<int> x = { 1,2,4 };
		vector<unordered_set<int>> s{ unordered_set<int>{4}, unordered_set<int>{1,4}, unordered_set<int>{2,4} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();
		
		fsmlib_assert("TC-HsTreeNode-NNNN",
			HsTreeNode::hSmallest.size() == 1,
			"HsTreeNode::expandNode() calculates all possible hitting sets of s, consisting of elements of x, "
			"and stores the minimal one in HsTreeNode::hSmallest.");
	}
}

// tests HsTreeNode::toDot()
void testHsTreeNodeToDot() {
	// no children (expandNode() is not invoked). x = {}
	{
		unordered_set<int> x = {  };
		vector<unordered_set<int>> s{ unordered_set<int>{1}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		cout << n.toDot() << endl;
		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.toDot().find("label=\"[]\"") != string::npos,
			"HsTreeNode::toDot() returns string that contains all expected node labels.");
	}

	// no children (expandNode() is not invoked). x = {1}
	{
		unordered_set<int> x = {1};
		vector<unordered_set<int>> s{ unordered_set<int>{1}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		cout << n.toDot() << endl;
		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.toDot().find("label=\"[1]\"") != string::npos,
			"HsTreeNode::toDot() returns string that contains all expected node labels.");
	}

	// no children (expandNode() is not invoked). x = {1,2}
	{
		unordered_set<int> x = { 1,2 };
		vector<unordered_set<int>> s{ unordered_set<int>{1}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		cout << n.toDot() << endl;
		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.toDot().find("label=\"[1, 2]\"") != string::npos,
			"HsTreeNode::toDot() returns string that contains all expected node labels.");
	}

	// with children (expandNode() is invoked). x = {1,2} s = <{1,2},{2}>
	{
		unordered_set<int> x = { 1,2 };
		vector<unordered_set<int>> s{ unordered_set<int>{1,2}, unordered_set<int>{2} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();
		cout << n.toDot() << endl;
		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.toDot().find("label=\"[1, 2]\"") != string::npos
			&& n.toDot().find("label=\"[2]\"") != string::npos,
			"HsTreeNode::toDot() returns string that contains all expected node labels.");
	}

	// with children (expandNode() is invoked). x = {1,2,4} s = <{4},{1,4},{2,4}>
	{
		unordered_set<int> x = { 1,2,4 };
		vector<unordered_set<int>> s{ unordered_set<int>{4}, unordered_set<int>{1,4}, unordered_set<int>{2,4} };
		HsTreeNode::hSmallest = x;
		HsTreeNode::maxNodeNum = 0;
		HsTreeNode n(x, s);
		n.expandNode();
		cout << n.toDot() << endl;
		fsmlib_assert("TC-HsTreeNode-NNNN",
			n.toDot().find("label=\"[1, 2, 4]\"") != string::npos
			&& n.toDot().find("label=\"[2, 4]\"") != string::npos
			&& n.toDot().find("label=\"[4]\"") != string::npos
			&& n.toDot().find("label=\"[1, 4]\"") != string::npos,
			"HsTreeNode::toDot() returns string that contains all expected node labels.");
	}

}

//===================================== Int2IntMap Tests ===================================================

// tests Constructor of Int2IntMap (Int2IntMap::Int2IntMap(const int maxInput))
void testInt2IntMapConstructor() {
	{
		Int2IntMap map(0);
		fsmlib_assert("TC-Int2IntMap-NNNN",
			map.size() == 1
			&& map.at(0) == -1,
			"Int2IntMap::Int2IntMap(const int maxInput) creates map containing maxInput + 1 keys mapped to -1");
	}
	{
		Int2IntMap map(1);
		fsmlib_assert("TC-Int2IntMap-NNNN",
			map.size() == 2
			&& map.at(0) == -1
			&& map.at(1) == -1,
			"Int2IntMap::Int2IntMap(const int maxInput) creates map containing maxInput + 1 keys mapped to -1");
	}
	{
		Int2IntMap map(2);
		fsmlib_assert("TC-Int2IntMap-NNNN",
			map.size() == 3
			&& map.at(0) == -1
			&& map.at(1) == -1
			&& map.at(2) == -1,
			"Int2IntMap::Int2IntMap(const int maxInput) creates map containing maxInput + 1 keys mapped to -1");
	}
}

//===================================== PkTableRow Tests ===================================================

// tests PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c)
// Positive Case
void testPkTableRowIsEquivalentPositive(){
	// ioMap1 = [0], i2pMap1 = [0], ioMap2 = [1], i2pMap2 = [0], s2cMap = [0]
	{
		int maxInput = 0;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 0;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0], i2pMap1 = [-1],
	// ioMap2 = [1], i2pMap2 = [-1], 
	// s2cMap = [0]
	{
		int maxInput = 0;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = -1;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = -1;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 0;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0], i2pMap1 = [0],
	// ioMap2 = [1], i2pMap2 = [0], 
	// s2cMap = [1]
	{
		int maxInput = 0;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0,1], i2pMap1 = [1,0],
	// ioMap2 = [0,1], i2pMap2 = [1,0], 
	// s2cMap = [1,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 1;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 1;
		i2pMap1[1] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 1;
		i2pMap2[1] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 1;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0,1], i2pMap1 = [-1,0],
	// ioMap2 = [0,1], i2pMap2 = [-1,0], 
	// s2cMap = [1,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 1;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = -1;
		i2pMap1[1] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = -1;
		i2pMap2[1] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 1;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0,1], i2pMap1 = [1,0],
	// ioMap2 = [1,1], i2pMap2 = [0,0], 
	// s2cMap = [0,0,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 1;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 1;
		i2pMap1[1] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		ioMap2[1] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		i2pMap2[1] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput+1);
		s2cMap[0] = 0;
		s2cMap[1] = 0;
		s2cMap[2] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [0,1,1], i2pMap1 = [1,2,2],
	// ioMap2 = [1,1,1], i2pMap2 = [0,2,2], 
	// s2cMap = [1,1,0]
	{
		int maxInput = 2;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 1;
		ioMap1[2] = 1;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 1;
		i2pMap1[1] = 2;
		i2pMap1[2] = 2;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		ioMap2[1] = 1;
		ioMap2[2] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		i2pMap2[1] = 2;
		i2pMap2[2] = 2;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 1;
		s2cMap[1] = 1;
		s2cMap[2] = 0;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}

	// ioMap1 = [1,0,1], i2pMap1 = [1,3,4],
	// ioMap2 = [1,0,1], i2pMap2 = [0,2,2], 
	// s2cMap = [1,1,2,2,2]
	{
		int maxInput = 2;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 1;
		ioMap1[1] = 0;
		ioMap1[2] = 1;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 1;
		i2pMap1[1] = 3;
		i2pMap1[2] = 4;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 1;
		ioMap2[1] = 0;
		ioMap2[2] = 1;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		i2pMap2[1] = 2;
		i2pMap2[2] = 2;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput + 2);
		s2cMap[0] = 1;
		s2cMap[1] = 1;
		s2cMap[2] = 2;
		s2cMap[3] = 2;
		s2cMap[4] = 2;

		fsmlib_assert("TC-PkTableRow-NNNN",
			row1.isEquivalent(row2, s2cMap)
			&& row2.isEquivalent(row1, s2cMap),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns true if both rows are equivalent wrt. s2c");
	}
}

// tests PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c)
// Negative Case
void testPkTableRowIsEquivalentNegative() {
	// ioMap1 = [0], i2pMap1 = [0],
	// ioMap2 = [0], i2pMap2 = [1], 
	// s2cMap = [0,1]
	{
		int maxInput = 0;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 1;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput + 1);
		s2cMap[0] = 0;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}

	// ioMap1 = [0,0], i2pMap1 = [0,1],
	// ioMap2 = [0,0], i2pMap2 = [1,0], 
	// s2cMap = [0,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		i2pMap1[1] = 1;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 1;
		i2pMap2[1] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 0;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}

	// ioMap1 = [0], i2pMap1 = [0],
	// ioMap2 = [0], i2pMap2 = [1], 
	// s2cMap = [1,0]
	{
		int maxInput = 0;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 1;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput + 1);
		s2cMap[0] = 1;
		s2cMap[1] = 0;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}

	// ioMap1 = [0,0,0], i2pMap1 = [0,1,2],
	// ioMap2 = [0,0,0], i2pMap2 = [0,1,1], 
	// s2cMap = [1,2,1]
	{
		int maxInput = 2;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 0;
		ioMap1[2] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		i2pMap1[1] = 1;
		i2pMap1[2] = 2;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 0;
		ioMap2[2] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		i2pMap2[1] = 1;
		i2pMap2[2] = 1;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 1;
		s2cMap[1] = 2;
		s2cMap[2] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}

	// ioMap1 = [0,0], i2pMap1 = [0,-1],
	// ioMap2 = [0,0], i2pMap2 = [0,1], 
	// s2cMap = [0,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		i2pMap1[1] = -1;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = 0;
		i2pMap2[1] = 1;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 0;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}

	// ioMap1 = [0,0], i2pMap1 = [0,1],
	// ioMap2 = [0,0], i2pMap2 = [-1,0], 
	// s2cMap = [0,1]
	{
		int maxInput = 1;

		IOMap ioMap1(maxInput);
		ioMap1[0] = 0;
		ioMap1[1] = 0;
		I2PMap i2pMap1(maxInput);
		i2pMap1[0] = 0;
		i2pMap1[1] = 1;
		PkTableRow row1(ioMap1, i2pMap1);

		IOMap ioMap2(maxInput);
		ioMap2[0] = 0;
		ioMap2[1] = 0;
		I2PMap i2pMap2(maxInput);
		i2pMap2[0] = -1;
		i2pMap2[1] = 0;
		PkTableRow row2(ioMap2, i2pMap2);

		S2CMap s2cMap(maxInput);
		s2cMap[0] = 0;
		s2cMap[1] = 1;

		fsmlib_assert("TC-PkTableRow-NNNN",
			not (row1.isEquivalent(row2, s2cMap))
			&& not (row2.isEquivalent(row1, s2cMap)),
			"PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c) returns false if both rows aren't equivalent wrt. s2c");
	}
}

//===================================== PkTable Tests ===================================================

// tests PkTable::maxClassId()
void testPkTableMaxClassId() {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	int maxInput = 4;
	// s2c = [-1] -> all s2c ids are smaller than 0
	{
		int numStates = 1;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, -1);
		fsmlib_assert("TC-PkTable-NNNN",
			pkTable.maxClassId() == 0,
			"PkTable::maxClassId() returns 0 if all elements from s2c are smaller than 0");
	}

	// s2c = [-1,0] -> not all s2c ids are smaller than 0
	{
		int numStates = 2;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, -1);
		pkTable.setClass(1, 0);
		fsmlib_assert("TC-PkTable-NNNN",
			pkTable.maxClassId() == 0,
			"PkTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,0,1]
	{
		int numStates = 3;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 0);
		pkTable.setClass(1, 0);
		pkTable.setClass(2, 1);
		fsmlib_assert("TC-PkTable-NNNN",
			pkTable.maxClassId() == 1,
			"PkTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,1,2,1]
	{
		int numStates = 4;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 0);
		pkTable.setClass(1, 1);
		pkTable.setClass(2, 2);
		pkTable.setClass(3, 1);
		fsmlib_assert("TC-PkTable-NNNN",
			pkTable.maxClassId() == 2,
			"PkTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,1,3,3]
	{
		int numStates = 4;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 0);
		pkTable.setClass(1, 1);
		pkTable.setClass(2, 3);
		pkTable.setClass(3, 3);
		fsmlib_assert("TC-PkTable-NNNN",
			pkTable.maxClassId() == 3,
			"PkTable::maxClassId() returns greatest value of s2c");
	}

}

// Applies inVec on row. Creates and returns pointer to vector representing the output sequence.
shared_ptr<vector<int>> applyToPkTable(const vector<int> &inVec, shared_ptr<PkTableRow> row, PkTable &pkTable) {
	vector<int> outVec;
	shared_ptr<PkTableRow> currentRow = row;
	for (const auto &i : inVec) {
		if (currentRow->get(i) != -1) {
			outVec.push_back(currentRow->getIOMap().at(i));
			currentRow = pkTable.getRow(currentRow->get(i));
		}
		else {
			break;
		}
	}
	return make_shared<vector<int>>(outVec);
}

// Checks if states represented by row1 and row2 are k-equivalent.
// Applies each input vector from iolc at both states (row1 and row2).
// Both states are k-equivalent if each input vector results in the same output vector.  
bool checkPkEqualityOfRows(PkTable &pkTable, shared_ptr<PkTableRow> row1, shared_ptr<PkTableRow> row2,
							IOListContainer &iolc, int k) {
	for (const vector<int> &inVec : *(iolc.getIOLists())) {
		shared_ptr<vector<int>> outVec1 = applyToPkTable(inVec, row1, pkTable);
		shared_ptr<vector<int>> outVec2 = applyToPkTable(inVec, row2, pkTable);
		if (*outVec1 != *outVec2) {
			return false;
		}
	}
	return true;
}

// checks if all rows from pkTable with the same class in s2c are k-equivalent
bool checkPkEqualityOfClasses(PkTable &pkTable, IOListContainer &iolc, int k, int numStates) {
	for (int i = 0; i < numStates; ++i) {
		for (int j = i + 1; j < numStates; ++j) {
			if (pkTable.getClass(i) == pkTable.getClass(j)) {
				if (not checkPkEqualityOfRows(pkTable, pkTable.getRow(i), pkTable.getRow(j), iolc, k)) {
					return false;
				}
			}
		}
	}
	return true;
}

// Checks if all rows from pkTable with different classes in s2c are not k-equivalent
bool checkPkUnequalityOfClasses(PkTable &pkTable, IOListContainer &iolc, int k, int numStates) {
	for (int i = 0; i < numStates; ++i) {
		for (int j = i + 1; j < numStates; ++j) {
			if (pkTable.getClass(i) != pkTable.getClass(j)) {
				if (checkPkEqualityOfRows(pkTable, pkTable.getRow(i), pkTable.getRow(j), iolc, k)) {
					return false;
				}
			}
		}
	}
	return true;
}

struct PkTableTestCase
{
	shared_ptr<FsmPresentationLayer> presentationLayer;
	shared_ptr<PkTable> pkTable;
	int maxInput;
	int maxOutput;
	int numStates;
};

PkTableTestCase getPkTableTestCase1() {
	std::shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
	int maxInput = 2;
	int maxOutput = 1;

	// create all PkTableRows
	// (maps are declared static because their references are stored in the created rows)
	static IOMap ioMap1(maxInput);
	ioMap1[0] = 1;
	ioMap1[1] = 0;
	ioMap1[2] = 0;
	static I2PMap i2pMap1(maxInput);
	i2pMap1[0] = 1;
	i2pMap1[1] = 1;
	i2pMap1[2] = 4;
	PkTableRow row1(ioMap1, i2pMap1);

	static IOMap ioMap2(maxInput);
	ioMap2[0] = 0;
	ioMap2[1] = 1;
	ioMap2[2] = 1;
	static I2PMap i2pMap2(maxInput);
	i2pMap2[0] = 0;
	i2pMap2[1] = 3;
	i2pMap2[2] = 3;
	PkTableRow row2(ioMap2, i2pMap2);

	static IOMap ioMap3(maxInput);
	ioMap3[0] = 1;
	ioMap3[1] = 0;
	ioMap3[2] = 0;
	static I2PMap i2pMap3(maxInput);
	i2pMap3[0] = 1;
	i2pMap3[1] = 1;
	i2pMap3[2] = 4;
	PkTableRow row3(ioMap3, i2pMap3);

	static IOMap ioMap4(maxInput);
	ioMap4[0] = 0;
	ioMap4[1] = 1;
	ioMap4[2] = 1;
	static I2PMap i2pMap4(maxInput);
	i2pMap4[0] = 2;
	i2pMap4[1] = 1;
	i2pMap4[2] = 1;
	PkTableRow row4(ioMap4, i2pMap4);

	static IOMap ioMap5(maxInput);
	ioMap5[0] = 1;
	ioMap5[1] = 0;
	ioMap5[2] = 0;
	static I2PMap i2pMap5(maxInput);
	i2pMap5[0] = 5;
	i2pMap5[1] = 3;
	i2pMap5[2] = 2;
	PkTableRow row5(ioMap5, i2pMap5);

	static IOMap ioMap6(maxInput);
	ioMap6[0] = 0;
	ioMap6[1] = 1;
	ioMap6[2] = 1;
	static I2PMap i2pMap6(maxInput);
	i2pMap6[0] = 7;
	i2pMap6[1] = 8;
	i2pMap6[2] = 5;
	PkTableRow row6(ioMap6, i2pMap6);

	static IOMap ioMap7(maxInput);
	ioMap7[0] = 1;
	ioMap7[1] = 0;
	ioMap7[2] = 0;
	static I2PMap i2pMap7(maxInput);
	i2pMap7[0] = 5;
	i2pMap7[1] = 1;
	i2pMap7[2] = 7;
	PkTableRow row7(ioMap7, i2pMap7);

	static IOMap ioMap8(maxInput);
	ioMap8[0] = 1;
	ioMap8[1] = 0;
	ioMap8[2] = 0;
	static I2PMap i2pMap8(maxInput);
	i2pMap8[0] = 3;
	i2pMap8[1] = 3;
	i2pMap8[2] = 6;
	PkTableRow row8(ioMap8, i2pMap8);

	static IOMap ioMap9(maxInput);
	ioMap9[0] = 0;
	ioMap9[1] = 1;
	ioMap9[2] = 1;
	static I2PMap i2pMap9(maxInput);
	i2pMap9[0] = 6;
	i2pMap9[1] = 8;
	i2pMap9[2] = 6;
	PkTableRow row9(ioMap9, i2pMap9);

	std::vector<std::shared_ptr<PkTableRow>> rows;
	rows.push_back(make_shared<PkTableRow>(row1));
	rows.push_back(make_shared<PkTableRow>(row2));
	rows.push_back(make_shared<PkTableRow>(row3));
	rows.push_back(make_shared<PkTableRow>(row4));
	rows.push_back(make_shared<PkTableRow>(row5));
	rows.push_back(make_shared<PkTableRow>(row6));
	rows.push_back(make_shared<PkTableRow>(row7));
	rows.push_back(make_shared<PkTableRow>(row8));
	rows.push_back(make_shared<PkTableRow>(row9));

	// create PkTable from PkTableRows
	int numStates = 9;
	PkTable pkTable(numStates, maxInput, rows, presentationLayer);

	// set classes. pkTable becomes P1 Table.
	pkTable.setClass(0, 0);
	pkTable.setClass(1, 1);
	pkTable.setClass(2, 0);
	pkTable.setClass(3, 1);
	pkTable.setClass(4, 0);
	pkTable.setClass(5, 1);
	pkTable.setClass(6, 0);
	pkTable.setClass(7, 0);
	pkTable.setClass(8, 1);

	PkTableTestCase testStructure;
	testStructure.presentationLayer = presentationLayer;
	testStructure.pkTable = make_shared<PkTable>(pkTable);
	testStructure.maxInput = maxInput;
	testStructure.maxOutput = maxOutput;
	testStructure.numStates = numStates;

	return testStructure;
}

PkTableTestCase getPkTableTestCase2() {
	std::shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
	int maxInput = 2;
	int maxOutput = 1;

	// create all PkTableRows
	// (maps are declared static because their references are stored in the created rows)
	static IOMap ioMap1(maxInput);
	ioMap1[0] = 0;
	ioMap1[1] = 0;
	static I2PMap i2pMap1(maxInput);
	i2pMap1[0] = 0;
	i2pMap1[1] = 1;
	PkTableRow row1(ioMap1, i2pMap1);

	static IOMap ioMap2(maxInput);
	ioMap2[2] = 1;
	static I2PMap i2pMap2(maxInput);
	i2pMap2[2] = 2;
	PkTableRow row2(ioMap2, i2pMap2);

	static IOMap ioMap3(maxInput);
	ioMap3[0] = 0;
	ioMap3[1] = 0;
	static I2PMap i2pMap3(maxInput);
	i2pMap3[0] = 3;
	i2pMap3[1] = 1;
	PkTableRow row3(ioMap3, i2pMap3);

	static IOMap ioMap4(maxInput);
	ioMap4[0] = 0;
	ioMap4[1] = 0;
	ioMap4[2] = 0;
	static I2PMap i2pMap4(maxInput);
	i2pMap4[0] = 3;
	i2pMap4[1] = 4;
	i2pMap4[2] = 5;
	PkTableRow row4(ioMap4, i2pMap4);

	static IOMap ioMap5(maxInput);
	ioMap5[0] = 0;
	ioMap5[1] = 0;
	static I2PMap i2pMap5(maxInput);
	i2pMap5[0] = 0;
	i2pMap5[1] = 2;
	PkTableRow row5(ioMap5, i2pMap5);

	static IOMap ioMap6(maxInput);
	ioMap6[2] = 1;
	static I2PMap i2pMap6(maxInput);
	i2pMap6[2] = 2;
	PkTableRow row6(ioMap6, i2pMap6);

	std::vector<std::shared_ptr<PkTableRow>> rows;
	rows.push_back(make_shared<PkTableRow>(row1));
	rows.push_back(make_shared<PkTableRow>(row2));
	rows.push_back(make_shared<PkTableRow>(row3));
	rows.push_back(make_shared<PkTableRow>(row4));
	rows.push_back(make_shared<PkTableRow>(row5));
	rows.push_back(make_shared<PkTableRow>(row6));

	// create PkTable from PkTableRows
	int numStates = rows.size();
	PkTable pkTable(numStates, maxInput, rows, presentationLayer);

	// set classes. pkTable becomes P1 Table.
	pkTable.setClass(0, 0);
	pkTable.setClass(1, 1);
	pkTable.setClass(2, 0);
	pkTable.setClass(3, 2);
	pkTable.setClass(4, 0);
	pkTable.setClass(5, 1);

	PkTableTestCase testStructure;
	testStructure.presentationLayer = presentationLayer;
	testStructure.pkTable = make_shared<PkTable>(pkTable);
	testStructure.maxInput = maxInput;
	testStructure.maxOutput = maxOutput;
	testStructure.numStates = numStates;

	return testStructure;
}

// tests PkTable::getPkPlusOneTable()
void testPkTableGetPkPlusOneTable() {

	{
		PkTableTestCase pkTableTestCase = getPkTableTestCase1();
		int k = 1;
		shared_ptr<PkTable> currentTable = pkTableTestCase.pkTable;
		do {
			IOListContainer iolc(pkTableTestCase.maxInput, k, k, pkTableTestCase.presentationLayer);
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkEqualityOfClasses(*currentTable, iolc, k, pkTableTestCase.numStates),
				"PkTableRows with the same class in s2c are k equivalent");
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkUnequalityOfClasses(*currentTable, iolc, k, pkTableTestCase.numStates),
				"PkTableRows with the different classes in s2c are not k equivalent");

			currentTable = currentTable->getPkPlusOneTable();
			++k;
		} while (currentTable != nullptr);

		fsmlib_assert("TC-PkTable-NNNN",
			k == 5,
			"PkTable::getPkPlusOneTable() returns nullptr if no new equivalence classes are generated.");
	}
 
	{
		PkTableTestCase pkTableTestCase = getPkTableTestCase2();
		int k = 1;
		shared_ptr<PkTable> currentTable = pkTableTestCase.pkTable;
		do {
			IOListContainer iolc(pkTableTestCase.maxInput, k, k, pkTableTestCase.presentationLayer);
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkEqualityOfClasses(*currentTable, iolc, k, pkTableTestCase.numStates),
				"PkTableRows with the same class in s2c are k equivalent");
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkUnequalityOfClasses(*currentTable, iolc, k, pkTableTestCase.numStates),
				"PkTableRows with the different classes in s2c are not k equivalent");

			currentTable = currentTable->getPkPlusOneTable();
			++k;
		} while (currentTable != nullptr);

		fsmlib_assert("TC-PkTable-NNNN",
			k == 3,
			"PkTable::getPkPlusOneTable() returns nullptr if no new equivalence classes are generated.");
	}
}

// tests PkTable::getMembers(const int c)
void testPkTableGetMembers() {
	// table contains only one row (state). This state is in class 0.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 1;
		int maxInput = 5;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 0);
		string result = pkTable.getMembers(0);
		fsmlib_assert("TC-PkTable-NNNN",
			result.find("0") != string::npos,
			"Result of PkTable::getMembers(const int c) contains all state identifiers of states in class c");
	}
	// Table contains two rows (states). Both are in different classes.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 2;
		int maxInput = 5;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 1);
		pkTable.setClass(1, 0);
		string membersOf0 = pkTable.getMembers(0);
		string membersOf1 = pkTable.getMembers(1);
		fsmlib_assert("TC-PkTable-NNNN",
			membersOf0.find("1") != string::npos
			&& membersOf1.find("0") != string::npos,
			"Result of PkTable::getMembers(const int c) contains all state identifiers of states in class c");

		fsmlib_assert("TC-PkTable-NNNN",
			membersOf0.find("0") == string::npos
			&& membersOf1.find("1") == string::npos,
			"Result of PkTable::getMembers(const int c) contains only state identifiers of states in class c");
	}
	// Table contains four rows (states). Some class contains more than one element.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 4;
		int maxInput = 5;
		PkTable pkTable(numStates, maxInput, pl);
		pkTable.setClass(0, 0);
		pkTable.setClass(1, 0);
		pkTable.setClass(2, 1);
		pkTable.setClass(3, 2);
		string membersOf0 = pkTable.getMembers(0);
		string membersOf1 = pkTable.getMembers(1);
		string membersOf2 = pkTable.getMembers(2);
		fsmlib_assert("TC-PkTable-NNNN",
			membersOf0.find("0") != string::npos && membersOf0.find("1") != string::npos
			&& membersOf1.find("2") != string::npos
			&& membersOf2.find("3") != string::npos,
			"Result of PkTable::getMembers(const int c) contains all state identifiers of states in class c");

		fsmlib_assert("TC-PkTable-NNNN",
			membersOf0.find("2") == string::npos && membersOf0.find("3") == string::npos
			&& membersOf1.find("0") == string::npos && membersOf1.find("1") == string::npos && membersOf1.find("3") == string::npos
			&& membersOf2.find("0") == string::npos && membersOf2.find("1") == string::npos && membersOf2.find("2") == string::npos,
			"Result of PkTable::getMembers(const int c) contains only state identifiers of states in class c");
	}
}

// checks if given node matches given row (same number of transitions, correct target states and correct outputs)
bool matchPkTableRowWithFsmNode(shared_ptr<FsmNode> node, shared_ptr<PkTableRow> row, const PkTable &pkTable,
	const Dfsm &dfsm, int maxInput, shared_ptr<FsmPresentationLayer> pl) {
	for (int i = 0; i <= maxInput; ++i) {
		if (row->get(i) == -1) {
			// node has transition for input but row doesn't
			if (not node->after(i).empty()) {
				return false;
			}
		}
		else {
			// node has more than one transitions for input i
			if (node->after(i).size() != 1) {
				return false;
			}
			// transition has different / false target node 
			if (node->after(i)[0] != dfsm.getNodes()[pkTable.getClass(row->get(i))]) {
				return false;
			}
			OutputTrace outTrc(pl);
			node->apply(i, outTrc);
			// transition has wrong output
			if (outTrc.get()[0] != row->getIOMap().at(i)) {
				return false;
			}
		}
	}
	return true;
}

// returns pointer to first row in pkTable with s2c-class c. 
// returns nullptr if no such row exists.
shared_ptr<PkTableRow> selectRowWithClass(int c, int numStates, PkTable &pkTable) {
	for (int i = 0; i < numStates; ++i) {
		if (pkTable.getClass(i) == c) {
			return pkTable.getRow(i);
		}
	}
	return nullptr;
}

// checks if each PkTable s2c-class is correctly represented in given dfsm. (is there a FsmNode matching a row with some class?)
bool checkPkTableToFsmProperty(PkTable &pkTable, const Dfsm &dfsm, int maxInput, 
	shared_ptr<FsmPresentationLayer> pl, int numStates) {
	
	for (int c = 0; c <= pkTable.maxClassId(); ++c) {
		shared_ptr<FsmNode> node = dfsm.getNodes().at(c);
		shared_ptr<PkTableRow> row = selectRowWithClass(c, numStates, pkTable);
		if (not matchPkTableRowWithFsmNode(node, row, pkTable, dfsm, maxInput, pl)) {
			return false;
		}
	}
	return true;
}

// tests PkTable::toFsm(std::string name, const int maxOutput)
void testPkTableToFsm() {
	{
		PkTableTestCase pkTableTestCase = getPkTableTestCase1();

		shared_ptr<PkTable> currentTable = pkTableTestCase.pkTable;
		do {
			Dfsm dfsm = currentTable->toFsm("", pkTableTestCase.maxOutput);
			fsmlib_assert("TC-PkTable-NNNN",
				dfsm.getNodes().size() == (currentTable->maxClassId() + 1),
				"Number of FsmNodes of PkTable::toFsm(std::string name, const int maxOutput) matches the number of k-equivalence classes.");
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkTableToFsmProperty(*currentTable, dfsm, pkTableTestCase.maxInput, pkTableTestCase.presentationLayer, pkTableTestCase.numStates),
				"Each k-equivalence class has a matching FsmNode");
			currentTable = currentTable->getPkPlusOneTable();
		} while (currentTable != nullptr);
	}

	{
		PkTableTestCase pkTableTestCase = getPkTableTestCase2();

		shared_ptr<PkTable> currentTable = pkTableTestCase.pkTable;
		do {
			Dfsm dfsm = currentTable->toFsm("", pkTableTestCase.maxOutput);
			fsmlib_assert("TC-PkTable-NNNN",
				dfsm.getNodes().size() == (currentTable->maxClassId() + 1),
				"Number of FsmNodes of PkTable::toFsm(std::string name, const int maxOutput) matches the number of k-equivalence classes.");
			fsmlib_assert("TC-PkTable-NNNN",
				checkPkTableToFsmProperty(*currentTable, dfsm, pkTableTestCase.maxInput, pkTableTestCase.presentationLayer, pkTableTestCase.numStates),
				"Each k-equivalence class has a matching FsmNode");

			currentTable = currentTable->getPkPlusOneTable();
		} while (currentTable != nullptr);
	}
	
}

//===================================== DFSMTable Tests ===================================================

// checks if classes in the given pkTable are set according to 1-equivalence.
// row1.class == row2.class <=> row1.ioSection == row2.ioSection
bool check1EquivalenceProperty(shared_ptr<PkTable> pkTable, int numStates) {
	for (int i = 0; i < numStates; ++i) {
		for (int j = i + 1; j < numStates; ++j) {
			if (pkTable->getClass(i) == pkTable->getClass(j)) {
				if (pkTable->getRow(i)->getIOMap() != pkTable->getRow(j)->getIOMap()) {
					return false;
				}
			} 
			// => pkTable->getClass(i) != pkTable->getClass(j)
			else {
				if (pkTable->getRow(i)->getIOMap() == pkTable->getRow(j)->getIOMap()) {
					return false;
				}
			}
		}
	}
	return true;
}

// tests DFSMTable::getP1Table()
void testDFSMTableGetP1Table() {
	// DFSMTable consists of one DFSMTableRow
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 1;
		int maxInput = 1;

		shared_ptr<DFSMTableRow> row0 = make_shared<DFSMTableRow>(0, maxInput);
		row0->getioSection()[0] = 0;
		row0->getioSection()[1] = 0;
		row0->geti2postSection()[0] = 1;
		row0->geti2postSection()[1] = 0;

		DFSMTable dfsmTable(numStates, maxInput, pl);
		dfsmTable.setRow(0, row0);

		shared_ptr<PkTable> p1Table = dfsmTable.getP1Table();
		fsmlib_assert("TC-DFSMTable-NNNN",
			dfsmTable.getRow(0)->getioSection() == p1Table->getRow(0)->getIOMap()
			&& dfsmTable.getRow(0)->geti2postSection() == p1Table->getRow(0)->getI2PMap(),
			"DFSMTable::getP1Table() doesn't change IOMap and I2PMap in result.");

		fsmlib_assert("TC-DFSMTable-NNNN",
			p1Table->getClass(0) != -1,
			"DFSMTable::getP1Table() sets classes in created PkTable according to 1-equivalence (matching IO-sections)");
	}

	// DFSMTable consists of two DFSMTableRows. These aren't 1-equivalent. Not completely specified.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 2;
		int maxInput = 1;

		shared_ptr<DFSMTableRow> row0 = make_shared<DFSMTableRow>(0, maxInput);
		row0->getioSection()[0] = 0;
		row0->getioSection()[1] = 0;
		row0->geti2postSection()[0] = 1;
		row0->geti2postSection()[1] = 1;

		shared_ptr<DFSMTableRow> row1 = make_shared<DFSMTableRow>(1, maxInput);
		row1->getioSection()[0] = 1;
		row1->getioSection()[1] = -1;
		row1->geti2postSection()[0] = 0;
		row1->geti2postSection()[1] = -1;

		DFSMTable dfsmTable(numStates, maxInput, pl);
		dfsmTable.setRow(0, row0);
		dfsmTable.setRow(1, row1);

		shared_ptr<PkTable> p1Table = dfsmTable.getP1Table();
		for (int id = 0; id < numStates; ++id) {
			fsmlib_assert("TC-DFSMTable-NNNN",
				dfsmTable.getRow(id)->getioSection() == p1Table->getRow(id)->getIOMap()
				&& dfsmTable.getRow(id)->geti2postSection() == p1Table->getRow(id)->getI2PMap(),
				"DFSMTable::getP1Table() doesn't change IOMap and I2PMap in result.");
		}

		fsmlib_assert("TC-DFSMTable-NNNN",
			check1EquivalenceProperty(p1Table, numStates),
			"DFSMTable::getP1Table() sets classes in created PkTable according to 1-equivalence (matching IO-sections)");
	}

	// DFSMTable consists of four DFSMTableRows. Only two different io-sections => two different 1-equivalence classes expected
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 4;
		int maxInput = 1;

		shared_ptr<DFSMTableRow> row0 = make_shared<DFSMTableRow>(0, maxInput);
		row0->getioSection()[0] = 0;
		row0->getioSection()[1] = 0;
		row0->geti2postSection()[0] = 1;
		row0->geti2postSection()[1] = 1;

		shared_ptr<DFSMTableRow> row1 = make_shared<DFSMTableRow>(1, maxInput);
		row1->getioSection()[0] = 1;
		row1->getioSection()[1] = -1;
		row1->geti2postSection()[0] = 0;
		row1->geti2postSection()[1] = -1;

		shared_ptr<DFSMTableRow> row2 = make_shared<DFSMTableRow>(2, maxInput);
		row2->getioSection()[0] = 1;
		row2->getioSection()[1] = -1;
		row2->geti2postSection()[0] = 0;
		row2->geti2postSection()[1] = -1;

		shared_ptr<DFSMTableRow> row3 = make_shared<DFSMTableRow>(3, maxInput);
		row3->getioSection()[0] = 0;
		row3->getioSection()[1] = 0;
		row3->geti2postSection()[0] = 0;
		row3->geti2postSection()[1] = 0;

		DFSMTable dfsmTable(numStates, maxInput, pl);
		dfsmTable.setRow(0, row0);
		dfsmTable.setRow(1, row1);
		dfsmTable.setRow(2, row2);
		dfsmTable.setRow(3, row3);

		shared_ptr<PkTable> p1Table = dfsmTable.getP1Table();
		for (int id = 0; id < numStates; ++id) {
			fsmlib_assert("TC-DFSMTable-NNNN",
				dfsmTable.getRow(id)->getioSection() == p1Table->getRow(id)->getIOMap()
				&& dfsmTable.getRow(id)->geti2postSection() == p1Table->getRow(id)->getI2PMap(),
				"DFSMTable::getP1Table() doesn't change IOMap and I2PMap in result.");
		}

		fsmlib_assert("TC-DFSMTable-NNNN",
			check1EquivalenceProperty(p1Table, numStates),
			"DFSMTable::getP1Table() sets classes in created PkTable according to 1-equivalence (matching IO-sections)");
	}

	// DFSMTable consists of four DFSMTableRows. Three different io-sections => three different 1-equivalence classes expected
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 4;
		int maxInput = 1;

		shared_ptr<DFSMTableRow> row0 = make_shared<DFSMTableRow>(0, maxInput);
		row0->getioSection()[0] = 0;
		row0->getioSection()[1] = 0;
		row0->geti2postSection()[0] = 1;
		row0->geti2postSection()[1] = 1;

		shared_ptr<DFSMTableRow> row1 = make_shared<DFSMTableRow>(1, maxInput);
		row1->getioSection()[0] = 1;
		row1->getioSection()[1] = -1;
		row1->geti2postSection()[0] = 0;
		row1->geti2postSection()[1] = -1;

		shared_ptr<DFSMTableRow> row2 = make_shared<DFSMTableRow>(2, maxInput);
		row2->getioSection()[0] = 0;
		row2->getioSection()[1] = 0;
		row2->geti2postSection()[0] = 0;
		row2->geti2postSection()[1] = 0;

		shared_ptr<DFSMTableRow> row3 = make_shared<DFSMTableRow>(3, maxInput);
		row3->getioSection()[0] = 1;
		row3->getioSection()[1] = 0;
		row3->geti2postSection()[0] = 1;
		row3->geti2postSection()[1] = 0;

		DFSMTable dfsmTable(numStates, maxInput, pl);
		dfsmTable.setRow(0, row0);
		dfsmTable.setRow(1, row1);
		dfsmTable.setRow(2, row2);
		dfsmTable.setRow(3, row3);

		shared_ptr<PkTable> p1Table = dfsmTable.getP1Table();
		for (int id = 0; id < numStates; ++id) {
			fsmlib_assert("TC-DFSMTable-NNNN",
				dfsmTable.getRow(id)->getioSection() == p1Table->getRow(id)->getIOMap()
				&& dfsmTable.getRow(id)->geti2postSection() == p1Table->getRow(id)->getI2PMap(),
				"DFSMTable::getP1Table() doesn't change IOMap and I2PMap in result.");
		}

		fsmlib_assert("TC-DFSMTable-NNNN",
			check1EquivalenceProperty(p1Table, numStates),
			"DFSMTable::getP1Table() sets classes in created PkTable according to 1-equivalence (matching IO-sections)");
	}

}

//===================================== OFSMTableRow Tests ===================================================

// tests OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput)
void testOFSMTableRowConstructor() {
	// maxInput = 0, maxOutput = 0
	{
		int maxInput = 0;
		int maxOutput = 0;
		OFSMTableRow row(maxInput, maxOutput);
		for (int input = 0; input <= maxInput; ++input) {
			for (int output = 0; output <= maxOutput; ++output) {
				fsmlib_assert("TC-OFSMTableRow-NNNN",
					row.get(input, output) == -1,
					"OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput) initializes each possible "
					"combination of input and output with -1");
			}
		}
	}

	// maxInput = 1, maxOutput = 0
	{
		int maxInput = 1;
		int maxOutput = 0;
		OFSMTableRow row(maxInput, maxOutput);
		for (int input = 0; input <= maxInput; ++input) {
			for (int output = 0; output <= maxOutput; ++output) {
				fsmlib_assert("TC-OFSMTableRow-NNNN",
					row.get(input, output) == -1,
					"OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput) initializes each possible "
					"combination of input and output with -1");
			}
		}
	}

	// maxInput = 0, maxOutput = 1
	{
		int maxInput = 0;
		int maxOutput = 1;
		OFSMTableRow row(maxInput, maxOutput);
		for (int input = 0; input <= maxInput; ++input) {
			for (int output = 0; output <= maxOutput; ++output) {
				fsmlib_assert("TC-OFSMTableRow-NNNN",
					row.get(input, output) == -1,
					"OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput) initializes each possible "
					"combination of input and output with -1");
			}
		}
	}

	// maxInput = 1, maxOutput = 1
	{
		int maxInput = 1;
		int maxOutput = 1;
		OFSMTableRow row(maxInput, maxOutput);
		for (int input = 0; input <= maxInput; ++input) {
			for (int output = 0; output <= maxOutput; ++output) {
				fsmlib_assert("TC-OFSMTableRow-NNNN",
					row.get(input, output) == -1,
					"OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput) initializes each possible "
					"combination of input and output with -1");
			}
		}
	}
}

// tests OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r)
// Positive case
void testOFSMTableIoEqualsPositive() {
	// this=[[-1]] other=[[-1]], maxInput = 0, maxOutput = 0
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->ioEquals(rowOther)
			&& rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns true if both rows contain -1 entries "
			"for the same indices.");
	}

	// this=[[1]] other=[[2]], maxInput = 0, maxOutput = 0
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 2);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->ioEquals(rowOther)
			&& rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns true if both rows contain -1 entries "
			"for the same indices.");
	}

	// this=[[1,-1], [-1,1]] other=[[2,-1], [-1,3]], maxInput = 1, maxOutput = 1
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 1);
		rowThis->set(0, 1, -1);
		rowThis->set(1, 0, -1);
		rowThis->set(1, 1, 1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 2);
		rowOther->set(0, 1, -1);
		rowOther->set(1, 0, -1);
		rowOther->set(1, 1, 3);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->ioEquals(rowOther)
			&& rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns true if both rows contain -1 entries "
			"for the same indices.");
	}

	// this=[[2, 3, 1], [1, -1, -1]] other=[[0, 0, 0], [2, -1, -1]], maxInput = 1, maxOutput = 2
	{
		int maxInput = 1;
		int maxOutput = 2;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 2);
		rowThis->set(0, 1, 3);
		rowThis->set(0, 2, 1);
		rowThis->set(1, 0, 1);
		rowThis->set(1, 1, -1);
		rowThis->set(1, 2, -1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		rowOther->set(0, 1, 0);
		rowOther->set(0, 2, 0);
		rowOther->set(1, 0, 2);
		rowOther->set(1, 1, -1);
		rowOther->set(1, 2, -1);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->ioEquals(rowOther)
			&& rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns true if both rows contain -1 entries "
			"for the same indices.");
	}
}

// tests OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r)
// Negative case
void testOFSMTableIoEqualsNegative() {
	// this=[[-1]] other=[[0]], maxInput = 0, maxOutput = 0
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->ioEquals(rowOther)
			&& not rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns false if one row contains a -1 entry "
			"for some io-pair, for which the other row doesn't contain a -1 entry.");
	}

	// this=[[0,-1], [1,1]] other=[[0,0], [1,1]], maxInput = 1, maxOutput = 1
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		rowThis->set(0, 1, -1);
		rowThis->set(1, 0, 1);
		rowThis->set(1, 1, 1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		rowOther->set(0, 1, 0);
		rowOther->set(1, 0, 1);
		rowOther->set(1, 1, 1);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->ioEquals(rowOther)
			&& not rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns false if one row contains a -1 entry "
			"for some io-pair, for which the other row doesn't contain a -1 entry.");
	}

	// this=[[-1,-1,-1], [-1,-1,-1]] other=[[-1,-1,-1], [-1,0,-1]], maxInput = 1, maxOutput = 2
	{
		int maxInput = 1;
		int maxOutput = 2;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, -1);
		rowThis->set(0, 1, -1);
		rowThis->set(0, 2, -1);
		rowThis->set(1, 0, -1);
		rowThis->set(1, 1, -1);
		rowThis->set(1, 2, -1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, -1);
		rowOther->set(0, 1, -1);
		rowOther->set(0, 2, -1);
		rowOther->set(1, 0, -1);
		rowOther->set(1, 1, 0);
		rowOther->set(1, 2, -1);
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->ioEquals(rowOther)
			&& not rowOther->ioEquals(rowThis),
			"OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) returns false if one row contains a -1 entry "
			"for some io-pair, for which the other row doesn't contain a -1 entry.");
	}
}

// tests OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r)
// Positive case
void testOFSMTableClassEqualsPositive() {
	// this = [[0]], other = [[0]], s2c = {{0->0}}, maxInput = 0, maxOutput = 0 
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		S2CMap s2c(0);
		s2c[0] = 0;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->classEquals(s2c, rowOther)
			&& rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns true"
			" if each io-pair results in poststates s, s' with s2c[s] == s2c[s']");
	}

	// this = [[-1]], other = [[-1]], s2c = {{0->0}}, maxInput = 0, maxOutput = 0 
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, -1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, -1);
		S2CMap s2c(0);
		s2c[0] = 0;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->classEquals(s2c, rowOther)
			&& rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns true"
			" if each io-pair results in poststates s, s' with s2c[s] == s2c[s']");
	}

	// this = [[0]], other = [[1]], s2c = {{0->1},{1->1}}, maxInput = 0, maxOutput = 0 
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 1);
		S2CMap s2c(1);
		s2c[0] = 1;
		s2c[1] = 1;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->classEquals(s2c, rowOther)
			&& rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns true"
			" if each io-pair results in poststates s, s' with s2c[s] == s2c[s']");
	}

	// this = [[0,-1], [-1,2]], other = [[1,-1], [-1,2]], s2c = {{0->1},{1->1},{2->2}}, maxInput = 1, maxOutput = 1 
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		rowThis->set(0, 1, -1);
		rowThis->set(1, 0, -1);
		rowThis->set(1, 1, 2);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 1);
		rowOther->set(0, 1, -1);
		rowOther->set(1, 0, -1);
		rowOther->set(1, 1, 2);
		S2CMap s2c(2);
		s2c[0] = 1;
		s2c[1] = 1;
		s2c[2] = 2;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->classEquals(s2c, rowOther)
			&& rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns true"
			" if each io-pair results in poststates s, s' with s2c[s] == s2c[s']");
	}

	// this = [[0,1], [2,3]], other = [[1,4], [4,5]], s2c = {{0->2},{1->2},{2->2},{3->4},{4->2},{5->4}}, maxInput = 1, maxOutput = 1 
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		rowThis->set(0, 1, 1);
		rowThis->set(1, 0, 2);
		rowThis->set(1, 1, 3);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 1);
		rowOther->set(0, 1, 4);
		rowOther->set(1, 0, 4);
		rowOther->set(1, 1, 5);
		S2CMap s2c(5);
		s2c[0] = 2;
		s2c[1] = 2;
		s2c[2] = 2;
		s2c[3] = 4;
		s2c[4] = 2;
		s2c[5] = 4;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			rowThis->classEquals(s2c, rowOther)
			&& rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns true"
			" if each io-pair results in poststates s, s' with s2c[s] == s2c[s']");
	}
}

// tests OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r)
// Negative case
void testOFSMTableClassEqualsNegative() {
	// this = [[-1]], other = [[0]], s2c = {{0->0}}, maxInput = 0, maxOutput = 0 
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, -1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		S2CMap s2c(0);
		s2c[0] = 0;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->classEquals(s2c, rowOther)
			&& not rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns false"
			" if both rows aren't io-equivalent (wrt. OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r))");
	}

	// this = [[0]], other = [[1]], s2c = {{0->0},{1->1}}, maxInput = 0, maxOutput = 0 
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 1);
		S2CMap s2c(1);
		s2c[0] = 0;
		s2c[1] = 1;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->classEquals(s2c, rowOther)
			&& not rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns false"
			" if there is some io-pair which results in poststates s, s' with s2c[s] != s2c[s']");
	}

	// this = [[0,1], [-1,2]], other = [[0,2], [-1,2]], s2c = {{0->0},{1->1},{2->2}}, maxInput = 1, maxOutput = 1 
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		rowThis->set(0, 1, 1);
		rowThis->set(1, 0, -1);
		rowThis->set(1, 1, 2);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 0);
		rowOther->set(0, 1, 2);
		rowOther->set(1, 0, -1);
		rowOther->set(1, 1, 2);
		S2CMap s2c(2);
		s2c[0] = 0;
		s2c[1] = 1;
		s2c[2] = 2;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->classEquals(s2c, rowOther)
			&& not rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns false"
			" if there is some io-pair which results in poststates s, s' with s2c[s] != s2c[s']");
	}

	// this = [[0,-1], [-1,2], [-1,-1]], other = [[1,-1], [-1,3], [-1,-1]], s2c = {{0->2},{1->3},{2->1},{3->2}}, maxInput = 2, maxOutput = 1 
	{
		int maxInput = 2;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> rowThis = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowThis->set(0, 0, 0);
		rowThis->set(0, 1, -1);
		rowThis->set(1, 0, -1);
		rowThis->set(1, 1, 2);
		rowThis->set(2, 0, -1);
		rowThis->set(2, 1, -1);
		shared_ptr<OFSMTableRow> rowOther = make_shared<OFSMTableRow>(maxInput, maxOutput);
		rowOther->set(0, 0, 1);
		rowOther->set(0, 1, -1);
		rowOther->set(1, 0, -1);
		rowOther->set(1, 1, 3);
		rowOther->set(2, 0, -1);
		rowOther->set(2, 1, -1);
		S2CMap s2c(3);
		s2c[0] = 2;
		s2c[1] = 3;
		s2c[2] = 1;
		s2c[3] = 2;
		fsmlib_assert("TC-OFSMTableRow-NNNN",
			not rowThis->classEquals(s2c, rowOther)
			&& not rowOther->classEquals(s2c, rowThis),
			"OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r) returns false"
			" if there is some io-pair which results in poststates s, s' with s2c[s] != s2c[s']");
	}
}

//===================================== OFSMTable Tests ===================================================

// tests OFSMTable::OFSMTable(const vector<shared_ptr<FsmNode>>& nodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer> presentationLayer)
void testOFSMTableConstructor() {

	// nodes contains only one FsmNode which hasn't any transition.
	{
		int maxInput = 0;
		int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		vector<shared_ptr<FsmNode>> nodes{ n0 };
		OFSMTable table(nodes, maxInput, maxOutput, pl);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.get(n0->getId(), 0, 0) == -1,
			"Constructed OFSMTable marks each non-existing transitions with -1");
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.getS2C().at(n0->getId()) == 0,
			"Constructed OFSMTable maps each state to equivalence class 0.");
	}

	// nodes contains two FsmNodes with transitions.
	// n0 --(0/1)--> n0; n0 --(1/1)--> n1; n1 --(0/0)--> n0; n1 --(1/0)--> n0;
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 1, pl)));
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(0, 0, pl)));
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(1, 0, pl)));
		vector<shared_ptr<FsmNode>> nodes{ n0, n1 };
		OFSMTable table(nodes, maxInput, maxOutput, pl);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.get(n0->getId(), 0, 0) == -1
			&& table.get(n0->getId(), 0, 1) == n0->getId()
			&& table.get(n0->getId(), 1, 0) == -1
			&& table.get(n0->getId(), 1, 1) == n1->getId()
			&& table.get(n1->getId(), 0, 0) == n0->getId()
			&& table.get(n1->getId(), 0, 1) == -1
			&& table.get(n1->getId(), 1, 0) == n0->getId()
			&& table.get(n1->getId(), 1, 1) == -1,
			"Each row from constructed OFSMTable represents correctly the corresponding FsmNode (nodes[i] ~ rows[i])");
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.getS2C().at(n0->getId()) == 0
			&& table.getS2C().at(n1->getId()) == 0,
			"Constructed OFSMTable maps each state to equivalence class 0.");
	}

	// nodes contains three FsmNodes with transitions.
	// n0 --(0/0)--> n1; n1 --(1/1)--> n2; n2 --(0/0)--> n2; n2 --(1/0)--> n0;
	{
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(1, 1, pl)));
		n2->addTransition(make_shared<FsmTransition>(n2, n2, make_shared<FsmLabel>(0, 0, pl)));
		n2->addTransition(make_shared<FsmTransition>(n2, n0, make_shared<FsmLabel>(1, 0, pl)));
		vector<shared_ptr<FsmNode>> nodes{ n0, n1, n2 };
		OFSMTable table(nodes, maxInput, maxOutput, pl);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.get(n0->getId(), 0, 0) == n1->getId()
			&& table.get(n0->getId(), 0, 1) == -1
			&& table.get(n0->getId(), 1, 0) == -1
			&& table.get(n0->getId(), 1, 1) == -1
			&& table.get(n1->getId(), 0, 0) == -1
			&& table.get(n1->getId(), 0, 1) == -1
			&& table.get(n1->getId(), 1, 0) == -1
			&& table.get(n1->getId(), 1, 1) == n2->getId()
			&& table.get(n2->getId(), 0, 0) == n2->getId()
			&& table.get(n2->getId(), 0, 1) == -1
			&& table.get(n2->getId(), 1, 0) == n0->getId()
			&& table.get(n2->getId(), 1, 1) == -1,
			"Each row from constructed OFSMTable represents correctly the corresponding FsmNode (nodes[i] ~ rows[i])");
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.getS2C().at(n0->getId()) == 0
			&& table.getS2C().at(n1->getId()) == 0
			&& table.getS2C().at(n2->getId()) == 0,
			"Constructed OFSMTable maps each state to equivalence class 0.");
	}
}

// tests OFSMTable::maxClassId()
void testOFSMTableMaxClassId() {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	int maxInput = 4;
	int maxOutput = 2;
	vector<std::shared_ptr<OFSMTableRow>> rows;
	// s2c = [-1] -> all s2c ids are smaller than 0
	{
		int numStates = 1;		
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		S2CMap s2c(0);
		s2c[0] = -1;
		table.setS2C(s2c);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.maxClassId() == 0,
			"OFSMTable::maxClassId() returns 0 if all elements from s2c are smaller than 0");
	}

	// s2c = [-1,0] -> not all s2c ids are smaller than 0
	{
		int numStates = 2;
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		S2CMap s2c(1);
		s2c[0] = -1;
		s2c[1] = 0;
		table.setS2C(s2c);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.maxClassId() == 0,
			"OFSMTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,0,1]
	{
		int numStates = 3;
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		S2CMap s2c(2);
		s2c[0] = 0;
		s2c[1] = 0;
		s2c[2] = 1;
		table.setS2C(s2c);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.maxClassId() == 1,
			"OFSMTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,1,2,1]
	{
		int numStates = 4;
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		S2CMap s2c(3);
		s2c[0] = 0;
		s2c[1] = 1;
		s2c[2] = 2;
		s2c[3] = 1;
		table.setS2C(s2c);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.maxClassId() == 2,
			"OFSMTable::maxClassId() returns greatest value of s2c");
	}

	// s2c = [0,1,3,3]
	{
		int numStates = 4;
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		S2CMap s2c(3);
		s2c[0] = 0;
		s2c[1] = 1;
		s2c[2] = 3;
		s2c[3] = 3;
		table.setS2C(s2c);
		fsmlib_assert("TC-OFSMTable-NNNN",
			table.maxClassId() == 3,
			"OFSMTable::maxClassId() returns greatest value of s2c");
	}
}

// Checks if stateID1 and stateID2 (only) have transitions for the same io-pairs
bool isIOEquivalent(int stateID1, int stateID2, OFSMTable &table, int maxInput, int maxOutput) {
	for (int input = 0; input <= maxInput; ++input) {
		for (int output = 0; output <= maxOutput; ++output) {
			if (table.get(stateID1, input, output) < 0 && table.get(stateID2, input, output) >= 0) return false;

			if (table.get(stateID1, input, output) >= 0 && table.get(stateID2, input, output) < 0) return false;
		}
	}
	return true;
}

// Checks if all states in the same class are io-equivalent and if all io-equivalent states are in the same class.
// Can be used to check if table is a correct OFSM 1 Table
bool checkOFSM1TableProperty(OFSMTable &table, int numStates, int maxInput, int maxOutput) {
	for (int i = 0; i < numStates; ++i) {
		S2CMap s2c = table.getS2C();
		for (int j = i + 1; j < numStates; ++j) {
			if (s2c.at(i) == s2c.at(j)) {
				// => non io-equivalent states in the same ofsm 1 class
				if (not isIOEquivalent(i, j, table, maxInput, maxOutput)) {
					return false;
				}
			}
			else {
				// => io-equivalent states in different ofsm 1 classes.
				if (isIOEquivalent(i, j, table, maxInput, maxOutput)) {
					return false;
				}
			}
		}
	}
	return true;
}

struct OFSMTableTestCase
{
	shared_ptr<FsmPresentationLayer> presentationLayer;
	shared_ptr<OFSMTable> ofsmTable;
	int maxInput;
	int maxOutput;
	int numStates;
};

OFSMTableTestCase getOFSMTableTestCase1() {
	std::shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
	int maxInput = 1;
	int maxOutput = 1;

	// create all OFSMTableRows
	shared_ptr<OFSMTableRow> r0 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r0->set(0, 0, 1);
	r0->set(0, 1, 3);
	r0->set(1, 0, 6);
	r0->set(1, 1, 0);
	shared_ptr<OFSMTableRow> r1 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r1->set(0, 0, 2);
	r1->set(0, 1, 2);
	r1->set(1, 0, 3);
	r1->set(1, 1, -1);
	shared_ptr<OFSMTableRow> r2 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r2->set(0, 0, -1);
	r2->set(0, 1, 0);
	r2->set(1, 0, -1);
	r2->set(1, 1, 2);
	shared_ptr<OFSMTableRow> r3 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r3->set(0, 0, 4);
	r3->set(0, 1, 5);
	r3->set(1, 0, 1);
	r3->set(1, 1, -1);
	shared_ptr<OFSMTableRow> r4 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r4->set(0, 0, -1);
	r4->set(0, 1, 0);
	r4->set(1, 0, -1);
	r4->set(1, 1, 4);
	shared_ptr<OFSMTableRow> r5 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r5->set(0, 0, -1);
	r5->set(0, 1, 0);
	r5->set(1, 0, -1);
	r5->set(1, 1, 2);
	shared_ptr<OFSMTableRow> r6 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r6->set(0, 0, 7);
	r6->set(0, 1, 7);
	r6->set(1, 0, 3);
	r6->set(1, 1, -1);
	shared_ptr<OFSMTableRow> r7 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r7->set(0, 0, 8);
	r7->set(0, 1, 8);
	r7->set(1, 0, -1);
	r7->set(1, 1, 7);
	shared_ptr<OFSMTableRow> r8 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r8->set(0, 0, 7);
	r8->set(0, 1, 4);
	r8->set(1, 0, 6);
	r8->set(1, 1, 8);

	vector<shared_ptr<OFSMTableRow>> rows;
	rows.push_back(r0);
	rows.push_back(r1);
	rows.push_back(r2);
	rows.push_back(r3);
	rows.push_back(r4);
	rows.push_back(r5);
	rows.push_back(r6);
	rows.push_back(r7);
	rows.push_back(r8);

	// create OFSMTable from OFSMTableRows
	int numStates = rows.size();
	OFSMTable ofsmTable(numStates, maxInput, maxOutput, rows, presentationLayer);

	OFSMTableTestCase testStructure;
	testStructure.presentationLayer = presentationLayer;
	testStructure.ofsmTable = make_shared<OFSMTable>(ofsmTable);
	testStructure.maxInput = maxInput;
	testStructure.maxOutput = maxOutput;
	testStructure.numStates = numStates;

	return testStructure;
}

OFSMTableTestCase getOFSMTableTestCase2() {
	std::shared_ptr<FsmPresentationLayer> presentationLayer = make_shared<FsmPresentationLayer>();
	int maxInput = 2;
	int maxOutput = 2;

	// create all OFSMTableRows
	// a = 0, b = 1, c = 2
	shared_ptr<OFSMTableRow> r0 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r0->set(0, 0, -1); r0->set(0, 2, 4);
	r0->set(1, 1, 1); r0->set(1, 2, 0);
	r0->set(2, 1, 2);
	shared_ptr<OFSMTableRow> r1 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r1->set(0, 0, -1); r1->set(0, 2, 3);
	r1->set(1, 1, 2); r1->set(1, 2, 0);
	r1->set(2, 1, 1);
	shared_ptr<OFSMTableRow> r2 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r2->set(0, 0, 2); r2->set(0, 2, -1);
	r2->set(1, 1, -1); r2->set(1, 2, 3);
	r2->set(2, 1, 0);
	shared_ptr<OFSMTableRow> r3 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r3->set(0, 0, 2); r3->set(0, 2, 0);
	r3->set(1, 1, 3); r3->set(1, 2, -1);
	r3->set(2, 1, 3);
	shared_ptr<OFSMTableRow> r4 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r4->set(0, 0, 2); r4->set(0, 2, 3); 
	r4->set(1, 1, 2); r4->set(1, 2, 5);
	r4->set(2, 1, 6);
	shared_ptr<OFSMTableRow> r5 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r5->set(0, 0, 2); r5->set(0, 2, 7); 
	r5->set(1, 1, 8); r5->set(1, 2, 0);
	r5->set(2, 1, 9);
	shared_ptr<OFSMTableRow> r6 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r6->set(0, 0, -1); r6->set(0, 2, 10); 
	r6->set(1, 1, 4); r6->set(1, 2, 0);
	r6->set(2, 1, 4);
	shared_ptr<OFSMTableRow> r7 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r7->set(0, 0, 2); r7->set(0, 2, 10); 
	r7->set(1, 1, 4); r7->set(1, 2, 5);
	r7->set(2, 1, 7);
	shared_ptr<OFSMTableRow> r8 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r8->set(0, 0, 2); r8->set(0, 2, 5); 
	r8->set(1, 1, 9); r8->set(1, 2, 0);
	r8->set(2, 1, 8);
	shared_ptr<OFSMTableRow> r9 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r9->set(0, 0, 2); r9->set(0, 2, 0);
	r9->set(1, 1, 3); r9->set(1, 2, 3);
	r9->set(2, 1, 5);
	shared_ptr<OFSMTableRow> r10 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r10->set(0, 0, 2); r10->set(0, 2, 5);
	r10->set(1, 1, 9); r10->set(1, 2, 5);
	r10->set(2, 1, 11);
	shared_ptr<OFSMTableRow> r11 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r11->set(0, 0, 2); r11->set(0, 2, 12); // 2 12 10 0 10 
	r11->set(1, 1, 10); r11->set(1, 2, 0);
	r11->set(2, 1, 10);
	shared_ptr<OFSMTableRow> r12 = make_shared<OFSMTableRow>(maxInput, maxOutput);
	r12->set(0, 0, 2); r12->set(0, 2, 12);
	r12->set(1, 1, 10); r12->set(1, 2, 5);
	r12->set(2, 1, 12);

	vector<shared_ptr<OFSMTableRow>> rows;
	rows.push_back(r0);
	rows.push_back(r1);
	rows.push_back(r2);
	rows.push_back(r3);
	rows.push_back(r4);
	rows.push_back(r5);
	rows.push_back(r6);
	rows.push_back(r7);
	rows.push_back(r8);
	rows.push_back(r9);
	rows.push_back(r10);
	rows.push_back(r11);
	rows.push_back(r12);

	// create OFSMTable from OFSMTableRows
	int numStates = rows.size();
	OFSMTable ofsmTable(numStates, maxInput, maxOutput, rows, presentationLayer);

	OFSMTableTestCase testStructure;
	testStructure.presentationLayer = presentationLayer;
	testStructure.ofsmTable = make_shared<OFSMTable>(ofsmTable);
	testStructure.maxInput = maxInput;
	testStructure.maxOutput = maxOutput;
	testStructure.numStates = numStates;

	return testStructure;
}

// Checks if s2c-class of stateID is equal to the classes of the state ids in compareIDs.
bool classEquals(int stateID, const vector<int> &compareIDs, S2CMap & s2c) {
	int cls = s2c.at(stateID);
	for (int i : compareIDs) {
		if (cls != s2c.at(i)) {
			return false;
		}
	}
	return true;
}

// Checks if s2c-class of stateID is different than all classes of the states in compareIDs.
bool classUnequals(int stateID, const vector<int> &compareIDs, S2CMap & s2c) {
	int cls = s2c.at(stateID);
	for (int i : compareIDs) {
		if (cls == s2c.at(i)) {
			return false;
		}
	}
	return true;
}

// test OFSMTable::next()
void testOFSMTableNext() {
	// test OFSMTable::nextAfterZero() (call next() on table with id 0)
	// table:
	// q | 0/0 | 0/1
	// 0 |  -1 | 0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 1;
		int maxInput = 0;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> r0 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r0->set(0, 1, 0);
		vector<std::shared_ptr<OFSMTableRow>> rows{ r0 };
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		shared_ptr<OFSMTable> nextAfterZero = table.next();
		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getS2C().at(0) != -1,
			"OFSMTable::nextAfterZero() maps two states to the same s2c-class iff they are io-equivalent (wrt. OFSMTableRow::ioEquals)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getId() == 1,
			"OFSMTable::nextAfterZero() sets tblId to 1");
	}

	// test OFSMTable::nextAfterZero() (call next() on table with id 0)
	// table:
	// q | 0/0 | 0/1
	// 0 |  -1 | 0
	// 1 |  0  | -1
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 2;
		int maxInput = 0;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> r0 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r0->set(0, 1, 0);
		shared_ptr<OFSMTableRow> r1 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r1->set(0, 0, 0);
		vector<std::shared_ptr<OFSMTableRow>> rows{ r0, r1 };
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		shared_ptr<OFSMTable> nextAfterZero = table.next();
		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getS2C().at(0) != -1
			&& nextAfterZero->getS2C().at(1) != -1
			&& nextAfterZero->getS2C().at(0) != nextAfterZero->getS2C().at(1)
			&& checkOFSM1TableProperty(*nextAfterZero, numStates, maxInput, maxOutput),
			"OFSMTable::nextAfterZero() maps two states to the same s2c-class iff they are io-equivalent (wrt. OFSMTableRow::ioEquals)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getId() == 1,
			"OFSMTable::nextAfterZero() sets tblId to 1");
	}

	// test OFSMTable::nextAfterZero() (call next() on table with id 0)
	// table:
	// q | 0/0 | 0/1 | 1/0 | 1/1
	// 0 |  0  |  0  |  0  |  0
	// 1 |  -1 |  0  |  0  |  0
	// 2 |   1 |  2  |  0  |  0
	// 3 |  -1 |  1  |  0  |  0
	// 4 |  -1 |  -1 |  0  |  0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		int numStates = 5;
		int maxInput = 1;
		int maxOutput = 1;
		shared_ptr<OFSMTableRow> r0 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r0->set(0, 0, 0);
		r0->set(0, 1, 0);
		r0->set(1, 0, 0);
		r0->set(1, 1, 0);
		shared_ptr<OFSMTableRow> r1 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r1->set(0, 1, 0);
		r1->set(1, 0, 0);
		r1->set(1, 1, 0);
		shared_ptr<OFSMTableRow> r2 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r2->set(0, 0, 1);
		r2->set(0, 1, 2);
		r2->set(1, 0, 0);
		r2->set(1, 1, 0);
		shared_ptr<OFSMTableRow> r3 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r3->set(0, 1, 1);
		r3->set(1, 0, 0);
		r3->set(1, 1, 0);
		shared_ptr<OFSMTableRow> r4 = make_shared<OFSMTableRow>(maxInput, maxOutput);
		r4->set(1, 0, 0);
		r4->set(1, 1, 0);
		vector<std::shared_ptr<OFSMTableRow>> rows{ r0, r1, r2, r3, r4 };
		OFSMTable table(numStates, maxInput, maxOutput, rows, pl);
		shared_ptr<OFSMTable> nextAfterZero = table.next();
		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getS2C().at(0) != -1
			&& nextAfterZero->getS2C().at(1) != -1
			&& nextAfterZero->getS2C().at(2) != -1
			&& nextAfterZero->getS2C().at(3) != -1
			&& nextAfterZero->getS2C().at(4) != -1
			&& checkOFSM1TableProperty(*nextAfterZero, numStates, maxInput, maxOutput),
			"OFSMTable::nextAfterZero() maps two states to the same s2c-class iff they are io-equivalent (wrt. OFSMTableRow::ioEquals)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			nextAfterZero->getId() == 1,
			"OFSMTable::nextAfterZero() sets tblId to 1");
	}

	// test of successive invokations (full minimisation process) - using getOFSMTableTestCase1()
	{
		OFSMTableTestCase ofsmTableTestCase = getOFSMTableTestCase1();
		shared_ptr<OFSMTable> currentTable = ofsmTableTestCase.ofsmTable;
		currentTable = currentTable->next();
		S2CMap s2c = currentTable->getS2C();
		// currentTable is now OFSM Table 1
		fsmlib_assert("TC-OFSMTable-NNNN",
			s2c.at(0) == s2c.at(8)
			&& s2c.at(1) == s2c.at(3)
			&& s2c.at(1) == s2c.at(6)
			&& s2c.at(2) == s2c.at(4)
			&& s2c.at(2) == s2c.at(5)
			&& s2c.at(0) != s2c.at(1)
			&& s2c.at(0) != s2c.at(2)
			&& s2c.at(0) != s2c.at(7)
			&& s2c.at(1) != s2c.at(2)
			&& s2c.at(1) != s2c.at(7)
			&& s2c.at(2) != s2c.at(7),
			"OFSMTable::next() maps each state to the correct equivalence class.");

		currentTable = currentTable->next();
		s2c = currentTable->getS2C();
		fsmlib_assert("TC-OFSMTable-NNNN",
			s2c.at(1) == s2c.at(3)
			&& s2c.at(2) == s2c.at(4)
			&& s2c.at(2) == s2c.at(5)
			&& s2c.at(0) != s2c.at(1) && s2c.at(0) != s2c.at(2) && s2c.at(0) != s2c.at(6) && s2c.at(0) != s2c.at(7) && s2c.at(0) != s2c.at(8)
			&& s2c.at(1) != s2c.at(2) && s2c.at(1) != s2c.at(6) && s2c.at(1) != s2c.at(7) && s2c.at(1) != s2c.at(8)
			&& s2c.at(2) != s2c.at(6) && s2c.at(2) != s2c.at(7) && s2c.at(2) != s2c.at(8)
			&& s2c.at(6) != s2c.at(7) && s2c.at(6) != s2c.at(8)
			&& s2c.at(7) != s2c.at(8),
			"OFSMTable::next() maps each state to the correct equivalence class.");

		currentTable = currentTable->next();
		fsmlib_assert("TC-OFSMTable-NNNN",
			currentTable == nullptr,
			"OFSMTable::next() returns nullptr if no new class can be generated");
	}

	// test of successive invokations (full minimisation process) - using getOFSMTableTestCase2()
	{
		OFSMTableTestCase ofsmTableTestCase = getOFSMTableTestCase2();
		shared_ptr<OFSMTable> currentTable = ofsmTableTestCase.ofsmTable;
		currentTable = currentTable->next();
		S2CMap s2c = currentTable->getS2C();

		// currentTable is now OFSM Table 1
		fsmlib_assert("TC-OFSMTable-NNNN",
			currentTable->maxClassId() == 3,
			"Result of OFSMTable::next() contains the right number of classes.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			classEquals(0, vector<int>{0,1,6}, s2c) // class 0
			&& classEquals(2, vector<int>{2}, s2c)  // class 1
			&& classEquals(3, vector<int>{3}, s2c)  // class 2
			&& classEquals(4, vector<int>{4,5,7,8,9,10,11,12}, s2c)  // class 3
			&& classUnequals(0, vector<int>{2,3,4,5, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(2, vector<int>{0, 1, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(3, vector<int>{0, 1, 2, 4, 6, 5, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(4, vector<int>{0, 1, 2, 3, 6,}, s2c),
			"OFSMTable::next() maps each state to the correct equivalence class.");

		currentTable = currentTable->next();
		s2c = currentTable->getS2C();

		fsmlib_assert("TC-OFSMTable-NNNN",
			currentTable->maxClassId() == 8,
			"Result of OFSMTable::next() contains the right number of classes.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			classEquals(0, vector<int>{0}, s2c) // class 0
			&& classEquals(1, vector<int>{1}, s2c) // class 4
			&& classEquals(2, vector<int>{2}, s2c) // class 1
			&& classEquals(3, vector<int>{3}, s2c) // class 2
			&& classEquals(4, vector<int>{4}, s2c) // class 3
			&& classEquals(5, vector<int>{5,8,11}, s2c) // class 6
			&& classEquals(6, vector<int>{6}, s2c) // class 5
			&& classEquals(7, vector<int>{7, 10, 12}, s2c) // class 7
			&& classEquals(9, vector<int>{9}, s2c) // class 8
			&& classUnequals(0, vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(1, vector<int>{0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(2, vector<int>{0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(3, vector<int>{0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(4, vector<int>{0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(5, vector<int>{0, 1, 2, 3, 4, 6, 7, 9, 10, 12}, s2c)
			&& classUnequals(6, vector<int>{0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(7, vector<int>{0, 1, 2, 3, 4, 5, 6, 8, 9, 11}, s2c),
			"OFSMTable::next() maps each state to the correct equivalence class.");

		currentTable = currentTable->next();
		s2c = currentTable->getS2C();

		fsmlib_assert("TC-OFSMTable-NNNN",
			currentTable->maxClassId() == 12,
			"Result of OFSMTable::next() contains the right number of classes.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			classUnequals(0, vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(1, vector<int>{0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(2, vector<int>{0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(3, vector<int>{0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(4, vector<int>{0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(5, vector<int>{0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(6, vector<int>{0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(7, vector<int>{0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12}, s2c)
			&& classUnequals(8, vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12}, s2c)
			&& classUnequals(9, vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12}, s2c)
			&& classUnequals(10, vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12}, s2c)
			&& classUnequals(11, vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12}, s2c)
			&& classUnequals(12, vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, s2c),
			"OFSMTable::next() maps each state to the correct equivalence class.");


		currentTable = currentTable->next();
		fsmlib_assert("TC-OFSMTable-NNNN",
			currentTable == nullptr,
			"OFSMTable::next() returns nullptr if no new class can be generated");
	}
}

// Returns the index of the first OFSMTableRow in Table with class c. 
// Returns -1 if no such row exists.
int selectRowIndexWithClass(int c, shared_ptr<OFSMTable> table, int numStates) {
	S2CMap s2c = table->getS2C();
	for (int i = 0; i < numStates; ++i) {
		if (s2c.at(i) == c) {
			return i;
		}
	}
	return -1;
}

// Calculates the number of transitions that are defined in the state represented by OFSMTableRow at index rowIndex in table.
int getNumberOfTransitionsOfOFSMTableRow(shared_ptr<OFSMTable> table, int rowIndex, int maxInput, int maxOutput) {
	int counter = 0;
	for (int i = 0; i <= maxInput; ++i) {
		for (int j = 0; j <= maxOutput; ++j) {
			if (table->get(rowIndex, i, j) >= 0) {
				++counter;
			}
		}
	}
	return counter;
}

// Check if each FsmTransition of n is also represented in row at rowIndex.
bool matchFsmNodeAgainstOFSMTableRow(shared_ptr<FsmNode> n, shared_ptr<OFSMTable> table, int rowIndex) {
	S2CMap s2c = table->getS2C();
	for (shared_ptr<FsmTransition> t : n->getTransitions()) {
		int postStateID = table->get(rowIndex, t->getLabel()->getInput(), t->getLabel()->getOutput());
		if (s2c.at(postStateID) != t->getTarget()->getId()) {
			return false;
		}
	}
	return true;
}

// Select and return all FsmTransitions of n labeled with input/output.
vector<shared_ptr<FsmTransition>> getFsmTransitionsForIO(shared_ptr<FsmNode> n, int input, int output) {
	vector<shared_ptr<FsmTransition>> transitions;
	for (shared_ptr<FsmTransition> t : n->getTransitions()) {
		if (t->getLabel()->getInput() == input && t->getLabel()->getOutput() == output) {
			transitions.push_back(t);
		}
	}
	return transitions;
}

// Check if each transition in the row at rowIndex is also represented as a FsmTransition of n.
bool matchOFSMTableRowAgainstFsmNode(shared_ptr<FsmNode> n, shared_ptr<OFSMTable> table, int rowIndex, int maxInput, int maxOutput) {
	S2CMap s2c = table->getS2C();
	for (int i = 0; i <= maxInput; ++i) {
		for (int o = 0; o <= maxOutput; ++o) {
			int postState = table->get(rowIndex, i, o);
			// => no transition for i,o
			if (postState < 0) continue;
			
			vector<shared_ptr<FsmTransition>> transitions = getFsmTransitionsForIO(n, i, o);
			if (transitions.size() != 1) return false;
			if (transitions.at(0)->getTarget()->getId() != s2c.at(postState)) return false;
		}
	}
	return true;
}

// Checks if given FsmNode n corresponds to OFSMTable row in table at index rowIndex.
// True if both have the same number of transitions, io-labels of transitions match and targets of transitions match.
// False otherwise.
bool matchFsmNodeWithOFSMTableRow(shared_ptr<FsmNode> n, shared_ptr<OFSMTable> table, int rowIndex, int maxInput, int maxOutput) {
	if (n->getTransitions().size() != getNumberOfTransitionsOfOFSMTableRow(table, rowIndex, maxInput, maxOutput)) {
		return false;
	}
	if (not matchFsmNodeAgainstOFSMTableRow(n, table, rowIndex)) {
		return false;
	}
	if (not matchOFSMTableRowAgainstFsmNode(n, table, rowIndex, maxInput, maxOutput)) {
		return false;
	}
	return true;
}

// Check if each FsmNode from Fsm has corresponding OFSMTableRow in table.
bool matchAllFsmNodes(Fsm &fsm, shared_ptr<OFSMTable> table, int numStates, int maxInput, int maxOutput) {
	for (shared_ptr<FsmNode> n : fsm.getNodes()) {
		int rowIndex = selectRowIndexWithClass(n->getId(), table, numStates);
		if (not matchFsmNodeWithOFSMTableRow(n, table, rowIndex, maxInput, maxOutput)) {
			return false;
		}
	}
	return true;
}

// Check if each OFSMTableRow in table has a corresponding FsmNode in Fsm.
bool matchAllOFSMRows(Fsm &fsm, shared_ptr<OFSMTable> table, int numStates, int maxInput, int maxOutput) {
	S2CMap s2c = table->getS2C();
	for (int i = 0; i < numStates; ++i) {
		int nodeID = s2c.at(i);
		shared_ptr<FsmNode> node = fsm.getNodes().at(nodeID);
		if (not matchFsmNodeWithOFSMTableRow(node, table, i, maxInput, maxOutput)) {
			return false;
		}
	}
	return true;
}

//tests OFSMTable::toFsm(const string & name)
void testOFSMTableToFsm() {
	// using OFSMTable from lecture (Minimisation of OFSM Table)
	{
		OFSMTableTestCase ofsmTableTestCase = getOFSMTableTestCase1();
		shared_ptr<OFSMTable> table = ofsmTableTestCase.ofsmTable;
		shared_ptr<OFSMTable> next = table->next();
		while (next != nullptr) {
			table = next;
			next = next->next();
		}
		
		Fsm fsm = table->toFsm("");
		fsmlib_assert("TC-OFSMTable-NNNN",
			fsm.getNodes().size() == table->maxClassId() + 1,
			"OFSMTable::toFsm(const string & name) creates Fsm with correct number of FsmNodes "
			"(equal to the number of generated classes)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllFsmNodes(fsm, table, ofsmTableTestCase.numStates, ofsmTableTestCase.maxInput, ofsmTableTestCase.maxOutput),
			"OFSMTable::toFsm(const string & name): Each FsmNode from constructed Fsm has corresponding OFSMTableRow in table.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllOFSMRows(fsm, table, ofsmTableTestCase.numStates, ofsmTableTestCase.maxInput, ofsmTableTestCase.maxOutput),
			"OFSMTable::toFsm(const string & name): Each OFSMTableRow from the OFSMTable has a corresponding FsmNode in the constructed Fsm");
		
	}

	// Table is already minimal
	{
		OFSMTableTestCase ofsmTableTestCase = getOFSMTableTestCase2();
		shared_ptr<OFSMTable> table = ofsmTableTestCase.ofsmTable;
		shared_ptr<OFSMTable> next = table->next();
		while (next != nullptr) {
			table = next;
			next = next->next();
		}

		Fsm fsm = table->toFsm("");
		fsmlib_assert("TC-OFSMTable-NNNN",
			fsm.getNodes().size() == table->maxClassId() + 1,
			"OFSMTable::toFsm(const string & name) creates Fsm with correct number of FsmNodes "
			"(equal to the number of generated classes)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllFsmNodes(fsm, table, ofsmTableTestCase.numStates, ofsmTableTestCase.maxInput, ofsmTableTestCase.maxOutput),
			"OFSMTable::toFsm(const string & name): Each FsmNode from constructed Fsm has corresponding OFSMTableRow in table.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllOFSMRows(fsm, table, ofsmTableTestCase.numStates, ofsmTableTestCase.maxInput, ofsmTableTestCase.maxOutput),
			"OFSMTable::toFsm(const string & name): Each OFSMTableRow from the OFSMTable has a corresponding FsmNode in the constructed Fsm");
	}

	// Generation of random FSMs. These are transformed to OFSMTables, which are minimised with OFSMTable::next(). 
	// The last table is transformed to a FSM with toFsm(). 
	{
		int maxInput = 3;
		int maxOutput = 3;
		int maxState = 10;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<Fsm> originalFsm = Fsm::createRandomFsm("original", maxInput, maxOutput, maxState, pl);
		Fsm ofsm = originalFsm->transformToObservableFSM();
		shared_ptr<OFSMTable> table = make_shared<OFSMTable>(OFSMTable(ofsm.getNodes(), maxInput, maxOutput, pl));		
		shared_ptr<OFSMTable> next = table->next();
		while (next != nullptr) {
			table = next;
			next = next->next();
		}

		Fsm fsm = table->toFsm("");
		int numStates = ofsm.getNodes().size();


		fsmlib_assert("TC-OFSMTable-NNNN",
			fsm.getNodes().size() == table->maxClassId() + 1,
			"OFSMTable::toFsm(const string & name) creates Fsm with correct number of FsmNodes "
			"(equal to the number of generated classes)");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllFsmNodes(fsm, table, numStates, maxInput, maxOutput),
			"OFSMTable::toFsm(const string & name): Each FsmNode from constructed Fsm has corresponding OFSMTableRow in table.");

		fsmlib_assert("TC-OFSMTable-NNNN",
			matchAllOFSMRows(fsm, table, numStates, maxInput, maxOutput),
			"OFSMTable::toFsm(const string & name): Each OFSMTableRow from the OFSMTable has a corresponding FsmNode in the constructed Fsm");
	}
}

//===================================== FsmLabel Tests ===================================================

void testFsmLabelOperatorLessThan() {
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(1, 0, pl) < FsmLabel(1, 0, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns false if label1 == label2");
		fsmlib_assert("TC-FsmLabel-NNNN",
			(FsmLabel(1, 0, pl) < FsmLabel(1, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns true if inputs "
			"are the same but label1.output is less than label2.output");
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(1, 0, pl) < FsmLabel(0, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns false if label1.input > label2.input");
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(1, 1, pl) < FsmLabel(1, 0, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns false if inputs are the same but "
			"label1.output >= label2.output");
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(1, 1, pl) < FsmLabel(1, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns false if label1 == label2");
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(1, 1, pl) < FsmLabel(0, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns false if label1.input > label2.input");
		fsmlib_assert("TC-FsmLabel-NNNN",
			(FsmLabel(0, 1, pl) < FsmLabel(1, 0, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns true if label1.input < label2.input");
		fsmlib_assert("TC-FsmLabel-NNNN",
			(FsmLabel(0, 1, pl) < FsmLabel(1, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns true if label1.input < label2.input");
		fsmlib_assert("TC-FsmLabel-NNNN",
			not (FsmLabel(0, 1, pl) < FsmLabel(0, 1, pl)),
			"operator<(FsmLabel const & label1, FsmLabel const & label2) returns true if label1 == label2");
	}
}

//===================================== FsmNode Tests ===================================================

// tests FsmNode::addTransition(std::shared_ptr<FsmTransition> transition)
void testFsmNodeAddTransition() {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	// At the beginning n0 has no transition. A transition from n0 to n0 labeled with 0/1 is added to n0.
	shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);

	// original state
	int originalColor = n0->getColor();
	int originalID = n0->getId();
	string originalName = n0->getName();
	auto originalPair = n0->getPair();
	vector<string> originalSatisfied = n0->getSatisfied();

	shared_ptr<FsmTransition> tr0 = make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 1, pl));
	n0->addTransition(tr0);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions().size() == 1		
		&& n0->getTransitions().at(0) == tr0,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) adds transition to the node "
		"if it doesn't have a transition with equal target and label already.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");

	vector<shared_ptr<FsmTransition>> oldTransitions = n0->getTransitions();

	// Add a new transition which has the same label and target as tr0 (target = n0, label=0/1)
	shared_ptr<FsmTransition> tr1 = make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 1, pl));
	n0->addTransition(tr1);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions().size() == 1
		&& n0->getTransitions() == oldTransitions,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) ignores transition "
		"if node already has a transition with equal target and label.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");

	// Add transition from n0 to n1 with label 1/1.
	shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
	shared_ptr<FsmTransition> tr2 = make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl));
	n0->addTransition(tr2);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions().size() == 2
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr0) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr2) != n0->getTransitions().cend(),
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) adds transition to the node "
		"if it doesn't have a transition with equal target and label already.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");

	// Add transition from n0 to n1 with label 0/1.	
	shared_ptr<FsmTransition> tr3 = make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl));
	n0->addTransition(tr3);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions().size() == 3
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr0) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr2) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr3) != n0->getTransitions().cend(),
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) adds transition to the node "
		"if it doesn't have a transition with equal target and label already.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");

	oldTransitions = n0->getTransitions();

	// Add transition from n0 to n1 with label 1/1.	(Already contained in n0.transitions)
	shared_ptr<FsmTransition> tr4 = make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl));
	n0->addTransition(tr4);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions() == oldTransitions,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) ignores transition "
		"if node already has a transition with equal target and label.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");

	// Add transition from n0 to n0 with label 1/1.	(new transition)
	shared_ptr<FsmTransition> tr5 = make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(1, 1, pl));
	n0->addTransition(tr5);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getTransitions().size() == 4
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr0) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr2) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr3) != n0->getTransitions().cend()
		&& find(n0->getTransitions().cbegin(), n0->getTransitions().cend(), tr5) != n0->getTransitions().cend(),
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) adds transition to the node "
		"if it doesn't have a transition with equal target and label already.");

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->getColor() == originalColor
		&& n0->getId() == originalID
		&& n0->getName() == originalName
		&& n0->getPair() == originalPair
		&& n0->getSatisfied() == originalSatisfied,
		"FsmNode::addTransition(std::shared_ptr<FsmTransition> transition) changes only transitions");
}

// tests FsmNode::apply(const InputTrace& itrc, bool markAsVisited)
void testFsmNodeApply() {
	// n0 has no transition. 
	// inTrc = [1,2]
	// => inTrc starts with input for which n0 has no transition.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		InputTrace inTrc(vector<int>{1, 2}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		OutputTree expected(make_shared<TreeNode>(), inTrc, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree == expected,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) is the empty OutputTree if n0 has "
			"no outgoing transitions.");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited(),
			"Each expected FsmNode has been visited.");
	}

	// n0 --1/1--> n1 
	// inTrc = [2]
	// => inTrc starts with input for which n0 has no transition.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{2}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		OutputTree expected(make_shared<TreeNode>(), inTrc, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree == expected,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) is the empty OutputTree if itrc starts with some input "
			"that differs from the input of every transition of n0.");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited(),
			"Each expected FsmNode has been visited.");

		fsmlib_assert("TC-FsmNode-NNNN",
			not n1->hasBeenVisited(),
			"Only expected FsmNodes have been visited.");
	}

	// n0 --1/1--> n1 
	// inTrc = []
	// => inTrc is empty.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		OutputTree expected(make_shared<TreeNode>(), inTrc, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree == expected,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) is the empty OutputTree if itrc is empty.");

		fsmlib_assert("TC-FsmNode-NNNN",
			not n0->hasBeenVisited()
			&& not n1->hasBeenVisited(),
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) visits no FsmNode if itrc is empty.");
	}

	// n0 --1/1--> n1; n0 --0/0--> n0
	// inTrc = [0,1]
	// => inTrc matches exactly one path starting at n0 (completely).
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		InputTrace inTrc(vector<int>{0, 1}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		vector<OutputTrace> outTrcs = outTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 1,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		OutputTrace expectedTrc0(vector<int>{0, 1}, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited(),
			"Each expected FsmNode has been visited.");
	}

	// n0 --0/1--> n1; n0 --0/0--> n0; n1 --1/1--> n1
	// inTrc = [0,1]
	// => inTrc matches more than one path (nondeterminism of n0) completely.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --0/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{0, 1}, pl);
		OutputTree outTree = n0->apply(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		vector<OutputTrace> outTrcs = outTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 2,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		OutputTrace expectedTrc0(vector<int>{1, 1}, pl);
		OutputTrace expectedTrc1(vector<int>{0}, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend()
			&& find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc1) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");

		fsmlib_assert("TC-FsmNode-NNNN",
			not n0->hasBeenVisited()
			&& not n1->hasBeenVisited(),
			"If markAsVisited is set to false no FsmNode is marked as visited by "
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited)");
	}

	// n0 --0/0--> n1; n1 --1/1--> n1; n0 --1/1--> n2;
	// inTrc = [0,1,1] / [0,1,1,2]
	// => inTrc matches only one path and doesn't visit each FsmNode. In the [0,1,1,2]-case only a prefix can be applied.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --0/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --1/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{0, 1, 1}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		vector<OutputTrace> outTrcs = outTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 1,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		const OutputTrace expectedTrc0(vector<int>{0, 1, 1}, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited(),
			"Each expected FsmNode has been visited.");

		fsmlib_assert("TC-FsmNode-NNNN",
			not n2->hasBeenVisited(),
			"Only expected FsmNodes have been visited.");

		// now inTrc = [0,1,1,2] is applied
		inTrc = InputTrace(vector<int>{0, 1, 1, 2}, pl);
		OutputTree newOutTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			newOutTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		outTrcs = newOutTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 1,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		// expectedTrc0 is still [0,1,1]
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited(),
			"Each expected FsmNode has been visited.");

		fsmlib_assert("TC-FsmNode-NNNN",
			not n2->hasBeenVisited(),
			"Only expected FsmNodes have been visited.");
	}

	// n0 --0/0--> n1; n1 --1/1--> n1; n0 --0/1--> n2; n2 --1/3--> n2; n2 --1/2--> n1;
	// inTrc = [0,1,1]
	// => inTrc matches several paths.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --0/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --0/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(0, 1, pl)));
		// n2 --1/3--> n2
		n2->addTransition(make_shared<FsmTransition>(n2, n2, make_shared<FsmLabel>(1, 3, pl)));
		// n2 --1/2--> n1
		n2->addTransition(make_shared<FsmTransition>(n2, n1, make_shared<FsmLabel>(1, 2, pl)));
		InputTrace inTrc(vector<int>{0, 1, 1}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		vector<OutputTrace> outTrcs = outTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 4,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		const OutputTrace expectedTrc0(vector<int>{0, 1, 1}, pl);
		const OutputTrace expectedTrc1(vector<int>{1, 2, 1}, pl);
		const OutputTrace expectedTrc2(vector<int>{1, 3, 2}, pl);
		const OutputTrace expectedTrc3(vector<int>{1, 3, 3}, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend()
			&& find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc1) != outTrcs.cend()
			&& find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc2) != outTrcs.cend()
			&& find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc3) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited()
			&& n2->hasBeenVisited(),
			"Each expected FsmNode has been visited.");
	}

	// n0 --1/0--> n1; n1 --2/5--> n4; n0 --1/1--> n2; n2 --1/2--> n3; n3 --2/4--> n5;
	// itrc: [1,1,2]
	// itrc can be fully applied on one path but only partially on the other path of the Fsm.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n3 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n4 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n5 = make_shared<FsmNode>(2, pl);
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --2/5--> n4
		n1->addTransition(make_shared<FsmTransition>(n1, n4, make_shared<FsmLabel>(2, 5, pl)));
		// n0 --1/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 1, pl)));
		// n2 --1/2--> n3
		n2->addTransition(make_shared<FsmTransition>(n2, n3, make_shared<FsmLabel>(1, 2, pl)));
		// n3 --2/4--> n5
		n3->addTransition(make_shared<FsmTransition>(n3, n5, make_shared<FsmLabel>(2, 4, pl)));
		InputTrace inTrc(vector<int>{1,1,2}, pl);
		OutputTree outTree = n0->apply(inTrc, true);
		fsmlib_assert("TC-FsmNode-NNNN",
			outTree.getInputTrace() == inTrc,
			"FsmNode::apply(const InputTrace& itrc, bool markAsVisited) returns OutputTree which has itrc set as InputTrace.");

		vector<OutputTrace> outTrcs = outTree.getOutputTraces();
		fsmlib_assert("TC-FsmNode-NNNN",
			outTrcs.size() == 2,
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains the expected number of OutputTraces");

		const OutputTrace expectedTrc0(vector<int>{0}, pl);
		const OutputTrace expectedTrc1(vector<int>{1,2,4}, pl);
		
		fsmlib_assert("TC-FsmNode-NNNN",
			find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc0) != outTrcs.cend()
			&& find(outTrcs.cbegin(), outTrcs.cend(), expectedTrc1) != outTrcs.cend(),
			"Result of FsmNode::apply(const InputTrace& itrc, bool markAsVisited) contains each expected OutputTrace");
	}
}

// tests FsmNode::after(const InputTrace& itrc) and FsmNode::after(const vector<int>& itrc)
void testFsmNodeAfter1() {
	// n0 has no transition. 
	// inTrc = [1,2]
	// => inTrc starts with input for which n0 has no transition.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		InputTrace inTrc(vector<int>{1, 2}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.empty(),
			"Result of FsmNode::after(const InputTrace& itrc) is empty if n0 has no outgoing transition "
			"and itrc can't be fully applied");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --1/1--> n1 
	// inTrc = [2]
	// => inTrc starts with input for which n0 has no transition.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{2}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.empty(),
			"Result of FsmNode::after(const InputTrace& itrc) is empty if n0 has no outgoing transition "
			"for the first input of itrc.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --1/1--> n1 
	// inTrc = []
	// => inTrc is empty.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		InputTrace inTrc(vector<int>{}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1
			&& reached.count(n0) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains only n0 if itrc is empty.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --1/1--> n1, n0 --0/0--> n0 
	// inTrc = [1]
	// => inTrc applied to n0 reaches only n1.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		InputTrace inTrc(vector<int>{1}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",			
			reached.count(n1) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");		
	}

	// n0 --1/1--> n1, n0 --0/0--> n0 
	// inTrc = [0,1]
	// => inTrc applied to n0 reaches only n1.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		InputTrace inTrc(vector<int>{0,1}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n1) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --1/1--> n1, n0 --0/0--> n0 
	// inTrc = [1,1]
	// => inTrc can't be completely applied to n0 so no FsmNode is reached.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		InputTrace inTrc(vector<int>{1, 1}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.empty(),
			"Result of FsmNode::after(const InputTrace& itrc) is empty if itrc can't be completely applied to n0.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --0/1--> n1, n1 --1/2--> n2 
	// inTrc = [0,0,1]
	// => prefix and suffix of inTrc can be applied, but inTrc can't be completely applied.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --0/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl)));
		// n1 --1/2--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(1, 2, pl)));
		InputTrace inTrc(vector<int>{0, 0, 1}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.empty(),
			"Result of FsmNode::after(const InputTrace& itrc) is empty if itrc can't be completely applied to n0.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --0/0--> n1, n0 --0/1--> n2 
	// inTrc = [0]
	// => More than one FsmNode can be reached with itrc.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --0/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --0/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(0, 1, pl)));
		InputTrace inTrc(vector<int>{0}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n1) == 1
			&& reached.count(n2) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 2,
			"Result of FsmNode::after(const InputTrace& itrc) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}

	// n0 --0/0--> n1, n0 --0/1--> n2; n0 --0/3--> n3; n1 --1/2--> n2; n1 --1/3--> n0;
	// inTrc = [0, 1]
	// => Two FsmNodes can be reached with itrc.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);
		// n0 --0/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --0/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(0, 1, pl)));
		// n0 --0/3--> n3
		n0->addTransition(make_shared<FsmTransition>(n0, n3, make_shared<FsmLabel>(0, 3, pl)));
		// n1 --1/2--> n2 
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(1, 2, pl)));
		// n1 --1/3--> n0
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(1, 3, pl)));
		InputTrace inTrc(vector<int>{0, 1}, pl);
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(inTrc);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n0) == 1
			&& reached.count(n2) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 2,
			"Result of FsmNode::after(const InputTrace& itrc) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");


		// now apply inTrc = [0,1,0]
		inTrc = InputTrace(vector<int>{0, 1, 0}, pl);
		reached = n0->after(inTrc);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n1) == 1
			&& reached.count(n2) == 1
			&& reached.count(n3) == 1,
			"Result of FsmNode::after(const InputTrace& itrc) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 3,
			"Result of FsmNode::after(const InputTrace& itrc) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(inTrc.get()),
			"FsmNode::after(const InputTrace& itrc) == FsmNode::after(const vector<int>& itrc)");
	}
}

// tests FsmNode::after(const std::shared_ptr<TraceSegment> seg)
void testFsmNodeAfter2() {
	// n0 has no outgoing transition 
	// seg = [1,2] , 
	// seg.prefix = std::string::npos => seg can't be completely applied
	// seg.prefix = 0 => only initial node is reached.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		//shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		//shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		//shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);

		// n0 --0/0--> n1
		//n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		
		shared_ptr<TraceSegment> seg = make_shared<TraceSegment>(make_shared<vector<int>>(vector<int>{1,2}));
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(seg);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.empty(),
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) is empty if seg can't be completely applied to n0.");

		seg = make_shared<TraceSegment>(make_shared<vector<int>>(vector<int>{1, 2}), 0);
		reached = n0->after(seg);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1
			&& reached.count(n0) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) is initial node n0 if seg.prefix == 0");
	}

	// n0 --0/1--> n1; n1 --1/2--> n2; n0 --0/2--> n2; n2 --1/1--> n2
	// seg = [0,1] , 
	// seg.prefix = std::string::npos => seg can be completely applied.
	// seg.prefix = 1 => only prefix of seg can be applied
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		//shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);

		// n0 --0/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl)));
		// n1 --1/2--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(1, 2, pl)));
		// n0 --0/2--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(0, 2, pl)));
		// n2 --1/1--> n2
		n2->addTransition(make_shared<FsmTransition>(n2, n2, make_shared<FsmLabel>(1, 1, pl)));

		shared_ptr<TraceSegment> seg = make_shared<TraceSegment>(make_shared<vector<int>>(vector<int>{0, 1}));
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(seg);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n2) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(*seg->get()),
			"FsmNode::after(const std::shared_ptr<TraceSegment> seg) == FsmNode::after(const vector<int>& itrc) "
			"if seg.prefix == string::npos");

		seg->setPrefix(1);
		reached = n0->after(seg);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n1) == 1
			&& reached.count(n2) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 2,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");
	}

	// n0 --1/0--> n1; n1 --2/1--> n2; n2 --1/1--> n3;
	// seg = [1,2,1] , 
	// seg.prefix = npos => seg can be completely applied.
	// seg.prefix = 3 => seg can be completely applied (seg.prefix equals size of seg)
	// seg.prefix = 2 => only prefix of seg can be applied
	// seg.prefix = 0 => nothing of seg is applied. Only initial node n0 is reached.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);

		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --2/1--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(2, 1, pl)));
		// n2 --1/1--> n3
		n2->addTransition(make_shared<FsmTransition>(n2, n3, make_shared<FsmLabel>(1, 1, pl)));

		// prefix = string::npos
		shared_ptr<TraceSegment> seg = make_shared<TraceSegment>(make_shared<vector<int>>(vector<int>{1, 2, 1}));
		unordered_set<shared_ptr<FsmNode>> reached = n0->after(seg);
		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n3) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached == n0->after(*seg->get()),
			"FsmNode::after(const std::shared_ptr<TraceSegment> seg) == FsmNode::after(const vector<int>& itrc) "
			"if seg.prefix == string::npos");

		// prefix = 3
		seg->setPrefix(3);
		reached = n0->after(seg);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n3) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");

		// prefix = 2
		seg->setPrefix(2);
		reached = n0->after(seg);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n2) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");

		// prefix = 0
		seg->setPrefix(0);
		reached = n0->after(seg);

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.count(n0) == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains each expected FsmNode.");

		fsmlib_assert("TC-FsmNode-NNNN",
			reached.size() == 1,
			"Result of FsmNode::after(const std::shared_ptr<TraceSegment> seg) contains only expected FsmNodes.");
	}
}

// Checks if row contains transitions for other inputs (between 0 and maxInput) than those contained in except.
bool hasUnexpectedTransitionInRow(shared_ptr<DFSMTableRow> row, int maxInput, unordered_set<int> & except) {
	for (int i = 0; i <= maxInput; ++i) {
		if (except.count(i) == 1) continue;
		if (row->geti2postSection().at(i) != -1 || row->getioSection().at(i) != -1) return false;
	}
	return true;
}

// tests FsmNode::getDFSMTableRow(const int maxInput)
void testFsmNodeGetDFSMTableRow() {
	// n0 has no outgoing transition.
	// maxInput = 3
	{
		const int maxInput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);

		shared_ptr<DFSMTableRow> row = n0->getDFSMTableRow(maxInput);
		unordered_set<int> except;
		fsmlib_assert("TC-FsmNode-NNNN",
			hasUnexpectedTransitionInRow(row, maxInput, except),
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains no transition if n0 has no outgoing transition.");
	}

	// n0 --1/1--> n0
	// maxInput = 3
	{
		const int maxInput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);

		// n0 --1/1--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(1, 1, pl)));

		shared_ptr<DFSMTableRow> row = n0->getDFSMTableRow(maxInput);
		unordered_set<int> except{1};

		fsmlib_assert("TC-FsmNode-NNNN",
			row->geti2postSection().at(1) == 0 && row->getioSection().at(1) == 1,
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains each expected transition with the right values.");

		fsmlib_assert("TC-FsmNode-NNNN",
			hasUnexpectedTransitionInRow(row, maxInput, except),
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains no transition if n0 has no outgoing transition.");
	}

	// n0 --0/0--> n0; n0 --1/2--> n1
	// maxInput = 2
	{
		const int maxInput = 2;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);

		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --1/2--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 2, pl)));

		shared_ptr<DFSMTableRow> row = n0->getDFSMTableRow(maxInput);
		unordered_set<int> except{ 0,1 };

		fsmlib_assert("TC-FsmNode-NNNN",
			row->geti2postSection().at(0) == 0 && row->getioSection().at(0) == 0
			&& row->geti2postSection().at(1) == 1 && row->getioSection().at(1) == 2,
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains each expected transition with the right values.");

		fsmlib_assert("TC-FsmNode-NNNN",
			hasUnexpectedTransitionInRow(row, maxInput, except),
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains no transition if n0 has no outgoing transition.");
	}

	// n0 --0/0--> n0; n0 --1/1--> n0; n0 --2/0--> n1; n0 --3/4--> n2;
	// maxInput = 3
	{
		const int maxInput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);

		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --1/1--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --2/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(2, 0, pl)));
		// n0 --3/4--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(3, 4, pl)));

		shared_ptr<DFSMTableRow> row = n0->getDFSMTableRow(maxInput);
		unordered_set<int> except{ 0,1,2,3 };

		fsmlib_assert("TC-FsmNode-NNNN",
			row->geti2postSection().at(0) == 0 && row->getioSection().at(0) == 0
			&& row->geti2postSection().at(1) == 0 && row->getioSection().at(1) == 1
			&& row->geti2postSection().at(2) == 1 && row->getioSection().at(2) == 0
			&& row->geti2postSection().at(3) == 2 && row->getioSection().at(3) == 4,
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains each expected transition with the right values.");

		fsmlib_assert("TC-FsmNode-NNNN",
			hasUnexpectedTransitionInRow(row, maxInput, except),
			"FsmNode::getDFSMTableRow(const int maxInput): Constructed row contains no transition if n0 has no outgoing transition.");
	}
}

// tests FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w)
// Positive Case (InputTrace returned)
void testFsmNodeDistinguishedPositive() {
	// w contains only one Trace ([1,2]). This Trace produces different sets of OutputTraces (distinguishes both FsmNodes.)
	// n0 --1/1--> n0; n0 --2/2--> n1
	// o_n0 --1/1--> o_n0; o_n0 --2/2--> o_n1; o_n0 --2/3--> o_n2
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		InputTrace itrc(vector<int>{1, 2}, pl);

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);
		w->addToRoot(itrc.get());

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --2/2--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(2, 2, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> o_n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> o_n2 = make_shared<FsmNode>(2, pl);
		// o_n0 --1/1--> o_n0
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n0, make_shared<FsmLabel>(1, 1, pl)));
		// o_n0 --2/2--> o_n1
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n1, make_shared<FsmLabel>(2, 2, pl)));
		// o_n0 --2/3--> o_n2
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n2, make_shared<FsmLabel>(2, 3, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			*n0->distinguished(o_n0, w) == itrc
			&& *o_n0->distinguished(n0, w) == itrc,
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns pointer to distinguishing trace.");
	}

	// w contains [1,1] and [1,2]. Only one of the Traces distinguishes the FsmNodes.
	// n0 --1/1--> n1; n1 --1/1--> n1; n1 --2/0--> n0
	// o_n0 --1/1--> o_n0; o_n0 --2/1--> o_n1 
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> trc0{ 1,1 };
		vector<int> trc1{ 1,2 };

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);
		w->addToRoot(trc0);
		w->addToRoot(trc1);

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --2/0--> n0
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(2, 0, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> o_n1 = make_shared<FsmNode>(1, pl);
		// o_n0 --1/1--> o_n0
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n0, make_shared<FsmLabel>(1, 1, pl)));
		// o_n0 --2/1--> o_n1
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n1, make_shared<FsmLabel>(2, 1, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->distinguished(o_n0, w)->get() == trc1
			&& o_n0->distinguished(n0, w)->get() == trc1,
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns pointer to distinguishing trace.");
	}
	// w contains [1,1], [1,2] and [3,4]. More than one Trace distinguishes the FsmNodes.
	// n0 --1/1--> n1; n1 --1/0--> n1; n1 --2/1--> n0; n0 --3/3--> n2
	// o_n0 --1/1--> o_n1; o_n1 --1/0--> o_n0; o_n0 --3/3--> o_n2; o_n2 --4/3--> o_n2
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> trc0{ 1,1 };
		vector<int> trc1{ 1,2 };
		vector<int> trc2{ 3,4 };

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);
		w->addToRoot(trc0);
		w->addToRoot(trc1);
		w->addToRoot(trc2);

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --1/0--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --2/1--> n0
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(2, 1, pl)));
		// n0 --3/3--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(3, 3, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> o_n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> o_n2 = make_shared<FsmNode>(2, pl);
		// o_n0 --1/1--> o_n1
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n1, make_shared<FsmLabel>(1, 1, pl)));
		// o_n1 --1/0--> o_n0
		o_n1->addTransition(make_shared<FsmTransition>(o_n1, o_n0, make_shared<FsmLabel>(1, 0, pl)));
		// o_n0 --3/3--> o_n2
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n2, make_shared<FsmLabel>(3, 3, pl)));
		// o_n2 --4/3--> o_n2
		o_n2->addTransition(make_shared<FsmTransition>(o_n2, o_n2, make_shared<FsmLabel>(4, 3, pl)));
		fsmlib_assert("TC-FsmNode-NNNN",
			(n0->distinguished(o_n0, w)->get() == trc1 || n0->distinguished(o_n0, w)->get() == trc2)
			&& (o_n0->distinguished(n0, w)->get() == trc1 || o_n0->distinguished(n0, w)->get() == trc2),
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns pointer to distinguishing trace.");
	}
}

// tests FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w)
// Negative Case (nullptr returned)
void testFsmNodeDistinguishedNegative() {
	// w is empty. 
	// n0 --0/0--> n0
	// o_n0 --0/1--> o_n0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		// o_n0 --0/1--> o_n0
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n0, make_shared<FsmLabel>(0, 1, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->distinguished(o_n0, w) == nullptr
			&& o_n0->distinguished(n0, w) == nullptr,
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns nullptr "
			"if both FsmNodes can't be distinguished by w.");
	}

	// w contains only one Trace ([1,2]). This Trace doesn't distinguish the FsmNodes.
	// n0 --1/1--> n0; n0 --2/0--> n1
	// o_n0 --1/1--> o_n0; o_n0 --2/0--> o_n0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		InputTrace itrc(vector<int>{1, 2}, pl);

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);
		w->addToRoot(itrc.get());

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --2/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(2, 0, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		// o_n0 --1/1--> o_n0
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n0, make_shared<FsmLabel>(1, 1, pl)));
		// o_n0 --2/0--> o_n0
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n0, make_shared<FsmLabel>(2, 0, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->distinguished(o_n0, w) == nullptr
			&& o_n0->distinguished(n0, w) == nullptr,
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns nullptr "
			"if both FsmNodes can't be distinguished by w.");
	}

	// w contains [1,1] and [1,2]. 
	// n0 --1/1--> n1; n0 --1/0--> n1; n1 --1/1--> n1; n1 --2/2--> n2
	// o_n0 --1/1--> o_n1; o_n0 --1/0--> o_n2; o_n1 --1/1--> o_n1; o_n1 --2/2--> o_n2; o_n2 --1/1--> o_n2; o_n2 --2/2--> o_n1
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<int> trc0{ 1,1 };
		vector<int> trc1{ 1,2 };

		// construction of w.
		shared_ptr<Tree> w = make_shared<Tree>(make_shared<TreeNode>(), pl);
		w->addToRoot(trc0);
		w->addToRoot(trc1);

		// construction of first Fsm
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --2/2--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(2, 2, pl)));

		// construction of other Fsm
		shared_ptr<FsmNode> o_n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> o_n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> o_n2 = make_shared<FsmNode>(2, pl);
		// o_n0 --1/1--> o_n1
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n1, make_shared<FsmLabel>(1, 1, pl)));
		// o_n0 --1/0--> o_n2
		o_n0->addTransition(make_shared<FsmTransition>(o_n0, o_n2, make_shared<FsmLabel>(1, 0, pl)));
		// o_n1 --1/1--> o_n1
		o_n1->addTransition(make_shared<FsmTransition>(o_n1, o_n1, make_shared<FsmLabel>(1, 1, pl)));
		// o_n1 --2/2--> o_n2
		o_n1->addTransition(make_shared<FsmTransition>(o_n1, o_n2, make_shared<FsmLabel>(2, 2, pl)));
		// o_n2 --1/1--> o_n2
		o_n2->addTransition(make_shared<FsmTransition>(o_n2, o_n2, make_shared<FsmLabel>(1, 1, pl)));
		// o_n2 --2/2--> o_n1
		o_n2->addTransition(make_shared<FsmTransition>(o_n2, o_n1, make_shared<FsmLabel>(2, 2, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			n0->distinguished(o_n0, w) == nullptr
			&& o_n0->distinguished(n0, w) == nullptr,
			"FsmNode::distinguished(const shared_ptr<FsmNode> otherNode, shared_ptr<Tree> w) returns nullptr "
			"if both FsmNodes can't be distinguished by w.");
	}
}

// Calculates the distinguishing trace for each pair of FsmNodes in nodes with 
// FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode, const vector<shared_ptr<PkTable>>& pktblLst, const int maxInput).
// Checks if the calculated trace really distinguishes the current pair.
// If the calculated trace is empty, both FsmNodes have to be equivalent.
bool checkCalcDistinguishingTraceForAllPairs(const vector<shared_ptr<FsmNode>> &nodes,
	const vector<shared_ptr<PkTable>> &pkTables, const int maxInput) {
	for (shared_ptr<FsmNode> n : nodes) {
		for (shared_ptr<FsmNode> o_n : nodes) {
			InputTrace distTrc = n->calcDistinguishingTrace(o_n, pkTables, maxInput);
			if (distTrc.size() > 0) {
				if (not n->distinguished(o_n, distTrc.get())) {
					return false;
				}					
			}
			// distTrace.size() == 0 => n ~ o_n
			else {
				if (pkTables.at(pkTables.size() - 1)->getClass(n->getId()) != pkTables.at(pkTables.size() - 1)->getClass(o_n->getId())) {
					return false;
				}
			}
		}
	}
	return true;
}

// tests FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode, const vector<shared_ptr<PkTable>>& pktblLst, const int maxInput)
void testFsmNodeCalcDistinguishingTrace1() {
	// using the Dfsm specified at "../../../resources/TC-FsmNode-calcDistinguishingTrace1.fsm" (already minimal, completely specified)
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Dfsm dfsm("../../../resources/TC-FsmNode-calcDistinguishingTrace1.fsm", pl, "m1");
		const int maxInput = 4;

		dfsm.calcPkTables();

		fsmlib_assert("TC-FsmNode-NNNN",
			checkCalcDistinguishingTraceForAllPairs(dfsm.getNodes(), dfsm.getPktblLst(), maxInput),
			"Result of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode, const vector<shared_ptr<PkTable>>& pktblLst, const int maxInput)"
			" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
	}

	// using the Dfsm specified at "../../../resources/TC-FsmNode-calcDistinguishingTrace1.fsm" (not minimal, completely specified)
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Dfsm dfsm("../../../resources/TC-FsmNode-calcDistinguishingTrace2.fsm", pl, "m1");
		const int maxInput = 2;

		dfsm.calcPkTables();

		fsmlib_assert("TC-FsmNode-NNNN",
			checkCalcDistinguishingTraceForAllPairs(dfsm.getNodes(), dfsm.getPktblLst(), maxInput),
			"Result of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode, const vector<shared_ptr<PkTable>>& pktblLst, const int maxInput)"
			" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
	}

	// calculating and using random DFSMs. Testing each pair of FsmNodes.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		for (int maxNodes = 1; maxNodes <= 13; ++maxNodes) {
			for (int io = 0; io <= 3; ++io) {
				Dfsm dfsm("m", maxNodes, io, io, pl);
				dfsm.calcPkTables();
				fsmlib_assert("TC-FsmNode-NNNN",
					checkCalcDistinguishingTraceForAllPairs(dfsm.getNodes(), dfsm.getPktblLst(), io),
					"Result of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode, const vector<shared_ptr<PkTable>>& pktblLst, const int maxInput)"
					" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
			}
		}
	}
}

// This is a class only used for tests. It can be used to access protected attributes and methods.
class FsmTest : public Fsm{
public:
	FsmTest(const std::string& fname, const std::shared_ptr<FsmPresentationLayer> presentationLayer,
		const std::string& fsmName): Fsm(fname, presentationLayer, fsmName){
	}

	FsmTest(const std::string & fsmName,
		const int maxInput,
		const int maxOutput,
		const std::vector<std::shared_ptr<FsmNode>> lst,
		const std::shared_ptr<FsmPresentationLayer> presentationLayer): Fsm(fsmName, maxInput, maxOutput, lst, presentationLayer){}

	std::vector<std::shared_ptr<OFSMTable>> getOfsmTableLst(){
		return ofsmTableLst;
	}

	void calcOFSMTables() {
		Fsm::calcOFSMTables();
	}

	static std::shared_ptr<FsmTest>
		createRandomFsmTest(const std::string & fsmName, const int maxInput, const int maxOutput, const int maxState,
			const std::shared_ptr<FsmPresentationLayer> presentationLayer, const unsigned seed = 0) {
		shared_ptr<Fsm> fsm = Fsm::createRandomFsm(fsmName, maxInput, maxOutput, maxState, presentationLayer, seed);
		return make_shared<FsmTest>(fsm->getName(), fsm->getMaxInput(), fsm->getMaxOutput(), fsm->getNodes(), fsm->getPresentationLayer());
	}
	//Fsm(name + "_O", maxInput, maxOutput, nodeLst, obsPl)
	FsmTest transformToObservableFSM() const{
		Fsm fsm = Fsm::transformToObservableFSM();
		return FsmTest(fsm.getName(), fsm.getMaxInput(), fsm.getMaxOutput(), fsm.getNodes(), fsm.getPresentationLayer());
	}
};

// Calculates the distinguishing trace for each pair of FsmNodes in nodes with 
// calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput).
// Checks if the calculated trace really distinguishes the current pair.
// If the calculated trace is empty, both FsmNodes have to be equivalent.
bool checkCalcDistinguishingTraceForAllPairs(const vector<shared_ptr<FsmNode>> &nodes,
	const vector<shared_ptr<OFSMTable>> &ofsmTables, const int maxInput, const int maxOutput) {
	for (shared_ptr<FsmNode> n : nodes) {
		for (shared_ptr<FsmNode> o_n : nodes) {
			InputTrace distTrc = n->calcDistinguishingTrace(o_n, ofsmTables, maxInput, maxOutput);
			if (distTrc.size() > 0) {
				if (not n->distinguished(o_n, distTrc.get())) {
					return false;
				}
			}
			// distTrace.size() == 0 => n ~ o_n
			else {
				S2CMap s2c = ofsmTables.at(ofsmTables.size() - 1)->getS2C();
				if (s2c.at(n->getId()) != s2c.at(o_n->getId())) {
					return false;
				}
			}
		}
	}
	return true;
}

// tests FsmNode::calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput);
void testFsmNodeCalcDistinguishingTrace2() {
	// using FSM specified at "../../../resources/TC-FsmNode-calcDistinguishingTrace3.fsm" (completely specified and observable) 
	{
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		FsmTest fsm_t("../../../resources/TC-FsmNode-calcDistinguishingTrace3.fsm", pl, "M");
		fsm_t.calcOFSMTables();
		
		fsmlib_assert("TC-FsmNode-NNNN",
			checkCalcDistinguishingTraceForAllPairs(fsm_t.getNodes(), fsm_t.getOfsmTableLst(), maxInput, maxOutput),
			"Result of "
			"calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput)"
			" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
	}

	// using test case from the script (page 54 / example 3) (specified at "../../../resources/TC-FsmNode-calcDistinguishingTrace4.fsm")
	// This FSM is nondeterministic, observable and minimal
	{
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		FsmTest fsm_t("../../../resources/TC-FsmNode-calcDistinguishingTrace4.fsm", pl, "M");
		fsm_t.calcOFSMTables();

		fsmlib_assert("TC-FsmNode-NNNN",
			checkCalcDistinguishingTraceForAllPairs(fsm_t.getNodes(), fsm_t.getOfsmTableLst(), maxInput, maxOutput),
			"Result of "
			"calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput)"
			" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
	}

	// calculating and using random DFSMs. Testing each pair of FsmNodes.
	{
		

		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		for (int maxState = 0; maxState <= 6; ++maxState) {
			for (int io = 0; io <= 3; ++io) {

				FsmTest o_fsm_t = FsmTest::createRandomFsmTest("M", io, io, maxState, pl)->transformToObservableFSM();
				o_fsm_t.calcOFSMTables();

				fsmlib_assert("TC-FsmNode-NNNN",
					checkCalcDistinguishingTraceForAllPairs(o_fsm_t.getNodes(), o_fsm_t.getOfsmTableLst(), io, io),
					"Result of "
					"calcDistinguishingTrace(const std::shared_ptr<FsmNode> otherNode, const std::vector<std::shared_ptr<OFSMTable>>& ofsmTblLst, const int maxInput, const int maxOutput)"
					" is empty if both FsmNodes can't be distinguished. Otherwise it's a non-empty trace that distinguishes both FsmNodes.");
			}
		}
	}
}

// tests FsmNode::isObservable()
void testFsmNodeIsObservable() {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	// positive case:
	{
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		fsmlib_assert("TC-FsmNode-NNNN",
			n0->isObservable(),
			"FsmNode::isObservable() returns true if n0 has no outgoing transition.");

		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		fsmlib_assert("TC-FsmNode-NNNN",
			n0->isObservable(),
			"FsmNode::isObservable() returns true if each outgoing transition of n0 has a unique label.");

		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --0/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl)));
		fsmlib_assert("TC-FsmNode-NNNN",
			n0->isObservable(),
			"FsmNode::isObservable() returns true if each outgoing transition of n0 has a unique label.");
	}

	// negative case:
	{
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);

		// n0 --0/1--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 1, pl)));
		// n0 --0/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 1, pl)));
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			not n0->isObservable(),
			"FsmNode::isObservable() returns false if n0 has more than one outgoing transition with the same label.");

		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/1--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 1, pl)));

		fsmlib_assert("TC-FsmNode-NNNN",
			not n0->isObservable(),
			"FsmNode::isObservable() returns false if n0 has more than one outgoing transition with the same label.");
	}
}

// tests FsmNode::isDeterministic()
void testFsmNodeIsDeterministic() {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->isDeterministic(),
		"FsmNode::isDeterministic() returns true if n0 has at most one outgoing transition for each element in the input alphabet.");

	// n0 --0/1--> n0
	n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 1, pl)));

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->isDeterministic(),
		"FsmNode::isDeterministic() returns true if n0 has at most one outgoing transition for each element in the input alphabet.");

	shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
	// n0 --1/2--> n1
	n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 2, pl)));

	fsmlib_assert("TC-FsmNode-NNNN",
		n0->isDeterministic(),
		"FsmNode::isDeterministic() returns true if n0 has at most one outgoing transition for each element in the input alphabet.");

	shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
	// n0 --1/1--> n2
	n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 1, pl)));
	
	fsmlib_assert("TC-FsmNode-NNNN",
		not n0->isDeterministic(),
		"FsmNode::isDeterministic() returns false if n0 has more than one outgoing transition for some element in the input alphabet.");
}

//===================================== Fsm Traversal Tests ===================================================

//Fsm(const std::string & fsmName,
//	const int maxInput,
//	const int maxOutput,
//	const std::vector<std::shared_ptr<FsmNode>> lst,
//	const std::shared_ptr<FsmPresentationLayer> presentationLayer);

// tests the traversal of a Fsm with help of the Fsm::accept(FsmVisitor& v) method.
void testFsmAccept() {
	// Fsm contains only FsmNode n0 which has no outgoing transition.
	{
		const int maxInput = 2;
		const int maxOutput = 2;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		vector<shared_ptr<FsmNode>> nodes;
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		nodes.push_back(n0);
		Fsm fsm("M", maxInput, maxOutput, nodes, pl);
		FsmPrintVisitor v;
		//FsmSimVisitor v;
		//FsmOraVisitor v;
		fsm.accept(v);
		cout << endl;
		fsmlib_assert("TC-Fsm-NNNN",
			n0->hasBeenVisited(),
			"Fsm::accept(FsmVisitor& v) causes a complete traversal of the Fsm. Each FsmNode reachable from the initial node is "
			"marked as visited.");
	}

	// Fsm contains only FsmNode n0 which has only one outgoing transition to itself.
	// n0 --0/0--> n0
	{
		const int maxInput = 0;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		vector<shared_ptr<FsmNode>> nodes;
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		nodes.push_back(n0);

		Fsm fsm("M", maxInput, maxOutput, nodes, pl);
		//FsmPrintVisitor v;
		FsmSimVisitor v;
		//FsmOraVisitor v;
		fsm.accept(v);
		cout << endl;

		fsmlib_assert("TC-Fsm-NNNN",
			n0->hasBeenVisited(),
			"Fsm::accept(FsmVisitor& v) causes a complete traversal of the Fsm. Each FsmNode reachable from the initial node is "
			"marked as visited.");

	}

	// Each FsmNode of the Fsm is reachable from the inition node n0. The Fsm is non-observable.
	// n0 --0/0--> n0; n0 --0/0--> n1; n1 --1/1--> n1; n1 --1/1--> n0
	{
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		vector<shared_ptr<FsmNode>> nodes;
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --0/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(0, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --1/1--> n0
		n1->addTransition(make_shared<FsmTransition>(n1, n0, make_shared<FsmLabel>(1, 1, pl)));
		nodes.push_back(n0);
		nodes.push_back(n1);

		Fsm fsm("M", maxInput, maxOutput, nodes, pl);
		FsmPrintVisitor v;
		//FsmSimVisitor v;
		//FsmOraVisitor v;
		fsm.accept(v);
		cout << endl;

		fsmlib_assert("TC-Fsm-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited(),
			"Fsm::accept(FsmVisitor& v) causes a complete traversal of the Fsm. Each FsmNode reachable from the initial node is "
			"marked as visited.");

	}

	// n4 is not reachable from initial node n0. The Fsm is deterministic.
	// n0 --1/2--> n2; n0 --3/4--> n1; n1 --0/0--> n2; n2 --5/5--> n3; n3 --1/2--> n1; n3 --2/3--> n0;    n4
	{
		const int maxInput = 5;
		const int maxOutput = 5;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

		vector<shared_ptr<FsmNode>> nodes;
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);
		shared_ptr<FsmNode> n4 = make_shared<FsmNode>(4, pl);
		// n0 --1/2--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 2, pl)));
		// n0 --3/4--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(3, 4, pl)));
		// n1 --0/0--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(0, 0, pl)));
		// n2 --5/5--> n3
		n2->addTransition(make_shared<FsmTransition>(n2, n3, make_shared<FsmLabel>(5, 5, pl)));
		// n3 --1/2--> n1
		n3->addTransition(make_shared<FsmTransition>(n3, n1, make_shared<FsmLabel>(1, 2, pl)));
		// n3 --2/3--> n0
		n3->addTransition(make_shared<FsmTransition>(n3, n0, make_shared<FsmLabel>(2, 3, pl)));
		nodes.push_back(n0);
		nodes.push_back(n1);
		nodes.push_back(n2);
		nodes.push_back(n3);
		nodes.push_back(n4);

		Fsm fsm("M", maxInput, maxOutput, nodes, pl);
		//FsmPrintVisitor v;
		//FsmSimVisitor v;
		FsmOraVisitor v;
		fsm.accept(v);
		cout << endl;

		fsmlib_assert("TC-Fsm-NNNN",
			n0->hasBeenVisited()
			&& n1->hasBeenVisited()
			&& n2->hasBeenVisited()
			&& n3->hasBeenVisited(),
			"Fsm::accept(FsmVisitor& v) causes a complete traversal of the Fsm. Each FsmNode reachable from the initial node is "
			"marked as visited.");

		fsmlib_assert("TC-Fsm-NNNN",
			not n4->hasBeenVisited(),
			"Unreachable FsmNodes are not visited.");

	}
}

//===================================== Fsm Tests ===================================================

// Checks if tr2 is a deep copy of tr1. True if both have the same label, the same target IDs, but both are different Objects and
// the their FsmLabels are different Objects.
// False otherwise.
bool isFsmTransitionDeepCopy(shared_ptr<FsmTransition> tr1, shared_ptr<FsmTransition> tr2) {
	// same object => copy can't be deep.
	if (tr1 == tr2) return false;
	if (tr1->getLabel() == tr2->getLabel()) return false;
	if (tr1->getSource() == tr2->getSource()) return false;
	if (tr1->getTarget() == tr2->getTarget()) return false;
	// => copy has to by deep. Now check if both have the same values.
	if (tr1->getSource()->getId() != tr2->getSource()->getId()) return false;
	if (tr1->getTarget()->getId() != tr2->getTarget()->getId()) return false;
	if (not(*tr1->getLabel() == *tr2->getLabel())) return false;

	return true;
}

// checks if n2 is a deep copy of n1. True if both have the same id, the same number of outgoing transitions and if
// for all 0 <= i < n1.transitions.size: n2.transitions[i] is a deep copy of n1.transitions[i]. 
// Both FsmNodes have to be different objects/identities.
// False otherwise.
bool isFsmNodeDeepCopy(shared_ptr<FsmNode> n1, shared_ptr<FsmNode> n2) {
	// same object => copy can't be deep
	if (n1 == n2) {
		return false;
	}
	if (n1->getId() != n2->getId()) {
		return false;
	}
	if (n1->getTransitions().size() != n2->getTransitions().size()) {
		return false;
	}
	for (int i = 0; i < n1->getTransitions().size(); ++i) {
		if (not isFsmTransitionDeepCopy(n1->getTransitions().at(i), n2->getTransitions().at(i))) return false;
	}
	return true;
}

// Checks if fsm2Nodes is a copy of fsm1Nodes. 
// True if both lists have the same size and fsm2Nodes[i] is a deep copy of fsm1Nodes[i], for all 0 <= i < fsm2Nodes.size.
// False otherwise.
bool isNodeLstDeepCopy(const vector<shared_ptr<FsmNode>> &fsm1Nodes, const vector<shared_ptr<FsmNode>> &fsm2Nodes) {	
	if (fsm1Nodes.size() != fsm2Nodes.size()) {
		return false;
	}
	for (int i = 0; i < fsm1Nodes.size(); ++i) {
		if (not isFsmNodeDeepCopy(fsm1Nodes.at(i), fsm2Nodes.at(i))) {
			return false;
		}
	}
	return true;
}

// tests Fsm::Fsm(const Fsm& other)
void testFsmDeepCopyConstructor() {
	// Fsm consists of one FsmNode without outgoing transitions.
	{
		const int maxInput = 5;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}

	// n0 --0/0--> n0
	{
		const int maxInput = 3;
		const int maxOutput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}

	// n0 --1/0--> n1
	{
		const int maxInput = 1;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}

	// n0 --0/0--> n0; n0 --1/0--> n1
	{
		const int maxInput = 1;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}

	// n0 --1/0--> n1; n1 --1/1--> n1; n0 --2/2--> n2; n2 --2/1--> n0; n2 --2/1--> n1
	// (non observable fsm)
	{
		const int maxInput = 2;
		const int maxOutput = 2;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --2/2--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(2, 2, pl)));
		// n2 --2/1--> n0
		n2->addTransition(make_shared<FsmTransition>(n2, n0, make_shared<FsmLabel>(2, 1, pl)));
		// n2 --2/1--> n1
		n2->addTransition(make_shared<FsmTransition>(n2, n1, make_shared<FsmLabel>(2, 1, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1, n2}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}

	// n0 --1/1--> n1; n0 --1/0--> n2; n1 --2/2--> n2; n1 --2/3--> n3; n4
	// (non observable fsm)
	{
		const int maxInput = 2;
		const int maxOutput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		shared_ptr<FsmNode> n3 = make_shared<FsmNode>(3, pl);
		shared_ptr<FsmNode> n4 = make_shared<FsmNode>(4, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n0 --1/0--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 0, pl)));
		// n1 --2/2--> n2
		n1->addTransition(make_shared<FsmTransition>(n1, n2, make_shared<FsmLabel>(2, 2, pl)));
		// n1 --2/3--> n3
		n1->addTransition(make_shared<FsmTransition>(n1, n3, make_shared<FsmLabel>(2, 3, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1, n2, n3, n4}, pl);
		Fsm copy(fsm);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == copy.getInitStateIdx()
			&& fsm.getName() == copy.getName()
			&& fsm.getMaxInput() == copy.getMaxInput()
			&& fsm.getMaxOutput() == copy.getMaxOutput()
			&& fsm.getMaxNodes() == copy.getMaxNodes()
			&& fsm.getPresentationLayer() == copy.getPresentationLayer()
			&& fsm.isMinimal() == copy.isMinimal(),
			"Fsm::Fsm(const Fsm& other): Created Fsm has the same values for the expected attributes as other.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(fsm.getNodes(), copy.getNodes()),
			"Fsm::Fsm(const Fsm& other): The node list of the created Fsm is a deep copy of the node list of other.");

		fsmlib_assert("TC-Fsm-NNNN",
			copy.getInitialState()->isInitial(),
			"Fsm::Fsm(const Fsm& other): The expected initial state is marked as 'initial'");
	}
}

// Returns FsmNode from node list of the given fsm with the given id if there is one.
// Returns nullptr otherwise.
shared_ptr<FsmNode> getFsmNodeWithId(Fsm &fsm, int id) {
	for (const auto &n : fsm.getNodes()) {
		if (n == nullptr) continue;
		if (n->getId() == id) {
			return n;
		}
	}
	return nullptr;
}

// Checks if given fsm has a transition from Node with id source to Node with id target, labeled with input, output.
bool hasTransition(int source, int input, int output, int target, Fsm &fsm) {
	shared_ptr<FsmNode> s = getFsmNodeWithId(fsm, source);
	if (s == nullptr) return false;
	for (const auto &tr : s->getTransitions()) {
		if (tr->getTarget()->getId() == target && tr->getLabel()->getInput() == input && tr->getLabel()->getOutput() == output) {
			return true;
		}
	}
	return false;
}

// tests Fsm::Fsm(const string & fname, const string & fsmName, const int maxNodes, const int maxInput, const int maxOutput, const shared_ptr<FsmPresentationLayer> presentationLayer)
void testFsmConstructor1() {
	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor1.fsm
	// Contains only one line: 0 0 0 0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxState = 0;
		const int maxInput = 0;
		const int maxOutput = 0;
		const string name = "M";
		Fsm fsm("../../../resources/TC-Fsm-Constructor1.fsm", name, maxState, maxInput, maxOutput, pl);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == (maxState + 1)
			&& fsm.getMaxInput() == maxInput
			&& fsm.getMaxOutput() == maxOutput
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm contains expected number of Nodes and has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 0, 0, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor2.fsm
	// Contains only one line: 0 1 2 2
	// State with ID 1 has no incoming or outgoing transitions.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxState = 2;
		const int maxInput = 1;
		const int maxOutput = 2;
		const string name = "M";
		Fsm fsm("../../../resources/TC-Fsm-Constructor2.fsm", name, maxState, maxInput, maxOutput, pl);
		fsmlib_assert("TC-Fsm-NNNN",
			/*fsm.getMaxNodes() == (maxState + 1)
			&& */fsm.getMaxInput() == maxInput
			&& fsm.getMaxOutput() == maxOutput
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 1, 2, 2, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 1) == nullptr,
			"States without incoming and outgoing Transitions are nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 0,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor3.fsm
	// Contains three lines / transitions. ID of the initial state is 1.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxState = 1;
		const int maxInput = 2;
		const int maxOutput = 4;
		const string name = "M";
		Fsm fsm("../../../resources/TC-Fsm-Constructor3.fsm", name, maxState, maxInput, maxOutput, pl);
		fsmlib_assert("TC-Fsm-NNNN",
			/*fsm.getMaxNodes() == (maxState + 1)
			&& */fsm.getMaxInput() == maxInput
			&& fsm.getMaxOutput() == maxOutput
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 1
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(1, 1, 1, 0, fsm)
			&& hasTransition(0, 2, 3, 0, fsm)
			&& hasTransition(1, 2, 4, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 2,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor4.fsm
	// Contains six lines / transitions. Specified Fsm is non-observable. 
	// State 1 is not reachable, but has outgoing transitions to reachable states.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxState = 3;
		const int maxInput = 4;
		const int maxOutput = 5;
		const string name = "M";
		Fsm fsm("../../../resources/TC-Fsm-Constructor4.fsm", name, maxState, maxInput, maxOutput, pl);
		fsmlib_assert("TC-Fsm-NNNN",
			/*fsm.getMaxNodes() == (maxState + 1)
			&& */fsm.getMaxInput() == maxInput
			&& fsm.getMaxOutput() == maxOutput
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 0, 0, 2, fsm)
			&& hasTransition(3, 4, 5, 2, fsm)
			&& hasTransition(1, 1, 1, 0, fsm)
			&& hasTransition(1, 1, 1, 1, fsm)
			&& hasTransition(0, 2, 0, 3, fsm)
			&& hasTransition(0, 0, 0, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0) != nullptr
			&& getFsmNodeWithId(fsm, 1) != nullptr
			&& getFsmNodeWithId(fsm, 2) != nullptr
			&& getFsmNodeWithId(fsm, 3) != nullptr,
			"States with incoming or outgoing Transitions aren't nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 3
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 2
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 0
			&& getFsmNodeWithId(fsm, 3)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor5.fsm
	// Contains four lines / transitions. 
	// File contains two identical lines (1 1 1 2). => It is expected that only one transition is created for both identical lines.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const int maxState = 2;
		const int maxInput = 2;
		const int maxOutput = 1;
		const string name = "M";
		Fsm fsm("../../../resources/TC-Fsm-Constructor5.fsm", name, maxState, maxInput, maxOutput, pl);
		fsmlib_assert("TC-Fsm-NNNN",
			/*fsm.getMaxNodes() == (maxState + 1)
			&& */fsm.getMaxInput() == maxInput
			&& fsm.getMaxOutput() == maxOutput
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 1, 0, 1, fsm)
			&& hasTransition(1, 1, 1, 2, fsm)
			&& hasTransition(2, 2, 1, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}
}

// tests Fsm::Fsm(const string& fname, const shared_ptr<FsmPresentationLayer> presentationLayer, const string& fsmName)
void testFsmConstructor2() {
	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor1.fsm
	// Contains only one line: 0 0 0 0
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const string name = "FSM";
		Fsm fsm("../../../resources/TC-Fsm-Constructor1.fsm", pl, name);

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == 1
			&& fsm.getMaxInput() == 0
			&& fsm.getMaxOutput() == 0
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm contains expected number of Nodes and has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 0, 0, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor2.fsm
	// Contains only one line: 0 1 2 2
	// State with ID 1 has no incoming or outgoing transitions.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const string name = "FSM";
		Fsm fsm("../../../resources/TC-Fsm-Constructor2.fsm", pl, name);
		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == 3
			&& fsm.getMaxInput() == 1
			&& fsm.getMaxOutput() == 2
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 1) == nullptr,
			"States without incoming and outgoing Transitions are nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0) != nullptr
			&& getFsmNodeWithId(fsm, 2) != nullptr,
			"States with incoming or outgoing Transitions aren't nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			not getFsmNodeWithId(fsm, 2)->isInitial(),
			"Only the expected FsmNode is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 1, 2, 2, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 0,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor3.fsm
	// Contains three lines / transitions. ID of the initial state is 1.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const string name = "FSM";
		Fsm fsm("../../../resources/TC-Fsm-Constructor3.fsm", pl, name);
		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == 2
			&& fsm.getMaxInput() == 2
			&& fsm.getMaxOutput() == 4
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0) != nullptr
			&& getFsmNodeWithId(fsm, 1) != nullptr,
			"States with incoming or outgoing Transitions aren't nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 1
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			not getFsmNodeWithId(fsm, 0)->isInitial(),
			"Only the expected FsmNode is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(1, 1, 1, 0, fsm)
			&& hasTransition(0, 2, 3, 0, fsm)
			&& hasTransition(1, 2, 4, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 2,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor4.fsm
	// Contains six lines / transitions. Specified Fsm is non-observable. 
	// State 1 is not reachable, but has outgoing transitions to reachable states.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const string name = "FSM";
		Fsm fsm("../../../resources/TC-Fsm-Constructor4.fsm", pl, name);
		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == 4
			&& fsm.getMaxInput() == 4
			&& fsm.getMaxOutput() == 5
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0) != nullptr
			&& getFsmNodeWithId(fsm, 1) != nullptr
			&& getFsmNodeWithId(fsm, 2) != nullptr
			&& getFsmNodeWithId(fsm, 3) != nullptr,
			"States with incoming or outgoing Transitions aren't nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			not getFsmNodeWithId(fsm, 1)->isInitial()
			&& not getFsmNodeWithId(fsm, 2)->isInitial()
			&& not getFsmNodeWithId(fsm, 3)->isInitial(),
			"Only the expected FsmNode is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 0, 0, 2, fsm)
			&& hasTransition(3, 4, 5, 2, fsm)
			&& hasTransition(1, 1, 1, 0, fsm)
			&& hasTransition(1, 1, 1, 1, fsm)
			&& hasTransition(0, 2, 0, 3, fsm)
			&& hasTransition(0, 0, 0, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 3
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 2
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 0
			&& getFsmNodeWithId(fsm, 3)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}

	// uses Fsm specified in file ../../../resources/TC-Fsm-Constructor5.fsm
	// Contains four lines / transitions. 
	// File contains two identical lines (1 1 1 2). => It is expected that only one transition is created for both identical lines.
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		const string name = "FSM";
		Fsm fsm("../../../resources/TC-Fsm-Constructor5.fsm", pl, name);
		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getMaxNodes() == 3
			&& fsm.getMaxInput() == 2
			&& fsm.getMaxOutput() == 1
			&& fsm.getName() == name
			&& fsm.isMinimal() == Minimal::Maybe,
			"Constructed Fsm has expected maxInput, expected maxOutput and expected name.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0) != nullptr
			&& getFsmNodeWithId(fsm, 1) != nullptr
			&& getFsmNodeWithId(fsm, 2) != nullptr,
			"States with incoming or outgoing Transitions aren't nullptr.");

		fsmlib_assert("TC-Fsm-NNNN",
			fsm.getInitStateIdx() == 0
			&& fsm.getInitialState()->isInitial(),
			"initStateIdx is set to the first value of the file and the initial state is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			not getFsmNodeWithId(fsm, 1)->isInitial()
			&& not getFsmNodeWithId(fsm, 2)->isInitial(),
			"Only the expected FsmNode is marked as initial.");

		fsmlib_assert("TC-Fsm-NNNN",
			hasTransition(0, 1, 0, 1, fsm)
			&& hasTransition(1, 1, 1, 2, fsm)
			&& hasTransition(2, 2, 1, 0, fsm),
			"Constructed Fsm contains each expected FsmTransition.");

		fsmlib_assert("TC-Fsm-NNNN",
			getFsmNodeWithId(fsm, 0)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 1)->getTransitions().size() == 1
			&& getFsmNodeWithId(fsm, 2)->getTransitions().size() == 1,
			"Constructed Fsm contains only the expected FsmTransitions.");
	}
}

// Checks if given fsm contains a transition with a label containing an input > maxInput or an output > maxOutput.
bool containsTransitionWithGreaterIO(const shared_ptr<const Fsm> fsm, const int maxInput, const int maxOutput) {
	for (shared_ptr<FsmNode> n : fsm->getNodes()) {
		for (shared_ptr<FsmTransition> tr : n->getTransitions()) {
			if (tr->getLabel()->getInput() > maxInput) return true;
			if (tr->getLabel()->getOutput() > maxOutput) return true;
		}
	}
	return false;
}

void testFsmCreateRandomFsm(const int maxInput, const int maxOutput, const int maxState, const string &name, const unsigned seed) {
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	shared_ptr<Fsm> fsm = Fsm::createRandomFsm(name, maxInput, maxOutput, maxState, pl, seed);

	fsmlib_assert("TC-Fsm-NNNN",
		fsm->getNodes().size() == maxState + 1,
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"Constructed Fsm contains maxState + 1 FsmNodes.");

	fsmlib_assert("TC-Fsm-NNNN",
		fsm->getMaxInput() == maxInput
		&& fsm->getMaxOutput() == maxOutput
		&& fsm->getName() == name
		&& fsm->isMinimal() == Minimal::Maybe,
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"Constructed Fsm has expected maxInput, maxOutput and name values.");

	fsmlib_assert("TC-Fsm-NNNN",
		not containsTransitionWithGreaterIO(fsm, maxInput, maxOutput),
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"Constructed Fsm contains no transition labeled with some input > maxInput or some output > maxOutput.");

	vector<shared_ptr<FsmNode>> unreachableNodes;
	fsmlib_assert("TC-Fsm-NNNN",
		not fsm->removeUnreachableNodes(unreachableNodes),
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"Constructed Fsm contains no unreachable FsmNode.");

	fsmlib_assert("TC-Fsm-NNNN",
		fsm->isCompletelyDefined(),
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"Constructed Fsm is completely specified.");

	int numberOfInitialNodes = 0;
	for (int i = 0; i <= maxState; ++i) {
		if (fsm->getNodes().at(i)->isInitial()) ++numberOfInitialNodes;
	}

	fsmlib_assert("TC-Fsm-NNNN",
		fsm->getInitStateIdx() == 0
		&& fsm->getInitialState()->isInitial()
		&& numberOfInitialNodes == 1,
		"Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed): "
		"initStateIdx is set to 0, the initial state is marked as initial and there is only one initial state.");
}

// tests Fsm::createRandomFsm(const string & fsmName, const int maxInput, const int maxOutput, const int maxState, const shared_ptr<FsmPresentationLayer> pl, const unsigned seed)
void testFsmCreateRandomFsm() {
	{		
		const int maxInput = 0;
		const int maxOutput = 0;
		const int maxState = 0;
		const unsigned seed = 0;
		const string name = "FSM";
		testFsmCreateRandomFsm(maxInput, maxOutput, maxState, name, seed);
	}

	{
		const int maxInput = 1;
		const int maxOutput = 2;
		const int maxState = 1;
		const unsigned seed = 0;
		const string name = "FSM";
		testFsmCreateRandomFsm(maxInput, maxOutput, maxState, name, seed);
	}

	{
		const int maxInput = 3;
		const int maxOutput = 3;
		const int maxState = 5;
		const unsigned seed = 1;
		const string name = "FSM";
		testFsmCreateRandomFsm(maxInput, maxOutput, maxState, name, seed);
	}

	{
		const int maxInput = 5;
		const int maxOutput = 7;
		const int maxState = 10;
		const unsigned seed = 0;
		const string name = "FSM";
		testFsmCreateRandomFsm(maxInput, maxOutput, maxState, name, seed);
	}
}

// Checks if all Nodes with the same index from original and mutant have the same number of outgoing transitions.
// It is assumed that original and mutant have the same number of nodes.
// (for all 0 <= i <= maxState: original.nodes[i].transitions.size() == mutant.nodes[i].transitions.size())
bool checkIfEachNodeHasSameNumberOfTransitions(const Fsm &original, const Fsm &mutant) {
	for (int i = 0; i < original.getNodes().size(); ++i) {
		if (original.getNodes().at(i)->getTransitions().size() != mutant.getNodes().at(i)->getTransitions().size())
		{
			return false;
		}
	}
	return true;
}

// Checks if a transition from mutant differs from the corresponding transition in original wrt. the input.
// It is assumed that original and mutant have the same number of nodes and that all corresponding nodes have the same number of transitions.
bool checkIfInputsOfTransitionsChanged(const Fsm &original, const Fsm &mutant) {
	for (int i = 0; i < original.getNodes().size(); ++i) {
		shared_ptr<FsmNode> originalNode = original.getNodes().at(i);
		shared_ptr<FsmNode> mutantNode = mutant.getNodes().at(i);
		for (int j = 0; j < originalNode->getTransitions().size(); ++j) {
			if (originalNode->getTransitions().at(j)->getLabel()->getInput() != mutantNode->getTransitions().at(j)->getLabel()->getInput()) {
				return true;
			}
		}
	}
	return false;
}

// Calculates the number of output faults in mutant.
// It is assumed that original and mutant have the same number of nodes and that all corresponding nodes have the same number of transitions.
size_t calcNumberOfOutputFaults(Fsm &original, Fsm &mutant) {
	int numOutputFaults = 0;
	for (int i = 0; i < original.getNodes().size(); ++i) {
		shared_ptr<FsmNode> originalNode = original.getNodes().at(i);
		shared_ptr<FsmNode> mutantNode = mutant.getNodes().at(i);
		for (int j = 0; j < originalNode->getTransitions().size(); ++j) {
			if (originalNode->getTransitions().at(j)->getLabel()->getOutput() != mutantNode->getTransitions().at(j)->getLabel()->getOutput()) {
				++numOutputFaults;
			}
		}
	}
	return numOutputFaults;
}

// Calculates the number of transition faults in mutant.
// It is assumed that original and mutant have the same number of nodes and that all corresponding nodes have the same number of transitions
// and the same id.
size_t calcNumberOfTransitionFaults(Fsm &original, Fsm &mutant) {
	int numTransitionFaults = 0;
	for (int i = 0; i < original.getNodes().size(); ++i) {
		shared_ptr<FsmNode> originalNode = original.getNodes().at(i);
		shared_ptr<FsmNode> mutantNode = mutant.getNodes().at(i);
		for (int j = 0; j < originalNode->getTransitions().size(); ++j) {
			if (originalNode->getTransitions().at(j)->getTarget()->getId() != mutantNode->getTransitions().at(j)->getTarget()->getId()) {
				++numTransitionFaults;
			}
		}
	}
	return numTransitionFaults;
}

// Checks if some Node in mutant has several identical (same target id, same label) transitions.
bool hasIdenticalTransitions(shared_ptr<Fsm> mutant) {
	for (shared_ptr<FsmNode> n : mutant->getNodes()) {
		for (int i = 0; i < n->getTransitions().size(); ++i) {
			for (int j = i + 1; j < n->getTransitions().size(); ++j) {
				shared_ptr<FsmTransition> tr1 = n->getTransitions().at(i);
				shared_ptr<FsmTransition> tr2 = n->getTransitions().at(j);
				if (*tr1->getLabel() == *tr2->getLabel() && tr1->getTarget()->getId() == tr2->getTarget()->getId()) {
					return true;
				}
			}
		}
	}
	return false;
}

// tests Fsm::createMutant(const std::string & fsmName, const size_t numOutputFaults, const size_t numTransitionFaults)
void testFsmCreateMutant(Fsm &original, const int numOutputFaults, const int numTransitionFaults, const int maxInput, const int maxOutput) {
	Fsm copyOfOriginal(original);
	shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);

	fsmlib_assert("TC-Fsm-NNNN",
		isNodeLstDeepCopy(original.getNodes(), copyOfOriginal.getNodes()),
		"Fsm::createMutant(...) doesn't change the nodes and transitions of original.");

	fsmlib_assert("TC-Fsm-NNNN",
		original.getNodes().size() == mutant->getNodes().size(),
		"Fsm::createMutant(...) creates mutant with the same number of nodes.");

	fsmlib_assert("TC-Fsm-NNNN",
		original.getMaxInput() == mutant->getMaxInput()
		&& original.getMaxOutput() == mutant->getMaxOutput()
		&& original.getInitStateIdx() == mutant->getInitStateIdx(),
		"Fsm::createMutant(...) generates mutant which has the same maxInput, maxOutput and the same initStateIdx.");

	fsmlib_assert("TC-Fsm-NNNN",
		checkIfEachNodeHasSameNumberOfTransitions(original, *mutant),
		"Each Node in the mutant has the same number of transitions as it had in the original fsm.");

	fsmlib_assert("TC-Fsm-NNNN",
		not checkIfInputsOfTransitionsChanged(original, *mutant),
		"Each transition in the mutant has the same input as it has in the original fsm.");

	fsmlib_assert("TC-Fsm-NNNN",
		not containsTransitionWithGreaterIO(mutant, maxInput, maxOutput),
		"Each input/output in the mutant is less than or equal to maxInput/maxOutput.");

	fsmlib_assert("TC-Fsm-NNNN",
		original.isCompletelyDefined() == mutant->isCompletelyDefined(),
		"The mutant is completely defined iff original is completely defined.");

	fsmlib_assert("TC-Fsm-NNNN",
		mutant->getInitialState()->isInitial(),
		"Initial Node of the mutant is marked as initial.");

	fsmlib_assert("TC-Fsm-NNNN",
		calcNumberOfOutputFaults(original, *mutant) <= numOutputFaults,
		"The mutant contains at most numOutputFaults many output faults.");

	fsmlib_assert("TC-Fsm-NNNN",
		calcNumberOfTransitionFaults(original, *mutant) <= numTransitionFaults,
		"The mutant contains at most numTransitionFaults many transition faults.");

	//fsmlib_assert("TC-Fsm-NNNN",
	//	not hasIdenticalTransitions(mutant),
	//	"Each transition in mutant is contained only once.");
}

// tests Fsm::createMutant(const std::string & fsmName, const size_t numOutputFaults, const size_t numTransitionFaults)
void testFsmCreateMutant() {
	// numOutputFault = 0 and numTransitionFaults = 0 => mutant should be identical to original
	// Fsm: TC-Fsm-CreateMutant1.fsm (non-observable, #states: 4, #transitions: 7, maxInput: 3, maxOutput: 3)
	{
		const int numOutputFaults = 0;
		const int numTransitionFaults = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Fsm original("../../../resources/TC-Fsm-CreateMutant1.fsm", pl, "Original");
		Fsm copyOfOriginal(original);
		shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(), copyOfOriginal.getNodes()),
			"Fsm::createMutant(...) doesn't change the nodes and transitions of original.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(), mutant->getNodes())
			&& original.getMaxInput() == mutant->getMaxInput()
			&& original.getMaxOutput() == mutant->getMaxOutput()
			&& original.getInitStateIdx() == mutant->getInitStateIdx(),
			"Fsm::createMutant(...) generates mutant which is identical to original if numOutputFaults = 0 and numTransitionFaults = 0.");
	}

	// numOutputFault = 1 and numTransitionFaults = 2
	// Fsm consists of only one Node with 0 outgoing transitions.
	//{
	//	const int numOutputFaults = 1;
	//	const int numTransitionFaults = 2;
	//	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	//	shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
	//	Fsm original("Original", 0, 0, vector<shared_ptr<FsmNode>>{n0}, pl);
	//	Fsm copyOfOriginal(original);
	//	shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);
	//}

	// numOutputFault = 1 and numTransitionFaults = 0 => Mutant differs from original in one output only.
	// Fsm: n0 --1/1--> n1; n1 --1/1--> n1
	{
		const int numOutputFaults = 1;
		const int numTransitionFaults = 0;
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		Fsm original("Original", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1}, pl);
		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}

	// numOutputFault = 0 and numTransitionFaults = 1 => Mutant differs from original in one transition-target only.
	// Fsm: n0 --1/1--> n1; n1 --1/1--> n1
	{
		const int numOutputFaults = 0;
		const int numTransitionFaults = 1;
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --1/1--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(1, 1, pl)));
		Fsm original("Original", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1}, pl);

		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}
	
	// Fsm: n0 --0/0--> n0;
	// numOutputFault = 0 and numTransitionFaults = 1 => The target state of the only transition can't change because there is only one state.
	{
		const int numOutputFaults = 0;
		const int numTransitionFaults = 1;
		const int maxInput = 0;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		Fsm original("Original", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		Fsm copyOfOriginal(original);
		shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(), copyOfOriginal.getNodes()),
			"Fsm::createMutant(...) doesn't change the nodes and transitions of original.");

		fsmlib_assert("TC-Fsm-NNNN",
			mutant->getNodes().size() == 1
			&& hasTransition(0, 0, 0, 0, *mutant),
			"Fsm::createMutant(...) doesn't insert transition faults if original has only one node.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.getNodes().size() == mutant->getNodes().size(),
			"Fsm::createMutant(...) creates mutant with the same number of nodes.");
		
		fsmlib_assert("TC-Fsm-NNNN",
			original.getMaxInput() == mutant->getMaxInput()
			&& original.getMaxOutput() == mutant->getMaxOutput()
			&& original.getInitStateIdx() == mutant->getInitStateIdx(),
			"Fsm::createMutant(...) generates mutant which has the same maxInput, maxOutput and the same initStateIdx.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachNodeHasSameNumberOfTransitions(original, *mutant),
			"Each Node in the mutant has the same number of transitions as it had in the original fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.isCompletelyDefined() == mutant->isCompletelyDefined(),
			"The mutant is completely defined iff original is completely defined.");

		fsmlib_assert("TC-Fsm-NNNN",
			mutant->getInitialState()->isInitial(),
			"Initial Node of the mutant is marked as initial.");
	}

	// Fsm: n0 --0/0--> n0;
	// numOutputFault = 1 and numTransitionFaults = 0 => The output of the only transition can't change because 0 is the maxOutput.
	{
		const int numOutputFaults = 1;
		const int numTransitionFaults = 0;
		const int maxInput = 0;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		Fsm original("Original", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		Fsm copyOfOriginal(original);
		shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(), copyOfOriginal.getNodes()),
			"Fsm::createMutant(...) doesn't change the nodes and transitions of original.");

		fsmlib_assert("TC-Fsm-NNNN",
			mutant->getNodes().size() == 1
			&& hasTransition(0, 0, 0, 0, *mutant),
			"Fsm::createMutant(...) doesn't insert output faults if maxOutput is 0.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.getNodes().size() == mutant->getNodes().size(),
			"Fsm::createMutant(...) creates mutant with the same number of nodes.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.getMaxInput() == mutant->getMaxInput()
			&& original.getMaxOutput() == mutant->getMaxOutput()
			&& original.getInitStateIdx() == mutant->getInitStateIdx(),
			"Fsm::createMutant(...) generates mutant which has the same maxInput, maxOutput and the same initStateIdx.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachNodeHasSameNumberOfTransitions(original, *mutant),
			"Each Node in the mutant has the same number of transitions as it had in the original fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.isCompletelyDefined() == mutant->isCompletelyDefined(),
			"The mutant is completely defined iff original is completely defined.");

		fsmlib_assert("TC-Fsm-NNNN",
			mutant->getInitialState()->isInitial(),
			"Initial Node of the mutant is marked as initial.");
	}

	// numOutputFault = 2 and numTransitionFaults = 1 => Both fault types != 0.
	// Fsm: n0 --1/1--> n1; n1 --0/2--> n1; n0 --2/2--> n2; n2 --1/0--> n1 (deterministic)
	{
		const int numOutputFaults = 2;
		const int numTransitionFaults = 1;
		const int maxInput = 2;
		const int maxOutput = 2;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/1--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 1, pl)));
		// n1 --0/2--> n1
		n1->addTransition(make_shared<FsmTransition>(n1, n1, make_shared<FsmLabel>(0, 2, pl)));
		// n0 --2/2--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(2, 2, pl)));
		// n2 --1/0--> n1
		n2->addTransition(make_shared<FsmTransition>(n2, n1, make_shared<FsmLabel>(1, 0, pl)));
		Fsm original("Original", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1, n2}, pl);

		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}

	// numOutputFault = 2 and numTransitionFaults = 2 => Both fault types != 0.
	// Fsm: TC-Fsm-CreateMutant1.fsm (non-observable, #states: 4, #transitions: 7, maxInput: 3, maxOutput: 3)
	{
		const int numOutputFaults = 2;
		const int numTransitionFaults = 2;
		const int maxInput = 3;
		const int maxOutput = 3;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Fsm original("../../../resources/TC-Fsm-CreateMutant1.fsm", pl, "Original");

		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}

	// numOutputFault = 10 and numTransitionFaults = 10 => Both fault types != 0.
	// Fsm: TC-Fsm-CreateMutant2.fsm (non-observable, #states: 3, #transitions: 6, maxInput: 1, maxOutput: 1)
	{
		const int numOutputFaults = 10;
		const int numTransitionFaults = 10;
		const int maxInput = 1;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Fsm original("../../../resources/TC-Fsm-CreateMutant2.fsm", pl, "Original");

		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}

	// numOutputFault = 1 and numTransitionFaults = 0
	// Fsm: TC-Fsm-CreateMutant3.fsm (#states: 2, #transitions: 4, maxInput: 0, maxOutput: 1)
	// => each output fault creates a transition, that is already contained.
	{
		const int numOutputFaults = 1;
		const int numTransitionFaults = 0;
		const int maxInput = 0;
		const int maxOutput = 1;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Fsm original("../../../resources/TC-Fsm-CreateMutant3.fsm", pl, "Original");
		Fsm copyOfOriginal(original);
		shared_ptr<Fsm> mutant = original.createMutant("Mutant", numOutputFaults, numTransitionFaults);

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(), copyOfOriginal.getNodes()),
			"Fsm::createMutant(...) doesn't change the nodes and transitions of original.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.getNodes().size() == mutant->getNodes().size(),
			"Fsm::createMutant(...) creates mutant with the same number of nodes.");

		fsmlib_assert("TC-Fsm-NNNN",
			original.getMaxInput() == mutant->getMaxInput()
			&& original.getMaxOutput() == mutant->getMaxOutput()
			&& original.getInitStateIdx() == mutant->getInitStateIdx(),
			"Fsm::createMutant(...) generates mutant which has the same maxInput, maxOutput and the same initStateIdx.");

		fsmlib_assert("TC-Fsm-NNNN",
			isNodeLstDeepCopy(original.getNodes(),mutant->getNodes()),
			"Mutant is identical to original because it's impossible to add output faults.");
	}

	// numOutputFault = 1 and numTransitionFaults = 1
	// Fsm: TC-Fsm-CreateMutant4.fsm (deterministic, #states: 5, #transitions: 7, maxInput: 2, maxOutput: 2)
	// not all nodes are reachable
	{
		const int numOutputFaults = 1;
		const int numTransitionFaults = 1;
		const int maxInput = 2;
		const int maxOutput = 2;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		Fsm original("../../../resources/TC-Fsm-CreateMutant4.fsm", pl, "Original");

		testFsmCreateMutant(original, numOutputFaults, numTransitionFaults, maxInput, maxOutput);
	}
}

struct FileFormatTransition
{
	size_t source;
	size_t input;
	size_t output;
	size_t target;
};

// Extracts all Transitions contained in file under fileName. ("pre-state input output post-state"-format)
// Returns vector containing those transitions.
// If file can't be opened a empty vector is returned.
shared_ptr<vector<FileFormatTransition>> extractTransitionsFromFile(const string fileName) {
	vector<FileFormatTransition> transitions;
	ifstream inputFile(fileName);
	if (inputFile.is_open())
	{
		string line;
		size_t source;
		size_t input;
		size_t output;
		size_t target;
		while (getline(inputFile, line)) {
			stringstream ss(line);
			ss >> source;
			ss >> input;
			ss >> output;
			ss >> target;
			transitions.push_back(FileFormatTransition{ source, input, output, target });
		}
		inputFile.close();
	}

	else
	{
		cout << "Unable to open input file" << endl;
		//exit(EXIT_FAILURE);
	}
	return make_shared<vector<FileFormatTransition>>(transitions);
}

// Checks if the given FsmTransition is contained in the list fileTransitions extracted from a .fsm file
bool checkIfTransitionContainedInFile(shared_ptr<vector<FileFormatTransition>> fileTransitions, shared_ptr<FsmTransition> tr) {
	for (const auto &fileTr : *fileTransitions) {
		if (fileTr.source == tr->getSource()->getId()
			&& fileTr.input == tr->getLabel()->getInput()
			&& fileTr.output == tr->getLabel()->getOutput()
			&& fileTr.target == tr->getTarget()->getId()) {
			return true;
		}
	}
	return false;
}

// Checks if each transition from the given fsm is contained in the list fileTransitions extracted from a .fsm file.
bool checkIfEachTransitionContainedInFile(shared_ptr<vector<FileFormatTransition>> fileTransitions, Fsm &fsm) {
	for (const auto n : fsm.getNodes()) {
		for (const auto tr : n->getTransitions()) {
			if (not checkIfTransitionContainedInFile(fileTransitions, tr)) return false;
		}
	}
	return true;
}

// Checks if each transition in the given list fileTransitions extracted from a .fsm file is contained in the given fsm.
bool checkIfEachTransitionContainedInFsm(Fsm &fsm, shared_ptr<vector<FileFormatTransition>> fileTransitions) {
	for (const auto &tr : *fileTransitions) {
		if (not hasTransition(tr.source, tr.input, tr.output, tr.target, fsm)) {
			return false;
		}
	}
	return true;
}

// Calculates the total number of FsmTransitions of the given fsm.
size_t getNumberOfTransitions(Fsm &fsm) {
	size_t numTransitions = 0;
	for (const auto n : fsm.getNodes()) {
		numTransitions += n->getTransitions().size();
	}
	return numTransitions;
}

// tests Fsm::dumpFsm(std::ofstream & outputFile)
void testFsmDumpFsm() {
	// Fsm: n0 --0/0--> n0 (contains only one state and no transition)
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		const int maxInput = 0;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->empty(),
			"The file contains no transition if the fsm contains no transition.");
	}

	// Fsm: n0 --0/0--> n0 (contains only one state and one transition)
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		const int maxInput = 0;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);		
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0}, pl);
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFile(transitions, fsm),
			"Each transition of the Fsm is written to the file.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFsm(fsm, transitions),
			"Each transition that was written to the file is a FsmTransition of the Fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			getNumberOfTransitions(fsm) == transitions->size(),
			"The Fsm and the file contain the same number of transitions.");

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->at(0).source == fsm.getInitStateIdx(),
			"The first value that was written to the file is the id of the initial state.");		
	}

	// Fsm: n0 --1/0--> n1 (contains two states and one transition)
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		const int maxInput = 1;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		// n0 --1/0--> n1
		n0->addTransition(make_shared<FsmTransition>(n0, n1, make_shared<FsmLabel>(1, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1}, pl);
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFile(transitions, fsm),
			"Each transition of the Fsm is written to the file.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFsm(fsm, transitions),
			"Each transition that was written to the file is a FsmTransition of the Fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			getNumberOfTransitions(fsm) == transitions->size(),
			"The Fsm and the file contain the same number of transitions.");

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->at(0).source == fsm.getInitStateIdx(),
			"The first value that was written to the file is the id of the initial state.");
	}

	// Fsm: n0 --1/0--> n2; n0 --0/0--> n0;
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		const int maxInput = 1;
		const int maxOutput = 0;
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		shared_ptr<FsmNode> n0 = make_shared<FsmNode>(0, pl);
		shared_ptr<FsmNode> n1 = make_shared<FsmNode>(1, pl);
		shared_ptr<FsmNode> n2 = make_shared<FsmNode>(2, pl);
		// n0 --1/0--> n2
		n0->addTransition(make_shared<FsmTransition>(n0, n2, make_shared<FsmLabel>(1, 0, pl)));
		// n0 --0/0--> n0
		n0->addTransition(make_shared<FsmTransition>(n0, n0, make_shared<FsmLabel>(0, 0, pl)));
		Fsm fsm("M", maxInput, maxOutput, vector<shared_ptr<FsmNode>>{n0, n1, n2}, pl);
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFile(transitions, fsm),
			"Each transition of the Fsm is written to the file.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFsm(fsm, transitions),
			"Each transition that was written to the file is a FsmTransition of the Fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			getNumberOfTransitions(fsm) == transitions->size(),
			"The Fsm and the file contain the same number of transitions.");

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->at(0).source == fsm.getInitStateIdx(),
			"The first value that was written to the file is the id of the initial state.");
	}

	// Fsm: TC-Fsm-Constructor3.fsm (First value in the file is != 0  =>  id of initial state is 1)
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		Fsm fsm("../../../resources/TC-Fsm-Constructor3.fsm", make_shared<FsmPresentationLayer>(), "M");
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFile(transitions, fsm),
			"Each transition of the Fsm is written to the file.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFsm(fsm, transitions),
			"Each transition that was written to the file is a FsmTransition of the Fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			getNumberOfTransitions(fsm) == transitions->size(),
			"The Fsm and the file contain the same number of transitions.");

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->at(0).source == fsm.getInitStateIdx(),
			"The first value that was written to the file is the id of the initial state.");
	}

	// Fsm: TC-Fsm-Constructor4.fsm (non-observable, #transitions = 6)
	{
		const string fileName = "../../../resources/TC-Fsm-DumpFsm.fsm";
		std::ofstream ofs(fileName);
		Fsm fsm("../../../resources/TC-Fsm-Constructor4.fsm", make_shared<FsmPresentationLayer>(), "M");
		fsm.dumpFsm(ofs);
		ofs.close();

		shared_ptr<vector<FileFormatTransition>>  transitions = extractTransitionsFromFile(fileName);
	
		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFile(transitions, fsm),
			"Each transition of the Fsm is written to the file.");

		fsmlib_assert("TC-Fsm-NNNN",
			checkIfEachTransitionContainedInFsm(fsm, transitions),
			"Each transition that was written to the file is a FsmTransition of the Fsm.");

		fsmlib_assert("TC-Fsm-NNNN",
			getNumberOfTransitions(fsm) == transitions->size(),
			"The Fsm and the file contain the same number of transitions.");

		fsmlib_assert("TC-Fsm-NNNN",
			transitions->at(0).source == fsm.getInitStateIdx(),
			"The first value that was written to the file is the id of the initial state.");
	}
}

// Checks if there is a transition in intersection from source to target labeled with input/output and
// if fsm1 and fsm2 both have a corresponding transition with the same label, between the corresponding source and target nodes.
// (specified in the derivedFromPair attribute)
//bool checkTransitionOfIntersection(shared_ptr<FsmNode> source, shared_ptr<FsmNode> target, int input, int output, Fsm & intersection, Fsm & fsm1, Fsm & fsm2) {
//	if (hasTransition(source->getId(), input, output, target->getId(), intersection)) {
//		if()
//	}
//}

// Returns list containing all FsmTransitions of the given fsm.
vector<shared_ptr<FsmTransition>> getAllTransitionsFromFsm(Fsm& fsm) {
	vector<shared_ptr<FsmTransition>> transitions;
	for (const auto n : fsm.getNodes()) {
		transitions.insert(transitions.end(), n->getTransitions().begin(), n->getTransitions().end());
	}
	return transitions;
}

// Checks for all transitions in intersection, if (q_i,s_j)--x/y-->(q_i',s_j') => q_i --x/y--> q_i' and s_j --x/y--> s_j'  
// (Checks if each transition in intersection has corresponding transitions in fsm1 and fsm2 )
bool checkEachTransitionOfIntersection(Fsm &intersection, Fsm &fsm1, Fsm &fsm2) {
	for (const auto tr : getAllTransitionsFromFsm(intersection)) {
		const auto sourcePair = tr->getSource()->getPair();
		const auto targetPair = tr->getTarget()->getPair();
		shared_ptr<FsmLabel> lbl = tr->getLabel();
		cout << "Fsm1: " << sourcePair->first->getId() << "," << lbl->getInput() << "," << lbl->getOutput() << "," << targetPair->first->getId() << endl;
		cout << "Fsm2: " << sourcePair->second->getId() << "," << lbl->getInput() << "," << lbl->getOutput() << "," << targetPair->second->getId() << endl;
		if (not hasTransition(sourcePair->first->getId(), lbl->getInput(), lbl->getOutput(), targetPair->first->getId(), fsm1)) {
			return false;
		}
		if (not hasTransition(sourcePair->second->getId(), lbl->getInput(), lbl->getOutput(), targetPair->second->getId(), fsm2)) {
			return false;
		}
	}
	return true;
}

// Gets FsmNode from intersection that is derived from pair n1,n2.
// Returns nullptr if no such FsmNode exists.
shared_ptr<FsmNode> getNodeWithPair(Fsm &intersection, shared_ptr<FsmNode> n1, shared_ptr<FsmNode> n2) {
	for (shared_ptr<FsmNode> n : intersection.getNodes()) {
		if (n->getPair()->first == n1 && n->getPair()->second == n2) {
			return n;
		}
	}
	return nullptr;
}

// Checks if nodeDerivedFromPair has a transition with the same label as tr1 and tr2 (it is assumed that tr1.label == tr2.label),
// and if the target of this transition is derived from tr1.target and tr2.target.
bool hasMatchingTransition(shared_ptr<FsmNode> nodeDerivedFromPair, shared_ptr<FsmTransition> tr1, shared_ptr<FsmTransition> tr2) {
	for (const auto tr : nodeDerivedFromPair->getTransitions()) {
		//cout << "tr from derivedNode: " << *tr << endl;
		if (not (*tr->getLabel() == *tr1->getLabel())) continue;
		if (tr->getTarget()->getPair()->first == tr1->getTarget() && tr->getTarget()->getPair()->second == tr2->getTarget()) {
			//cout << "matching targets: " << tr1->getTarget()->getId() << ", " << tr2->getTarget()->getId() << endl;
			return true;
		}
	}
	return false;
}

// Checks if each pair of transitions of n1 and n2 which have the same labels has a corresponding transition in nodeDerivedFromPair.
bool checkEachTransitionOfPair(shared_ptr<FsmNode> nodeDerivedFromPair, shared_ptr<FsmNode> n1, shared_ptr<FsmNode> n2) {
	for (const auto tr1 : n1->getTransitions()) {
		for (const auto tr2 : n2->getTransitions()) {
			if (not (*tr1->getLabel() == *tr2->getLabel())) continue;
			// same label => nodeDerivedFromPair must have a transition with same label
			//cout << *tr1 << " and " << *tr2 << endl;
			if (not hasMatchingTransition(nodeDerivedFromPair, tr1, tr2)) {
				return false;
			}
		}
	}
	return true;
}

// Checks if (q,x,y,q') in fsm1 and (s,x,y,s') in fsm2 implies ((q,s),x,y,(q',s')) in intersection,
// but only for (q,s), that are contained in intersection 
// (intersection contains only (q,s) that are reachable from the initial state (q0,s0))
// Note that we separatly check that (q0,s0) is contained as the initial node in the intersection, so that all 
// Nodes that can be independently reached in fsm1 and fsm2 by the same InputTrace form a tuple/state in the intersection.
bool checkEachTransitionPair(Fsm &fsm1, Fsm &fsm2, Fsm &intersection) {
	for (const auto n1 : fsm1.getNodes()) {
		for (const auto n2 : fsm2.getNodes()) {
			shared_ptr<FsmNode> nodeDerivedFromPair = getNodeWithPair(intersection, n1, n2);
			if (nodeDerivedFromPair == nullptr) continue;
			if (not checkEachTransitionOfPair(nodeDerivedFromPair, n1, n2)) return false;			
		}
	}
	return true;
}


// tests Fsm::intersect(const Fsm & f)
//void testFsmIntersect() {
//	// fsm1 = fsm2 = TC-Fsm-Intersect1.fsm 
//	{
//		Fsm fsm1("../../../resources/TC-Fsm-Intersect1.fsm", make_shared<FsmPresentationLayer>(), "M");
//		Fsm fsm2("../../../resources/TC-Fsm-Intersect1.fsm", make_shared<FsmPresentationLayer>(), "M");
//
//		cout << fsm1.getNodes().size() << endl;
//
//		Fsm intersection = fsm1.intersect(fsm2);
//		cout << intersection.getNodes().size() << endl;
//		IOListContainer iolc(4, intersection.getNodes().size(), intersection.getNodes().size(), make_shared<FsmPresentationLayer>());
//		//for (const auto & v : *iolc.getIOLists()) {
//			//cout << "------------------------------" << endl;
//			//InputTrace itrc{ v, make_shared<FsmPresentationLayer>() };
//			//OutputTree ot1 = intersection.apply(itrc);
//			//cout << ot1 << endl;
//			//OutputTree ot2 = fsm1.apply(itrc);
//			//cout << ot2 << endl;
//			//cout << "------------------------------" << endl;
//		//}
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			intersection.getInitialState()->getPair()->first == fsm1.getInitialState()
//			&& intersection.getInitialState()->getPair()->second == fsm2.getInitialState(),
//			"The initial state (q,s) of the intersection is derived from the initial states q and s.");
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			checkEachTransitionOfIntersection(intersection,fsm1,fsm2),
//			"((q,s),x,y,(q',s')) in intersection => (q,x,y,q') in fsm1 and (s,x,y,s') in fsm2.");
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			checkEachTransitionPair(fsm1, fsm2, intersection),
//			"(q,x,y,q') in fsm1 and (s,x,y,s') in fsm2 => ((q,s),x,y,(q',s')) in intersection.");
//
//		//cout << fsm1 << endl;
//
//		//for (auto tr : getAllTransitionsFromFsm(intersection)) {
//		//	cout << *tr << endl;
//		//	cout << tr->getSource()->getPair()->first->getId() << " ";
//		//	cout << tr->getTarget()->getPair()->first->getId() << endl;
//		//}
//	}
//
//	// fsm1 = TC-Fsm-Intersect2.fsm , TC-Fsm-Intersect3.fsm (example from the lecture)
//	{
//		Fsm fsm1("../../../resources/TC-Fsm-Intersect2.fsm", make_shared<FsmPresentationLayer>(), "M1");
//		Fsm fsm2("../../../resources/TC-Fsm-Intersect3.fsm", make_shared<FsmPresentationLayer>(), "M2");
//
//		cout << fsm1.getNodes().size() << endl;
//
//		Fsm intersection = fsm1.intersect(fsm2);
//		cout << "intersection size: " << intersection.getNodes().size() << endl;
//
//		for (const auto n : intersection.getNodes()) {
//			cout << n->getPair()->first->getId() << "," << n->getPair()->second->getId() << endl;
//		}
//		IOListContainer iolc(1, intersection.getNodes().size(), intersection.getNodes().size(), make_shared<FsmPresentationLayer>());
//		for (const auto & v : *iolc.getIOLists()) {
//			cout << "------------------------------" << endl;
//			InputTrace itrc{ v, make_shared<FsmPresentationLayer>() };
//			OutputTree ot1 = intersection.apply(itrc); 
//			for (auto trc : ot1.getOutputTraces()) cout << trc << endl;
//			//cout << ot1 << endl;
//			cout << "===========" << endl;
//			OutputTree ot2 = fsm1.apply(itrc);
//			for (auto trc : ot2.getOutputTraces()) cout << trc << endl;
//			//cout << ot2 << endl;
//			cout << "------------------------------" << endl;
//		}
//
//		cout << "deterministic: " << intersection.isDeterministic() << endl;
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			intersection.getInitialState()->getPair()->first == fsm1.getInitialState()
//			&& intersection.getInitialState()->getPair()->second == fsm2.getInitialState(),
//			"The initial state (q,s) of the intersection is derived from the initial states q and s.");
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			checkEachTransitionOfIntersection(intersection, fsm1, fsm2),
//			"((q,s),x,y,(q',s')) in intersection => (q,x,y,q') in fsm1 and (s,x,y,s') in fsm2.");
//
//		fsmlib_assert("TC-Fsm-NNNN",
//			checkEachTransitionPair(fsm1, fsm2, intersection),
//			"(q,x,y,q') in fsm1 and (s,x,y,s') in fsm2 => ((q,s),x,y,(q',s')) in intersection.");
//
//		//cout << fsm1 << endl;
//
//		//for (auto tr : getAllTransitionsFromFsm(intersection)) {
//		//	cout << *tr << endl;
//		//	cout << tr->getSource()->getPair()->first->getId() << " ";
//		//	cout << tr->getTarget()->getPair()->first->getId() << endl;
//		//}
//	}
//}

// Nearly an exact copy of Fsm::apply() but here no element of the inputtrace can be ignored while constructing the outputtraces
// => only outputtraces generated for the whole inputtrace or a prefix of it is added to the returned outputtree.
OutputTree applyCompleteTraceWithPrefixes(shared_ptr<FsmNode> fsm_root, shared_ptr<FsmPresentationLayer> presentationLayer, const InputTrace& itrc, bool markAsVisited)
{
	deque<shared_ptr<TreeNode>> tnl;
	unordered_map<shared_ptr<TreeNode>, shared_ptr<FsmNode>> t2f;

	shared_ptr<TreeNode> root = make_shared<TreeNode>();
	OutputTree ot = OutputTree(root, itrc, presentationLayer);

	if (itrc.get().size() == 0)
	{
		return ot;
	}

	t2f[root] = fsm_root;

	for (auto it = itrc.cbegin(); it != itrc.cend(); ++it)
	{
		int x = *it;

		vector< shared_ptr<TreeNode> > vaux = ot.getLeaves();
		//cout << "it - itrc.cbegin()" << it - itrc.cbegin() << endl;
		for (auto n : vaux) {
			//cout << "n->getPath().size()" << n->getPath().size() << endl;
			if (n->getPath().size() < it - itrc.cbegin()) continue;
			tnl.push_back(n);
		}

		while (!tnl.empty())
		{
			shared_ptr<TreeNode> thisTreeNode = tnl.front();
			tnl.pop_front();

			shared_ptr<FsmNode> thisState = t2f.at(thisTreeNode);
			if (markAsVisited) thisState->setVisited();

			for (shared_ptr<FsmTransition> tr : thisState->getTransitions())
			{
				if (tr->getLabel()->getInput() == x)
				{
					int y = tr->getLabel()->getOutput();
					shared_ptr<FsmNode> tgtState = tr->getTarget();
					shared_ptr<TreeNode> tgtNode = make_shared<TreeNode>();
					shared_ptr<TreeEdge> te = make_shared<TreeEdge>(y, tgtNode);
					thisTreeNode->add(te);
					t2f[tgtNode] = tgtState;
					if (markAsVisited) tgtState->setVisited();
				}
			}
		}
	}
	return ot;
}

int main(int argc, char** argv)
{
    
    
    
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
	//testFsmIntersect();


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



