/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <string>
#include <chrono>
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
#include <fsm/FsmTransition.h>
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
#include <deque>

using namespace std;
using namespace Json;

using std::chrono::system_clock;

static const string csvSep = ";";

static string nowText;

static const string ascTestDirectory = "../../../resources/asc-tests/";
static const string ascTestResultDirectory = "../../../resources/asc-test-results/";
static const string ascCsvDirectory = "../../../csv/";

static shared_ptr<FsmPresentationLayer> plTestSpec = make_shared<FsmPresentationLayer>(
            ascTestDirectory + "adaptive-test-in.txt",
            ascTestDirectory + "adaptive-test-out.txt",
            ascTestDirectory + "adaptive-test-state-spec.txt");

static shared_ptr<FsmPresentationLayer> plTestIut = make_shared<FsmPresentationLayer>(
            ascTestDirectory + "adaptive-test-in.txt",
            ascTestDirectory + "adaptive-test-out.txt",
            ascTestDirectory + "adaptive-test-state-iut.txt");


static const string testSepLine = "###################################################################";

enum class TestIteration
{
    INPUT,
    OUTPUT,
    STATE,
    OUTPUT_FAULT,
    TRANSITION_FAULT,
    DEGREE_COMPLETENESS,
    END
};

struct TestIterationHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

enum class CsvField
{
    TEST_NAME,
    NUM_STATES,
    NUM_INPUTS,
    NUM_OUTPUTS,
    NUM_D_REACHABLE_STATES,
    NUM_SETS_OF_MAXIMAL_R_DIST_STATES,
    NUM_OUT_FAULTS,
    NUM_TRANS_FAULTS,
    DEGREE_OF_COMPLETENESS,
    DEGREE_OF_NON_DETERMINISM,
    IUT_IS_REDUCTION,
    REMOVED_TRANSITIONS,
    FAIL_TRACE_FOUND,
    FAIL_TRACE_FOUND_SIZE,
    OBSERVED_TRACES_SIZE,
    LONGEST_OBSERVED_TRACE,
    LONGEST_OBSERVED_TRACE_SIZE,
    ADAPTIVE_STATE_COUNTING_RESULT,
    CREATE_RANDOM_FSM_SEED,
    CREATE_MUTANT_SEED,
    ITERATIONS,
    DURATION_MS,
    DURATION_M,
    PASS,
    END,
    BEGIN = TEST_NAME
};

struct CsvFieldHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

static map<CsvField, string> csvHeaders;
static map<TestIteration, string> csvIterations;
static unordered_map<TestIteration, unordered_map<CsvField, deque<string>, CsvFieldHash>, TestIterationHash> csvContent;

void initCsvHeaders()
{
    csvHeaders.insert(make_pair(CsvField::TEST_NAME, "testName"));
    csvHeaders.insert(make_pair(CsvField::NUM_STATES, "numStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_INPUTS, "numInputs"));
    csvHeaders.insert(make_pair(CsvField::NUM_OUTPUTS, "numOutputs"));
    csvHeaders.insert(make_pair(CsvField::NUM_D_REACHABLE_STATES, "numDReachableStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_SETS_OF_MAXIMAL_R_DIST_STATES, "numSetsOfMaximalRDistStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_OUT_FAULTS, "numOutFaults"));
    csvHeaders.insert(make_pair(CsvField::NUM_TRANS_FAULTS, "numTransFaults"));
    csvHeaders.insert(make_pair(CsvField::DEGREE_OF_COMPLETENESS, "degreeOfCompleteness"));
    csvHeaders.insert(make_pair(CsvField::DEGREE_OF_NON_DETERMINISM, "degreeOfNonDeterminism"));
    csvHeaders.insert(make_pair(CsvField::IUT_IS_REDUCTION, "iutIsReduction"));
    csvHeaders.insert(make_pair(CsvField::REMOVED_TRANSITIONS, "removedTransitions"));
    csvHeaders.insert(make_pair(CsvField::FAIL_TRACE_FOUND, "failTraceFound"));
    csvHeaders.insert(make_pair(CsvField::FAIL_TRACE_FOUND_SIZE, "failTraceFoundSize"));
    csvHeaders.insert(make_pair(CsvField::OBSERVED_TRACES_SIZE, "observedTracesSize"));
    csvHeaders.insert(make_pair(CsvField::LONGEST_OBSERVED_TRACE, "longestObservedTrace"));
    csvHeaders.insert(make_pair(CsvField::LONGEST_OBSERVED_TRACE_SIZE, "longestObservedTraceSize"));
    csvHeaders.insert(make_pair(CsvField::ADAPTIVE_STATE_COUNTING_RESULT, "adaptiveStateCountingResult"));
    csvHeaders.insert(make_pair(CsvField::CREATE_RANDOM_FSM_SEED, "createRandomFsmSeed"));
    csvHeaders.insert(make_pair(CsvField::CREATE_MUTANT_SEED, "createMutantSeed"));
    csvHeaders.insert(make_pair(CsvField::ITERATIONS, "iterations"));
    csvHeaders.insert(make_pair(CsvField::DURATION_MS, "durationMS"));
    csvHeaders.insert(make_pair(CsvField::DURATION_M, "durationM"));
    csvHeaders.insert(make_pair(CsvField::PASS, "pass"));
}

void initCsvIterations()
{
    csvIterations.insert(make_pair(TestIteration::INPUT, "Inputs"));
    csvIterations.insert(make_pair(TestIteration::OUTPUT, "Outputs"));
    csvIterations.insert(make_pair(TestIteration::STATE, "States"));
    csvIterations.insert(make_pair(TestIteration::OUTPUT_FAULT, "Outout Faults"));
    csvIterations.insert(make_pair(TestIteration::TRANSITION_FAULT, "Transition Faults"));
    csvIterations.insert(make_pair(TestIteration::DEGREE_COMPLETENESS, "Degree of Completeness"));
}

string initialize()
{
    initCsvHeaders();
    initCsvIterations();

    const system_clock::time_point now = system_clock::now();
    std::time_t tNow = system_clock::to_time_t(now);
    char nowTextRaw[21];
    struct tm buf;
    strftime(nowTextRaw, 21, "%Y-%m-%d--%H-%M-%S", localtime_r(&tNow, &buf));
    return string(nowTextRaw);
}

struct CsvConfig
{
    bool logEveryIteration = false;
    TestIteration context = TestIteration::END;
    vector<CsvField> fieldsContext;
    vector<CsvField> fieldsEveryIteration;
};

struct LoggingConfig
{
    bool toDot = false;
    bool toFsm = false;
    bool printSetsOfMaximalRDistStates = false;
    bool printObservedTraces = false;
};

struct AdaptiveTestConfigDebug
{
    string prefix = "DEBUG";
    int numStates;
    int numInput;
    int numOutput;
    int numOutFaults;
    int numTransFaults;
    unsigned int createRandomFsmSeed;
    unsigned int createMutantSeed;
    bool createReduction;
    float degreeOfCompleteness;
    float maxDegreeOfNonDeterminism;
    CsvConfig csvConfig;
    LoggingConfig loggingConfig;
};

struct AdaptiveTestConfig
{

    // Required
    string testName;
    int numFsm = -1;
    int maxInput = -1;
    int maxOutput = -1;
    int maxStates = -1;
    int maxOutFaults = -1;
    int maxTransFaults = -1;

    bool createReduction = false;

    float minDegreeOfCompleteness = 1.0f;
    float maxDegreeOfCompleteness = 1.0f;


    float maxDegreeOfNonDeterminism = 1.0f;

    // Optional
    CsvConfig csvConfig;
    int minStates = 2;
    int minInput = 2;
    int minOutput = 2;
    int minOutFaults = 1;
    int minTransFaults = 1;
    unsigned int seed = 0;
    bool dontTestReductions = false;
    // The mutant has to comply with the given parameters.
    bool forceTestParameters = true;
    LoggingConfig loggingConfig;
};

struct AdaptiveTestResult
{
    string testName;
    int numStates = -2;
    int numInputs = -2;
    int numOutputs = -2;
    int numDReachableStates = -1;
    vector<vector<shared_ptr<FsmNode>>> setsOfMaximalRDistStates;
    int numSetsOfMaximalRDistStates;
    int numOutFaults = -1;
    int numTransFaults = -1;
    float degreeOfCompleteness = -1;
    float degreeOfNonDeterminism = -1;
    bool iutIsReduction;
    int removedTransitions = -1;
    shared_ptr<IOTrace> failTraceFound;
    IOTraceContainer observedTraces;
    int numObservedTraces;
    shared_ptr<IOTrace> longestObservedTrace;
    bool adaptiveStateCountingResult;
    unsigned createRandomFsmSeed = 0;
    unsigned createMutantSeed = 0;
    int iterations = -1;
    long durationMS = -1;
    long durationM = -1;
    bool pass = 0;
};

void assertInconclusive(string tc, string comment = "") {
    
    string sVerdict("INCONCLUSIVE");
    CLOG(INFO, logging::testParameters) << sVerdict << ": " << tc << " : " << comment <<  endl;
}

void fsmlib_assert(string tc, bool verdict, string comment = "") {
    
    string out = (verdict) ? "PASS" : "FAIL";
    out += ": " + tc;
    if (!comment.empty())
    {
        out += ": " + comment;
    }

    CLOG(INFO, logging::testParameters) << out;
}

void assertOnFail(string tc, bool verdict, string comment = "") {

    string sVerdict = (verdict) ? "PASS: " + tc : "FAIL: " + tc + ": " + comment;
    if (verdict)
    {
        CLOG(INFO, logging::testParameters) << sVerdict;
    }
    else
    {
        CLOG(ERROR, logging::testParameters) << sVerdict;
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

string getFieldFromResult(const AdaptiveTestResult& result, const CsvField& field)
{
    std::stringstream out;
    switch (field) {
    case CsvField::TEST_NAME:
        out << result.testName;
        break;
    case CsvField::NUM_STATES:
        out << result.numStates;
        break;
    case CsvField::NUM_INPUTS:
        out << result.numInputs;
        break;
    case CsvField::NUM_OUTPUTS:
        out << result.numOutputs;
        break;
    case CsvField::NUM_D_REACHABLE_STATES:
        out << result.numDReachableStates;
        break;
    case CsvField::NUM_SETS_OF_MAXIMAL_R_DIST_STATES:
        out << result.numSetsOfMaximalRDistStates;
        break;
    case CsvField::NUM_OUT_FAULTS:
        out << result.numOutFaults;
        break;
    case CsvField::NUM_TRANS_FAULTS:
        out << result.numTransFaults;
        break;
    case CsvField::DEGREE_OF_COMPLETENESS:
        out << result.degreeOfCompleteness;
        break;
    case CsvField::DEGREE_OF_NON_DETERMINISM:
        out << result.degreeOfNonDeterminism;
        break;
    case CsvField::IUT_IS_REDUCTION:
        out << result.iutIsReduction;
        break;
    case CsvField::REMOVED_TRANSITIONS:
        out << result.removedTransitions;
        break;
    case CsvField::FAIL_TRACE_FOUND:
        if (result.failTraceFound)
        {
            out << *result.failTraceFound;
        }
        else
        {
            out << "";
        }
        break;
    case CsvField::FAIL_TRACE_FOUND_SIZE:
        if (result.failTraceFound)
        {
            out << result.failTraceFound->size();
        }
        else
        {
            out << "";
        }
        break;
    case CsvField::OBSERVED_TRACES_SIZE:
        out << result.numObservedTraces;
        break;
    case CsvField::LONGEST_OBSERVED_TRACE:
        if (result.longestObservedTrace)
        {
            out << *result.longestObservedTrace;
        }
        else
        {
            out << "";
        }
        break;
    case CsvField::LONGEST_OBSERVED_TRACE_SIZE:
        if (result.longestObservedTrace)
        {
            out << result.longestObservedTrace->size();
        }
        else
        {
            out << "";
        }
        break;
    case CsvField::ADAPTIVE_STATE_COUNTING_RESULT:
        out << result.adaptiveStateCountingResult;
        break;
    case CsvField::CREATE_RANDOM_FSM_SEED:
        out << result.createRandomFsmSeed;
        break;
    case CsvField::CREATE_MUTANT_SEED:
        out << result.createMutantSeed;
        break;
    case CsvField::ITERATIONS:
        out << result.iterations;
        break;
    case CsvField::DURATION_MS:
        out << result.durationMS;
        break;
    case CsvField::DURATION_M:
        out << result.durationM;
        break;
    case CsvField::PASS:
        out << result.pass;
        break;
    default:
        CLOG(FATAL, logging::globalLogger) << "Unhandled case in CSV Logging.";
    }
    return out.str();
}




void logToCsv(ofstream& csvOut, const AdaptiveTestResult& result, const CsvConfig& config)
{
    string output;


    size_t size = config.fieldsEveryIteration.size();
    if (size == 0)
    {
        size = csvHeaders.size() - 1;
        size_t i = 0;
        for (const auto& pair : csvHeaders)
        {
            const CsvField& f = pair.first;
            output += getFieldFromResult(result, f);
            if (i++ != size)
            {
                output += csvSep;
            }
        }
    }
    else
    {
        --size;
        for (size_t i = 0; i <= size; ++i)
        {
            CsvField field = config.fieldsEveryIteration.at(i);
            output += getFieldFromResult(result, field);
            if (i != size)
            {
                output += csvSep;
            }
        }
    }
    csvOut << output << endl;
}

void initQueue(const TestIteration& context, const CsvField& field, const int& rows)
{
    deque<string> queue;
    queue.push_back(csvIterations.at(context));

    for (size_t r = 0; r <= static_cast<size_t>(rows); ++r)
    {
        queue.push_back("");
    }
    csvContent[context][field] = queue;
}

void writeContextHeaderToQueue(const TestIteration& context, const CsvField& field, vector<string> values)
{
    deque<string>& queue = csvContent[context][field];

    string header;
    for (const string& s : values)
    {
        header += csvSep + s;
    }
    queue.at(0) += header;

}

void nextColumnInQueue(const TestIteration& context, const vector<CsvField>& fields)
{
    for (const CsvField& f : fields)
    {
        deque<string>& queue = csvContent[context][f];

        for (size_t i = 1; i < queue.size(); ++i)
        {
            queue.at(i) += csvSep;
        }
    }
}

void writeToQueue(const TestIteration& context, const CsvField& field,
                  const int& row, const string& value)
{
    csvContent[context][field].at(static_cast<size_t>(row) + 1) += value;
}

void writeResultToQueue(const AdaptiveTestResult& result, const TestIteration& context,
                        const vector<CsvField>& fields, const int& row)
{
    for (const CsvField& f : fields)
    {
        writeToQueue(context, f, row, getFieldFromResult(result, f));
    }
}

void flushQueues(string testName)
{
    for (auto contextIt : csvContent)
    {
        for (auto fieldIt : contextIt.second)
        {
            string fieldName = csvHeaders.at(fieldIt.first);
            deque<string>& queue = fieldIt.second;
            ofstream out(ascCsvDirectory + "Results-" + nowText + "-" + testName + "---" + fieldName + ".csv");

            for (const string& l : queue)
            {
                out << l << endl;
            }

            out.close();
        }
    }
}

void writeCsvHeader(ofstream& csvOut, const CsvConfig& config)
{
    string header = "";
    size_t size = config.fieldsEveryIteration.size();
    if (size == 0)
    {
        size = csvHeaders.size() - 1;
        size_t i = 0;
        for (const auto& pair : csvHeaders)
        {
            const string& h = pair.second;
            header += h;
            if (i++ != size)
            {
                header += csvSep;
            }
        }
    }
    else
    {
        --size;
        for (size_t i = 0; i <= size; ++i)
        {
            CsvField field = config.fieldsEveryIteration.at(i);
            header += csvHeaders.at(field);
            if (i != size)
            {
                header += csvSep;
            }
        }
    }
    csvOut << header << endl;
}

ofstream newCsvFile(string testName, const CsvConfig& csvConfig)
{
    ofstream csvOut(ascCsvDirectory + "Results-" + nowText + "-" + testName + ".csv");
    writeCsvHeader(csvOut, csvConfig);
    return csvOut;
}

void printTestBegin(string name)
{
    CLOG(INFO, logging::testParameters) << "#################### Test " << name << " ####################";
    CLOG(INFO, logging::globalLogger) << "-------------------- Test " << name << " --------------------";
}

void printSummary(const string& testName,
                  const int& executed,
                  const int& passed,
                  const int& notExecuted,
                  const long& durationS,
                  const long& durationM)
{
    CLOG(INFO, logging::globalLogger) << "#################### SUMMARY ####################";
    CLOG(INFO, logging::globalLogger) << "# Test name    : " << testName;
    CLOG(INFO, logging::globalLogger) << "# Total tests  : " << executed;
    CLOG(INFO, logging::globalLogger) << "# Passed       : " << passed;
    CLOG(INFO, logging::globalLogger) << "# Failed       : " << executed - passed;
    CLOG(INFO, logging::globalLogger) << "# Not executed : " << notExecuted;
    CLOG(INFO, logging::globalLogger) << "# Duration     : " << durationS << " s (" << durationM << " min).";
    CLOG(INFO, logging::globalLogger) << "#################################################" << endl;
}

void printTestConfig(const AdaptiveTestConfig& config)
{
    CLOG(INFO, logging::testParameters) << "-------------------- Test Config --------------------";
    CLOG(INFO, logging::testParameters) << "numFsm: " << config.numFsm;
    CLOG(INFO, logging::testParameters) << "minInput: " << config.minInput + 1;
    CLOG(INFO, logging::testParameters) << "maxInput: " << config.maxInput + 1;
    CLOG(INFO, logging::testParameters) << "minOutput: " << config.minOutput + 1;
    CLOG(INFO, logging::testParameters) << "maxOutput: " << config.maxOutput + 1;
    CLOG(INFO, logging::testParameters) << "minStates: " << config.minStates + 1;
    CLOG(INFO, logging::testParameters) << "maxStates: " << config.maxStates + 1;

    CLOG(INFO, logging::testParameters) << "minOutFaults: " << config.minOutFaults;
    CLOG(INFO, logging::testParameters) << "maxOutFaults: " << config.maxOutFaults;
    CLOG(INFO, logging::testParameters) << "minTransFaults: " << config.minTransFaults;
    CLOG(INFO, logging::testParameters) << "maxTransFaults: " << config.maxTransFaults;

    CLOG(INFO, logging::testParameters) << "minDegreeOfCompleteness: " << config.minDegreeOfCompleteness;
    CLOG(INFO, logging::testParameters) << "maxDegreeOfCompleteness: " << config.maxDegreeOfCompleteness;
    CLOG(INFO, logging::testParameters) << "maxDegreeOfNonDeterminism: " << config.maxDegreeOfNonDeterminism;

    CLOG(INFO, logging::testParameters) << "dontTestReductions: " << std::boolalpha << config.dontTestReductions;
    CLOG(INFO, logging::testParameters) << "forceTestParameters: " << std::boolalpha << config.forceTestParameters;
    CLOG(INFO, logging::testParameters) << "seed: " << config.seed;
    CLOG(INFO, logging::testParameters) << "--------------------------------------------------------";
}

void printTestResult(AdaptiveTestResult& result, const CsvConfig& csvConfig,
                     const LoggingConfig& loggingConfig, ofstream& csvOut)
{

    CLOG(INFO, logging::testParameters) << "Test                       : " << result.testName;
    CLOG(INFO, logging::testParameters) << "numStates                  : " << result.numStates;
    CLOG(INFO, logging::testParameters) << "numInputs                  : " << result.numInputs;
    CLOG(INFO, logging::testParameters) << "numOutputs                 : " << result.numOutputs;
    CLOG(INFO, logging::testParameters) << "numDReachableStates        : " << result.numDReachableStates;


    if (loggingConfig.printSetsOfMaximalRDistStates)
    {
        CLOG(INFO, logging::testParameters) << "setsOfMaximalRDistStates   : ";
        for (auto v : result.setsOfMaximalRDistStates)
        {
            stringstream ss;
            ss << "{";
            for (auto e : v)
            {
                ss << e->getName() << ", ";
            }
            ss << "}";
            CLOG(INFO, logging::testParameters) << ss.str();
        }
    }

    CLOG(INFO, logging::testParameters) << "numSetsOfMaximalRDistStates: " << result.numSetsOfMaximalRDistStates;
    CLOG(INFO, logging::testParameters) << "numOutFaults               : " << result.numOutFaults;
    CLOG(INFO, logging::testParameters) << "numTransFaults             : " << result.numTransFaults;
    CLOG(INFO, logging::testParameters) << "degreeOfCompleteness       : " << result.degreeOfCompleteness;
    CLOG(INFO, logging::testParameters) << "degreeOfNonDeterminism     : " << result.degreeOfNonDeterminism;
    CLOG(INFO, logging::testParameters) << "iutIsReduction             : " << std::boolalpha << result.iutIsReduction;
    CLOG(INFO, logging::testParameters) << "removedTransitions         : " << result.removedTransitions;
    if (result.failTraceFound)
    {
        CLOG(INFO, logging::testParameters) << "failTraceFound             : " << *result.failTraceFound;
    }
    else
    {
        CLOG(INFO, logging::testParameters) << "failTraceFound             : None";
    }
    CLOG(INFO, logging::testParameters) << "adaptiveStateCountingResult: " << std::boolalpha
                                      << result.adaptiveStateCountingResult;
    if (result.createRandomFsmSeed != 0)
    {
        CLOG(INFO, logging::testParameters) << "createRandomFsmSeed        : " << result.createRandomFsmSeed;
    }
    if (result.createMutantSeed != 0)
    {
        CLOG(INFO, logging::testParameters) << "createMutantSeed           : " << result.createMutantSeed;
    }
    CLOG(INFO, logging::testParameters) << "longestObservedTrace       : " << *result.longestObservedTrace;
    CLOG(INFO, logging::testParameters) << "length longestObservedTrace: " << result.longestObservedTrace->size();
    CLOG(INFO, logging::testParameters) << "observedTraces size        : " << result.numObservedTraces;
    if (loggingConfig.printObservedTraces)
    {
        CLOG(INFO, logging::testParameters) << "observedTraces             : " << result.observedTraces;
    }
    CLOG(INFO, logging::testParameters) << "iterations                 : " << result.iterations;
    CLOG(INFO, logging::testParameters) << "Calculation took " << result.durationMS << " ms ("
                                      << result.durationM << " minutes).";
    CLOG(INFO, logging::testParameters) << "-------------------------------------------";

    if (csvConfig.logEveryIteration)
    {
        logToCsv(csvOut, result, csvConfig);
    }
}

void printTestResult(AdaptiveTestResult& result, const CsvConfig& csvConfig,
                     const LoggingConfig& loggingConfig)
{
    ofstream dummyout = ofstream();
    printTestResult(result, csvConfig, loggingConfig, dummyout);
}

bool isReduction(Fsm& spec, Fsm& iut, string intersectionName, shared_ptr<Fsm>& intersection)
{
    Fsm inter = spec.intersect(iut, intersectionName);
    intersection = make_shared<Fsm>(inter);
    return !inter.hasFailure();
}

void executeAdaptiveTest(const string& testName, Fsm& spec, Fsm& iut, size_t m, string intersectionName,
                         const bool& toDot, const bool& toFsm, const bool& dontTestReductions, AdaptiveTestResult& result)
{
    Fsm specMin = spec.minimise("", "", false);
    Fsm iutMin = iut.minimise("", "", false);

    result.degreeOfCompleteness = specMin.getDegreeOfCompleteness();
    result.degreeOfNonDeterminism = specMin.getDegreeOfNonDeterminism();

    result.numStates = specMin.getMaxNodes();
    result.numInputs = specMin.getMaxInput() + 1;
    result.numOutputs = specMin.getMaxOutput() + 1;

    shared_ptr<Fsm> intersection;
    result.iutIsReduction = isReduction(specMin, iutMin, intersectionName, intersection);


    if (toDot)
    {
        spec.toDot(ascTestResultDirectory + testName + "-" + spec.getName());
        iut.toDot(ascTestResultDirectory + testName + "-"  + iut.getName());
        specMin.toDot(ascTestResultDirectory + testName + "-"  + spec.getName() + "-min");
        iutMin.toDot(ascTestResultDirectory + testName + "-"  + iut.getName() + "-min");
        intersection->toDot(ascTestResultDirectory + testName + "-"  + intersection->getName());
    }
    if (toFsm)
    {
        ofstream out(ascTestResultDirectory + spec.getName() + ".fsm");
        spec.dumpFsm(out);
        out.close();

        out.open(ascTestResultDirectory + iut.getName() + ".fsm");
        iut.dumpFsm(out);
        out.close();
    }

    if (result.iutIsReduction && dontTestReductions) {
        CLOG_IF(VLOG_IS_ON(1), INFO, logging::globalLogger) << "Won't test this one, since it is a reduction.";
        throw unexpected_reduction("Interrupting testing, since IUT is an unexcpected reduction of the specification.");
    }

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    result.adaptiveStateCountingResult = Fsm::adaptiveStateCounting(specMin, iutMin, m,
                                                                    result.observedTraces, result.failTraceFound, result.iterations);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    result.durationMS = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    result.durationM = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();

    if (result.failTraceFound)
    {
        result.failTraceFound = make_shared<IOTrace>(result.failTraceFound->removeLeadingEpsilons());
    }

    result.pass = (result.iutIsReduction == result.adaptiveStateCountingResult);
    result.setsOfMaximalRDistStates = specMin.getMaximalSetsOfRDistinguishableStates();
    result.numSetsOfMaximalRDistStates = static_cast<int>(result.setsOfMaximalRDistStates.size());
    result.numDReachableStates = static_cast<int>(specMin.getDReachableStates().size());
    result.longestObservedTrace = result.observedTraces.getLongestTrace();
    result.numObservedTraces = static_cast<int>(result.observedTraces.size());
}


void createAndExecuteAdaptiveTest(
        ofstream& csvOut,
        const string& prefix,
        const int numStates,
        const int numInputs,
        const int numOutputs,
        const int numOutFaults,
        const int numTransFaults,
        const float degreeOfCompleteness,
        const float maxDegreeOfNonDeterminism,
        const unsigned int createRandomFsmSeed,
        const unsigned int createMutantSeed,
        const shared_ptr<FsmPresentationLayer>& plSpec,
        const shared_ptr<FsmPresentationLayer>& plIut,
        const bool createReduction,
        const bool dontTestReductions,
        const CsvConfig& csvConfig,
        const LoggingConfig loggingConfig,
        AdaptiveTestResult& result)
{

    CLOG(INFO, logging::testParameters) << testSepLine;
    CLOG(INFO, logging::testParameters) << "createAndExecuteAdaptiveTest()";
    CLOG(INFO, logging::testParameters) << "Name                     : " << result.testName;
    CLOG(INFO, logging::testParameters) << "numStates                : " << numStates + 1;
    CLOG(INFO, logging::testParameters) << "numInputs                : " << numInputs + 1;
    CLOG(INFO, logging::testParameters) << "numOutputs               : " << numOutputs + 1;
    CLOG(INFO, logging::testParameters) << "numOutFaults             : " << numOutFaults;
    CLOG(INFO, logging::testParameters) << "numTransFaults           : " << numTransFaults;
    CLOG(INFO, logging::testParameters) << "degreeOfCompleteness     : " << degreeOfCompleteness;
    CLOG(INFO, logging::testParameters) << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism;
    CLOG(INFO, logging::testParameters) << "createReduction          : " << boolalpha << createReduction;
    CLOG(INFO, logging::testParameters) << "createRandomFsmSeed      : " << createRandomFsmSeed;
    CLOG(INFO, logging::testParameters) << "createMutantSeed         : " << createMutantSeed;
    CLOG(INFO, logging::testParameters) << "-------------------------------------------";

    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating FSM.";
    shared_ptr<Fsm> spec = Fsm::createRandomFsm(prefix + "-spec",
                                                numInputs,
                                                numOutputs,
                                                numStates,
                                                plSpec,
                                                degreeOfCompleteness,
                                                maxDegreeOfNonDeterminism,
                                                createReduction,
                                                true,
                                                true,
                                                createRandomFsmSeed);

    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating mutant.";

    shared_ptr<Fsm> iut;
    if (createReduction)
    {
        iut = spec->createReduction(prefix + "-iut", true, result.removedTransitions, createMutantSeed, plIut);
    }
    else
    {
        iut = spec->createMutant(prefix + "-iut",
                                 numOutFaults,
                                 numTransFaults,
                                 true,
                                 createMutantSeed,
                                 plIut);
        result.removedTransitions = 0;
    }

    result.numStates = numStates + 1;
    result.numInputs = numInputs + 1;
    result.numOutputs = numOutputs + 1;
    result.numOutFaults = numOutFaults;
    result.numTransFaults = numTransFaults;
    result.createRandomFsmSeed = createRandomFsmSeed;
    result.createMutantSeed = createMutantSeed;


    CLOG(INFO, logging::testParameters) << "Starting adaptive state counting.";
    CLOG(INFO, logging::testParameters) << "Name                     : " << result.testName;
    CLOG(INFO, logging::testParameters) << "numStates                : " << result.numStates;
    CLOG(INFO, logging::testParameters) << "numInputs                : " << result.numInputs;
    CLOG(INFO, logging::testParameters) << "numOutputs               : " << result.numOutputs;
    CLOG(INFO, logging::testParameters) << "numOutFaults             : " << result.numOutFaults;
    CLOG(INFO, logging::testParameters) << "numTransFaults           : " << result.numTransFaults;
    CLOG(INFO, logging::testParameters) << "degreeOfCompleteness tgt : " << degreeOfCompleteness;
    CLOG(INFO, logging::testParameters) << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism;
    CLOG(INFO, logging::testParameters) << "createReduction          : " << boolalpha << createReduction;
    CLOG(INFO, logging::testParameters) << "removedTransitions       : " << result.removedTransitions;
    CLOG(INFO, logging::testParameters) << "createRandomFsmSeed      : " << result.createRandomFsmSeed;
    CLOG(INFO, logging::testParameters) << "createMutantSeed         : " << result.createMutantSeed;
    CLOG(INFO, logging::testParameters) << "-------------------------------------------";

    executeAdaptiveTest(result.testName, *spec, *iut, static_cast<size_t>(iut->getMaxNodes()),
                        prefix + "-intersect", loggingConfig.toDot, loggingConfig.toFsm, dontTestReductions, result);

    printTestResult(result, csvConfig, loggingConfig, csvOut);
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

void adaptiveTestRandom(AdaptiveTestConfig& config)
{

    logging::setLogfileSuffix(config.testName, logging::testParameters);
    printTestBegin(config.testName);
    ofstream csvOut = newCsvFile(config.testName, config.csvConfig);

    std::chrono::steady_clock::time_point totalStart = std::chrono::steady_clock::now();

    if (config.numFsm < 0 ||
            config.maxInput < 0 ||
            config.maxOutput < 0 ||
            config.maxStates < 0 ||
            config.maxOutFaults < 0 ||
            config.maxTransFaults < 0 ||
            config.maxDegreeOfNonDeterminism < 0)
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

    const int numberDigits = ((config.numFsm <= 1)? 1 : static_cast<int>(log10(config.numFsm)) + 1);


    const int diffInput = config.maxInput - config.minInput + 1;
    const int diffOutput = config.maxOutput - config.minOutput + 1;
    const int diffStates = config.maxStates - config.minStates + 1;
    const int diffOutFaults = config.maxOutFaults - config.minOutFaults + 1;
    const int diffTransFaults = config.maxTransFaults - config.minTransFaults + 1;
    int degreeOfCompletenessIterations = 1;

    if (config.minDegreeOfCompleteness > 0 && config.maxDegreeOfCompleteness <= 1)
    {
        degreeOfCompletenessIterations = static_cast<int>(config.maxDegreeOfCompleteness - config.minDegreeOfCompleteness) * 10 + 1;
    }

    if (diffInput <= 0 || diffOutput <= 0 || diffStates <= 0 || diffOutFaults <= 0 || diffTransFaults <= 0 || degreeOfCompletenessIterations <= 0)
    {
        CLOG(FATAL, logging::globalLogger) << "Please check the test parameters.";
    }

    if (config.maxDegreeOfNonDeterminism <= 0 && config.createReduction)
    {
        CLOG(FATAL, logging::globalLogger) << "Can't create reduction of a deterministic FSM.";
    }

    float divisor = (diffInput * diffOutput * diffStates * diffOutFaults * diffTransFaults * degreeOfCompletenessIterations);
    int innerIterations = static_cast<int>(ceil(static_cast<float>(config.numFsm) / divisor));
    int totalIterations = static_cast<int>(innerIterations * divisor);

    printTestConfig(config);

    CLOG(INFO, logging::testParameters) << "Seed           : " << config.seed;
    CLOG(INFO, logging::testParameters) << "divisor        : " << divisor;
    CLOG(INFO, logging::testParameters) << "innerIterations: " << innerIterations;
    CLOG(INFO, logging::testParameters) << "totalIterations: " << totalIterations;
    CLOG(INFO, logging::testParameters) << "";
    CLOG(INFO, logging::testParameters) << "diffInput      : " << diffInput;
    CLOG(INFO, logging::testParameters) << "diffOutput     : " << diffOutput;
    CLOG(INFO, logging::testParameters) << "diffStates     : " << diffStates;
    CLOG(INFO, logging::testParameters) << "diffOutFaults  : " << diffOutFaults;
    CLOG(INFO, logging::testParameters) << "diffTransFaults: " << diffTransFaults;
    CLOG(INFO, logging::testParameters) << "degreeOfCompletenessIterations: " << degreeOfCompletenessIterations;
    CLOG(INFO, logging::testParameters) << "";

    const float degStep = 0.1f;

    unordered_map<CsvField, vector<string>, CsvFieldHash> csvHeaderValues;
    if (config.csvConfig.context != TestIteration::END)
    {
        for (const CsvField& f : config.csvConfig.fieldsContext)
        {
            initQueue(config.csvConfig.context, f, innerIterations);
        }
    }

    std::mt19937 gen(config.seed);

    int executed = 0;
    int passed = 0;
    int i = 0;

    const TestIteration& loggingContext = config.csvConfig.context;

    for (int ctInput = config.minInput; ctInput <= config.maxInput; ++ctInput)
    {
        if (loggingContext == TestIteration::INPUT)
        {
            nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
        }
        for (int ctOutput = config.minOutput; ctOutput <= config.maxOutput; ++ctOutput)
        {
            if (loggingContext == TestIteration::OUTPUT)
            {
                nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
            }
            for (int ctState = config.minStates; ctState <= config.maxStates ; ++ctState)
            {
                if (loggingContext == TestIteration::STATE)
                {
                    nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
                }
                for (int ctOutFault = config.minOutFaults; ctOutFault <= config.maxOutFaults; ++ctOutFault)
                {
                    if (loggingContext == TestIteration::OUTPUT_FAULT)
                    {
                        nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
                    }
                    for (int ctTransFault = config.minTransFaults; ctTransFault <= config.maxTransFaults; ++ctTransFault)
                    {
                        if (loggingContext == TestIteration::TRANSITION_FAULT)
                        {
                            nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
                        }
                        for (int degreeCount = 0; degreeCount < degreeOfCompletenessIterations; ++degreeCount)
                        {
                            if (loggingContext == TestIteration::DEGREE_COMPLETENESS)
                            {
                                nextColumnInQueue(config.csvConfig.context, config.csvConfig.fieldsContext);
                            }
                            float degreeOfCompleteness;

                            if (config.minDegreeOfCompleteness <= 0 && config.maxDegreeOfCompleteness <= 0)
                            {
                                degreeOfCompleteness = -1;
                                CLOG(INFO, logging::testParameters) << "Degree of completeness does not matter.";
                            }
                            else if (config.minDegreeOfCompleteness > 0 && config.maxDegreeOfCompleteness <= 0)
                            {
                                degreeOfCompleteness = config.minDegreeOfCompleteness +
                                        static_cast<float>(rand()) / (static_cast <float> (RAND_MAX/(1.0f - config.minDegreeOfCompleteness)));
                                CLOG(INFO, logging::testParameters) << "Selected random degree of completeness with a minimal value of " <<
                                                                     config.minDegreeOfCompleteness << ": " << degreeOfCompleteness;
                            }
                            else if (config.minDegreeOfCompleteness <= 0 && config.maxDegreeOfCompleteness <=1)
                            {
                                degreeOfCompleteness = 0.1f +
                                        static_cast<float>(rand()) / (static_cast <float> (RAND_MAX/(1.0f - 0.1f)));
                                CLOG(INFO, logging::testParameters) << "Selected random degree of completeness with a maximal value of " <<
                                                                     config.maxDegreeOfCompleteness << ": " << degreeOfCompleteness;
                            }
                            else
                            {

                                degreeOfCompleteness = config.minDegreeOfCompleteness + (degStep * degreeCount);
                                CLOG(INFO, logging::testParameters) << "Degree of completeness: " << degreeOfCompleteness;
                            }

                            for (int ctInnerIt = 0; ctInnerIt < innerIterations; ++ctInnerIt)
                            {
                                CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger)
                                        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                                CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << ctInput << " "
                                                                                    << ctOutput << " "
                                                                                    << ctState << " "
                                                                                    << ctOutFault << " "
                                                                                    << ctTransFault << " "
                                                                                    << ctInnerIt;

                                shared_ptr<FsmPresentationLayer> plTestSpecCopy = make_shared<FsmPresentationLayer>(*plTestSpec);
                                shared_ptr<FsmPresentationLayer> plTestIutCopy = make_shared<FsmPresentationLayer>(*plTestIut);

                                stringstream ss;
                                ss << setw(numberDigits) << setfill('0') << i;
                                string iteration = ss.str();

                    #ifdef ENABLE_DEBUG_MACRO
                                logging::setLogfileSuffix(config.testName + "-" + iteration);
                    #endif


                                int numInputs = ctInput;
                                int numOutputs = ctOutput;
                                int numStates = ctState;
                                int numOutFaults = ctOutFault;
                                int numTransFaults = ctTransFault;


                                if (numOutputs < 1 && numOutFaults > 0)
                                {
                                    CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
                                    CLOG(INFO, logging::globalLogger) << "numInput: " << numInputs + 1;
                                    CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutputs + 1;
                                    CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
                                    CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
                                    CLOG(WARNING, logging::globalLogger) << "Too little outputs. Can not create requested number of "
                                                                         << "output faults. Could not create mutant. Skipping.";
                                    ++i;
                                    continue;
                                }

                                unsigned int createRandomFsmSeed = static_cast<unsigned int>(getRandom(gen));
                                unsigned int createMutantSeed = static_cast<unsigned int>(getRandom(gen));

                                float maxDegNonDet = config.maxDegreeOfNonDeterminism * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                                if (config.createReduction && maxDegNonDet < 0.5f)
                                {
                                    maxDegNonDet = 0.5f + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1.0f - 0.5f)));
                                    if (config.maxDegreeOfNonDeterminism < 0.5f || maxDegNonDet > config.maxDegreeOfNonDeterminism)
                                    {
                                        CLOG_IF(VLOG_IS_ON(2), WARNING, logging::globalLogger)
                                                << "Chosen maximal degree of non-determinism very low. Adjusting.";
                                    }
                                }

                                TIMED_SCOPE(timerBlkObj, "heavy-iter");
                                AdaptiveTestResult result;
                                result.testName = config.testName + "-" + iteration;

                                CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "testName: " << result.testName;

                                bool couldCreateIut = false;
                                bool abort = false;
                                bool newSeeds = false;
                                do
                                {
                                    if (newSeeds)
                                    {
                                        CLOG_IF(VLOG_IS_ON(1), INFO, logging::globalLogger) << "Generating new seeds.";
                                        createRandomFsmSeed = static_cast<unsigned int>(getRandom(gen));
                                        createMutantSeed = static_cast<unsigned int>(getRandom(gen));
                                        newSeeds = false;
                                    }
                                    try {
                                        createAndExecuteAdaptiveTest(
                                                    csvOut,
                                                    iteration,
                                                    numStates,
                                                    numInputs,
                                                    numOutputs,
                                                    numOutFaults,
                                                    numTransFaults,
                                                    degreeOfCompleteness,
                                                    maxDegNonDet,
                                                    createRandomFsmSeed,
                                                    createMutantSeed,
                                                    plTestSpecCopy,
                                                    plTestIutCopy,
                                                    config.createReduction,
                                                    config.dontTestReductions,
                                                    config.csvConfig,
                                                    config.loggingConfig,
                                                    result);
                                        couldCreateIut = true;

                                        if (config.csvConfig.context != TestIteration::END)
                                        {
                                            writeResultToQueue(result, config.csvConfig.context, config.csvConfig.fieldsContext, ctInnerIt);

                                        }
                                        result.setsOfMaximalRDistStates.clear();
                                        result.observedTraces.clear();

                                    }
                                    catch (unexpected_reduction& e)
                                    {
                                        CLOG_IF(VLOG_IS_ON(1), INFO, logging::globalLogger) << "IUT is a reduction of the specification.";
                                        if (config.forceTestParameters)
                                        {
                                            newSeeds = true;
                                        }
                                        else
                                        {
                                            abort = true;
                                        }
                                    }
                                    catch (too_many_transition_faults& e)
                                    {
                                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create mutant.";
                                        if (!config.forceTestParameters
                                                && numTransFaults - 1 >= config.minTransFaults && numTransFaults - 1 > 0) {
                                            --numTransFaults;
                                            CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing transition faults.";
                                            continue;
                                        }
                                        else if (config.forceTestParameters)
                                        {
                                            newSeeds = true;

                                        }
                                        else
                                        {
                                            abort = true;
                                        }
                                    }
                                    catch (too_many_output_faults& e)
                                    {
                                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create mutant.";
                                        if (!config.forceTestParameters
                                                && numOutFaults - 1 >= config.minOutFaults && numOutFaults - 1 > 0) {
                                            --numOutFaults;
                                            CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing output faults.";
                                            continue;
                                        }
                                        else if (!config.forceTestParameters
                                                 && numTransFaults - 1 >= config.minTransFaults && numTransFaults - 1 > 0)
                                        {
                                            --numTransFaults;
                                            CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing transition faults.";
                                            continue;
                                        }
                                        else if (config.forceTestParameters)
                                        {
                                            newSeeds = true;

                                        }
                                        else
                                        {
                                            abort = true;
                                        }
                                    }
                                    catch (reduction_not_possible& e)
                                    {
                                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create reduction.";
                                        newSeeds = true;
                                    }
                                } while (!couldCreateIut && !abort);
                                if (!couldCreateIut)
                                {
                                    ++i;
                                    CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
                                    CLOG(INFO, logging::globalLogger) << "numInput: " << numInputs + 1;
                                    CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutputs + 1;
                                    CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
                                    CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
                                    CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed: " << createRandomFsmSeed;
                                    CLOG(INFO, logging::globalLogger) << "createMutantSeed: " << createMutantSeed;
                                    CLOG(WARNING, logging::globalLogger) << "Could not create requested mutant. Skipping.";
                                    continue;
                                }

                                if (result.pass)
                                {
                                    ++passed;
                                }

                                fsmlib_assert(result.testName, result.pass);
                                ++executed;
                                ++i;
                            }  // End inner loop

                            if (loggingContext == TestIteration::DEGREE_COMPLETENESS)
                            {
                                for (const CsvField& f : config.csvConfig.fieldsContext)
                                {
                                    csvHeaderValues[f].push_back(to_string(static_cast<int>((degreeOfCompleteness * 10)) * 10));
                                }
                            }

                        } // End degree of completeness loop

                        if (loggingContext == TestIteration::TRANSITION_FAULT)
                        {
                            for (const CsvField& f : config.csvConfig.fieldsContext)
                            {
                                csvHeaderValues[f].push_back(to_string(ctTransFault));
                            }
                        }
                    }  // End trans faults

                    if (loggingContext == TestIteration::OUTPUT_FAULT)
                    {
                        for (const CsvField& f : config.csvConfig.fieldsContext)
                        {
                            csvHeaderValues[f].push_back(to_string(ctOutFault));
                        }
                    }
                }  // End out faults

                if (loggingContext == TestIteration::STATE)
                {
                    for (const CsvField& f : config.csvConfig.fieldsContext)
                    {
                        csvHeaderValues[f].push_back(to_string(ctState + 1));
                    }
                }
            } // End states

            if (loggingContext == TestIteration::OUTPUT)
            {
                for (const CsvField& f : config.csvConfig.fieldsContext)
                {
                    csvHeaderValues[f].push_back(to_string(ctOutput + 1));
                }
            }
        } // End outputs

        if (loggingContext == TestIteration::INPUT)
        {
            for (const CsvField& f : config.csvConfig.fieldsContext)
            {
                csvHeaderValues[f].push_back(to_string(ctInput + 1));
            }
        }
    } // End inputs


    std::chrono::steady_clock::time_point totalEnd = std::chrono::steady_clock::now();
    long durationS = std::chrono::duration_cast<std::chrono::seconds>(totalEnd - totalStart).count();
    long durationM = std::chrono::duration_cast<std::chrono::minutes>(totalEnd - totalStart).count();

    for (const CsvField& f : config.csvConfig.fieldsContext)
    {
        writeContextHeaderToQueue(config.csvConfig.context, f, csvHeaderValues[f]);
    }

    flushQueues(config.testName);
    printSummary(config.testName, executed, passed, i - executed, durationS, durationM);
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

void trial(bool debug)
{
    AdaptiveTestConfig config;
    config.testName = "TRIAL";
    config.numFsm = 1000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 10;
    config.maxStates = 10;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 1;
    config.maxOutFaults = 1;

    config.maxDegreeOfNonDeterminism = 1.0f;

    config.dontTestReductions = true;
    config.forceTestParameters = true;

    config.seed = 1337;

    config.csvConfig.logEveryIteration = true;

    config.loggingConfig.toDot = false;
    config.loggingConfig.printSetsOfMaximalRDistStates = false;

    if (debug)
    {
        CLOG(INFO, logging::globalLogger) << "############## Debugging ##############";

        AdaptiveTestConfigDebug debugConfig;
        debugConfig.numStates = 4;
        debugConfig.numInput = 4;
        debugConfig.numOutput = 1;
        debugConfig.numOutFaults = 0;
        debugConfig.numTransFaults = 0;
        debugConfig.createRandomFsmSeed = 98824417;
        debugConfig.createMutantSeed = 1642605748;
        debugConfig.createReduction = false;
        debugConfig.degreeOfCompleteness = 1.0f;
        debugConfig.maxDegreeOfNonDeterminism = 0.508369f;
        debugConfig.csvConfig.logEveryIteration = false;
        debugConfig.loggingConfig.toDot = false;

        shared_ptr<FsmPresentationLayer> plTestSpecCopy = make_shared<FsmPresentationLayer>(*plTestSpec);
        shared_ptr<FsmPresentationLayer> plTestIutCopy = make_shared<FsmPresentationLayer>(*plTestIut);

        AdaptiveTestResult result;

        ofstream dummyStream = ofstream();

        createAndExecuteAdaptiveTest(
                    dummyStream,
                    debugConfig.prefix,
                    debugConfig.numStates - 1,
                    debugConfig.numInput - 1,
                    debugConfig.numOutput - 1,
                    debugConfig.numOutFaults,
                    debugConfig.numTransFaults,
                    debugConfig.degreeOfCompleteness,
                    debugConfig.maxDegreeOfNonDeterminism,
                    debugConfig.createRandomFsmSeed,
                    debugConfig.createMutantSeed,
                    plTestSpecCopy,
                    plTestIutCopy,
                    debugConfig.createReduction,
                    false,
                    debugConfig.csvConfig,
                    debugConfig.loggingConfig,
                    result);

        assertOnFail(debugConfig.prefix, result.pass);
    }
    else
    {
        adaptiveTestRandom(config);
    }
}

void test00_00()
{
    CsvConfig csvConfig;
    csvConfig.logEveryIteration = false;

    LoggingConfig loggingConfig;

    AdaptiveTestResult result;
    result.testName = "00-00";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/00-spec.dot", "00-00-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/00-iut.dot", "00-00-iut");
    result.numOutFaults = 0;
    result.numTransFaults = 0;
    executeAdaptiveTest(result.testName, spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-00-inter", true, false, false, result);
    printTestResult(result, csvConfig, loggingConfig);
    fsmlib_assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

void test00_01()
{
    CsvConfig csvConfig;
    csvConfig.logEveryIteration = false;

    LoggingConfig loggingConfig;

    AdaptiveTestResult result;
    result.testName = "00-01";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/01-spec.dot", "00-01-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/01-iut.dot", "00-01-iut");

    result.numOutFaults = 1;
    result.numTransFaults = 0;
    executeAdaptiveTest(result.testName, spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-01-inter", true, false, false, result);
    printTestResult(result, csvConfig, loggingConfig);
    fsmlib_assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

/**
 * Testing FSMs with number of states varying
 * and always one output fault.
 */
void OUTF_STA()
{
    AdaptiveTestConfig config;
    config.testName = "OUTF-STA";
    config.numFsm = 20000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 1;
    config.maxStates = 20;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 1;
    config.maxOutFaults = 1;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 1337;

    adaptiveTestRandom(config);
}

/**
 * Testing FSMs with number of inputs varying
 * and always one output fault.
 */
void OUTF_INP()
{
    AdaptiveTestConfig config;
    config.testName = "OUTF-INP";
    config.numFsm = 35000;

    config.minInput = 1;
    config.maxInput = 35;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 5;
    config.maxStates = 5;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 1;
    config.maxOutFaults = 1;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 7364746;

    adaptiveTestRandom(config);
}

void OUTF_OUTP()
{
    AdaptiveTestConfig config;
    config.testName = "OUTF-OUTP";
    config.numFsm = 34000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 2;
    config.maxOutput = 35;

    config.minStates = 5;
    config.maxStates = 5;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 1;
    config.maxOutFaults = 1;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 895;

    adaptiveTestRandom(config);
}

void TRANSF_STA()
{
    AdaptiveTestConfig config;
    config.testName = "TRANSF-STA";
    config.numFsm = 13000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 2;
    config.maxStates = 14;

    config.minTransFaults = 1;
    config.maxTransFaults = 1;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 7331;

    adaptiveTestRandom(config);
}

void TRANSF_INP()
{
    AdaptiveTestConfig config;
    config.testName = "TRANSF-INP";
    config.numFsm = 35000;

    config.minInput = 1;
    config.maxInput = 35;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 5;
    config.maxStates = 5;

    config.minTransFaults = 1;
    config.maxTransFaults = 1;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 89786;

    adaptiveTestRandom(config);
}

void TRANSF_OUTP()
{
    AdaptiveTestConfig config;
    config.testName = "TRANSF-OUTP";
    config.numFsm = 34000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 2;
    config.maxOutput = 35;

    config.minStates = 5;
    config.maxStates = 5;

    config.minTransFaults = 1;
    config.maxTransFaults = 1;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = true;

    config.csvConfig.logEveryIteration = true;

    config.seed = 101;

    adaptiveTestRandom(config);
}

void AEQ_STA()
{
    AdaptiveTestConfig config;
    config.testName = "AEQ-STA";
    config.numFsm = 3000;

    config.minInput = 2;
    config.maxInput = 2;

    config.minOutput = 2;
    config.maxOutput = 2;

    config.minStates = 1;
    config.maxStates = 3;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 785676;

    adaptiveTestRandom(config);
}

void AEQ_INP()
{
    AdaptiveTestConfig config;
    config.testName = "AEQ-INP";
    config.numFsm = 25000;

    config.minInput = 1;
    config.maxInput = 25;

    config.minOutput = 2;
    config.maxOutput = 2;

    config.minStates = 2;
    config.maxStates = 2;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 546457;

    adaptiveTestRandom(config);
}

void AEQ_OUTP()
{
    AdaptiveTestConfig config;
    config.testName = "AEQ-OUTP";
    config.numFsm = 24000;

    config.minInput = 2;
    config.maxInput = 2;

    config.minOutput = 2;
    config.maxOutput = 25;

    config.minStates = 2;
    config.maxStates = 2;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 805645;

    adaptiveTestRandom(config);
}

void RED_STA()
{
    AdaptiveTestConfig config;
    config.testName = "RED-STA";
    config.numFsm = 3000;

    config.minInput = 2;
    config.maxInput = 2;

    config.minOutput = 2;
    config.maxOutput = 2;

    config.minStates = 1;
    config.maxStates = 3;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.createReduction = true;
    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 4334;

    adaptiveTestRandom(config);
}

void RED_INP()
{
    AdaptiveTestConfig config;
    config.testName = "RED-INP";
    config.numFsm = 25000;

    config.minInput = 1;
    config.maxInput = 25;

    config.minOutput = 2;
    config.maxOutput = 2;

    config.minStates = 2;
    config.maxStates = 2;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.createReduction = true;
    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 56457;

    adaptiveTestRandom(config);
}

void RED_OUTP()
{
    AdaptiveTestConfig config;
    config.testName = "RED-OUTP";
    config.numFsm = 24000;

    config.minInput = 2;
    config.maxInput = 2;

    config.minOutput = 2;
    config.maxOutput = 25;

    config.minStates = 2;
    config.maxStates = 2;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.createReduction = true;
    config.dontTestReductions = false;

    config.csvConfig.logEveryIteration = true;

    config.seed = 7896774;

    adaptiveTestRandom(config);
}


void runAdaptiveStateCountingTests()
{
    // Input faults
    OUTF_STA();
    OUTF_INP();
    OUTF_OUTP();

    // Transition faults
    TRANSF_STA();
    TRANSF_INP();
    TRANSF_OUTP();

    // Equivalence
    AEQ_STA();
    AEQ_INP();
    AEQ_OUTP();

    // Reductions
    RED_STA();
    RED_INP();
    RED_OUTP();
}

int main(int argc, char** argv)
{

    START_EASYLOGGINGPP(argc, argv);
    nowText = initialize();
    logging::initLogging(nowText);


#ifdef ENABLE_DEBUG_MACRO
    CLOG(INFO, logging::globalLogger) << "This is a debug build!";
#endif

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

    /** Uncomment to run Adaptive State Counting tests **/
    // runAdaptiveStateCountingTests();

    exit(0);

}

