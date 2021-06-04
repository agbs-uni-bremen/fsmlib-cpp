/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <unordered_map>

#include "interface/FsmPresentationLayer.h"
#include "fsm/Dfsm.h"
#include "fsm/PkTable.h"
#include "fsm/FsmNode.h"
#include "fsm/IOTrace.h"
#include "fsm/SegmentedTrace.h"
#include "fsm/StrongReductionTestSuiteGenerator.h"

#include "trees/IOListContainer.h"
#include "trees/InputTree.h"
#include "trees/OutputTree.h"
#include "trees/TestSuite.h"
#include "trees/TreeNode.h"

#include "fsm/generalized-h-method.hpp"

#include "json/json.h"

#define DBG 0
using namespace std;
using namespace Json;


/**
 *   Program execution parameters and associated types
 */
typedef enum {
    FSM_CSV,
    FSM_JSON,
    FSM_BASIC
} model_type_t;

typedef enum {
    WMETHOD,
    WPMETHOD,
    SAFE_WMETHOD,
    SAFE_WPMETHOD,
    SAFE_HMETHOD,
    HMETHOD,
    HSIMETHOD,
    STRONG_REDUCTION_METHOD
} generation_method_t;


/** File containing the reference model */
static model_type_t modelType;
static string modelFile;

/** only for generation method SAFE_WPMETHOD */
static model_type_t modelAbstractionType;
static string modelAbstractionFile;
static string plStateFile;
static string plInputFile;
static string plOutputFile;
static string fsmName;
static string testSuiteFileName;
static string tcFilePrefix;
static generation_method_t genMethod;
static unsigned int numAddStates;

static shared_ptr<FsmPresentationLayer> pl = nullptr;
static shared_ptr<Dfsm> dfsm = nullptr;
static shared_ptr<Dfsm> dfsmAbstraction = nullptr;
static shared_ptr<Fsm> fsm = nullptr;
static shared_ptr<Fsm> fsmAbstraction = nullptr;

static bool isDeterministic = false;
static bool rttMbtStyle = false;


/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    cerr << "usage: " << name
    << " [-w|-wp|-h|-hsi|-sr] [-s] [-n fsmname] [-p infile outfile statefile] "
    << "[-a additionalstates] [-t testsuitename] [-rtt <prefix>] modelfile "
    << "[model abstraction file]" << endl;
}

/**
 *  Determine the model type of a model specified in an *.fsm file.
 *
 *  @param modelFile Name of the file containing the model
 *
 *  @return FSM_BASIC, if the model file has extension .fsm and
 *          contains the low-level encoding.
 *
 *  @return FSM_JSON, if the model file has extension .fsm and
 *                    contains the JSON encoding.
 *
 *  @return FSM_CSV, if the model file has extension .csv
 *
 */
static model_type_t getModelType(const string& mf) {
    
    if ( mf.find(".csv") != string::npos ) {
        return FSM_CSV;
    }
    
    model_type_t t = FSM_BASIC;
    
    ifstream inputFile(mf);
    string line;
    getline(inputFile,line);
    // Basic encoding does not contain any { or [
    if ( line.find("{") != string::npos or
        line.find("[") != string::npos) {
        t = FSM_JSON;
    }
    
    inputFile.close();
    
    return t;
    
}

/**
 * Parse parameters, stop execution if parameters are illegal.
 *
 * @param argc parameter 1 from main() invocation
 * @param argv parameter 2 from main() invocation
 */
static void parseParameters(int argc, char* argv[]) {
    
    // Set default parameters
    genMethod = WPMETHOD;
    fsmName = string("FSM");
    testSuiteFileName = string("testsuite.txt");
    numAddStates = 0;
    
    bool haveModelFileName = false;
    
    for ( int p = 1; p < argc; p++ ) {
        
        if ( strcmp(argv[p],"-w") == 0 ) {
            switch (genMethod) {
                case WPMETHOD: genMethod = WMETHOD;
                    break;
                case SAFE_WPMETHOD: genMethod = SAFE_WMETHOD;
                    break;
                default:
                    break;
            }
        }
        else if ( strcmp(argv[p],"-wp") == 0 ) {
            if ( genMethod == SAFE_WMETHOD or
                genMethod == SAFE_WPMETHOD ) {
                genMethod= SAFE_WPMETHOD;
            }
            else {
                genMethod = WPMETHOD;
            }
        }
        else if ( strcmp(argv[p],"-h") == 0 ) {
            if ( genMethod == SAFE_WMETHOD or
                genMethod == SAFE_WPMETHOD or
                genMethod == SAFE_HMETHOD ) {
                genMethod= SAFE_HMETHOD;
            }
            else {
                genMethod = HMETHOD;
            }
        }
        else if ( strcmp(argv[p],"-hsi") == 0 ) {
            genMethod = HSIMETHOD;
        }
        else if ( strcmp(argv[p],"-sr") == 0 ) {
            genMethod = STRONG_REDUCTION_METHOD;
        }
        else if ( strcmp(argv[p],"-s") == 0 ) {
            switch (genMethod) {
                case WPMETHOD: genMethod = SAFE_WPMETHOD;
                    break;
                case WMETHOD: genMethod = SAFE_WMETHOD;
                    break;
                case HMETHOD: genMethod = SAFE_HMETHOD;
                    break;
                default:
                    break;
            }
        }
        else if ( strcmp(argv[p],"-n") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing FSM name" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                fsmName = string(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-t") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing test suite name" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                testSuiteFileName = string(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-a") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing number of additional states" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                numAddStates = atoi(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-rtt") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing prefix for RTT-MBT test suite files" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                rttMbtStyle = true;
                tcFilePrefix = string(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-p") == 0 ) {
            if ( argc < p+4 ) {
                cerr << argv[0] << ": missing presentation layer files" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                plInputFile = string(argv[++p]);
                plOutputFile = string(argv[++p]);
                plStateFile = string(argv[++p]);
            }
        }
        else if ( strstr(argv[p],".csv")  ) {
            haveModelFileName = true;
            modelFile = string(argv[p]);
            modelType = getModelType(modelFile);
        }
        else if ( strstr(argv[p],".fsm")  ) {
            haveModelFileName = true;
            modelFile = string(argv[p]);
            modelType = getModelType(modelFile);
        }
        else {
            cerr << argv[0] << ": illegal parameter `" << argv[p] << "'" << endl;
            printUsage(argv[0]);
            exit(1);
        }
        
        if ( haveModelFileName and
            (genMethod == SAFE_WPMETHOD or
             genMethod == SAFE_WMETHOD or
             genMethod == SAFE_HMETHOD) ) {
                p++;
                if ( p >= argc ) {
                    cerr << argv[0] << ": missing model abstraction file" << endl;
                    printUsage(argv[0]);
                    exit(1);
                }
                modelAbstractionFile = string(argv[p]);
                modelAbstractionType = getModelType(modelAbstractionFile); 
            }
        
    }
    
    if ( modelFile.empty() ) {
        cerr << argv[0] << ": missing model file" << endl;
        printUsage(argv[0]);
        exit(1);
    }
    
}


/**
 *   Instantiate DFSM or FSM from input file according to
 *   the different input formats which are supported.
 */
static void readModel(model_type_t mtp,
                      string thisFileName,
                      string thisFsmName,
                      shared_ptr<Fsm>& myFsm,
                      shared_ptr<Dfsm>& myDfsm) {
    
    myDfsm = nullptr;
    myFsm = nullptr;
    
    
    switch ( mtp ) {
        case FSM_CSV:
            isDeterministic = true;
            myDfsm = make_shared<Dfsm>(thisFileName,thisFsmName);
            pl = myDfsm->getPresentationLayer();
            
            break;
            
        case FSM_JSON:
        {
            Reader jReader;
            Value root;
            stringstream document;
            ifstream inputFile(modelFile);
            document << inputFile.rdbuf();
            inputFile.close();
            
            if ( jReader.parse(document.str(),root) ) {
                myDfsm = make_shared<Dfsm>(root);
                pl = myDfsm->getPresentationLayer();
            }
            else {
                cerr << "Could not parse JSON model - exit." << endl;
                exit(1);
            }
        }
            break;
            
        case FSM_BASIC:
            if ( plStateFile.empty() ) {
                pl = make_shared<FsmPresentationLayer>();
            }
            else {
                pl = make_shared<FsmPresentationLayer>(plInputFile,plOutputFile,plStateFile);
            }
            myFsm = make_shared<Fsm>(modelFile,pl,fsmName);
            if ( myFsm->isDeterministic() ) {
                isDeterministic = true;
                myDfsm = make_shared<Dfsm>(modelFile,pl,fsmName);
                myFsm = nullptr;
            }
            break;
    }
    
    if ( myFsm != nullptr ) {
        myFsm->toDot(fsmName);
    }
    else if ( myDfsm != nullptr ) {
        myDfsm->toDot(fsmName);
        myDfsm->toCsv(fsmName);
    }
    
}


static void readModelAbstraction(model_type_t mtp,
                                 string thisFileName,
                                 string thisFsmName,
                                 shared_ptr<Dfsm>& myDfsm,
                                 shared_ptr<FsmPresentationLayer> plRef) {
    
    myDfsm = nullptr;
    
    
    switch ( mtp ) {
        case FSM_CSV:
            isDeterministic = true;
            myDfsm = make_shared<Dfsm>(thisFileName,thisFsmName,plRef);
            break;
            
        case FSM_JSON:
        {
            Reader jReader;
            Value root;
            stringstream document;
            ifstream inputFile(thisFileName);
            document << inputFile.rdbuf();
            inputFile.close();
            
            if ( jReader.parse(document.str(),root) ) {
                myDfsm = make_shared<Dfsm>(root,plRef);
            }
            else {
                cerr << "Could not parse JSON model - exit." << endl;
                exit(1);
            }
        }
            break;
            
        case FSM_BASIC:
            cerr << "ERROR. Model abstraction for SAFE W/WP/H METHOD may only be specified in CSV or JSON format - exit." << endl;
            exit(1);
            break;
    }
    
    if ( myDfsm != nullptr ) {
        myDfsm->toDot("ABS_"+fsmName);
        myDfsm->toCsv("ABS_"+fsmName);
    }

}

typedef vector<int> TCTrace;
typedef pair < TCTrace, TCTrace > TracePair;

struct TCTraceHash
{
    size_t operator() (const vector<int> & itrc) const
    {
        size_t seed = itrc.size();
        for(auto& i : itrc) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct TracePairHash
{
    size_t operator() (const TracePair & tracePair) const
    {
        auto trc1 = tracePair.first;
        auto trc2 = tracePair.second;

        size_t seed1 = trc1.size();
        for(auto& i : trc1) {
            seed1 ^= i + 0x9e3779b9 + (seed1 << 6) + (seed1 >> 2);
        }
        size_t seed2 = trc2.size();
        for(auto& i : trc1) {
            seed2 ^= i + 0x9e3779b9 + (seed2 << 6) + (seed2 >> 2);
        }
        return seed1 ^ seed2;
    }
};

typedef pair < TCTrace, TracePair > TC2TracePair;
typedef pair < TracePair, TCTrace > TracePair2Trace;
typedef unordered_map < TracePair, TCTrace, TracePairHash > Traces2GammaMap;
typedef unordered_multimap < TCTrace, TracePair, TCTraceHash > TC2TracesMap;

bool inPrefixRelation(vector<int> aPath, vector<int> bPath)
{
    if (aPath.empty() || bPath.empty())
        return false;
    for (unsigned i = 0; i<aPath.size() && i < bPath.size(); i++)
    {
        if (aPath[i] != bPath[i])
        {
            return false;
        }
    }
    return true;
}

shared_ptr<Tree> getPrefixRelationTreeWithoutTrace(const shared_ptr<Tree> & a, const shared_ptr<Tree> & b, const vector<int> & trc)
{
    IOListContainer aIOlst = a->getIOLists();
    IOListContainer bIOlst = b->getIOLists();

    shared_ptr<vector<vector<int>>> aPrefixes = aIOlst.getIOLists();
    shared_ptr<vector<vector<int>>> bPrefixes = bIOlst.getIOLists();

    shared_ptr<TreeNode> r = make_shared<TreeNode>();
    shared_ptr<Tree> tree = make_shared<Tree>(r, pl);

    if (aPrefixes->at(0).empty() && bPrefixes->at(0).empty())
    {
        return tree;
    }

    if (aPrefixes->at(0).empty())
    {
        return b;
    }
    if (bPrefixes->at(0).empty())
    {
        return a;
    }


    for (const auto &aPrefix : *aPrefixes)
    {
        for (const auto &bPrefix : *bPrefixes)
        {
            if (inPrefixRelation(aPrefix, bPrefix))
            {
                if (trc != aPrefix && trc != bPrefix)
                {
                    tree->addToRoot(aPrefix);
                    tree->addToRoot(bPrefix);
                }
            }
        }
    }
    return tree;

}

static int costMatrix[3][3];

static void initCostMatrix() {
    costMatrix[0][0] = 0;
    costMatrix[0][1] = 1;
    costMatrix[0][2] = 3;
    costMatrix[1][0] = 1;
    costMatrix[1][1] = 2;
    costMatrix[1][2] = 4;
    costMatrix[2][0] = 3;
    costMatrix[2][1] = 4;
    costMatrix[2][2] = 5;
}

static int insertionCosts(int trc1Costs, int trc2Costs) {
    return costMatrix[trc1Costs][trc2Costs];
}

static void addSHTraces(deque<pair<shared_ptr<SegmentedTrace>,shared_ptr<SegmentedTrace>>>  X,
                        Dfsm& refDfsm,
                        Dfsm& distDfsm,
                        Tree &testSuiteTree,
                        vector<int>* dfsmMinNodes2dfsmNodes = nullptr) {
    
    
    // Add all traces from X-pairs, extended by distinguishing traces
    for ( const auto p : X ) {
        const shared_ptr<SegmentedTrace> tr1 = p.first;
        const shared_ptr<SegmentedTrace> tr2 = p.second;
        shared_ptr<FsmNode> s1 = tr1->getTgtNode();
        shared_ptr<FsmNode> s2 = tr2->getTgtNode();
        
#if DBG
        cout << "============================================ " << endl;
        cout << "CHECK: " << *tr1 << "  NODE " << tr1->getTgtNode()->getId() << endl;
        cout << "CHECK: " << *tr2 << "  NODE " << tr2->getTgtNode()->getId() << endl;
#endif
        
        // Only handle trace pairs leading to
        // distinguishable target nodes.
        // If dfsmMinNodes2dfsmNodes is not null,
        // we have to translate the node numbers from
        // those of the minimised reference machine to
        // the original numbers of the unminimised machine.
        // The latter are identical to the node numbers of
        // the unminimised abstraction used here as distDfsm.
        shared_ptr<FsmNode> d1;
        shared_ptr<FsmNode> d2;
        if ( dfsmMinNodes2dfsmNodes == nullptr) {
            d1 = s1;
            d2 = s2;
        }
        else {
            d1 = distDfsm.getNodes().at(dfsmMinNodes2dfsmNodes->at(s1->getId()));
            d2 = distDfsm.getNodes().at(dfsmMinNodes2dfsmNodes->at(s2->getId()));
#if DBG
            cout << "ORIGINAL NODES " << d1->getId() << " " << d2->getId() << endl;
#endif
        }
        if ( not distDfsm.distinguishable(*d1, *d2) ) {
#if DBG
            cout << "Indistinguishable." << endl;
#endif
            continue;
        }
        
        // Now process this trace pair, extended by
        // a distinguishing trace
        
        // Which are the candidates to distinguish s1 and s2 ?
        vector< shared_ptr< vector<int> > > v01 = refDfsm.getDistTraces(*s1,*s2);
        shared_ptr<SegmentedTrace> tr1Ext = make_shared<SegmentedTrace>(*tr1);
        shared_ptr<SegmentedTrace> tr2Ext = make_shared<SegmentedTrace>(*tr2);
        
        shared_ptr< vector<int> > vBest = v01[0];
        auto seg = make_shared<TraceSegment>(vBest);
        tr1Ext->add(seg);
        tr2Ext->add(seg);
        
        int bestEffect1 = testSuiteTree.tentativeAddToRoot(*tr1Ext);
        int bestEffect2 = testSuiteTree.tentativeAddToRoot(*tr2Ext);
        
        for ( size_t i = 1; i < v01.size() and bestEffect1 + bestEffect2 > 0; i++ ) {
            shared_ptr< vector<int> > vAux = v01[i];
            seg = make_shared<TraceSegment>(vAux);
            shared_ptr<SegmentedTrace> tr1Aux = make_shared<SegmentedTrace>(*tr1);
            shared_ptr<SegmentedTrace> tr2Aux = make_shared<SegmentedTrace>(*tr2);
            tr1Aux->add(seg);
            tr2Aux->add(seg);
            
            int effAux1 = testSuiteTree.tentativeAddToRoot(*tr1Aux);
            int effAux2 = testSuiteTree.tentativeAddToRoot(*tr2Aux);
            
            if ( insertionCosts(effAux1, effAux2) < insertionCosts(bestEffect1, bestEffect2) ) {
                vBest = vAux;
                tr1Ext = tr1Aux;
                tr2Ext = tr2Aux;
                bestEffect1 = effAux1;
                bestEffect2 = effAux2;
            }
        }
        
        vector<int> v1 = tr1Ext->getCopy();
        vector<int> v2 = tr2Ext->getCopy();
        
        if ( bestEffect1 > 0 ) {
#if DBG
            cout << "ADD: " << *tr1Ext << endl;
#endif
            testSuiteTree.addToRoot(v1);
        }
        if ( bestEffect2 > 0 ) {
#if DBG
            cout << "ADD: " << *tr2Ext << endl;
#endif
            testSuiteTree.addToRoot(v2);
        }
        
    }
    
#if DBG
    IOListContainer testCasesSH = testSuiteTree.getIOLists();
    cout << "TEST SUITE STATUS:" << endl<< refDfsm.createTestSuite(testCasesSH) << endl;
#endif
}

#if 0

static void safeHMethod(const shared_ptr<TestSuite> &testSuite) {
    
    // Minimise original reference DFSM
    Dfsm dfsmRefMin = dfsm->minimise();

    // Minimise abstracted Dfsm
    Dfsm dfsmAbstractionMin = dfsmAbstraction->minimise();
    
    string fsmNameMinimal(fsmName + "_MINIMAL");
    string absFsmNameMinimal("ABS_" + fsmName + "_MINIMAL");
    
    dfsmRefMin.toDot(fsmNameMinimal);
    dfsmAbstractionMin.toDot(absFsmNameMinimal);
    dfsmAbstractionMin.toCsv(absFsmNameMinimal);
    
    cout << "REF    size = " << dfsm->size() << endl;
    cout << "REFMIN size = " << dfsmRefMin.size() << endl;
    cout << "ABSMIN size = " << dfsmAbstractionMin.size() << endl;
    
    shared_ptr<FsmNode> s0 = dfsmRefMin.getInitialState();
    shared_ptr<FsmPresentationLayer> pl = dfsmRefMin.getPresentationLayer();

    shared_ptr<Tree> iTreeH = dfsmRefMin.getStateCover();
    shared_ptr<Tree> iTreeSH = dfsmRefMin.getStateCover();

    // Let A be a set consisting of alpha.w, beta.w for any alpha != beta in V
    // and a distinguishing trace w. q0-after-alpha !~ q0-after-beta
    IOListContainer iolcV = iTreeH->getIOListsWithPrefixes();
    shared_ptr<vector<vector<int>>> iolV = iolcV.getIOLists();
    
    // Let B = V.(union_(i=1)^(m-n+1) Sigma_I)
    shared_ptr<Tree> B = dfsmRefMin.getStateCover();
    IOListContainer inputEnum = IOListContainer(dfsmRefMin.getMaxInput(),
                                                1,
                                                numAddStates + 1,
                                                pl);
    B->add(inputEnum);
    iTreeH->unionTree(B);
    iTreeSH->unionTree(B);

    Traces2GammaMap tracePair2gamma;
    TC2TracesMap tc2traces;

    vector<TracePair> tracesToCompare;

    // A
    for (unsigned i = 0; i < iolV->size(); i++)
    {
        shared_ptr<InputTrace> alpha = make_shared<InputTrace>(iolV->at(i), pl);
        for (unsigned j = i + 1; j < iolV->size(); j++)
        {
            shared_ptr<InputTrace> beta = make_shared<InputTrace>(iolV->at(j), pl);

            shared_ptr<Tree> alphaTree = iTreeH->getSubTree(alpha);
            shared_ptr<Tree> betaTree = iTreeH->getSubTree(beta);

            shared_ptr<Tree> prefixRelationTree = alphaTree->getPrefixRelationTree(betaTree);

            InputTrace gamma = dfsmRefMin.calcDistinguishingTrace(alpha, beta, prefixRelationTree);

            shared_ptr<InputTrace> iAlphaGamma = make_shared<InputTrace>(alpha->get(), pl);
            iAlphaGamma->append(gamma.get());

            shared_ptr<InputTrace> iBetaGamma = make_shared<InputTrace>(beta->get(), pl);
            iBetaGamma->append(gamma.get());

            iTreeH->addToRoot(iAlphaGamma->get());
            iTreeH->addToRoot(iBetaGamma->get());

            iTreeSH->addToRoot(iAlphaGamma->get());
            iTreeSH->addToRoot(iBetaGamma->get());

            TracePair2Trace tracePair2Gamma1(make_pair(alpha->get(), beta->get()), gamma.get());
            tracePair2gamma.insert(tracePair2Gamma1);
        }
    }


    // B
    shared_ptr<vector<vector<int>>> iolB = inputEnum.getIOLists();
    shared_ptr<FsmNode> abs_s0 = dfsmAbstractionMin.getInitialState();

    for (const auto &beta : *iolB)
    {
        for (auto alpha : *iolV)
        {
            shared_ptr<InputTrace> iAlphaBeta = make_shared<InputTrace>(alpha, pl);
            iAlphaBeta->append(beta);

            for (auto o : *iolV)
            {
                shared_ptr<InputTrace> iOmega = make_shared<InputTrace>(o, pl);

                shared_ptr<FsmNode> s0AfterAlphaBeta = *s0->after(*iAlphaBeta).begin();
                shared_ptr<FsmNode> s0AfterOmega = *s0->after(*iOmega).begin();
                if (s0AfterAlphaBeta == s0AfterOmega) continue;
                tracesToCompare.emplace_back(iAlphaBeta->get(),iOmega->get());
            }
        }
    }

    // C
    for (const auto &v : *iolV)
    {
        for (const auto &g1 : *iolB)
        {
            shared_ptr<InputTrace> iVG1 = make_shared<InputTrace>(v, pl);
            iVG1->append(g1);

            for (unsigned i = 1; i < g1.size(); ++i)
            {
                vector<int> g2(begin(g1), begin(g1) + i);
                shared_ptr<InputTrace> iVG2 = make_shared<InputTrace>(v, pl);
                iVG2->append(g2);

                shared_ptr<FsmNode> s0AfterVG1 = *s0->after(*iVG1).begin();
                shared_ptr<FsmNode> s0AfterVG2 = *s0->after(*iVG2).begin();
                if (s0AfterVG1 == s0AfterVG2) continue;
                tracesToCompare.emplace_back(iVG1->get(),iVG2->get());
            }
        }
    }

    /* h method */
    for (TracePair tracePair : tracesToCompare)  {
        shared_ptr<InputTrace> alpha = make_shared<InputTrace>(tracePair.first, pl);
        shared_ptr<InputTrace> beta = make_shared<InputTrace>(tracePair.second, pl);

        shared_ptr<Tree> alphaTree = iTreeH->getSubTree(alpha);
        shared_ptr<Tree> betaTree = iTreeH->getSubTree(beta);
        shared_ptr<Tree> prefixRelationTree = alphaTree->getPrefixRelationTree(betaTree);
        InputTrace gamma = dfsmRefMin.calcDistinguishingTrace(alpha, beta, prefixRelationTree);

        shared_ptr<InputTrace> iAlphaGamma = make_shared<InputTrace>(alpha->get(), pl);
        iAlphaGamma->append(gamma.get());
        shared_ptr<InputTrace> iBetaGamma = make_shared<InputTrace>(beta->get(), pl);
        iBetaGamma->append(gamma.get());

        iTreeH->addToRoot(iAlphaGamma->get());
        iTreeH->addToRoot(iBetaGamma->get());

        shared_ptr<FsmNode> afterAlpha = *abs_s0->after(*alpha).begin();
        shared_ptr<FsmNode> afterBeta = *abs_s0->after(*beta).begin();
        if (afterAlpha == afterBeta) continue;

        TC2TracePair tc2trace1(iAlphaGamma->get(), tracePair);
        TC2TracePair tc2trace2(iBetaGamma->get(), tracePair);
        tc2traces.insert(tc2trace2);
        tc2traces.insert(tc2trace1);

        TracePair2Trace tracePair2Gamma1(tracePair, gamma.get());
        tracePair2gamma.insert(tracePair2Gamma1);
    }

    /* safety h */
    for (TracePair tracePair : tracesToCompare)  {
        shared_ptr<InputTrace> alpha = make_shared<InputTrace>(tracePair.first, pl);
        shared_ptr<InputTrace> beta = make_shared<InputTrace>(tracePair.second, pl);

        shared_ptr<FsmNode> afterAlpha = *abs_s0->after(*alpha).begin();
        shared_ptr<FsmNode> afterBeta = *abs_s0->after(*beta).begin();
        if (afterAlpha == afterBeta) continue;

        // Each TracePair (alpha, beta) generates two test cases in h method.
        // alpha.gamma and beta.gamma
        // Add those two test cases to safety-H-Tree
        auto gamma = tracePair2gamma[tracePair];
        auto iAlphaGamma = alpha->get();
        iAlphaGamma.insert(iAlphaGamma.end(), gamma.begin(), gamma.end());
        auto iBetaGamma = beta->get();
        iBetaGamma.insert(iBetaGamma.end(), gamma.begin(), gamma.end());

        iTreeSH->addToRoot(iAlphaGamma);
        iTreeSH->addToRoot(iBetaGamma);
    }

    /**
     * This loop finds and removes redundant test cases
     *
       - iterate over all safety-h test cases
         - iterate over all trace pairs (a, b) needing this test case tc
           - note: tc1 = a.g . There is also a partner test case tc2 = b.g
           - find a better gamma g' for (a,b)
         - if a better gamma can be found for all trace pairs:
           - remove this test case (a.g) and its partner (b.g)
           - add a.g' and b.g' to iTreeSH

    *  By adding a.g' and b.g' some existend test cases can be lengthened
    *  but no new test cases are added
    *
    *  NOTE: Going through this loop multiple times improves the result
    *        One could loop until there is no change.
    */

    /* test cases that should not be changed */
    vector< vector<int> > perfectTestCases;
    for (unsigned int i = 0; i < 3; i++) {
    vector< vector<int> > tests = *iTreeSH->getIOLists().getIOLists();
    for (vector<int> testCase : tests)
    {
        if (find(begin(perfectTestCases), end(perfectTestCases), testCase) != end(perfectTestCases))
        {
            continue;
        }

        bool betterGammaForAllPairs = true;
        Traces2GammaMap pair2NewGamma;

        auto range = tc2traces.equal_range(testCase);
        for (auto it = range.first; it != range.second; ++it )
        {
            TracePair tracePair = it->second;
            vector<int> gamma = tracePair2gamma[tracePair];

            auto alpha = make_shared<InputTrace>(tracePair.first, pl);
            auto beta =  make_shared<InputTrace>(tracePair.second, pl);

            shared_ptr<Tree> alphaTree = iTreeSH->getSubTree(make_shared<InputTrace>(alpha->get(),pl));
            shared_ptr<Tree> betaTree = iTreeSH->getSubTree(make_shared<InputTrace>(beta->get(),pl));
            shared_ptr<Tree> prefixRelationTree = getPrefixRelationTreeWithoutTrace(alphaTree, betaTree, gamma);

            if (prefixRelationTree->size() == 1)
            {
                // no chance to find a better gamma
                betterGammaForAllPairs = false;
                break;
            }

            InputTrace newGamma = dfsmRefMin.calcDistinguishingTraceInTree(alpha, beta, prefixRelationTree);
            if (newGamma.size() > 0)
            {
                // store newGamma for this pair
                TracePair2Trace p2g(tracePair, newGamma.get());
                pair2NewGamma.insert(p2g);
            } else {
                betterGammaForAllPairs = false;
                break;
            }
        }
        // if for this test case all trace pairs can have better gammas
        if (betterGammaForAllPairs)
        {

            //            shared_ptr<InputTrace> iTestCase = make_shared<InputTrace>(testCase, pl);
            //            shared_ptr<TreeNode> afterTC = iTreeSH->getRoot()->after(iTestCase->cbegin(), iTestCase->cend());
            //            if (afterTC && afterTC->isLeaf())
            //                afterTC->deleteSingleNode();

            // append new gammas to all TracePairs (a,b)
            auto range = tc2traces.equal_range(testCase);
            vector<TC2TracePair> newTestCases;
            for (auto it = range.first; it != range.second; ++it)
            {
                TracePair tracePair = it->second;

                // delete this test case and its partner if possible
                auto oldGamma = tracePair2gamma[tracePair];
                auto iAlphaGammaOld = make_shared<InputTrace>(tracePair.first, pl);
                iAlphaGammaOld->append(oldGamma);
                auto iBetaGammaOld = make_shared<InputTrace>(tracePair.first, pl);
                iBetaGammaOld->append(oldGamma);
                shared_ptr<TreeNode> afterAGOld = iTreeSH->getRoot()->after(iAlphaGammaOld->cbegin(), iAlphaGammaOld->cend());
                if (afterAGOld && afterAGOld->isLeaf())
                    afterAGOld->deleteNode();
                shared_ptr<TreeNode> afterBGOld = iTreeSH->getRoot()->after(iBetaGammaOld->cbegin(), iBetaGammaOld->cend());
                if (afterBGOld && afterBGOld->isLeaf())
                    afterBGOld->deleteNode();

                // insert new two new test cases
                // actual insertion is done after iteration
                auto newGamma = pair2NewGamma[tracePair];
                auto iAlphaGamma = make_shared<InputTrace>(tracePair.first, pl);
                iAlphaGamma->append(newGamma);
                auto iBetaGamma =  make_shared<InputTrace>(tracePair.second, pl);
                iBetaGamma->append(newGamma);

                newTestCases.emplace_back(iAlphaGamma->get(), tracePair);
                newTestCases.emplace_back(iBetaGamma->get(), tracePair);
                perfectTestCases.push_back(iAlphaGamma->get());
                perfectTestCases.push_back(iBetaGamma->get());

                tracePair2gamma.erase(tracePair);
                TracePair2Trace tracePair2NewGamma(tracePair, newGamma);
                tracePair2gamma.insert(tracePair2NewGamma);

                iTreeSH->addToRoot(iAlphaGamma->get());
                iTreeSH->addToRoot(iBetaGamma->get());

            }
            for (const TC2TracePair &tc2trace : newTestCases){
                tc2traces.insert(tc2trace);
            }

            tc2traces.erase(testCase);
        }
    }
    }

    IOListContainer testCasesSH = iTreeSH->getIOLists();
    *testSuite = dfsmRefMin.createTestSuite(testCasesSH);
}

#else


static void safeHMethod(const shared_ptr<TestSuite> &testSuite) {
    
    initCostMatrix();
    Dfsm dfsmRefMin = dfsm->minimise();
    dfsmRefMin.calculateDistMatrix();
    
    // Map from node numbers in dfsmRefMin to
    // representatives in dfsm
    vector<int> dfsmMinNodes2dfsmNodes;
    for ( size_t n = 0; n < dfsmRefMin.size(); n++ ) {
        dfsmMinNodes2dfsmNodes.push_back(0);
    }
    shared_ptr<PkTable> pDfsm = dfsm->getPktblLst().back();
    for ( size_t n = 0; n < dfsm->size(); n++ ) {
        dfsmMinNodes2dfsmNodes[pDfsm->getClass(n)] = n;
    }
    
#if DBG
    for ( int n = 0; n < dfsmRefMin.size(); n++ ) {
        for ( int m = n+1; m < dfsmRefMin.size(); m++ ) {
            vector< shared_ptr< vector<int> > > v01 =
            dfsmRefMin.getDistTraces(*dfsmRefMin.getNodes()[n],
                                     *dfsmRefMin.getNodes()[m]);
            
            cout << endl << "DIST TRC: " << n << " " << m << endl;
            for ( const auto trc : v01 ) {
                for ( const auto x : *trc ) {
                    cout << x << ", ";
                }
                cout << endl;
            }
        }
    }
#endif
    
    
    // Minimise abstracted Dfsm
    Dfsm dfsmAbstractionMin = dfsmAbstraction->minimise();
    
    string fsmNameMinimal(fsmName + "_MINIMAL");
    string absFsmNameMinimal("ABS_" + fsmName + "_MINIMAL");
    
    dfsmRefMin.toDot(fsmNameMinimal);
    dfsmAbstractionMin.toDot(absFsmNameMinimal);
    dfsmAbstractionMin.toCsv(absFsmNameMinimal);
    cout << "REF    size = " << dfsm->size() << endl;
    cout << "REFMIN size = " << dfsmRefMin.size() << endl;
    cout << "ABSMIN size = " << dfsmAbstractionMin.size() << endl;
    
    shared_ptr<FsmNode> s0 = dfsmRefMin.getInitialState();
    shared_ptr<FsmPresentationLayer> pl = dfsmRefMin.getPresentationLayer();
    
    // Create an empty test suite as a tree of input traces
    shared_ptr<Tree> testSuiteTree = make_shared<Tree>(make_shared<TreeNode>(),pl);
    
    shared_ptr<Tree> V = dfsmRefMin.getStateCover();
    
    // Complete state cover as IOListContainer
    IOListContainer Vcontainer = V->getIOListsWithPrefixes();
    
    // State cover as vector of vectors
    shared_ptr< vector< vector<int> > > Vvectors = Vcontainer.getIOLists();
    
    // Empty deque of pointers to segmented traces
    deque< shared_ptr<SegmentedTrace> > Vtraces;
    
    // Fill Vtraces with new instances of segmented traces,
    // each trace consisting of a single segment from V, together
    // with its target node.
    for ( const auto v : *Vvectors ) {
        shared_ptr<FsmNode> tgtNode = *(s0->after(v).begin());
        shared_ptr< vector<int> > vPtr =
            make_shared< vector<int> >(v.begin(),v.end());
        shared_ptr<TraceSegment> seg = make_shared<TraceSegment>(vPtr,
                                                                 string::npos,
                                                                 tgtNode);
        deque< shared_ptr<TraceSegment> > d;
        d.push_back(seg);
        Vtraces.push_back(make_shared<SegmentedTrace>(d));
    }
    
    // Fill deque A with all pairs of Vtraces leading to distinct
    // target nodes
    deque< pair< shared_ptr<SegmentedTrace>,shared_ptr<SegmentedTrace> > > A;
    for ( size_t i = 0; i < Vtraces.size(); i++ ) {
        auto st1 = Vtraces.at(i);
        auto tgtNode1 = st1->getTgtNode();
        for ( size_t j = i+1; j < Vtraces.size(); j++ ) {
            auto st2 = Vtraces.at(j);
            if ( st2->getTgtNode() != tgtNode1 ) {
                A.push_back(make_pair(st1, st2));
            }
        }
    }
    
    // Construct all Sigma_I traces of
    // length 1..(m-n+1) (input enumeration).
    // The IOListContainer containing the input enumerations
    // is transformed into a deque of trace segments
    IOListContainer inputEnum = IOListContainer(dfsmRefMin.getMaxInput(),
                                                1,
                                                numAddStates + 1,
                                                pl);
    shared_ptr< vector< vector<int> > > inputEnumVec = inputEnum.getIOLists();
    deque< shared_ptr<TraceSegment> > inputEnumDeq;
    for ( const auto v : *inputEnumVec ) {
        shared_ptr< vector<int> > vPtr =
            make_shared< vector<int> >(v.begin(),v.end());
        shared_ptr<TraceSegment> seg = make_shared<TraceSegment>(vPtr,
                                                                 string::npos);
        inputEnumDeq.push_back(seg);
    }
    
    // Create deque of all Vtraces extended by suffixes from inpuEnumDeque.
    // Each of these extended traces needs to be in the test suite. As a consequence,
    // all traces of the state cover are also contained in the test suite as prefixes.
    deque< shared_ptr<SegmentedTrace> > V_inputEnumTraces;
    for ( const auto v : Vtraces ) {
        for ( const auto seg : inputEnumDeq ) {
            // Calculate the target node reached via v.seg
            shared_ptr<FsmNode> tgtNode =
                *(v->getTgtNode()->after(*seg->get()).begin());
            shared_ptr<TraceSegment> s = make_shared<TraceSegment>(*seg);
            s->setTgtNode(tgtNode);
            // Append s to copy of v
            shared_ptr<SegmentedTrace> u = make_shared<SegmentedTrace>(*v);
            u->add(s);
            
            // Do not add u if it is already contained in Vtraces
            bool containedInVtraces = false;
            for ( const auto w : Vtraces ) {
                if ( *w == *u ) {
                    containedInVtraces = true;
                    break;
                }
            }
            
            if ( not containedInVtraces ) {
                V_inputEnumTraces.push_back(u);
                // Put u into the test suite
                testSuiteTree->addToRoot(u->getCopy());
            }
        }
    }
    
#if DBG
    IOListContainer xxx = testSuiteTree->getIOLists();
    cout << "TEST SUITE STATUS (initial):" << endl<< dfsmRefMin.createTestSuite(xxx) << endl;
#endif
    
    // Create deque B containing all pairs of elements of Vtraces
    // V_inputEnumTraces, such that their target nodes differ
    deque< pair< shared_ptr<SegmentedTrace>,shared_ptr<SegmentedTrace> > > B;
    for ( const auto v : Vtraces ) {
        for ( const auto u : V_inputEnumTraces ) {
            if ( v->getTgtNode() != u->getTgtNode() ) {
                B.push_back(make_pair(v,u));
            }
        }
    }
    
    // Create deque C containing all pairs <alpha.gamma',alpha.gamma>
    // of elements from V_inputEnumTraces, such that alpha is in Vtraces,
    // and gamma' is a non-empty prefix of gamma, and the target nodes
    // differ. Recall that the alpha part is always contained in the
    // first segment of a segmented trace, while the input enumeration
    // is in the second segment.
    deque< pair< shared_ptr<SegmentedTrace>,shared_ptr<SegmentedTrace> > > C;
    for ( const auto v : V_inputEnumTraces ) {
        shared_ptr<TraceSegment> seg = v->back();
        // Loop over all true gamma-prefixes of v = alpha.gamma
        for ( size_t prefix = seg->size() - 1; prefix > 0; prefix-- ) {
            shared_ptr<FsmNode> tgtNode1 = v->getTgtNode();

            shared_ptr<TraceSegment> seg1 = v->front();
            shared_ptr<TraceSegment> s1 = make_shared<TraceSegment>(*seg1);
            shared_ptr<TraceSegment> s2 = make_shared<TraceSegment>(*seg);
            s2->setPrefix(prefix);
            
            shared_ptr<FsmNode> tgtNode2 =
                *(s1->getTgtNode()->after(s2).begin());
            
            if ( tgtNode1 != tgtNode2 ) {
                deque< std::shared_ptr<TraceSegment> > segs;
                segs.push_back(s1);
                segs.push_back(s2);
                shared_ptr<SegmentedTrace> vPref = make_shared<SegmentedTrace>(segs);
                C.push_back(make_pair(vPref,v));
            }
        }
    }
    
    addSHTraces(A,dfsmRefMin,dfsmRefMin,*testSuiteTree);
    addSHTraces(B,dfsmRefMin,*dfsmAbstraction,*testSuiteTree,&dfsmMinNodes2dfsmNodes);
    addSHTraces(C,dfsmRefMin,*dfsmAbstraction,*testSuiteTree,&dfsmMinNodes2dfsmNodes);
    
    IOListContainer testCasesSH = testSuiteTree->getIOLists();
    *testSuite = dfsmRefMin.createTestSuite(testCasesSH);
    
}

#endif


static void safeWpMethod(const shared_ptr<TestSuite> &testSuite) {
    
    // Minimise original reference DFSM
    // Dfsm dfsmRefMin = dfsm->minimise();
    Fsm dfsmRefMin = dfsm->minimiseObservableFSM();
    
    dfsmRefMin.toDot("REFMIN");
    cout << "REF    size = " << dfsm->size() << endl;
    cout << "REFMIN size = " << dfsmRefMin.size() << endl;
    
    // Get state cover of original model
    shared_ptr<Tree> scov = dfsmRefMin.getStateCover();
    
    // Get transition cover of the original model
    shared_ptr<Tree> tcov = dfsmRefMin.getTransitionCover();
    
    // Get characterisation set of original model
    IOListContainer w = dfsmRefMin.getCharacterisationSet();
    
    cout << "W = " << w << endl;
    
    // Minimise the abstracted reference model
    Dfsm dfsmAbstractionMin = dfsmAbstraction->minimise();
    //Fsm dfsmAbstractionMin = dfsmAbstraction->minimiseObservableFSM();
    
    
    dfsmAbstractionMin.toDot("ABSMIN");
    cout << "ABSMIN size = " << dfsmAbstractionMin.size() << endl;
    
    
    // Get W_s, the characterisation set of dfsmAbstractionMin
    IOListContainer wSafe = dfsmAbstractionMin.getCharacterisationSet();
    
    cout << "wSafe = " << wSafe << endl;
    
    
    // Get W_sq, the state identification sets of dfsmAbstractionMin
    dfsmAbstractionMin.calcStateIdentificationSets();
    
    // Calc W1 = V.W, W from original model
    shared_ptr<Tree> W1 = dfsmRefMin.getStateCover();
    
    W1->add(w);
    
    // Calc W21 = V.wSafe
    shared_ptr<Tree> W2 = dfsmRefMin.getStateCover();
    W2->add(wSafe);
    
    // Calc W22 = V.(union_(i=1)^(m-n) Sigma_I).wSafe)
    shared_ptr<Tree> W22;
    if ( numAddStates > 0 ) {
        W22 = dfsmRefMin.getStateCover();
        IOListContainer inputEnum = IOListContainer(dfsm->getMaxInput(),
                                                    1,
                                                    numAddStates,
                                                    pl);
        W22->add(inputEnum);
        W22->add(wSafe);
        W2->unionTree(W22);
    }
    
    // Calc W3 = V.Sigma_I^(m - n + 1) oplus
    //           {Wis | Wis is state identification set of csmAbsMin}
    shared_ptr<Tree> W3 = dfsmRefMin.getStateCover();
    IOListContainer inputEnum2 = IOListContainer(dfsm->getMaxInput(),
                                                 (numAddStates+1),
                                                 (numAddStates+1),
                                                 pl);
    W3->add(inputEnum2);
    
    dfsmAbstractionMin.appendStateIdentificationSets(W3);
    
    // Union of all test cases: W1 union W2 union W3
    // Collected again in W1
    W1->unionTree(W2);
    W1->unionTree(W3);
    
    IOListContainer iolc = W1->getTestCases();
    *testSuite = dfsm->createTestSuite(iolc);
    
}

static void safeWMethod(const shared_ptr<TestSuite> &testSuite) {
    
    // Minimise original reference DFSM
    Dfsm dfsmRefMin = dfsm->minimise();
    
    cout << "REF    size = " << dfsm->size() << endl;
    cout << "REFMIN size = " << dfsmRefMin.size() << endl;
    
    // Get state cover of original model
    shared_ptr<Tree> scov = dfsmRefMin.getStateCover();
    
    // Get characterisation set of original model
    IOListContainer w = dfsmRefMin.getCharacterisationSet();
    
    cout << "W = " << w << endl;
    
    // Minimise the abstracted reference model
    Dfsm dfsmAbstractionMin = dfsmAbstraction->minimise();
    
    cout << "ABSMIN size = " << dfsmAbstractionMin.size() << endl;
    
    
    // Get W_s, the characterisation set of dfsmAbstractionMin
    IOListContainer wSafe = dfsmAbstractionMin.getCharacterisationSet();
    
    cout << "wSafe = " << wSafe << endl;
    
    // Calc W1 = V.W, W from original model
    shared_ptr<Tree> W1 = dfsmRefMin.getStateCover();
    
    W1->add(w);
    
    // Calc W21 = V.W_s
    shared_ptr<Tree> W21 = dfsmRefMin.getStateCover();
    W21->add(wSafe);
    
    // Calc W22 = V.(union_(i=1)^(m-n+1) Sigma_I).wSafe)
    shared_ptr<Tree> W22 = dfsmRefMin.getStateCover();
    
    IOListContainer inputEnum = IOListContainer(dfsm->getMaxInput(),
                                                1,
                                                numAddStates+1,
                                                pl);
    W22->add(inputEnum);
    
    W22->add(wSafe);
    
    // Union of all test cases: W1 union W2 union W3
    // Collected again in W1
    W1->unionTree(W21);
    W1->unionTree(W22);
    
    IOListContainer iolc = W1->getTestCases();
    *testSuite = dfsm->createTestSuite(iolc);
    
}


static void generateStrongReductionTestSuite() {

    StrongReductionTestSuiteGenerator gen(fsm,true);
    InputTree testSuite = gen.generateTestSuite(fsm->getNodes().size() + numAddStates);

    ofstream out(testSuiteFileName);
    out << testSuite;
    out.close();
    
    if ( rttMbtStyle ) {
        cout << "RTT-MBT style is not supported for strong reduction testing " << endl;        
    }    
    
    cout << "Number of test cases (input sequences): " << testSuite.getNumberOfSequences() << endl;
    cout << "Total length (inputs)                 : " << testSuite.getTotalLengthOfSequences() << endl;
}



static void generateTestSuite() {

    // test suites for strong reduction are not represented using type
    // TestSuite but instead are only represented as lists of input
    // sequences
    if ( genMethod == STRONG_REDUCTION_METHOD ) {
        generateStrongReductionTestSuite();
        return;
    }


    shared_ptr<TestSuite> testSuite =
    make_shared<TestSuite>();
    
    switch ( genMethod ) {
        case WMETHOD:
            if ( dfsm != nullptr ) {
                IOListContainer iolc = dfsm->wMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
            }
            else {
                IOListContainer iolc = fsm->wMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(fsm->apply(*itrc));
                }
            }
            break;
            
        case WPMETHOD:
            if ( dfsm != nullptr ) {
                IOListContainer iolc = dfsm->wpMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
            }
            else {
                IOListContainer iolc = fsm->wpMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(fsm->apply(*itrc));
                }
            }
            break;
            
        case HMETHOD:
            if ( dfsm != nullptr ) {
                Dfsm dfsmMin = dfsm->minimise();
                IOListContainer iolc =
                dfsmMin.hMethodOnMinimisedDfsm(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
                //if(isApplicable(dfsmMin)) {
                //    auto hMethodTestSuite = generateHMethodTestSuite(dfsmMin, numAddStates);
                //    for ( auto inVec : hMethodTestSuite ) {
                //        shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                //        testSuite->push_back(dfsm->apply(*itrc));
                //    }
                //} else {
                //    exit(-1);
                //    std::cout << "Invalid Dfsm for H-Method." << std::endl;
                //}
            } else {
                Fsm fsmMin = fsm->minimise();
                if(isApplicable(fsmMin)) {
                    auto hMethodTestSuite = generateHMethodTestSuite(fsmMin, numAddStates);
                    for ( auto inVec : hMethodTestSuite ) {
                        shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                        testSuite->push_back(fsm->apply(*itrc));
                    }
                } else {
                    std::cout << "Invalid Fsm for H-Method." << std::endl;
                }
            }
            break;
            
        case HSIMETHOD:
            if ( dfsm != nullptr ) {
                IOListContainer iolc = dfsm->hsiMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(dfsm->apply(*itrc));
                }
            }
            else {
                IOListContainer iolc = fsm->hsiMethod(numAddStates);
                for ( auto inVec : *iolc.getIOLists() ) {
                    shared_ptr<InputTrace> itrc = make_shared<InputTrace>(inVec,pl);
                    testSuite->push_back(fsm->apply(*itrc));
                }
            }
            break;
            
        case SAFE_HMETHOD:
            safeHMethod(testSuite);
            break;
        case SAFE_WPMETHOD:
            safeWpMethod(testSuite);
            break;
            
        case SAFE_WMETHOD:
            safeWMethod(testSuite);
            break;

        default: 
            cout << "unsupported test method" << endl;
            return;
    }
    
    testSuite->save(testSuiteFileName);
    
    if ( rttMbtStyle ) {
        int numTc = 0;
        for ( size_t tIdx = 0; tIdx < testSuite->size(); tIdx++ ) {
            
            OutputTree ot = testSuite->at(tIdx);
            vector<IOTrace> iotrcVec;
            ot.toIOTrace(iotrcVec);
            
            for ( size_t iIdx = 0; iIdx < iotrcVec.size(); iIdx++ ) {
                ostringstream tcFileName;
                tcFileName << tcFilePrefix << tIdx << "_" << iIdx << ".log";
                ofstream outFile(tcFileName.str());
                outFile << iotrcVec[iIdx].toRttString();
                outFile.close();
                numTc++;
            }
            
        }
        
    }
    
    cout << "Number of test cases: " << testSuite->size() << endl;
    cout << "        total length: " << testSuite->totalLength() << endl;
    
}

int main(int argc, char* argv[])
{
    
    parseParameters(argc,argv);
    readModel(modelType,modelFile,fsmName,fsm,dfsm);
    
    if ( genMethod == SAFE_WPMETHOD or
        genMethod == SAFE_WMETHOD or
        genMethod == SAFE_HMETHOD) {
        if ( dfsm == nullptr ) {
            cerr << "SAFE W/WP METHOD only operates on deterministic FSMs - exit."
            << endl;
            exit(1);
        }
        
        shared_ptr<FsmPresentationLayer> plRef = dfsm->getPresentationLayer();
        
        readModelAbstraction(modelAbstractionType,
                             modelAbstractionFile,
                             "ABS_"+fsmName,
                             dfsmAbstraction,
                             plRef);
    }
    
    generateTestSuite();
    
    exit(0);
    
}



