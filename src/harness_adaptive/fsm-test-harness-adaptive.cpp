#include "sut_wrapper_adaptive.h"
#include "vPrimeEnumerator.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <utility>


#include "../interface/FsmPresentationLayer.h"
#include "../fsm/Dfsm.h"
#include "../fsm/PkTable.h"
#include "../fsm/FsmNode.h"
#include "../fsm/IOTrace.h"
#include "../fsm/SegmentedTrace.h"
#include "../fsm/IOTraceContainer.h"

#include "../trees/AdaptiveTreeNode.h"
#include "../trees/IOListContainer.h"
#include "../trees/IOTreeContainer.h"
#include "../trees/OutputTree.h"
#include "../trees/TestSuite.h"

#include "../utils/Logger.hpp"

#define DBG 0
using namespace std;
using namespace Json;


/* Required inputs:
 *  - SUT-handling (init, reset, apply)
 *     --> via sut_wrapper
 *  - k (for applying tests k-times)
 *     --> via additional input
 *  - Spec
 *     --> via additional input
 *  - m
 *     --> via additional input
 */


/**
 *   Program execution parameters and associated types
 */
typedef enum {
    FSM_CSV,
    FSM_JSON,
    FSM_BASIC
} model_type_t;


/** File containing the reference model */
static model_type_t modelType;
static string modelFile; 

static string plStateFile;
static string plInputFile;
static string plOutputFile;

static shared_ptr<FsmPresentationLayer> pl = nullptr;
static shared_ptr<Fsm> fsm = nullptr;


static int k; // number of times to repeat application of test
static int m; // upper bound on the number of states of the FSM-representation of the SUT
static int numAddStates; // m - (size of reference model)




/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    cerr << "usage: " << name
    << "[-a additionalstates] [-k input application repetitions] modelfile " << endl;
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
    
    
    
    // set default values
    k = 10;
    numAddStates = 0;
    
    
    for ( int p = 1; p < argc; p++ ) {
        
        if ( strcmp(argv[p],"-a") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing number of additional states" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                numAddStates = atoi(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-k") == 0 ) {
            if ( argc < p+2 ) {
                cerr << argv[0] << ": missing number of input application repetitions" << endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                k = atoi(argv[++p]);
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
            modelFile = string(argv[p]);
            modelType = getModelType(modelFile);
        }
        else if ( strstr(argv[p],".fsm")  ) {
            modelFile = string(argv[p]);
            modelType = getModelType(modelFile);
        }
        else {
            cerr << argv[0] << ": illegal parameter `" << argv[p] << "'" << endl;
            printUsage(argv[0]);
            exit(1);
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
                      shared_ptr<Fsm>& myFsm) {
    
    myFsm = nullptr;
    
    
    switch ( mtp ) {
        case FSM_CSV:
            myFsm = make_shared<Dfsm>(thisFileName,"ReferenceModel");

            pl = myFsm->getPresentationLayer();
            
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
                myFsm = make_shared<Dfsm>(root);
                pl = myFsm->getPresentationLayer();
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
            myFsm = make_shared<Fsm>(modelFile,pl,"ReferenceModel");
            if ( myFsm->isDeterministic() ) {
                myFsm = make_shared<Dfsm>(modelFile,pl,"ReferenceModel");
                myFsm = nullptr;
            }
            break;
    }
    
    /*
    if ( myFsm != nullptr ) {
        myFsm->toDot(fsmName);
    }
    else if ( myFsm != nullptr ) {
        myFsm->toDot(fsmName);
        myFsm->toCsv(fsmName);
    }
    */


    // set m
    m = myFsm->size() + numAddStates;
    
}




bool exceedsBound(  const size_t m,
                    const IOTrace& base,
                    const IOTrace& suffix,
                    const vector<shared_ptr<FsmNode>>& states,
                    const IOTreeContainer& adaptiveTestCases,
                    unordered_set<IOTraceContainer> bOmegaT,
                    const IOTraceContainer& vDoublePrime,
                    const vector<shared_ptr<FsmNode>>& dReachableStates,
                    const Fsm& spec)
{
    size_t lB = Fsm::lowerBound(base, suffix, states, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
    LOG("VERBOSE_1") << "lB: " << lB << std::endl;
    return lB > m;
}

size_t lowerBound(const IOTrace& base,
                       const IOTrace& suffix,
                       const vector<shared_ptr<FsmNode>>& states,
                       //const IOTreeContainer& adaptiveTestCases,
                       //unordered_set<IOTraceContainer> bOmegaT,
                       unordered_set<shared_ptr<IOTraceContainer>> responseSets,
                       unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseMaps,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec)
{
    LOG("VERBOSE_1") << "lowerBound()" << std::endl;
    LOG("VERBOSE_1") << "base: " << base << std::endl;
    LOG("VERBOSE_1") << "suffix: " << suffix << std::endl;
    LOG("VERBOSE_1") << "states:" << std::endl;
    for (auto s : states)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    //LOG("VERBOSE_1") << "adaptiveTestCases: " << adaptiveTestCases << std::endl;
    LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;
    LOG("VERBOSE_1") << "dReachableStates: " << std::endl;
    for (auto s : dReachableStates)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    size_t result = 0;
    LOG("VERBOSE_1") << "lb result: " << result << std::endl;

    LOG("VERBOSE_1") << "bOmegaT:" << std::endl;
    for (const auto& cont : bOmegaT)
    {
        LOG("VERBOSE_1") << "  " << cont << std::endl;
    }

    for (shared_ptr<FsmNode> state : states)
    {
        const IOTraceContainer& rResult = spec.r(state, base, suffix);
        LOG("VERBOSE_1") << "--- state: " << state->getName() << std::endl;
        LOG("VERBOSE_1") << "rResult(" << state->getName() << ", " << base << ", " << suffix << "): " << rResult << std::endl;
        result += rResult.size();
        LOG("VERBOSE_1") << "lb result: " << result << std::endl;
        if(find(dReachableStates.begin(), dReachableStates.end(), state) != dReachableStates.end()) {
            ++result;
            LOG("VERBOSE_1") << "State " << state->getName() << " is d-reachable. Incrementing." << std::endl;
            LOG("VERBOSE_1") << "lb result: " << result << std::endl;
        }

        IOTraceContainer rPlusResult = spec.rPlus(state, base, suffix, vDoublePrime, true);
        rPlusResult.add(rResult);
        LOG("VERBOSE_1") << "rPlusResult: " << rPlusResult << std::endl;
        for (auto traceIt = rPlusResult.cbegin(); traceIt != rPlusResult.cend(); ++traceIt)
        {
            const shared_ptr<const IOTrace>& trace = *traceIt;
            IOTraceContainer traces = iut.bOmega(adaptiveTestCases, *trace);
            LOG("VERBOSE_1") << "Removing " << traces << " from testTraces." << std::endl;

            IOTraceContainer::remove(bOmegaT, traces);

            LOG("VERBOSE_1") << "testTraces:" << std::endl;
            for (const auto& cont : bOmegaT)
            {
                LOG("VERBOSE_1") << "  " << cont << std::endl;
            }
        }
    }
    LOG("VERBOSE_1") << "bOmegaT size: " << bOmegaT.size() << std::endl;
    LOG("VERBOSE_1") << "bOmegaT:" << std::endl;
    for (const auto& cont : bOmegaT)
    {
        LOG("VERBOSE_1") << "  " << cont << std::endl;
    }
    result += bOmegaT.size();
    LOG("VERBOSE_1") << "lowerBound() result: " << result << std::endl;
    return result;
}



// static void performAdaptiveTesting(Fsm& spec, IOTraceContainer& observedTraces, const IOListContainer& rCharacterisationSet, const IOTreeContainer& adaptiveTestCases, const IOListContainer& adaptiveList, const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates, InputTraceSet detStateCover, const vector<shared_ptr<FsmNode>>& dReachableStates) 

    
// {
    


//     {
//         ++iterations;
//         stringstream ss;
// #ifdef ENABLE_DEBUG_MACRO
//         ss << "tC: ";
//         for (auto w : tC)
//         {
//             ss << *w << ", ";
//         }
//         LOG("VERBOSE_1") << ss.str() << std::endl;
//         ss.str(std::string());
//         ss << "t: ";
//         for (auto w : t)
//         {
//             ss << *w << ", ";
//         }
//         LOG("VERBOSE_1") << ss.str() << std::endl;
//         ss.str(std::string());
// #endif
//         LOG("VERBOSE_1") << "adaptiveTestCases as input traces:" << std::endl;
//         LOG("VERBOSE_1") << adaptiveList << std::endl;
//         map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;
//         size_t numberInputTraces = tC.size();
//         size_t inputTraceCount = 0;
//         // Applying all input traces from T_c to this FSM.
//         // All observed outputs are bein recorded.
//         // If the FSM observes a failure, adaptive state counting terminates.
//         for (const shared_ptr<InputTrace>& inputTrace : tC)
//         {
//             LOG("VERBOSE_1") << "############################################################" << std::endl;
//             LOG("VERBOSE_1") << "  Applying inputTrace " << ++inputTraceCount << " of " << numberInputTraces << ": " << *inputTrace << std::endl;
//             /**
//              * Hold the produced output traces for the current input trace.
//              */
//             vector<shared_ptr<OutputTrace>> producedOutputsSpec;
//             vector<shared_ptr<OutputTrace>> producedOutputsIut;
//             /**
//              * Hold the reached nodes for the current input trace.
//              */
//             vector<shared_ptr<FsmNode>> reachedNodesSpec;
//             vector<shared_ptr<FsmNode>> reachedNodesIut;

//             spec.apply(*inputTrace, producedOutputsSpec, reachedNodesSpec);
//             iut.apply(*inputTrace, producedOutputsIut, reachedNodesIut);
// #ifdef ENABLE_DEBUG_MACRO
//             ss << "    producedOutputs spec: ";
//             for (size_t i = 0; i < producedOutputsSpec.size(); ++i)
//             {
//                 ss << *producedOutputsSpec.at(i);
//                 if (i != producedOutputsSpec.size() - 1)
//                 {
//                     ss << ", ";
//                 }
//             }
//             LOG("VERBOSE_1") << ss.str() << std::endl;
//             ss.str(std::string());
//             ss << "    producedOutputs IUT: ";
//             for (size_t i = 0; i < producedOutputsIut.size(); ++i)
//             {
//                 ss << *producedOutputsIut.at(i);
//                 if (i != producedOutputsIut.size() - 1)
//                 {
//                     ss << ", ";
//                 }
//             }
//             LOG("VERBOSE_1") << ss.str() << std::endl;
//             ss.str(std::string());
//             ss << "    reachedNodes spec: ";
//             for (size_t i = 0; i < reachedNodesSpec.size(); ++i)
//             {
//                 ss << reachedNodesSpec.at(i)->getName();
//                 if (i != reachedNodesSpec.size() - 1)
//                 {
//                     ss << ", ";
//                 }
//             }
//             LOG("VERBOSE_1") << ss.str() << std::endl;
//             ss.str(std::string());
//             ss << "    reachedNodes IUT: ";
//             for (size_t i = 0; i < reachedNodesIut.size(); ++i)
//             {
//                 ss << reachedNodesIut.at(i)->getName();
//                 if (i != reachedNodesIut.size() - 1)
//                 {
//                     ss << ", ";
//                 }
//             }
//             LOG("VERBOSE_1") << ss.str() << std::endl;
//             ss.str(std::string());
// #endif
//             observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputsIut));

//             for (const shared_ptr<OutputTrace>& oTrace : producedOutputsIut)
//             {
//                 observedTraces.add(make_shared<const IOTrace>(*inputTrace, *oTrace));
//             }

//             LOG("VERBOSE_1") << "Checking produced outputs for failures" << std::endl;
//             //Chek if the IUT has produced any output that can not be produced by the specification.
//             for (size_t i = 0; i < producedOutputsIut.size(); ++i)
//             {
//                 const shared_ptr<OutputTrace>& outIut = producedOutputsIut.at(i);
//                 bool allowed = false;
//                 for (size_t j = 0; j < producedOutputsSpec.size(); ++j)
//                 {
//                     const shared_ptr<OutputTrace>& outSpec = producedOutputsSpec.at(j);
//                     if (*outIut == *outSpec)
//                     {
//                         allowed = true;
//                         // No need to apply adaptive test cases, if there are no adaptive test cases.
//                         if (adaptiveTestCases.size() > 0)
//                         {
//                             // Applying adaptive test cases to every node reached by the current input/output trace.
//                             LOG("VERBOSE_1") << "----------------- Getting adaptive traces -----------------" << std::endl;
//                             IOTraceContainer observedAdaptiveTracesIut;
//                             IOTraceContainer observedAdaptiveTracesSpec;
//                             const shared_ptr<FsmNode>& nodeIut = reachedNodesIut.at(i);
//                             const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

//                             iut.addPossibleIOTraces(nodeIut, adaptiveTestCases, observedAdaptiveTracesIut);
//                             spec.addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

//                             LOG("VERBOSE_1") << "  observedAdaptiveTracesIut (" << nodeIut->getName() << "): " << observedAdaptiveTracesIut << std::endl;
//                             LOG("VERBOSE_1") << "  observedAdaptiveTracesSpec (" << nodeSpec->getName() << "): " << observedAdaptiveTracesSpec << std::endl;

//                             bool failure = false;
//                             for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
//                             {
//                                 const shared_ptr<const IOTrace>& trace = *traceIt;
//                                 if (!observedAdaptiveTracesSpec.contains(trace))
//                                 {
//                                     LOG("INFO") << "  Specification does not contain " << *trace << std::endl;
//                                     failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
//                                     IOTrace traceCopy = IOTrace(*trace);
//                                     failTrace->append(traceCopy);
//                                     LOG("INFO") << "failTrace: " << *failTrace << std::endl;
//                                     failure = true;
//                                     break;
//                                 }
//                             }
//             //                PERFORMANCE_CHECKPOINT_WITH_ID(timerBlkObj, "after observedAdaptiveTracesIut loop");
//                             LOG("VERBOSE_1") << "  concatenating: " << *inputTrace << "/" << *outIut << std::endl;
//                             observedAdaptiveTracesIut.concatenateToFront(inputTrace, outIut);
//                             LOG("VERBOSE_1") << "  observedAdaptiveTraces after concatenation to front: " << observedAdaptiveTracesIut << std::endl;
//                             observedTraces.add(observedAdaptiveTracesIut);
//                             if (failure)
//                             {
//                                 // IUT produced an output that can not be produced by the specification.
//                                 LOG("INFO") << "  Failure observed:" << std::endl;
//                                 LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
//                                 LOG("INFO") << "    Observed adaptive traces:" << std::endl;
//                                 LOG("INFO") << observedAdaptiveTracesIut << std::endl;
//                                 LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
//                                 return false;
//                             }
//                         }
//                         // No failure observed, IUT output is allowed by specification.
//                         // No need to search through remaining specification outputs.
//                         break;
//                     }
//                 }
//                 if (!allowed)
//                 {
//                     // IUT produced an output that can not be produced by the specification.
//                     LOG("INFO") << "  Failure observed:" << std::endl;
//                     LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
//                     ss << "    Produced Outputs Iut: ";
//                     for (size_t i = 0; i < producedOutputsIut.size(); ++i)
//                     {
//                         ss << *producedOutputsIut.at(i);
//                         if (i != producedOutputsIut.size() - 1)
//                         {
//                             ss << ", ";
//                         }
//                     }
//                     LOG("INFO") << ss.str() << std::endl;
//                     ss.str(std::string());
// #ifdef ENABLE_DEBUG_MACRO
//                     ss << "    Produced Outputs Spec: ";
//                     for (size_t i = 0; i < producedOutputsSpec.size(); ++i)
//                     {
//                         ss << *producedOutputsSpec.at(i);
//                         if (i != producedOutputsSpec.size() - 1)
//                         {
//                             ss << ", ";
//                         }
//                     }
//                     ss << "    Reached nodes: ";
//                     for (size_t i = 0; i < reachedNodesIut.size(); ++i)
//                     {
//                         ss << reachedNodesIut.at(i)->getName();
//                         if (i != reachedNodesIut.size() - 1)
//                         {
//                             ss << ", ";
//                         }
//                     }
//                     LOG("VERBOSE_1") << ss.str() << std::endl;
//                     ss.str(std::string());
// #endif
//                     LOG("VERBOSE_1") << "Specification does not produce output " << *outIut << "." << std::endl;
//                     LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
//                     failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
//                     LOG("INFO") << "failTrace: " << *failTrace << std::endl;
//                     return false;
//                 }
//             }

//             if (producedOutputsIut.size() != reachedNodesIut.size())
//             {
//                 cerr << "Number of produced outputs and number of reached nodes do not match.";
//                 exit(EXIT_FAILURE);
//             }
//         }

//         long numberToCheck = 0;
//         LOG("VERBOSE_1") << "observedOutputsTCElements:" << std::endl;
//         for (auto e : observedOutputsTCElements)
//         {
//             LOG("VERBOSE_1") << "  " << *e.first << ":" << std::endl;
//             for (auto o : e.second)
//             {
//                 LOG("VERBOSE_1") << "    " << *o << std::endl;
//                 ++numberToCheck;
//             }
//         }
//         LOG("VERBOSE_1") << "Number of input/output combinations: " << numberToCheck << std::endl;
//         InputTraceSet newT = t;
//         InputTraceSet newTC;
//         inputTraceCount = 0;
//         for (shared_ptr<InputTrace> inputTrace : tC)
//         {
//             bool inputTraceMeetsCriteria = true;
//             LOG("INFO") << "check inputTrace: " << *inputTrace << " (" << ++inputTraceCount << " of " << numberInputTraces << ")" << std::endl;
//             vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
//             LOG("VERBOSE_1") << "producedOutputs:" << std::endl;
//             for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
//             {
//                 LOG("VERBOSE_1") << "  " << *outputTrace << std::endl;
//             }
//             long outputTraceCount = 0;
//             size_t numberOutputTraces = producedOutputs.size();

//             shared_ptr<const InputTrace> maxInputPrefixInV = nullptr;
//             for (const shared_ptr<InputTrace>& detStateTransition : detStateCover)
//             {
//                 if (inputTrace->isPrefix(*detStateTransition, false, true) &&
//                         ( !maxInputPrefixInV || maxInputPrefixInV->isEmptyTrace() || detStateTransition->size() > maxInputPrefixInV->size()))
//                 {
//                     maxInputPrefixInV = detStateTransition;
//                 }
//             }
//             if (!maxInputPrefixInV)
//             {
//                 stringstream ss;
// ss << "No prefix for input trace " << *inputTrace << " found in V. This should not happen.";
// std::cerr << ss.str();
// throw ss.str();
//             }
//             LOG("VERBOSE_1") << "maxInputPrefixInV: " << *maxInputPrefixInV << std::endl;

//             for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
//             {
//                 if (!inputTraceMeetsCriteria)
//                 {
//                     break;
//                 }
//                 LOG("INFO") << "outputTrace: " << *outputTrace << " (" << ++outputTraceCount << " of " << numberOutputTraces << ")" << std::endl;
//                 IOTrace currentTrace(*inputTrace, *outputTrace);
//                 LOG("VERBOSE_1") << "currentTrace (x_1/y_1): " << currentTrace << std::endl;
//                 bool outputTraceMeetsCriteria = false;
//                 vPrimeLazy.reset();

//                 LOG("VERBOSE_1") << "maxInputPrefixInV.size(): " << maxInputPrefixInV->size() << std::endl;
//                 shared_ptr<const IOTrace> maxIOPrefixInV = make_shared<const IOTrace>(*static_pointer_cast<const Trace>(maxInputPrefixInV),
//                                                                                       *outputTrace->getPrefix(maxInputPrefixInV->size(), true));
//                 LOG("VERBOSE_1") << "maxIOPrefixInV (v/v'): " << *maxIOPrefixInV << std::endl;
//                 IOTrace suffix(InputTrace(spec.presentationLayer), OutputTrace(spec.presentationLayer));
//                 suffix = currentTrace.getSuffix(*maxIOPrefixInV);
//                 LOG("VERBOSE_1") << "suffix (x/y): " << suffix << std::endl;

//                 LOG("VERBOSE_1") << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext() << std::endl;
//                 while (vPrimeLazy.hasNext())
//                 {
//                     const IOTraceContainer& vDoublePrime = vPrimeLazy.getNext();
//                     if (outputTraceMeetsCriteria)
//                     {
//                         break;
//                     }
//                     LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;

//                     if (!vDoublePrime.contains(maxIOPrefixInV))
//                     {
//                         LOG("VERBOSE_1") << "vDoublePrime does not contain prefix " << *maxIOPrefixInV << ". Skipping." << std::endl;
//                         LOG("VERBOSE_1") << "vPrimeLazy.hasNext(): " << vPrimeLazy.hasNext() << std::endl;
//                         continue;
//                     }
//                     for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
//                     {
//                         LOG("VERBOSE_1") << "rDistStates:" << std::endl;
//                         for (auto r : rDistStates)
//                         {
//                             LOG("VERBOSE_1") << "  " << r->getName() << std::endl;
//                         }
//                          //size_t lB = Fsm::lowerBound(*maxPrefix, suffix, t, rDistStates, adaptiveTestCases, vDoublePrime, dReachableStates, spec, iut);
//                         //LOG("VERBOSE_1") << "lB: " << lB << std::endl;
//                         bool exceedsBound = Fsm::exceedsBound(m, *maxIOPrefixInV, suffix, rDistStates, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
//                         LOG("VERBOSE_1") << "exceedsBound: " << exceedsBound << std::endl;
//                         if (exceedsBound)
//                         {
//                             LOG("VERBOSE_1") << "Exceeded lower bound. Output trace " << *outputTrace << " meets criteria." << std::endl;
//                             outputTraceMeetsCriteria = true;
//                             break;
//                         }
//                     }
//                 }
//                 if (outputTraceMeetsCriteria == false)
//                 {
//                     inputTraceMeetsCriteria = false;
//                 }
//             }

//             if (!inputTraceMeetsCriteria)
//             {
//                 // Keeping current input trace in T_C
//                 LOG("VERBOSE_1") << "Keeping " << *inputTrace << " in T_C." << std::endl;
//                 newTC.insert(inputTrace);
//                 // Next input trace.
//                 continue;
//             }
//             else
//             {
//                 LOG("VERBOSE_1") << "Removing " << *inputTrace << " from T_C." << std::endl;
//             }
//         }
//         ss << "newTC: ";
//         for (auto w : newTC)
//         {
//             ss << *w << ", ";
//         }
//         LOG("VERBOSE_1") << ss.str() << endl << std::endl;
//         ss.str(std::string());
//         // Expanding sequences.
//         InputTraceSet expandedTC;
//         InputTraceSet tracesAddedToT;
//         LOG("INFO") << "Expanding input sequences." << std::endl;
//         for (int x = 0; x <= spec.maxInput; ++x)
//         {
//             for (const shared_ptr<InputTrace>& inputTrace : newTC)
//             {

//                 shared_ptr<InputTrace> concat;
//                 if (inputTrace->isEmptyTrace())
//                 {
//                     concat = make_shared<InputTrace>(inputTrace->getPresentationLayer());
//                 }
//                 else
//                 {
//                     concat = make_shared<InputTrace>(*inputTrace);
//                 }

//                 concat->add(x);
//                 if (!InputTrace::contains(t, concat))
//                 {
//                     expandedTC.insert(concat);
//                 }

//                 if (newT.insert(concat).second)
//                 {
//                     tracesAddedToT.insert(concat);
//                 }
//             }
//         }
//         LOG("INFO") << "Finished expansion." << std::endl;
//         iut.bOmega(adaptiveTestCases, tracesAddedToT, bOmegaT);
//         LOG("INFO") << "Finished calculating bOmega." << std::endl;

//         ss << "expandedTC: ";
//         for (auto w : expandedTC)
//         {
//             ss << *w << ", ";
//         }
//         LOG("VERBOSE_1") << ss.str() << endl << std::endl;
//         ss.str(std::string());
//         ss << "newT: ";
//         for (auto w : newT)
//         {
//             ss << *w << ", ";
//         }
//         LOG("VERBOSE_1") << ss.str() << endl << std::endl;
//         ss.str(std::string());
//         tC = expandedTC;
//         t = newT;
//     }
//     LOG("VERBOSE_1") << "  RESULT: " << observedTraces << std::endl;
//     LOG("VERBOSE_1") << "IUT is a reduction of the specification." << std::endl;
    
    
//     return;
// }


int applyInputToSUT(int input) {
    
    // TODO: handle invalid input
    string inputStr = pl->getInId(input);
        
    // TODO: handle invalid output
    string outputStr = sut(inputStr);
    int output = pl->out2Num(outputStr);
    
    return output;
}

int applyInputToSUT(int input, string& outputStr) {
    
    // TODO: handle invalid input
    string inputStr = pl->getInId(input);
        
    // TODO: handle invalid output
    outputStr = sut(inputStr);
    int output = pl->out2Num(outputStr);
    
    return output;
}


IOTrace getSingleResponseToInputTrace(const InputTrace& inputTrace) {
    sut_reset();

    vector<int> outputs;
    
    for (auto traceIt = inputTrace.cbegin(); traceIt != inputTrace.cend(); ++traceIt)
    {
        outputs.push_back(applyInputToSUT(*traceIt));
    }

    OutputTrace outputTrace(outputs, inputTrace.getPresentationLayer());
    IOTrace ioTrace(inputTrace, outputTrace);
    return ioTrace;
}

IOTraceContainer getResponsesToInputTrace(const InputTrace& inputTrace, int repetitions) {
    IOTraceContainer cont;
    for (int i = 0; i < repetitions; ++i) {
        shared_ptr<const IOTrace> response = make_shared<const IOTrace>(getSingleResponseToInputTrace(inputTrace));
        cont.add(response);
    }
    return cont;
}

vector<shared_ptr<OutputTrace>> getOutputTracesToInputTrace(const InputTrace& inputTrace, int repetitions) {
    vector<shared_ptr<OutputTrace>> cont;
    for (int i = 0; i < repetitions; ++i) {
        IOTrace response = getSingleResponseToInputTrace(inputTrace);
        cont.push_back(make_shared<OutputTrace>(response.getOutputTrace()));
    }
    return cont;
}

vector<IOTraceContainer> calculateVPrime(const InputTraceSet& detStateCover) {
    vector<IOTraceContainer> vPrime;

    for(auto trace : detStateCover) {
        vPrime.push_back(getResponsesToInputTrace(*trace,k));
    }

    return vPrime;
}

void applyAdaptiveTestCaseAfterInputTrace(const InputOutputTree& adaptiveTestCase, const InputTrace& inputTrace, IOTrace& traceResponse, IOTrace& testCaseResponse) {
    sut_reset();

    vector<int> outputs;
    vector<int> responseInputs;
    vector<int> responseOutputs;
    
    for (auto traceIt = inputTrace.cbegin(); traceIt != inputTrace.cend(); ++traceIt)
    {
        outputs.push_back(applyInputToSUT(*traceIt));
    }

    OutputTrace outputTrace(outputs, inputTrace.getPresentationLayer());
    traceResponse = IOTrace(inputTrace,outputTrace);
    

    // shortcut for empty adaptive test case
    if (adaptiveTestCase.isEmpty()) {
        return;
    }

    
    shared_ptr<TreeNode> rootNode = adaptiveTestCase.getRoot();
    shared_ptr<AdaptiveTreeNode> node = static_pointer_cast<AdaptiveTreeNode>(rootNode);
    bool edgeExists = true; 
    while (!node->isLeaf() && edgeExists) {

        int input = node->getInput();
        responseInputs.push_back(input);

        int output = applyInputToSUT(input);
        responseOutputs.push_back(output);


        edgeExists = false;

        // search outgoing edges for the observed output
        for (const std::shared_ptr<TreeEdge> edge : *node->getChildren())
        {
            if (edge->getIO() == output)
            {
                node = static_pointer_cast<AdaptiveTreeNode>(edge->getTarget());
                edgeExists = true;
                break;
            }
        }

    }

    testCaseResponse = IOTrace(InputTrace(responseInputs, pl),OutputTrace(responseOutputs, pl));
}


/*
IOTraceContainer applyAdaptiveTestCasesAfterInputSequence(const IOTreeContainer& adaptiveTestCases, const InputTrace& trace, int maxInput, int repetitions) {
    IOTraceContainer cont;

    for (const auto tree : *adaptiveTestCases.getList()) {
        IOTraceContainer subCont;
        for (int i = 0; i < repetitions; ++i) {
            shared_ptr<const IOTrace> observedTrace = make_shared<const IOTrace>(applyAdaptiveTestCaseAfterInputTrace(*tree, trace, maxInput));
            subCont.add(observedTrace);
        }
    }

    return cont;
}

IOTraceContainer applyAdaptiveTestCasesAfterInputSequences(const IOTreeContainer& adaptiveTestCases, const InputTraceSet& traces, int maxInput, int repetitions) {
    IOTraceContainer cont;

    for (const auto trace : traces) {
        cont.add(applyAdaptiveTestCasesAfterInputSequence(adaptiveTestCases,*trace,maxInput,repetitions));
    }

    return cont;
}
*/


unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> collectResponseMapForInputTrace(const IOTreeContainer& adaptiveTestCases, const InputTrace& trace, int repetitions) {

    unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseSetsForIOTraces;

    for (const auto tree : *adaptiveTestCases.getList()) {
        for (int i = 0; i < repetitions; ++i) {
            IOTrace traceResponse(pl);
            IOTrace testCaseResponse(pl);
            applyAdaptiveTestCaseAfterInputTrace(*tree, trace, traceResponse, testCaseResponse);

            shared_ptr<const IOTrace> testCaseResponsePtr = make_shared<const IOTrace>(testCaseResponse);
            responseSetsForIOTraces[make_shared<IOTrace>(traceResponse)]->add(testCaseResponsePtr);
        }
    }
}


unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> collectResponseMapForInputTraces(const IOTreeContainer& adaptiveTestCases, const InputTraceSet& traces, int repetitions) {
    unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseMap;

    for (const auto trace : traces) {
        unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> traceMap = collectResponseMapForInputTrace(adaptiveTestCases,*trace,repetitions);
        responseMap.insert(traceMap.begin(), traceMap.end());
    }

    return responseMap;
}


unordered_set<shared_ptr<IOTraceContainer>> collectResponseSetsForInputTrace(const IOTreeContainer& adaptiveTestCases, const InputTrace& trace, int repetitions) {

    /*unordered_map<IOTrace,IOTraceContainer> responseSetsForIOTraces;

    for (const auto tree : *adaptiveTestCases.getList()) {
        for (int i = 0; i < repetitions; ++i) {
            IOTrace traceResponse(pl);
            IOTrace testCaseResponse(pl);
            applyAdaptiveTestCaseAfterInputTrace(*tree, trace, traceResponse, testCaseResponse);

            shared_ptr<const IOTrace> testCaseResponsePtr = make_shared<const IOTrace>(testCaseResponse);
            responseSetsForIOTraces[traceResponse].add(testCaseResponsePtr);
        }
    }*/

    unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseMap = collectResponseMapForInputTrace(adaptiveTestCases,trace,repetitions);

    unordered_set<shared_ptr<IOTraceContainer>> responseSets;
    for (auto kv : responseMap) {
        responseSets.insert(kv.second);
    }
    return responseSets;

}



unordered_set<shared_ptr<IOTraceContainer>> collectResponseSetsForInputTraces(const IOTreeContainer& adaptiveTestCases, const InputTraceSet& traces, int repetitions) {
    unordered_set<shared_ptr<IOTraceContainer>> responseSets;

    for (const auto trace : traces) {
        unordered_set<shared_ptr<IOTraceContainer>> traceSets = collectResponseSetsForInputTrace(adaptiveTestCases,*trace,repetitions);
        responseSets.insert(traceSets.begin(), traceSets.end());
    }

    return responseSets;
}



unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> collectResponseMapAndSetsForInputTraces(const IOTreeContainer& adaptiveTestCases, const InputTraceSet& traces, int repetitions, unordered_set<shared_ptr<IOTraceContainer>>& responseSets) {
    unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseMap;

    responseSets = unordered_set<shared_ptr<IOTraceContainer>>();

    for (const auto trace : traces) {
        unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> traceMap = collectResponseMapForInputTrace(adaptiveTestCases,*trace,repetitions);
        responseMap.insert(traceMap.begin(), traceMap.end());
    }

    for (auto kv : responseMap) {
        responseSets.insert(kv.second);
    }

    return responseMap;
}



int main(int argc, char* argv[])
{
    shared_ptr<IOTrace> failTrace;

    sut_init();

    parseParameters(argc,argv);
    readModel(modelType,modelFile,fsm);
    fsm->minimise();

    int maxInput = fsm->getMaxInput();
    int maxOutput = fsm->getMaxOutput();


    /** calculation of basic structures */
    fsm->calcRDistinguishableStates();
    IOListContainer rCharacterisationSet = fsm->getRCharacterisationSet();
    IOTreeContainer adaptiveTestCases = fsm->getAdaptiveRCharacterisationSet();
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = fsm->getMaximalSetsOfRDistinguishableStates();
    InputTraceSet detStateCover;
    const vector<shared_ptr<FsmNode>>& dReachableStates = fsm->calcDReachableStates(detStateCover);
    vector<IOTraceContainer> vPrime = calculateVPrime(detStateCover);

    VPrimeEnumerator vPrimeEnumerator(vPrime);

    
    IOTraceContainer observedTraces;

    /**
     * T - set of input sequences that have been followed by Ω.
     */
    InputTraceSet t = detStateCover;
    /**
     * Holds all B_Ω(T) for the current t.
     */
    //unordered_set<shared_ptr<IOTraceContainer>> bOmegaT;
    unordered_set<shared_ptr<IOTraceContainer>> responseSets;
    unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> responseMaps;

    //// TODO: Execute first iteration later
    responseMaps = collectResponseMapAndSetsForInputTraces(adaptiveTestCases,t,k,responseSets);
    LOG("VERBOSE_1") << "Response sets observed for:" << std::endl;
    for (const auto& kv : responseMaps)
    {
        LOG("VERBOSE_1") << "\t" << kv.first << std::endl;
        for (const auto& response : *kv.second) {
            LOG("VERBOSE_1") << "\t\t" << response << std::endl;
        }
    }

    /**
     * T_c - set of current elements of T: those that are being considered in the search
     * through state space. The elements in T_c are the maximal sequences considered that
     * do not meet the termination criterion.
     */
    InputTraceSet tC = detStateCover;
    int iterations = 0;

    while (tC.size() != 0) 
        {
        ++iterations;
        stringstream ss;

        map<shared_ptr<InputTrace>, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;
        size_t numberInputTraces = tC.size();
        size_t inputTraceCount = 0;

        unordered_map<shared_ptr<IOTrace>, shared_ptr<IOTraceContainer>> observedIUTResponsesToATCs = collectResponseMapForInputTraces(adaptiveTestCases,tC,k);
        
        // Applying all input traces from T_c to this FSM.
        // All observed outputs are being recorded.
        // If the FSM observes a failure, adaptive state counting terminates.
        for (const shared_ptr<InputTrace>& inputTrace : tC)
        {
            /**
             * Hold the produced output traces for the current input trace.
             */
            vector<shared_ptr<OutputTrace>> producedOutputsSpec;
            vector<shared_ptr<OutputTrace>> producedOutputsIut;
            /**
             * Hold the reached nodes for the current input trace.
             */
            vector<shared_ptr<FsmNode>> reachedNodesSpec;
            ////vector<shared_ptr<FsmNode>> reachedNodesIut;

            fsm->apply(*inputTrace, producedOutputsSpec, reachedNodesSpec);

            ////iut.apply(*inputTrace, producedOutputsIut, reachedNodesIut);
            producedOutputsIut = getOutputTracesToInputTrace(*inputTrace,k);

            observedOutputsTCElements.insert(make_pair(inputTrace, producedOutputsIut));

            for (const shared_ptr<OutputTrace>& oTrace : producedOutputsIut)
            {
                observedTraces.add(make_shared<const IOTrace>(*inputTrace, *oTrace));
            }

            LOG("VERBOSE_1") << "Checking produced outputs for failures" << std::endl;
            //Check if the IUT has produced any output that can not be produced by the specification.
            for (size_t i = 0; i < producedOutputsIut.size(); ++i)
            {
                const shared_ptr<OutputTrace>& outIut = producedOutputsIut.at(i);
                bool allowed = false;
                for (size_t j = 0; j < producedOutputsSpec.size(); ++j)
                {
                    const shared_ptr<OutputTrace>& outSpec = producedOutputsSpec.at(j);
                    if (*outIut == *outSpec)
                    {
                        allowed = true;
                        // No need to apply adaptive test cases, if there are no adaptive test cases.
                        if (adaptiveTestCases.size() > 0)
                        {
                            // get IOTraces for the two output traces combined with their shared input trace
                            shared_ptr<IOTrace> observedTraceSpec = make_shared<IOTrace>(IOTrace(*inputTrace,*outSpec));
                            shared_ptr<IOTrace> observedTraceIut = make_shared<IOTrace>(IOTrace(*inputTrace,*outIut));

                            // Applying adaptive test cases to every node reached by the current input/output trace.
                            LOG("VERBOSE_1") << "----------------- Getting adaptive traces -----------------" << std::endl;
                            IOTraceContainer observedAdaptiveTracesIut;
                            IOTraceContainer observedAdaptiveTracesSpec;
                            //const shared_ptr<FsmNode>& nodeIut = reachedNodesIut.at(i);
                            const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

                            ////iut.addPossibleIOTraces(nodeIut, adaptiveTestCases, observedAdaptiveTracesIut);
                            observedAdaptiveTracesIut = *observedIUTResponsesToATCs[observedTraceIut];
                            fsm->addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

                            
                            bool failure = false;
                            for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
                            {
                                const shared_ptr<const IOTrace>& trace = *traceIt;
                                if (!observedAdaptiveTracesSpec.contains(trace))
                                {
                                    LOG("INFO") << "  Specification does not contain " << *trace << std::endl;
                                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                                    IOTrace traceCopy = IOTrace(*trace);
                                    failTrace->append(traceCopy);
                                    LOG("INFO") << "failTrace: " << *failTrace << std::endl;
                                    failure = true;
                                    break;
                                }
                            }
                            LOG("VERBOSE_1") << "  concatenating: " << *inputTrace << "/" << *outIut << std::endl;
                            observedAdaptiveTracesIut.concatenateToFront(inputTrace, outIut);
                            LOG("VERBOSE_1") << "  observedAdaptiveTraces after concatenation to front: " << observedAdaptiveTracesIut << std::endl;
                            observedTraces.add(observedAdaptiveTracesIut);
                            if (failure)
                            {
                                // IUT produced an output that can not be produced by the specification.
                                LOG("INFO") << "  Failure observed:" << std::endl;
                                LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
                                LOG("INFO") << "    Observed adaptive traces:" << std::endl;
                                LOG("INFO") << observedAdaptiveTracesIut << std::endl;
                                LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
                                return false;
                            }
                        }
                        // No failure observed, IUT output is allowed by specification.
                        // No need to search through remaining specification outputs.
                        break;
                    }
                }
                if (!allowed)
                {
                    // IUT produced an output that can not be produced by the specification.
                    LOG("INFO") << "  Failure observed:" << std::endl;
                    LOG("INFO") << "    Input Trace: " << *inputTrace << std::endl;
                    ss << "    Produced Outputs Iut: ";
                    for (size_t i = 0; i < producedOutputsIut.size(); ++i)
                    {
                        ss << *producedOutputsIut.at(i);
                        if (i != producedOutputsIut.size() - 1)
                        {
                            ss << ", ";
                        }
                    }
                    LOG("INFO") << ss.str() << std::endl;
                    ss.str(std::string());
                    LOG("VERBOSE_1") << "Specification does not produce output " << *outIut << "." << std::endl;
                    LOG("VERBOSE_1") << "IUT is not a reduction of the specification." << std::endl;
                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                    LOG("INFO") << "failTrace: " << *failTrace << std::endl;
                    return false;
                }
            }


            // TODO: remove
            // if (producedOutputsIut.size() != reachedNodesIut.size())
            // {
            //     cerr << "Number of produced outputs and number of reached nodes do not match.";
            //     exit(EXIT_FAILURE);
            // }
        }

        long numberToCheck = 0;
        LOG("VERBOSE_1") << "observedOutputsTCElements:" << std::endl;
        for (auto e : observedOutputsTCElements)
        {
            LOG("VERBOSE_1") << "  " << *e.first << ":" << std::endl;
            for (auto o : e.second)
            {
                LOG("VERBOSE_1") << "    " << *o << std::endl;
                ++numberToCheck;
            }
        }
        LOG("VERBOSE_1") << "Number of input/output combinations: " << numberToCheck << std::endl;
        InputTraceSet newT = t;
        InputTraceSet newTC;
        inputTraceCount = 0;
        for (shared_ptr<InputTrace> inputTrace : tC)
        {
            bool inputTraceMeetsCriteria = true;
            LOG("INFO") << "check inputTrace: " << *inputTrace << " (" << ++inputTraceCount << " of " << numberInputTraces << ")" << std::endl;
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(inputTrace);
            LOG("VERBOSE_1") << "producedOutputs:" << std::endl;
            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                LOG("VERBOSE_1") << "  " << *outputTrace << std::endl;
            }
            long outputTraceCount = 0;
            size_t numberOutputTraces = producedOutputs.size();

            shared_ptr<const InputTrace> maxInputPrefixInV = nullptr;
            for (const shared_ptr<InputTrace>& detStateTransition : detStateCover)
            {
                if (inputTrace->isPrefix(*detStateTransition, false, true) &&
                        ( !maxInputPrefixInV || maxInputPrefixInV->isEmptyTrace() || detStateTransition->size() > maxInputPrefixInV->size()))
                {
                    maxInputPrefixInV = detStateTransition;
                }
            }
            if (!maxInputPrefixInV)
            {
                stringstream ss;
                ss << "No prefix for input trace " << *inputTrace << " found in V. This should not happen.";
                std::cerr << ss.str();
                throw ss.str();
            }
            LOG("VERBOSE_1") << "maxInputPrefixInV: " << *maxInputPrefixInV << std::endl;

            for (shared_ptr<OutputTrace> outputTrace : producedOutputs)
            {
                if (!inputTraceMeetsCriteria)
                {
                    break;
                }
                LOG("INFO") << "outputTrace: " << *outputTrace << " (" << ++outputTraceCount << " of " << numberOutputTraces << ")" << std::endl;
                IOTrace currentTrace(*inputTrace, *outputTrace);
                LOG("VERBOSE_1") << "currentTrace (x_1/y_1): " << currentTrace << std::endl;
                bool outputTraceMeetsCriteria = false;
                
                ////vPrimeLazy.reset();
                vPrimeEnumerator.reset();

                LOG("VERBOSE_1") << "maxInputPrefixInV.size(): " << maxInputPrefixInV->size() << std::endl;
                shared_ptr<const IOTrace> maxIOPrefixInV = make_shared<const IOTrace>(*static_pointer_cast<const Trace>(maxInputPrefixInV),
                                                                                      *outputTrace->getPrefix(maxInputPrefixInV->size(), true));
                LOG("VERBOSE_1") << "maxIOPrefixInV (v/v'): " << *maxIOPrefixInV << std::endl;
                
                
                ////IOTrace suffix(InputTrace(spec.presentationLayer), OutputTrace(spec.presentationLayer));
                IOTrace suffix(InputTrace( fsm->getPresentationLayer() ), OutputTrace( fsm->getPresentationLayer() ));

                
                suffix = currentTrace.getSuffix(*maxIOPrefixInV);
                LOG("VERBOSE_1") << "suffix (x/y): " << suffix << std::endl;



                // TODO: modify method of calculating next bOmegaT

                LOG("VERBOSE_1") << "vPrimeEnumerator.hasNext(): " << vPrimeEnumerator.hasNext() << std::endl;
                while (vPrimeEnumerator.hasNext())
                {
                    const IOTraceContainer& vDoublePrime = vPrimeEnumerator.getNext();
                    if (outputTraceMeetsCriteria)
                    {
                        break;
                    }
                    LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;

                    if (!vDoublePrime.contains(maxIOPrefixInV))
                    {
                        LOG("VERBOSE_1") << "vDoublePrime does not contain prefix " << *maxIOPrefixInV << ". Skipping." << std::endl;
                        LOG("VERBOSE_1") << "vPrimeEnumerator.hasNext(): " << vPrimeEnumerator.hasNext() << std::endl;
                        continue;
                    }
                    for (const vector<shared_ptr<FsmNode>>& rDistStates : maximalSetsOfRDistinguishableStates)
                    {
                        LOG("VERBOSE_1") << "rDistStates:" << std::endl;
                        for (auto r : rDistStates)
                        {
                            LOG("VERBOSE_1") << "  " << r->getName() << std::endl;
                        }
                        
                        bool exceedsBound = Fsm::exceedsBound(m, *maxIOPrefixInV, suffix, rDistStates, adaptiveTestCases, bOmegaT, vDoublePrime, dReachableStates, spec, iut);
                        LOG("VERBOSE_1") << "exceedsBound: " << exceedsBound << std::endl;
                        if (exceedsBound)
                        {
                            LOG("VERBOSE_1") << "Exceeded lower bound. Output trace " << *outputTrace << " meets criteria." << std::endl;
                            outputTraceMeetsCriteria = true;
                            break;
                        }
                    }
                }
                if (outputTraceMeetsCriteria == false)
                {
                    inputTraceMeetsCriteria = false;
                }
            }

            if (!inputTraceMeetsCriteria)
            {
                // Keeping current input trace in T_C
                LOG("VERBOSE_1") << "Keeping " << *inputTrace << " in T_C." << std::endl;
                newTC.insert(inputTrace);
                // Next input trace.
                continue;
            }
            else
            {
                LOG("VERBOSE_1") << "Removing " << *inputTrace << " from T_C." << std::endl;
            }
        }
        ss << "newTC: ";
        for (auto w : newTC)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        // Expanding sequences.
        InputTraceSet expandedTC;
        InputTraceSet tracesAddedToT;
        LOG("INFO") << "Expanding input sequences." << std::endl;
        for (int x = 0; x <= fsm->getMaxInput(); ++x)
        {
            for (const shared_ptr<InputTrace>& inputTrace : newTC)
            {

                shared_ptr<InputTrace> concat;
                if (inputTrace->isEmptyTrace())
                {
                    concat = make_shared<InputTrace>(inputTrace->getPresentationLayer());
                }
                else
                {
                    concat = make_shared<InputTrace>(*inputTrace);
                }

                concat->add(x);
                if (!InputTrace::contains(t, concat))
                {
                    expandedTC.insert(concat);
                }

                if (newT.insert(concat).second)
                {
                    tracesAddedToT.insert(concat);
                }
            }
        }
        LOG("INFO") << "Finished expansion." << std::endl;
        iut.bOmega(adaptiveTestCases, tracesAddedToT, bOmegaT);
        LOG("INFO") << "Finished calculating bOmega." << std::endl;

        ss << "expandedTC: ";
        for (auto w : expandedTC)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        ss << "newT: ";
        for (auto w : newT)
        {
            ss << *w << ", ";
        }
        LOG("VERBOSE_1") << ss.str() << endl << std::endl;
        ss.str(std::string());
        tC = expandedTC;
        t = newT;
    }
    LOG("VERBOSE_1") << "  RESULT: " << observedTraces << std::endl;
    LOG("VERBOSE_1") << "IUT is a reduction of the specification." << std::endl;


    
    
    exit(0);
    
}

