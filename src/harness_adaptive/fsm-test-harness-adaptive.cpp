/**
 * This is the main class for performing the online adaptive state counting desribed in
 *      R. M. Hierons, "Testing from a nondeterministic finite state machine using adaptive state counting," 
 *      in IEEE Transactions on Computers, vol. 53, no. 10, pp. 1330-1342, Oct. 2004, doi: 10.1109/TC.2004.85.
 * 
 * The specification is read from an FSM supplied in the arguments to main.
 * The SUT is to be accessed via a wrapper - see sut_wrapper_adaptive.h.
 * 
 * Example usage from the resources/harness_adative_example directory
 *  <path to fsm-harness-adaptive> -a 1 -k 100 -p example.in example.out example.state example.fsm 
 */

#include "sut_wrapper_adaptive.h"
// replace the include for sut_wrapper_adaptive.h with the following include for testing purposes.
//#include "sut_wrapper_adaptive_test.h"

#include "vPrimeEnumerator.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


#include "interface/FsmPresentationLayer.h"
#include "fsm/Dfsm.h"
#include "fsm/PkTable.h"
#include "fsm/FsmNode.h"
#include "fsm/IOTrace.h"
#include "fsm/InputTrace.h"
#include "fsm/SegmentedTrace.h"
#include "fsm/IOTraceContainer.h"

#include "trees/AdaptiveTreeNode.h"
#include "trees/IOListContainer.h"
#include "trees/IOTreeContainer.h"
#include "trees/OutputTree.h"
#include "trees/InputOutputTree.h"
#include "trees/TestSuite.h"
#include "trees/TreeEdge.h"

#include "utils/Logger.hpp"

#include "json/json.h"

#define DBG 0
using namespace std;
using namespace Json;


/**
 *  Hashing support for unordered sets of IOTrace-elements
 */
namespace std {
    template <> struct hash<unordered_set<IOTrace>>
    {
        size_t operator()(const unordered_set<IOTrace>& traces) const
        {

            size_t hashValue = 0;
            std::hash<IOTrace> hasher;

            for (const auto& t : traces)
            {
                // probably inefficient / non-uniform distribution
                hashValue += hasher(t);
            }

            return hashValue;
        }
    };
}


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


static size_t k; // number of times to repeat application of test
static size_t m; // upper bound on the number of states of the FSM-representation of the SUT
static size_t numAddStates; // m - (size of reference model)




/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    cerr << "usage: " << name
    << "[-a additionalstates] [-k input application repetitions] [-p presentation layer files (.in, .out, .states)] modelfile " << endl;
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
    
    // set m as <size of minimised and observable reference model> + <number of additional states>
    *myFsm = myFsm->minimise();
    m = myFsm->size() + numAddStates;
    
}







/**
 * Calculates the lower bound that may be placed on the number of states
 * of the FSM if there has been no repetition in the states of the product
 * machine for a given input/output sequence.
 *
 * Let x/y be the given input/output sequence, being observed when applying
 * a determinisitc state cover input sequence.
 * @param base Input/output sequence with the input sequence being an element
 * of the deterministic state cover and the output sequence being a possible
 * output to said input sequence. `base` is a prefix of x/y.
 * @param suffix Suffix x/y with base.suffix = x/y
 * @param states A maximum set of r-distinguishable states, that has been used
 * during the calculation of `base` and `suffix`.
 * @param responseSets Response sets observed during application of adaptive
 * test cases to the SUT.
 * @param responseMap Map from IO-traces observed in the SUT to the response 
 * set observed after that trace.
 * @param vDoublePrime A set of input/output sequences that represents one
 * possible combination of all input traces from the deterministic state cover
 * and theirs corresponding output sequences.
 * @param dReachableStates All d-reachable states.
 * @param spec The reference model.
 * @return The lower bound that may be placed on the number of states
 * of the FSM if there has been no repetition in the states of the product
 * machine for a given input/output sequence.
 */
size_t lowerBound(const IOTrace& base,
                       const IOTrace& suffix,
                       const vector<shared_ptr<FsmNode>>& states,
                       const unordered_set<unordered_set<IOTrace>>& responseSets,
                       const unordered_map<IOTrace, const shared_ptr<const unordered_set<IOTrace>>>& responseMap,
                       const IOTraceContainer& vDoublePrime,
                       const vector<shared_ptr<FsmNode>>& dReachableStates,
                       const Fsm& spec)
{
    // response sets observed along suffix applied after base
    // TODO: efficient as non-pointer set?
    unordered_set<unordered_set<IOTrace>> observedResponseSetsAlongSuffix;

    LOG("VERBOSE_1") << "lowerBound()" << std::endl;
    LOG("VERBOSE_1") << "base: " << base << std::endl;
    LOG("VERBOSE_1") << "suffix: " << suffix << std::endl;
    LOG("VERBOSE_1") << "states:" << std::endl;
    for (auto s : states)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    LOG("VERBOSE_1") << "vDoublePrime: " << vDoublePrime << std::endl;
    LOG("VERBOSE_1") << "dReachableStates: " << std::endl;
    for (auto s : dReachableStates)
    {
        LOG("VERBOSE_1") << "  " << s->getName() << std::endl;
    }
    size_t result = 0;
    LOG("VERBOSE_1") << "lb result: " << result << std::endl;

    // TODO: reduce overall logging?
    // LOG("VERBOSE_1") << "bOmegaT:" << std::endl;
    // for (const auto& cont : bOmegaT)
    // {
    //     LOG("VERBOSE_1") << "  " << cont << std::endl;
    // }

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
            
            auto findResult = responseMap.find(*trace);
            

            if (findResult == responseMap.end()) {
                LOG("VERBOSE_1") << "Trace " << *trace << " not observed in IUT " << std::endl;    
                continue;
            }

            const unordered_set<IOTrace>& traces = *findResult->second;

            LOG("VERBOSE_1") << "Removing from testTraces: {" << std::endl;
            for (const auto& t : traces) {
                LOG("VERBOSE_1") << "\t\t" << t  << std::endl;
            }
            LOG("VERBOSE_1") << "\t}" << std::endl;

            observedResponseSetsAlongSuffix.insert(traces);
        }


    }
    

    auto nonObservedResponseSets = responseSets.size() - observedResponseSetsAlongSuffix.size();
    LOG("VERBOSE_1") << "ResponseSets#           : " << responseSets.size()  << std::endl;
    if (responseSets.size() > m) {
        LOG("ERROR") << "WARNING: observed more response sets (" << responseSets.size() << ") than the IUT is assumed to have states (" << m << ")" << std::endl;
    }
    LOG("VERBOSE_1") << "ObservedResponseSets#   : " << observedResponseSetsAlongSuffix.size() << std::endl;
    LOG("VERBOSE_1") << "nonObservedResponseSets#: " << nonObservedResponseSets << std::endl;


    unsigned int response_set_counter = 0;
    LOG("VERBOSE_1") << "All observed response sets:" << std::endl;
    for (const auto& cont : responseSets)
    {
        LOG("VERBOSE_1") << "\tResponse set #" << response_set_counter++ << std::endl;
        
        for (const auto& t : cont) {
            LOG("VERBOSE_1") << "\t\t" << t  << std::endl;
        }
    }

    

    unsigned int observed_response_set_counter = 0;
    LOG("VERBOSE_1") << "Observed response sets along suffix:" << std::endl;
    for (const auto& cont : observedResponseSetsAlongSuffix)
    {
        LOG("VERBOSE_1") << "\tObserved response set #" << observed_response_set_counter++ << std::endl;
        for (const auto& t : cont) {
            LOG("VERBOSE_1") << "\t\t" << t  << std::endl;
        }
    }
    result += nonObservedResponseSets;

    LOG("VERBOSE_1") << "lowerBound() result: " << result << std::endl;
    return result;
}

/**
 * Calculates the lower bound that may be placed on the number of states
 * of the FSM if there has been no repetition in the states of the product
 * machine for a given input/output sequence and checks if it exceeds some
 * upper bound.
 *
 * Let x/y be the given input/output sequence, being observed when applying
 * a determinisitc state cover input sequence.
 * @param upperBound The upper bound to compare the lower bound against.
 * @param base Input/output sequence with the input sequence being an element
 * of the deterministic state cover and the output sequence being a possible
 * output to said input sequence. `base` is a prefix of x/y.
 * @param suffix Suffix x/y with base.suffix = x/y
 * @param states A maximum set of r-distinguishable states, that has been used
 * during the calculation of `base` and `suffix`.
 * @param responseSets Response sets observed during application of adaptive
 * test cases to the SUT.
 * @param responseMap Map from IO-traces observed in the SUT to the response 
 * set observed after that trace.
 * @param vDoublePrime A set of input/output sequences that represents one
 * possible combination of all input traces from the deterministic state cover
 * and theirs corresponding output sequences.
 * @param dReachableStates All d-reachable states.
 * @param spec The reference model.
 * @return The lower bound that may be placed on the number of states
 * of the FSM if there has been no repetition in the states of the product
 * machine for a given input/output sequence.
 */
bool exceedsBound(  const size_t upperBound,
                    const IOTrace& base,
                    const IOTrace& suffix,
                    const vector<shared_ptr<FsmNode>>& states,
                    const unordered_set<unordered_set<IOTrace>>& responseSets,
                    const unordered_map<IOTrace, const shared_ptr<const unordered_set<IOTrace>>>& responseMaps,
                    const IOTraceContainer& vDoublePrime,
                    const vector<shared_ptr<FsmNode>>& dReachableStates,
                    const Fsm& spec)
{
    size_t lB = lowerBound(base, suffix, states, responseSets, responseMaps, vDoublePrime, dReachableStates, spec);
    LOG("VERBOSE_1") << "lB: " << lB << std::endl;
    return lB > upperBound;
}







/**
 * Apply a single input to the SUT.
 * @param input The input to apply.
 * @return The index of the observed output. -1 if the string returned by the SUT is not a known output.
 */
int applyInputToSUT(int input) {
    
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "Applying input: " << input << std::endl;

    string inputStr = pl->getInId(input);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tInput string is: " << inputStr << std::endl;        

    string outputStr = sut(inputStr);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tOutput string is: " << outputStr << std::endl;
    
    int output = pl->out2Num(outputStr);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tOutput is: " << output << std::endl;
    
    return output;
}

/**
 * Apply a single input to the SUT.
 * @param input The input to apply.
 * @param outputStr Output parameter to store the string returned by the SUT.
 * @return The index of the observed output. -1 if the string returned by the SUT is not a known output.
 */
int applyInputToSUT(int input, string& outputStr) {
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "Applying input: " << input << std::endl;
   
    // TODO: handle invalid input
    string inputStr = pl->getInId(input);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tInput string is: " << inputStr << std::endl;
        
    // TODO: handle invalid output
    outputStr = sut(inputStr);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tOutput string is: " << outputStr << std::endl;

    int output = pl->out2Num(outputStr);
    LOG("VERBOSE_SUT_APPLICATIONS_2") << "\tOutput is: " << output << std::endl;
    
    return output;
}


/**
 * Apply a single input trace to the SUT.
 * @param inputTrace The input trace to apply.
 * @param trace Output parameter to store the observed IO-trace.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool getSingleResponseToInputTrace(const InputTrace& inputTrace, IOTrace& trace) {

    LOG("VERBOSE_SUT_APPLICATIONS_1") << "Applying input trace: " << inputTrace << std::endl;

    sut_reset();

    vector<int> outputs;
    
    for (auto traceIt = inputTrace.cbegin(); traceIt != inputTrace.cend(); ++traceIt)
    {
        string outputStr;
        int output = applyInputToSUT(*traceIt, outputStr);
        
        if (output == -1) {
            
            LOG("ERROR") << "Failure detected: Observed invalid output \"" << outputStr << "\" in trace "; 

            if (outputs.size() == 0) {
                LOG("ERROR") << "(" << pl->getInId(inputTrace.get().at(0)) << "/" << outputStr << ")";
            } else {
                LOG("ERROR") << "(" << pl->getInId(inputTrace.get().at(0)) << "/" << pl->getOutId(outputs.at(0)) << ")";
                for (unsigned int i = 1; i < outputs.size(); ++i) { 
                    LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(i)) << "/" << pl->getOutId(outputs.at(i)) << ")";
                }
                LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(0)) << "/" << outputStr << ")";
            }
            LOG("ERROR") << std::endl;

            return false;
        }
        
        outputs.push_back(output);
    }

    OutputTrace outputTrace(outputs, inputTrace.getPresentationLayer());
    trace = IOTrace(inputTrace, outputTrace);
    return true;
}

/**
 * Apply a single input trace to the SUT multiple times.
 * @param inputTrace The input trace to apply.
 * @param repetition The number of times to apply the input trace.
 * @param cont Output parameter to store the observed IO-traces.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool getResponsesToInputTrace(const InputTrace& inputTrace, int repetitions, IOTraceContainer& cont) {
    for (int i = 0; i < repetitions; ++i) {
        IOTrace t(pl);
        if (!getSingleResponseToInputTrace(inputTrace, t)) {
            return false;
        }
        shared_ptr<const IOTrace> response = make_shared<const IOTrace>(t);
        cont.add(response);
    }
    return true;
}

/**
 * Calculate V''.
 * @param detStateCover The deterministic state cover to apply.
 * @param vPrime Output parameter to store the observed IO-traces for each trace in detStateCover.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool calculateVPrime(const InputTraceSet& detStateCover, vector<IOTraceContainer>& vPrime) {
    
    for(auto trace : detStateCover) {
        IOTraceContainer vResponses;
        if (!getResponsesToInputTrace(*trace,k,vResponses)) {
            return false;
        }
        vPrime.push_back(vResponses);
    }

    return true;
}

/**
 * Apply a single adaptive test case after applying a single input trace.
 * @param adaptiveTestCase The adaptive test case to apply.
 * @param inputTrace Input trace to apply.
 * @param traceResponse Output parameter to store the observed IO-trace to inputTrace.
 * @param testCaseResponse Output parameter to store the response to the adaptive test case.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool applyAdaptiveTestCaseAfterInputTrace(const InputOutputTree& adaptiveTestCase, const InputTrace& inputTrace, IOTrace& traceResponse, IOTrace& testCaseResponse) {
    sut_reset();

    vector<int> outputs;
    vector<int> responseInputs;
    vector<int> responseOutputs;

    LOG("VERBOSE_SUT_APPLICATIONS_1") << "Applying input trace: " << inputTrace << " to be followed by an ATC" << std::endl;
    
    for (auto traceIt = inputTrace.cbegin(); traceIt != inputTrace.cend(); ++traceIt)
    {
        //outputs.push_back(applyInputToSUT(*traceIt));
        string outputStr;
        int output = applyInputToSUT(*traceIt, outputStr);
        
        if (output == -1) {
            // TODO: terminate execution
            LOG("ERROR") << "Failure detected: Observed invalid output \"" << outputStr << "\" in trace "; 

            if (outputs.size() == 0) {
                LOG("ERROR") << "(" << pl->getInId(inputTrace.get().at(0)) << "/" << outputStr << ")";
            } else {
                LOG("ERROR") << "(" << pl->getInId(inputTrace.get().at(0)) << "/" << pl->getOutId(outputs.at(0)) << ")";
                for (unsigned int i = 1; i < outputs.size(); ++i) { 
                    LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(i)) << "/" << pl->getOutId(outputs.at(i)) << ")";
                }
                LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(0)) << "/" << outputStr << ")";
            }

            LOG("ERROR") << std::endl;
            return false;
            
        }
        
        outputs.push_back(output);
    }

    OutputTrace outputTrace(outputs, inputTrace.getPresentationLayer());
    traceResponse = IOTrace(inputTrace,outputTrace);
    

    // shortcut for empty adaptive test case
    if (adaptiveTestCase.isEmpty()) {
        return true;
    }

    
    shared_ptr<TreeNode> rootNode = adaptiveTestCase.getRoot();
    shared_ptr<AdaptiveTreeNode> node = static_pointer_cast<AdaptiveTreeNode>(rootNode);
    bool edgeExists = true; 
    while (!node->isLeaf() && edgeExists) {

        int input = node->getInput();
        responseInputs.push_back(input);

        
        string outputStr;
        int output = applyInputToSUT(input, outputStr);
        
        if (output == -1) {
            
            LOG("ERROR") << "Failure detected: Observed invalid output \"" << outputStr << "\" in trace "; 

            if (outputs.size() > 0) {
                LOG("ERROR") << "(" << pl->getInId(inputTrace.get().at(0)) << "/" << pl->getOutId(outputs.at(0)) << ")";
                for (unsigned int i = 1; i < outputs.size(); ++i) { 
                    LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(i)) << "/" << pl->getOutId(outputs.at(i)) << ")";
                }
                LOG("ERROR") << ".(" << pl->getInId(inputTrace.get().at(0)) << "/" << outputStr << ")";
                LOG("ERROR") << ".";
            }

            if (responseOutputs.size() == 0) {
                LOG("ERROR") << "(" << pl->getInId(input) << "/" << outputStr << ")";
            } else {
                LOG("ERROR") << "(" << pl->getInId(responseInputs.at(0)) << "/" << pl->getOutId(responseOutputs.at(0)) << ")";
                for (unsigned int i = 1; i < responseOutputs.size(); ++i) { 
                    LOG("ERROR") << ".(" << pl->getInId(responseInputs.at(i)) << "/" << pl->getOutId(responseOutputs.at(i)) << ")";
                }
                LOG("ERROR") << ".(" << pl->getInId(input) << "/" << outputStr << ")";
            }

            LOG("ERROR") << std::endl;
            return false;
            
        }


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

    return true;
}



/**
 * Apply multiple adaptive test case, each after applying a single input trace.
 * @param adaptiveTestCases The adaptive test cases to apply.
 * @param trace Input trace to apply.
 * @param repetition How often the input trace is to be followed by each adaptive test case.
 * @param responseMap Output parameter to store a map that uses observed IO-traces to the 
 *                    input trace as keys mapping to the response set of the state reached 
 *                    by that IO-trace to the adaptive test cases.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool collectResponseMapForInputTrace(const IOTreeContainer& adaptiveTestCases, const InputTrace& trace, int repetitions, unordered_map<IOTrace, shared_ptr<const unordered_set<IOTrace>>>& responseMap) {

    unordered_map<IOTrace, shared_ptr<unordered_set<IOTrace>>> responseSetsForIOTraces;

    unsigned int testCaseNo = 0;
    for (const auto tree : *adaptiveTestCases.getList()) {
        ++testCaseNo;
        for (int i = 0; i < repetitions; ++i) {

            LOG("VERBOSE_SUT_APPLICATIONS_1") << "Applying input trace " << trace << " followed by ATC #" << testCaseNo << " - repetition #" << i << std::endl;

            IOTrace traceResponse(pl);
            IOTrace testCaseResponse(pl);
            bool applicationSuccessful = applyAdaptiveTestCaseAfterInputTrace(*tree, trace, traceResponse, testCaseResponse);

            if (!applicationSuccessful) {
                return false;
            }

            LOG("VERBOSE_SUT_APPLICATIONS_1") << "\tTrace response   : " << traceResponse << std::endl;
            LOG("VERBOSE_SUT_APPLICATIONS_1") << "\tTestcase response: " << testCaseResponse << std::endl;

            auto findResult = responseSetsForIOTraces.find(traceResponse);
            if (findResult == responseSetsForIOTraces.end()) {
                responseSetsForIOTraces.emplace(traceResponse, make_shared<unordered_set<IOTrace>>(unordered_set<IOTrace>()));
            }

            responseSetsForIOTraces.at(traceResponse)->insert(testCaseResponse);
        }
    }

    // TODO: no better way?
    for (const auto& kv : responseSetsForIOTraces) {

        responseMap.emplace(kv.first,kv.second);
    }

    return true;
}


/**
 * Apply multiple adaptive test case, each after applying a single input trace.
 * Performed for each input trace in a given set of traces.
 * @param adaptiveTestCases The adaptive test cases to apply.
 * @param trace Input trace to apply.
 * @param repetition How often the input trace is to be followed by each adaptive test case.
 * @param responseSet Set into which the values of responseMap are inserted.
 * @param responseMap Output parameter to store a map that uses observed IO-traces to the 
 *                    input traces as keys mapping to the response set of the state reached 
 *                    by that IO-trace to the adaptive test cases.
 * @return True if no invalid output has been observed during application. False otherwise.
 */
bool collectResponseMapAndSetsForInputTraces(const IOTreeContainer& adaptiveTestCases, const InputTraceSet& traces, int repetitions, unordered_set<unordered_set<IOTrace>>& responseSets, unordered_map<IOTrace, const shared_ptr<const unordered_set<IOTrace>>>& responseMap) {
    for (const auto trace : traces) {
        unordered_map<IOTrace, shared_ptr<const unordered_set<IOTrace>>> traceMap;
        bool applicationSuccessful = collectResponseMapForInputTrace(adaptiveTestCases,*trace,repetitions, traceMap);
        if (!applicationSuccessful) {
            return false;
        }
        responseMap.insert(traceMap.begin(), traceMap.end());
    }

    for (auto kv : responseMap) {
        responseSets.insert(*kv.second);
    }

    return true;
}



/**
 *  Log an adaptive test case
 * @param node Root of the test case to log
 * @param depth Indentation depth
 */
void printATC(AdaptiveTreeNode& node, unsigned int depth) {
    
    if (node.isLeaf()) {
        LOG("VERBOSE_SPEC") << string(depth, '\t') << "INPUT: " << pl->getInId(node.getInput()) << " (" << node.getInput() << ") (LEAF)" << std::endl;
        return;
    }

    LOG("VERBOSE_SPEC") << string(depth, '\t') << "INPUT: " << pl->getInId(node.getInput()) << " (" << node.getInput() << ")" << std::endl;
    for (const std::shared_ptr<TreeEdge>& edge : *node.getChildren())
    {
        LOG("VERBOSE_SPEC") << string(depth+1, '\t') << "OUTPUT: " << pl->getOutId(edge->getIO()) << " (" << edge->getIO() << ")" << std::endl;
        printATC(*static_pointer_cast<AdaptiveTreeNode>(edge->getTarget()), depth+2);
    }

}




int main(int argc, char* argv[])
{

    // set logging levels
    LogCoordinator::getStandardLogger().bindAllToDevNull();
    
    LogCoordinator::getStandardLogger().createLogTargetAndBind("ERROR", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("WARNING", std::cout);    
    LogCoordinator::getStandardLogger().createLogTargetAndBind("FATAL", std::cout);
    
    /*
    LogCoordinator::getStandardLogger().createLogTargetAndBind("INFO", std::cout);

    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_SPEC", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_FAILURE_CHECK", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_SUT_APPLICATIONS_1", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_SUT_APPLICATIONS_2", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_SUT_INTERNAL", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_1", std::cout);
    LogCoordinator::getStandardLogger().createLogTargetAndBind("VERBOSE_2", std::cout);
    */

    shared_ptr<IOTrace> failTrace;

    

    parseParameters(argc,argv);
    readModel(modelType,modelFile,fsm);

    sut_init();    
    // replace the call to sut_init with the following call if sut_wrapper_adaptive_test.h is used
    //sut_init(fsm,false);


    /** calculation of basic structures */
    fsm->calcRDistinguishableStates();
    IOListContainer rCharacterisationSet = fsm->getRCharacterisationSet();
    IOTreeContainer adaptiveTestCases = fsm->getAdaptiveRCharacterisationSet();
    IOListContainer adaptiveList = adaptiveTestCases.toIOList();
    const vector<vector<shared_ptr<FsmNode>>>& maximalSetsOfRDistinguishableStates = fsm->getMaximalSetsOfRDistinguishableStates();
    InputTraceSet detStateCover;
    const vector<shared_ptr<FsmNode>>& dReachableStates = fsm->calcDReachableStates(detStateCover);
    
    // remove epsilons
    for (auto& x : detStateCover) {
        *x = x->removeEpsilon();
    }

    vector<IOTraceContainer> vPrime;
    if (!calculateVPrime(detStateCover,vPrime)) {
        return false;
    }
    VPrimeEnumerator vPrimeEnumerator(vPrime);

    unsigned int vDoublePrime_set_counter = 0;
    LOG("VERBOSE_SPEC") << "V''s: " << std::endl;
    while (vPrimeEnumerator.hasNext()) {
        LOG("VERBOSE_SPEC") << "\tV'' #" << vDoublePrime_set_counter++ << std::endl;
        IOTraceContainer vDoublePrime = vPrimeEnumerator.getNext();
        for (const auto& v : vDoublePrime) {
            LOG("VERBOSE_SPEC") << "\t\t" << *v << std::endl;
        }
    }

    unsigned int rd_set_counter = 0;
    LOG("VERBOSE_SPEC") << "Maximal r-d state sets: " << std::endl;
    for (const auto& d : maximalSetsOfRDistinguishableStates) {
        LOG("VERBOSE_SPEC") << "\t" <<  "Maximal r-d set #" << rd_set_counter++ << std::endl;
        for (const auto& s : d) {
            LOG("VERBOSE_SPEC") << "\t\t" <<  s->getId() << std::endl;
        }
    }

    LOG("VERBOSE_SPEC") << "Det. state cover: " << std::endl;
    for (const auto& x : detStateCover) {
        LOG("VERBOSE_SPEC") << "\t" <<  *x << std::endl;
    }
    
    LOG("VERBOSE_SPEC") << "D-r states: " << std::endl;
    for (const auto& s : dReachableStates) {
        LOG("VERBOSE_SPEC") << "\t" <<  s->getId() << std::endl;
    }

    unsigned int vprime_set_counter = 0;
    LOG("VERBOSE_SPEC") << "V': " << std::endl;
    for (const auto& d : vPrime) {
        LOG("VERBOSE_SPEC") << "\t" <<  "V' subset #" << vprime_set_counter++ << std::endl;
        for (auto traceIt = d.cbegin(); traceIt != d.cend(); ++traceIt) {
            const shared_ptr<const IOTrace>& trace = *traceIt;
            LOG("VERBOSE_SPEC") << "\t\t" << *trace << std::endl;
        }
    }

    // print ATCs
    LOG("VERBOSE_SPEC") << "ATCs: " << std::endl;
    unsigned int atcNo = 0;
    for (const auto& atcIt : *adaptiveTestCases.getList()) {
        LOG("VERBOSE_SPEC") << "\tATC #" << ++atcNo << std::endl;
        
        const InputOutputTree& atc = *atcIt;
        shared_ptr<AdaptiveTreeNode> root = static_pointer_cast<AdaptiveTreeNode>(atc.getRoot());
        printATC(*root,2);
    }

    
    unordered_set<IOTrace> observedTraces;

    /**
     * T - set of input sequences that have been followed by Ω.
     */
    InputTraceSet t = detStateCover;
    /**
     * Holds all B_Ω(T) for the current t.
     */
    

    // contains all observed response sets, not only the current ones
    unordered_set<unordered_set<IOTrace>> responseSets;
    
    // contains only current response sets
    unordered_map<IOTrace, const shared_ptr<const unordered_set<IOTrace>>> responseMap;

    

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

        unordered_map<InputTrace, vector<shared_ptr<OutputTrace>>> observedOutputsTCElements;
        size_t numberInputTraces = tC.size();
        size_t inputTraceCount = 0;

        
        unordered_map<IOTrace, const shared_ptr<const unordered_set<IOTrace>>> currentResponseMap;
        if (!collectResponseMapAndSetsForInputTraces(adaptiveTestCases,tC,k,responseSets,currentResponseMap)) {
            return false;
        }
        LOG("VERBOSE_1") << "Response sets observed for:" << std::endl;
        for (const auto& kv : currentResponseMap)
        {
            LOG("VERBOSE_1") << "\t" << kv.first << std::endl;
            for (auto traceIt = kv.second->cbegin(); traceIt != kv.second->cend(); ++traceIt) {
                
                const IOTrace& trace = *traceIt;
                LOG("VERBOSE_1") << "\t\t" << trace << std::endl;
            }

            // add currently observed response sets to previously observed response sets
            responseMap.insert(kv);
        }

        for (const auto& kv : responseMap)
        {
            auto findResult = observedOutputsTCElements.find(kv.first.getInputTrace());
            if (findResult == observedOutputsTCElements.end()) {
                observedOutputsTCElements.emplace(kv.first.getInputTrace(), vector<shared_ptr<OutputTrace>>());
            }
            observedOutputsTCElements.at(kv.first.getInputTrace()).push_back(make_shared<OutputTrace>(kv.first.getOutputTrace()));
        }

        LOG("VERBOSE_SUT_APPLICATIONS_1") << "Observed responses to tC" << std::endl;
        for (const auto& kv : observedOutputsTCElements) {
            LOG("VERBOSE_SUT_APPLICATIONS_1") << "\tObserved responses to " << kv.first << std::endl;
            for (const auto& resp : kv.second) {
                LOG("VERBOSE_SUT_APPLICATIONS_1") << "\t\t " << *resp << std::endl;
            }
        }



        
        // Applying all input traces from T_c to this FSM.
        // All observed outputs are being recorded.
        // If the FSM observes a failure, adaptive state counting terminates.
        for (const shared_ptr<InputTrace>& inputTrace : tC)
        {
            
            /**
             * Hold the produced output traces for the current input trace.
             */
            vector<shared_ptr<OutputTrace>> producedOutputsSpec;
            
            

            // TODO: use observedOutputsTCElements
            vector<shared_ptr<OutputTrace>>& producedOutputsIut = observedOutputsTCElements[*inputTrace];

            /**
             * Hold the reached nodes for the current input trace.
             */
            vector<shared_ptr<FsmNode>> reachedNodesSpec;
            
            fsm->apply(*inputTrace, producedOutputsSpec, reachedNodesSpec);

            // add empty output trace in case the empty trace is applied
            if (inputTrace->get().size() == 0) {
                vector<int> outputs; // empty output vector
                producedOutputsSpec.push_back(make_shared<OutputTrace>(outputs,inputTrace->getPresentationLayer()));
            }

            

            for (const shared_ptr<OutputTrace>& oTrace : producedOutputsIut)
            {
                observedTraces.insert(IOTrace(*inputTrace, *oTrace));
            }

            

            LOG("VERBOSE_1") << "Checking produced outputs for failures" << std::endl;
            LOG("VERBOSE_1") << "\tinputTrace            : " << *inputTrace << std::endl;
            //Check if the IUT has produced any output that can not be produced by the specification.
            
            for (size_t i = 0; i < producedOutputsIut.size(); ++i)
            {
                LOG("VERBOSE_1") << "\t\toutputTrace (IUT) : " << *producedOutputsIut.at(i) << std::endl;

                LOG("VERBOSE_1") << "\t\toutputTraces (REF): " << std::endl;
                for (const auto& t : producedOutputsSpec) {
                    LOG("VERBOSE_1") << "\t\t\t" << *t << std::endl;
                }

                const shared_ptr<OutputTrace>& outIut = producedOutputsIut.at(i);
                bool allowed = false;
                for (size_t j = 0; j < producedOutputsSpec.size(); ++j)
                {
                    
                    const shared_ptr<OutputTrace>& outSpec = producedOutputsSpec.at(j);

                    LOG("VERBOSE_1") << "\tinputTrace       : " << *inputTrace << std::endl;
                    LOG("VERBOSE_1") << "\t\tchecking outputTrace (IUT): " << *outIut << std::endl;
                    LOG("VERBOSE_1") << "\t\tagainst  outputTrace (REF): " << *outSpec << std::endl;

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
                            
                            unordered_set<IOTrace> observedAdaptiveTracesIut = *responseMap[*observedTraceIut];
                            unordered_set<IOTrace> observedAdaptiveTracesSpec;

                            const shared_ptr<FsmNode>& nodeSpec = reachedNodesSpec.at(j);

                            
                            fsm->addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpec);

                            
                            bool failure = false;
                            for (auto traceIt = observedAdaptiveTracesIut.cbegin(); traceIt != observedAdaptiveTracesIut.cend(); ++traceIt)
                            
                            {
                                const IOTrace& trace = *traceIt;

                                //if (!observedAdaptiveTracesSpec.contains(trace))
                                if (observedAdaptiveTracesSpec.find(trace) == observedAdaptiveTracesSpec.end())
                                {
                                    LOG("INFO") << "\t\tinputTrace       : " << *inputTrace << std::endl;
                                    LOG("INFO") << "\t\toutputTrace (IUT): " << *outIut << std::endl;
                                    LOG("INFO") << "\t\toutputTrace (REF): " << *outSpec << std::endl;

                                    LOG("INFO") << "\tobserved responses IUT:" << std::endl;
                                    for (const auto& t : observedAdaptiveTracesIut) {
                                        LOG("INFO") << "\t\t" << t << std::endl;
                                    }

                                    LOG("INFO") << "\tobserved responses REF:" << std::endl;
                                    for (const auto& t : observedAdaptiveTracesSpec) {
                                        LOG("INFO") << "\t\t" << t << std::endl;
                                    }

                                    IOTraceContainer observedAdaptiveTracesSpecC;
                                    fsm->addPossibleIOTraces(nodeSpec, adaptiveTestCases, observedAdaptiveTracesSpecC);
                                    LOG("INFO") << "\tobserved responses REF (C):" << std::endl;
                                    for (const auto& t : observedAdaptiveTracesSpecC) {
                                        LOG("INFO") << "\t\t" << *t << std::endl;
                                    }

                                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                                    LOG("ERROR") << "  Specification does not contain " << trace << " after trace " << *failTrace << std::endl;
                                    IOTrace traceCopy = IOTrace(trace);
                                    failTrace->append(traceCopy);
                                    LOG("ERROR") << "  failTrace: " << *failTrace << std::endl;
                                    failure = true;
                                    break;
                                }
                            }
                            LOG("VERBOSE_1") << "  concatenating: " << *inputTrace << "/" << *outIut << std::endl;
                            
                            IOTrace newIOTrace = IOTrace(*inputTrace, *outIut);
                            unordered_set<IOTrace> newSet;
                            for (const auto& oldTrace : observedAdaptiveTracesIut)
                            {
                                const IOTrace newTrace = IOTrace(oldTrace, newIOTrace, true);
                                newSet.insert(newTrace);
                            }
                            observedAdaptiveTracesIut = newSet;


                            LOG("VERBOSE_1") << "  observedAdaptiveTraces after concatenation to front: " << std::endl;
                            for (const auto& t : observedAdaptiveTracesIut) {
                                LOG("VERBOSE_1") << "\t\t" << t << std::endl;
                            }
                            
                            
                            for (const auto& t : observedAdaptiveTracesIut) {
                                observedTraces.insert(t);
                            }
                            if (failure)
                            {
                                // IUT produced an output that can not be produced by the specification.
                                LOG("ERROR") << "  Failure observed after input trace: " << *inputTrace << std::endl;
                                LOG("ERROR") << "  Observed responses of the SUT on " << *inputTrace << " and extensions: " << std::endl;
                                for (const auto& t : observedAdaptiveTracesIut) {
                                    LOG("ERROR") << "\t" << t << std::endl;
                                }
                                LOG("ERROR") << "IUT is not a reduction of the specification." << std::endl;
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
                    LOG("ERROR") << "  Failure observed:" << std::endl;
                    LOG("ERROR") << "    Input Trace: " << *inputTrace << std::endl;
                    ss << "    Produced Outputs Iut: ";
                    for (size_t i = 0; i < producedOutputsIut.size(); ++i)
                    {
                        ss << *producedOutputsIut.at(i);
                        if (i != producedOutputsIut.size() - 1)
                        {
                            ss << ", ";
                        }
                    }
                    LOG("ERROR") << ss.str() << std::endl;
                    ss.str(std::string());
                    LOG("ERROR") << "Specification does not produce output " << *outIut << " in response to " << *inputTrace << std::endl;
                    LOG("ERROR") << "IUT is not a reduction of the specification." << std::endl;
                    failTrace = make_shared<IOTrace>(*inputTrace, *outIut);
                    LOG("ERROR") << "failTrace: " << *failTrace << std::endl;
                    return false;
                }
            }
        }

        long numberToCheck = 0;
        LOG("VERBOSE_1") << "observedOutputsTCElements:" << std::endl;
        for (auto e : observedOutputsTCElements)
        {
            LOG("VERBOSE_1") << "  " << e.first << ":" << std::endl;
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
            vector<shared_ptr<OutputTrace>>& producedOutputs = observedOutputsTCElements.at(*inputTrace);
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
                LOG("ERROR") << ss.str();
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
                
                vPrimeEnumerator.reset();

                LOG("VERBOSE_1") << "maxInputPrefixInV.size(): " << maxInputPrefixInV->size() << std::endl;
                
                shared_ptr<const IOTrace> maxIOPrefixInV = make_shared<const IOTrace>(*static_pointer_cast<const Trace>(maxInputPrefixInV), outputTrace->getPrefix(maxInputPrefixInV->size(), true)->removeEpsilon());

                LOG("VERBOSE_1") << "maxIOPrefixInV (v/v'): " << *maxIOPrefixInV << std::endl;
                
                
                
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
                        

                        bool doesExceedBound = exceedsBound(m, *maxIOPrefixInV, suffix, rDistStates, responseSets, responseMap, vDoublePrime, dReachableStates, *fsm);
                        LOG("VERBOSE_1") << "exceedsBound: " << doesExceedBound << std::endl;
                        if (doesExceedBound)
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
    
    LOG("VERBOSE_1") << "  RESULT: " << std::endl;
    for (const auto& t : observedTraces) {
        LOG("VERBOSE_1") << "\t\t" << t << std::endl;
    }
    LOG("VERBOSE_1") << "IUT is a reduction of the specification." << std::endl;

    cout << "IUT is a reduction of the specification." << std::endl;
    cout << "finished" << endl;
    
    exit(0);
    
}

