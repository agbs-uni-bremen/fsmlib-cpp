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

static const string ascTestDirectory = "../../../resources/asc-tests/";
static const string ascTestResultDirectory = "../../../resources/asc-test-results/";

static shared_ptr<FsmPresentationLayer> plTestSpec = make_shared<FsmPresentationLayer>(
            ascTestDirectory + "adaptive-test-in.txt",
            ascTestDirectory + "adaptive-test-out.txt",
            ascTestDirectory + "adaptive-test-state-spec.txt");

static shared_ptr<FsmPresentationLayer> plTestIut = make_shared<FsmPresentationLayer>(
            ascTestDirectory + "adaptive-test-in.txt",
            ascTestDirectory + "adaptive-test-out.txt",
            ascTestDirectory + "adaptive-test-state-iut.txt");


static const string testSepLine = "------------------------------------------------------------------";


static const vector<string> csvHeaders = {
   "testName",
    "numStates",
    "numInputs",
    "numOutputs",
    "numDReachableStates",
    "numSetsOfMaximalRDistStates",
    "numOutFaults",
    "numTransFaults",
    "iutIsReduction",
    "failTraceFound",
    "failTraceFoundSize",
    "observedTracesSize",
    "longestObservedTrace",
    "adaptiveStateCountingResult",
    "createRandomFsmSeed",
    "createMutantSeed",
    "durationMS",
    "durationM",
    "pass"
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

    // Optional
    int minStates = 2;
    int minInput = 2;
    int minOutput = 2;
    int minOutFaults = 1;
    int minTransFaults = 1;
    unsigned int seed = 0;
    bool dontTestReductions = false;
};

struct AdaptiveTestResult
{
    string testName;
    int numStates = -2;
    int numInputs = -2;
    int numOutputs = -2;
    int numDReachableStates = -1;
    int numSetsOfMaximalRDistStates = -1;
    int numOutFaults = -1;
    int numTransFaults = -1;
    bool iutIsReduction;
    shared_ptr<IOTrace> failTraceFound;
    IOTraceContainer observedTraces;
    shared_ptr<IOTrace> longestObservedTrace;
    bool adaptiveStateCountingResult;
    unsigned createRandomFsmSeed = 0;
    unsigned createMutantSeed = 0;
    long durationMS = -1;
    long durationM = -1;
    bool pass = 0;
    shared_ptr<Fsm> intersection;
};

void assertInconclusive(string tc, string comment = "") {
    
    string sVerdict("INCONCLUSIVE");
    CLOG(INFO, logging::globalLogger) << sVerdict << ": " << tc << " : " << comment <<  endl;    
}

void assert(string tc, bool verdict, string comment = "") {
    
    string out = (verdict) ? "PASS" : "FAIL";
    out += ": " + tc;
    if (!comment.empty())
    {
        out += ": " + comment;
    }

    CLOG(INFO, logging::globalLogger) << out <<  endl;
    
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

void writeCsvHeader()
{
    string header = "";
    size_t size = csvHeaders.size() - 1;
    for (size_t i = 0; i <= size; ++i)
    {
        header += csvHeaders.at(i);
        if (i != size)
        {
            header += ",";
        }
    }
    CLOG(INFO, logging::csvLogger) << header;
}

void printTestBegin(string name)
{
    CLOG(INFO, logging::globalLogger) << "#################### Test " << name << "####################";
}

void printSummary(const int& executed,
                  const int& passed,
                  const int& notExecuted,
                  const long& durationS,
                  const long& durationM)
{
    CLOG(INFO, logging::globalLogger) << "";
    CLOG(INFO, logging::globalLogger) << "#################### SUMMARY ####################";
    CLOG(INFO, logging::globalLogger) << "# Total tests  : " << executed;
    CLOG(INFO, logging::globalLogger) << "# Passed       : " << passed;
    CLOG(INFO, logging::globalLogger) << "# Failed       : " << executed - passed;
    CLOG(INFO, logging::globalLogger) << "# Not executed : " << notExecuted;
    CLOG(INFO, logging::globalLogger) << "# Duration     : " << durationS << " s (" << durationM << " min).";
    CLOG(INFO, logging::globalLogger) << "#################################################";
}

void printTestConfig(const AdaptiveTestConfig& config, const int& diffStates, const int& subLoopIterations)
{
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
    CLOG(INFO, logging::globalLogger) << "dontTestReductions: " << std::boolalpha << config.dontTestReductions;
    CLOG(INFO, logging::globalLogger) << "seed: " << config.seed;

    CLOG(INFO, logging::globalLogger) << "diffStates: " << diffStates;
    CLOG(INFO, logging::globalLogger) << "subLoopIterations: " << subLoopIterations;
}

void printTestResult(AdaptiveTestResult& result, bool log, bool csv, bool printTraces)
{
    if (log)
    {
        CLOG(INFO, logging::globalLogger) << testSepLine;
        CLOG(INFO, logging::globalLogger) << "Test                       : " << result.testName;
        CLOG(INFO, logging::globalLogger) << "numStates                  : " << result.numStates;
        CLOG(INFO, logging::globalLogger) << "numInputs                  : " << result.numInputs;
        CLOG(INFO, logging::globalLogger) << "numOutputs                 : " << result.numOutputs;
        CLOG(INFO, logging::globalLogger) << "numDReachableStates        : " << result.numDReachableStates;
        CLOG(INFO, logging::globalLogger) << "numSetsOfMaximalRDistStates: " << result.numSetsOfMaximalRDistStates;
        CLOG(INFO, logging::globalLogger) << "numOutFaults               : " << result.numOutFaults;
        CLOG(INFO, logging::globalLogger) << "numTransFaults             : " << result.numTransFaults;
        CLOG(INFO, logging::globalLogger) << "iutIsReduction             : " << std::boolalpha << result.iutIsReduction;
        if (result.failTraceFound)
        {
            CLOG(INFO, logging::globalLogger) << "failTraceFound             : " << *result.failTraceFound;
        }
        else
        {
            CLOG(INFO, logging::globalLogger) << "failTraceFound             : None";
        }
        CLOG(INFO, logging::globalLogger) << "adaptiveStateCountingResult: " << std::boolalpha
                                          << result.adaptiveStateCountingResult;
        if (result.createRandomFsmSeed != 0)
        {
            CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed        : " << result.createRandomFsmSeed;
        }
        if (result.createMutantSeed != 0)
        {
            CLOG(INFO, logging::globalLogger) << "createMutantSeed           : " << result.createMutantSeed;
        }
        CLOG(INFO, logging::globalLogger) << "longestObservedTrace       : " << *result.longestObservedTrace;
        CLOG(INFO, logging::globalLogger) << "length longestObservedTrace: " << result.longestObservedTrace->size();
        CLOG(INFO, logging::globalLogger) << "observedTraces size        : " << result.observedTraces.size();
        if (printTraces)
        {
            CLOG(INFO, logging::globalLogger) << "observedTraces             : " << result.observedTraces;
        }
        CLOG(INFO, logging::globalLogger) << "Calculation took " << result.durationMS << " ms ("
                                          << result.durationM << " minutes).";
    }

    if (csv)
    {
        std::stringstream csvOutput;
        csvOutput << result.testName;
        csvOutput << "," << result.numStates;
        csvOutput << "," << result.numInputs;
        csvOutput << "," << result.numOutputs;
        csvOutput << "," << result.numDReachableStates;
        csvOutput << "," << result.numSetsOfMaximalRDistStates;
        csvOutput << "," << result.numOutFaults;
        csvOutput << "," << result.numTransFaults;
        csvOutput << "," << result.iutIsReduction;
        csvOutput << "," << result.adaptiveStateCountingResult;
        csvOutput << "," << result.durationMS;
        csvOutput << "," << result.pass;
        csvOutput << "," << result.createRandomFsmSeed;
        csvOutput << "," << result.createMutantSeed;
        CLOG(INFO, logging::csvLogger) << csvOutput.str();
    }
}

bool isReduction(Fsm& spec, Fsm& iut, string intersectionName, shared_ptr<Fsm>& intersection)
{
    Fsm inter = spec.intersect(iut, intersectionName);
    intersection = make_shared<Fsm>(inter);
    return !inter.hasFailure();
}

void executeAdaptiveTest(Fsm& spec, Fsm& iut, size_t m, string intersectionName, const bool& toDot, const bool& dontTestReductions, AdaptiveTestResult& result)
{
    Fsm specMin = spec.minimise(false, "", "", false);
    Fsm iutMin = iut.minimise(false, "", "", false);

    result.numStates = specMin.getMaxNodes();
    result.numInputs = specMin.getMaxInput() + 1;
    result.numOutputs = specMin.getMaxOutput() + 1;

    result.iutIsReduction = isReduction(specMin, iutMin, intersectionName, result.intersection);


    if (toDot)
    {
        spec.toDot(ascTestResultDirectory + spec.getName());
        iut.toDot(ascTestResultDirectory + iut.getName());
        specMin.toDot(ascTestResultDirectory + spec.getName() + "-min");
        iutMin.toDot(ascTestResultDirectory + iut.getName() + "-min");
        result.intersection->toDot(ascTestResultDirectory + result.intersection->getName());
    }

    if (result.iutIsReduction && dontTestReductions) {
        CLOG(INFO, logging::globalLogger) << "Won't test this one, since it is a reduction.";
        throw unexpected_reduction("Interrupting testing, since IUT is an unexcpected reduction of the specification.");
    }

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    result.adaptiveStateCountingResult = Fsm::adaptiveStateCounting(specMin, iutMin, m, result.observedTraces, result.failTraceFound);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    result.durationMS = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    result.durationM = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();

    result.pass = (result.iutIsReduction == result.adaptiveStateCountingResult);
    result.numSetsOfMaximalRDistStates = static_cast<int>(specMin.getMaximalSetsOfRDistinguishableStates().size());
    result.numDReachableStates = static_cast<int>(specMin.getDReachableStates().size());
    result.longestObservedTrace = make_shared<IOTrace>(result.observedTraces.getLongestTrace());
}


void createAndExecuteAdaptiveTest(
        const string& prefix,
        const int numStates,
        const int numInputs,
        const int numOutputs,
        const int numOutFaults,
        const int numTransFaults,
        const unsigned int createRandomFsmSeed,
        const unsigned int createMutantSeed,
        const shared_ptr<FsmPresentationLayer>& plSpec,
        const shared_ptr<FsmPresentationLayer>& plIut,
        const bool dontTestReductions,
        AdaptiveTestResult& result)
{

    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating FSM.";
    shared_ptr<Fsm> spec = Fsm::createRandomFsm(prefix + "-spec",
                                                numInputs,
                                                numOutputs,
                                                numStates,
                                                plSpec,
                                                true,
                                                createRandomFsmSeed);
    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating mutant.";


    shared_ptr<Fsm> iut = spec->createMutant(prefix + "-iut",
                                             numOutFaults,
                                             numTransFaults,
                                             true,
                                             createMutantSeed,
                                             plIut);

    result.numStates = numStates + 1;
    result.numInputs = numInputs + 1;
    result.numOutputs = numOutputs + 1;
    result.numOutFaults = numOutFaults;
    result.numTransFaults = numTransFaults;
    result.createRandomFsmSeed = createRandomFsmSeed;
    result.createMutantSeed = createMutantSeed;


    executeAdaptiveTest(*spec, *iut, static_cast<size_t>(iut->getMaxNodes()),
                        prefix + "-intersect", false, dontTestReductions, result);
    printTestResult(result, true, true, false);
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
    printTestBegin(config.testName);
    std::chrono::steady_clock::time_point totalStart = std::chrono::steady_clock::now();

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

    printTestConfig(config, diffStates, subLoopIterations);

    int executed = 0;
    int passed = 0;
    int i = 0;
    for (int numStates = config.minStates; numStates <= config.maxStates; ++numStates)
    {
        for (int j = 0; j < subLoopIterations; ++j)
        {

            shared_ptr<FsmPresentationLayer> plTestSpecCopy = make_shared<FsmPresentationLayer>(*plTestSpec);
            shared_ptr<FsmPresentationLayer> plTestIutCopy = make_shared<FsmPresentationLayer>(*plTestIut);

            stringstream ss;
            ss << setw(numberDigits) << setfill('0') << i;
            string iteration = ss.str();

#ifdef ENABLE_DEBUG_MACRO
            logging::setLogfileSuffix(iteration);
#endif

            int numInput = getRandom(config.minInput, config.maxInput, gen);
            int numOutput = getRandom(config.minOutput, config.maxOutput, gen);
            int numOutFaults = getRandom(config.minOutFaults, config.maxOutFaults, gen);
            int numTransFaults = getRandom(config.minTransFaults, config.maxTransFaults, gen);

            if (numOutput <= 1 && numOutFaults > 0)
            {
                if (config.maxOutput < 2) {
                    CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
                    CLOG(INFO, logging::globalLogger) << "numInput: " << numInput + 1;
                    CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutput + 1;
                    CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
                    CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
                    CLOG(WARNING, logging::globalLogger) << "Too little outputs. Can not create requested umber of " <<
                                                            "output faults. Could not create mutant. Skipping.";
                    continue;
                }
                else
                {
                    numOutput = getRandom(2, config.maxOutput, gen);
                }
            }

            const unsigned int createRandomFsmSeed = static_cast<unsigned int>(getRandom(gen));
            const unsigned int createMutantSeed = static_cast<unsigned int>(getRandom(gen));

            TIMED_SCOPE(timerBlkObj, "heavy-iter");
            AdaptiveTestResult result;
            result.testName = config.testName + "-" + iteration;

            bool couldCreateMutant = false;
            while (!couldCreateMutant)
            {
                try {
                    createAndExecuteAdaptiveTest(
                                iteration,
                                numStates,
                                numInput,
                                numOutput,
                                numOutFaults,
                                numTransFaults,
                                createRandomFsmSeed,
                                createMutantSeed,
                                plTestSpecCopy,
                                plTestIutCopy,
                                config.dontTestReductions,
                                result);
                    couldCreateMutant = true;
                }
                catch (unexpected_reduction& e)
                {
                    CLOG(WARNING, logging::globalLogger) << "IUT is a reduction of the specification. Skipping.";
                    break;
                }
                catch (too_many_transition_faults& e)
                {
                    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create mutant.";
                    if (numTransFaults - 1 >= config.minTransFaults && numTransFaults - 1 > 0) {
                        --numTransFaults;
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing transition faults.";
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
                catch (too_many_output_faults& e)
                {
                    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Could not create mutant.";
                    if (numOutFaults - 1 >= config.minOutFaults && numOutFaults - 1 > 0) {
                        --numOutFaults;
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing output faults.";
                        continue;
                    }
                    else if (numTransFaults - 1 >= config.minTransFaults && numTransFaults - 1 > 0)
                    {
                        --numTransFaults;
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing transition faults.";
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

            if (result.pass)
            {
                ++passed;
            }

            assertOnFail(result.testName, result.pass);
            ++executed;
            ++i;
        }
    }

    std::chrono::steady_clock::time_point totalEnd = std::chrono::steady_clock::now();
    long durationS = std::chrono::duration_cast<std::chrono::seconds>(totalEnd - totalStart).count();
    long durationM = std::chrono::duration_cast<std::chrono::minutes>(totalEnd - totalStart).count();

    printSummary(executed, passed, i - executed, durationS, durationM);
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

void test00_00()
{
    AdaptiveTestResult result;
    result.testName = "00-00";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/00-spec.dot", "00-00-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/00-iut.dot", "00-00-iut");
    result.numOutFaults = 0;
    result.numTransFaults = 0;
    executeAdaptiveTest(spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-00-inter", true, false, result);
    printTestResult(result, true, false, true);
    assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

void test00_01()
{
    AdaptiveTestResult result;
    result.testName = "00-01";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/01-spec.dot", "00-01-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/01-iut.dot", "00-01-iut");

    result.numOutFaults = 1;
    result.numTransFaults = 0;
    executeAdaptiveTest(spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-01-inter", true, false, result);
    printTestResult(result, true, false, true);
    assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

void runAdaptiveStateCountingTests()
{
    test00_00();
    test00_01();
}


int main(int argc, char* argv[])
{
    START_EASYLOGGINGPP(argc, argv);
    logging::initLogging();

    CLOG(INFO, logging::globalLogger) << "############## Starting Application ##############";
#ifdef ENABLE_DEBUG_MACRO
    CLOG(INFO, logging::globalLogger) << "This is a debug build!";
#else
    CLOG(INFO, logging::globalLogger) << "This is a release build!";
#endif

    /*
    runAdaptiveStateCountingTests();
    return 0;

    */

/*
    shared_ptr<FsmPresentationLayer> plSpec =
    make_shared<FsmPresentationLayer>("../../../resources/example-master-fehler.in",
                                      "../../../resources/example-master-fehler.out",
                                      "../../../resources/example-master-fehler.state");

    shared_ptr<FsmPresentationLayer> plIut =
    make_shared<FsmPresentationLayer>("../../../resources/example-master-fehler.in",
                                      "../../../resources/example-master-fehler.out",
                                      "../../../resources/example-master-fehler-iut.state");


    Fsm spec = Fsm("../../../resources/example-master-fehler.fsm", plSpec, "spec");
    Fsm iut = Fsm("../../../resources/example-master-fehler-iut.fsm", plIut, "iut");

    if (!spec.isCompletelyDefined()) {
        LOG(FATAL) << "Specification is not completely defined.";
    }
    if (!iut.isCompletelyDefined()) {
        LOG(FATAL) << "IUT is not completely defined.";
    }

    Fsm specMin = spec.minimise();
    Fsm iutMin = iut.minimise();

    const string dotPrefix = "../../../resources/example-master-fehler/";
    spec.toDot(dotPrefix + "spec");
    iut.toDot(dotPrefix + "iut");
    specMin.toDot(dotPrefix + "specMin");
    iutMin.toDot(dotPrefix + "iutMin");

    Fsm intersect = spec.intersect(iut);
    intersect.toDot(dotPrefix + "intersect");

    bool isReduction = !intersect.hasFailure();

    IOTraceContainer observedTraces;
    bool result = Fsm::adaptiveStateCounting(specMin, iutMin, static_cast<size_t>(iutMin.getMaxNodes()), observedTraces);

    CLOG(INFO, logging::globalLogger) << "isReduction: " << isReduction;
    if (failTrace)
    {
        CLOG(INFO, logging::globalLogger) << "failTrace: " << *failTrace;
    }
    CLOG(INFO, logging::globalLogger) << "result: " << result;

    return 0;
*/
    //TODO Analyze increasing memory usage with valgrind
    AdaptiveTestConfig config;
    config.testName = "TRIAL";
    config.numFsm = 10;

    config.minInput = 3;
    config.maxInput = 3;

    config.minOutput = 3;
    config.maxOutput = 3;

    config.minStates = 20;
    config.maxStates = 23;

    config.minTransFaults = 1;
    config.maxTransFaults = 1;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = true;

    config.seed = 1337;

    writeCsvHeader();

    bool debug = false;
    if (debug)
    {
        CLOG(INFO, logging::globalLogger) << "############## Debugging ##############";

        AdaptiveTestConfigDebug debugConfig;
        debugConfig.numStates = 4;
        debugConfig.numInput = 5;
        debugConfig.numOutput = 2;
        debugConfig.numOutFaults = 0;
        debugConfig.numTransFaults = 1;
        debugConfig.createRandomFsmSeed = 653887049;
        debugConfig.createMutantSeed = 1588831269;

        shared_ptr<FsmPresentationLayer> plTestSpecCopy = make_shared<FsmPresentationLayer>(*plTestSpec);
        shared_ptr<FsmPresentationLayer> plTestIutCopy = make_shared<FsmPresentationLayer>(*plTestIut);

        AdaptiveTestResult result;
        createAndExecuteAdaptiveTest(
                    debugConfig.prefix,
                    debugConfig.numStates - 1,
                    debugConfig.numInput - 1,
                    debugConfig.numOutput - 1,
                    debugConfig.numOutFaults,
                    debugConfig.numTransFaults,
                    debugConfig.createRandomFsmSeed,
                    debugConfig.createMutantSeed,
                    plTestSpecCopy,
                    plTestIutCopy,
                    false,
                    result);
        assertOnFail(debugConfig.prefix, result.pass);
    }
    else
    {
        adaptiveTestRandom(config);
    }

	cout << endl << endl;
}



