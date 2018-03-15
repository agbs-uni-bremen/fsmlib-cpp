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


bool executeAdaptiveTest(
        const string& prefix,
        const int numStates,
        const int numInput,
        const int numOutput,
        const size_t numOutFaults,
        const size_t numTransFaults,
        const unsigned int createRandomFsmSeed,
        const unsigned int createMutantSeed,
        const shared_ptr<FsmPresentationLayer>& pl,
        const bool dontTestReductions,
        bool& isReduction,
        const string iteration = "")
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


    shared_ptr<Fsm> iut = spec->createMutant("mutant" + prefix,
                                             numOutFaults,
                                             numTransFaults,
                                             true,
                                             createMutantSeed);

    CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
    CLOG(INFO, logging::globalLogger) << "numInput: " << numInput + 1;
    CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutput + 1;
    CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
    CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
    CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed: " << createRandomFsmSeed;
    CLOG(INFO, logging::globalLogger) << "createMutantSeed: " << createMutantSeed;

    Fsm specMin = spec->minimise(false);
    Fsm iutMin = iut->minimise(false);
    //Should not be in release:
    Fsm intersect = spec->intersect(*iut);

#ifdef ENABLE_DEBUG_MACRO
    const string dotPrefix = "../../../resources/adaptive-test/" + spec->getName() + "-";
    spec->toDot(dotPrefix + "spec");
    iut->toDot(dotPrefix + "iut");
    specMin.toDot(dotPrefix + "specMin");
    iutMin.toDot(dotPrefix + "iutMin");

    intersect.toDot(dotPrefix + "intersect");
#endif

    shared_ptr<IOTrace> failTrace;
    isReduction = !intersect.hasFailure(failTrace);

    if (failTrace)
    {
        failTrace = make_shared<IOTrace>(failTrace->removeLeadingEpsilons());
    }

    CLOG(INFO, logging::globalLogger) << "IUT is " + string((isReduction) ? "" : "NOT ") +
                                         "a reduction of the specification.";

    if (failTrace)
    {
        CLOG(INFO, logging::globalLogger) << "failTrace: " << *failTrace;
        CLOG(INFO, logging::globalLogger) << "failTrace length: " << failTrace->size();
    }

    if (isReduction && dontTestReductions) {
        CLOG(INFO, logging::globalLogger) << "Won't test this one, since it is a reduction.";
        throw unexpected_reduction("Interrupting testing, since IUT is an unexcpected reduction of the specification.");
    }
    else
    {
        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Testing.";
    }

    IOTraceContainer observedTraces;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    bool result = Fsm::adaptiveStateCounting(specMin, iutMin, static_cast<size_t>(iutMin.getMaxNodes()), observedTraces);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long durationMS = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    long durationMin = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();



    CLOG(INFO, logging::globalLogger) << "Calculation took " << durationMS << " ms (" << durationMin << " minutes).";
    LOG(INFO) << "observedTraces: " << observedTraces;


    std::stringstream csvOutput;
    csvOutput << iteration;
    csvOutput << "," << numStates + 1;
    csvOutput << "," << numInput + 1;
    csvOutput << "," << numOutput + 1;
    csvOutput << "," << specMin.getDReachableStates().size();
    csvOutput << "," << specMin.getMaximalSetsOfRDistinguishableStates().size();
    csvOutput << "," << numOutFaults;
    csvOutput << "," << numTransFaults;
    csvOutput << "," << isReduction;
    if (failTrace)
    {
        csvOutput << "," << failTrace->size();
    }
    else
    {
        csvOutput << ",-1";
    }
    csvOutput << "," << result;
    csvOutput << "," << durationMS;
    csvOutput << "," << (isReduction == result);
    csvOutput << "," << createRandomFsmSeed;
    csvOutput << "," << createMutantSeed;

    CLOG(INFO, logging::csvLogger) << csvOutput.str();

    return isReduction == result;
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
    bool dontTestReductions = false;
};

void adaptiveTest01(AdaptiveTestConfig& config)
{
    CLOG(INFO, logging::globalLogger) << "############## Adaptive Test 01 ##############";
    std::chrono::steady_clock::time_point totalStart = std::chrono::steady_clock::now();

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
    CLOG(INFO, logging::globalLogger) << "dontTestReductions: " << std::boolalpha << config.dontTestReductions;
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

#ifdef ENABLE_DEBUG_MACRO
            logging::setLogfileSuffix(iteration);
#endif

            int numInput = getRandom(config.minInput, config.maxInput, gen);
            int numOutput = getRandom(config.minOutput, config.maxOutput, gen);
            size_t numOutFaults = static_cast<size_t>(getRandom(config.minOutFaults, config.maxOutFaults, gen));
            size_t numTransFaults = static_cast<size_t>(getRandom(config.minTransFaults, config.maxTransFaults, gen));

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
            CLOG(INFO, logging::globalLogger) << "------------------------------------------------------------------";
            CLOG(INFO, logging::globalLogger) << "i: " << iteration;

            bool isReduction;
            bool result;

            bool couldCreateMutant = false;
            while (!couldCreateMutant)
            {
                try {
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
                                config.dontTestReductions,
                                isReduction,
                                iteration);
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
                    if (numTransFaults - 1 >= static_cast<size_t>(config.minTransFaults) && numTransFaults - 1 > 0) {
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
                    if (numOutFaults - 1 >= static_cast<size_t>(config.minOutFaults) && numOutFaults - 1 > 0) {
                        --numOutFaults;
                        CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Decreasing output faults.";
                        continue;
                    }
                    else if (numTransFaults - 1 >= static_cast<size_t>(config.minTransFaults) && numTransFaults - 1 > 0)
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

            if (result)
            {
                ++passed;
            }

            assertOnFail("TC-AT-01-" + iteration, result);

            ++executed;
            ++i;
        }
    }

    std::chrono::steady_clock::time_point totalEnd = std::chrono::steady_clock::now();
    long durationS = std::chrono::duration_cast<std::chrono::seconds>(totalEnd - totalStart).count();
    long durationMin = std::chrono::duration_cast<std::chrono::minutes>(totalEnd - totalStart).count();

    CLOG(INFO, logging::globalLogger) << "";
    CLOG(INFO, logging::globalLogger) << "#################### SUMMARY ####################";
    CLOG(INFO, logging::globalLogger) << "# Total tests  : " << executed;
    CLOG(INFO, logging::globalLogger) << "# Passed       : " << passed;
    CLOG(INFO, logging::globalLogger) << "# Failed       : " << executed - passed;
    CLOG(INFO, logging::globalLogger) << "# Not executed : " << i - executed;
    CLOG(INFO, logging::globalLogger) << "# Duration     : " << durationS << " s (" << durationMin << " min).";
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
#ifdef ENABLE_DEBUG_MACRO
    CLOG(INFO, logging::globalLogger) << "This is a debug build!";
#else
    CLOG(INFO, logging::globalLogger) << "This is a release build!";
#endif
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

    shared_ptr<IOTrace> failTrace;
    bool isReduction = !intersect.hasFailure(failTrace);

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
    config.numFsm = 1000000;

    config.minInput = 1;
    config.maxInput = 10;

    config.minOutput = 1;
    config.maxOutput = 10;

    config.minStates = 1;
    config.maxStates = 10;

    config.minTransFaults = 1;
    config.maxTransFaults = 1;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.dontTestReductions = true;

    config.seed = 1337;

    CLOG(INFO, logging::csvLogger) << "i,numStates,numInput,numOutput,numDReachable,numMaximalSetsOfRDistStates,numOutFaults,numTransFault,"
                                      "isReduction,failTraceSize,result,durationMS,pass,createRandomFsmSeed,createMutantSeed";

    bool debug = true;
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

        shared_ptr<FsmPresentationLayer> plTest =
        make_shared<FsmPresentationLayer>("../../../resources/adaptive-test-in.txt",
                                          "../../../resources/adaptive-test-out.txt",
                                          "../../../resources/adaptive-test-state.txt");
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
                    false,
                    isReduction);
        assertOnFail("TC-AT-DEBUG", result);
    }
    else
    {
        adaptiveTest01(config);
    }

	cout << endl << endl;
}



