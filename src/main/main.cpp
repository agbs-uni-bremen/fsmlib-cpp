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

using namespace std;

static const string csvSep = ";";

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


enum class TestIteration
{
    INPUT,
    OUTPUT,
    STATE,
    OUTPUT_FAULT,
    TRANSITION_FAULT,
    END
};

enum class CsvField
{
    TEST_NAME,
    NUM_STATES,
    NUM_INPUTS,
    NUM_D_REACHABLE_STATES,
    NUM_SETS_OF_MAXIMAL_R_DIST_STATES,
    NUM_OUT_FAULTS,
    NUM_TRANS_FAULTS,
    DEGREE_OF_COMPLETENESS,
    DEGREE_OF_NON_DETERMINISM,
    IUT_IS_REDUCTION,
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

static map<CsvField, string> csvHeaders;
static map<TestIteration, string> csvIterations;

void initCsvHeaders()
{
    csvHeaders.insert(make_pair(CsvField::TEST_NAME, "testName"));
    csvHeaders.insert(make_pair(CsvField::NUM_STATES, "numStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_INPUTS, "numInputs"));
    csvHeaders.insert(make_pair(CsvField::NUM_D_REACHABLE_STATES, "numDReachableStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_SETS_OF_MAXIMAL_R_DIST_STATES, "numSetsOfMaximalRDistStates"));
    csvHeaders.insert(make_pair(CsvField::NUM_OUT_FAULTS, "numOutFaults"));
    csvHeaders.insert(make_pair(CsvField::NUM_TRANS_FAULTS, "numTransFaults"));
    csvHeaders.insert(make_pair(CsvField::DEGREE_OF_COMPLETENESS, "degreeOfCompleteness"));
    csvHeaders.insert(make_pair(CsvField::DEGREE_OF_NON_DETERMINISM, "degreeOfNonDeterminism"));
    csvHeaders.insert(make_pair(CsvField::IUT_IS_REDUCTION, "iutIsReduction"));
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
}

void initialize()
{
    initCsvHeaders();
    initCsvIterations();
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
    bool logTestDetails = false;
    bool logTestResultsSummary = true;
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

    float degreeOfCompleteness = -1;
    float maxDegreeOfNonDeterminism = -1;

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
    int numOutFaults = -1;
    int numTransFaults = -1;
    float degreeOfCompleteness = -1;
    float degreeOfNonDeterminism = -1;
    bool iutIsReduction;
    shared_ptr<IOTrace> failTraceFound;
    IOTraceContainer observedTraces;
    shared_ptr<IOTrace> longestObservedTrace;
    bool adaptiveStateCountingResult;
    unsigned createRandomFsmSeed = 0;
    unsigned createMutantSeed = 0;
    int iterations = -1;
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
    case CsvField::NUM_D_REACHABLE_STATES:
        out << result.numDReachableStates;
        break;
    case CsvField::NUM_SETS_OF_MAXIMAL_R_DIST_STATES:
        out << result.setsOfMaximalRDistStates.size();
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
        out << result.observedTraces.size();
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

void writeCsvHeader(const CsvConfig& config, const string& logger)
{
    string header = "";
    if (logger == logging::csvLoggerEveryIteration)
    {
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
        CLOG(INFO, logging::csvLoggerEveryIteration) << header;
    }
}

void logToCsv(const AdaptiveTestResult& result, const CsvConfig& config)
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

    CLOG(INFO, logging::csvLoggerEveryIteration) << output;
}


void logToCsv(const vector<vector<AdaptiveTestResult>>& resultsMatrix, const CsvConfig& config)
{
    if (config.context == TestIteration::END || config.fieldsContext.size() == 0)
    {
        return;
    }

    string contextHeader = csvIterations.at(config.context) + csvSep;
    string fieldHeader = csvSep;

    int numStates = 0;


    size_t matrixHeight = 0;
    for (size_t i = 0; i < resultsMatrix.size(); ++i)
    {
        size_t h = resultsMatrix.at(i).size();
        if (i == 0)
        {
            matrixHeight = h;
        }
        else if (matrixHeight != h)
        {
            CLOG(WARNING, logging::globalLogger) << "The results matrix is not of consitent height.";
            if (h > matrixHeight)
            {
                matrixHeight = h;
            }
        }
    }

    vector<string> out(matrixHeight, "");

    for (size_t fieldIdx = 0; fieldIdx < config.fieldsContext.size(); ++fieldIdx)
    {
        const CsvField& field = config.fieldsContext.at(fieldIdx);
        const string& fieldName = csvHeaders.at(field);

        for (size_t col = 0; col < resultsMatrix.size(); ++col)
        {
            for (size_t row = 0; row < resultsMatrix.at(col).size(); ++row)
            {
                const AdaptiveTestResult& result = resultsMatrix.at(col).at(row);
                if (row == 0)
                {
                    numStates = result.numStates;
                }
                else if (result.numStates != numStates)
                {
                    CLOG(FATAL, logging::globalLogger)
                            << "Trying to log context, but number of states differ.";
                }
                out.at(row) += getFieldFromResult(result, field) + csvSep;
            }
            contextHeader += to_string(numStates) + csvSep;
            fieldHeader += fieldName + csvSep;
        }
        for (size_t l = 0; l < out.size(); ++l)
        {
            out.at(l) += csvSep;
        }
        contextHeader += csvSep;
        fieldHeader += csvSep;
    }


    CLOG(INFO, logging::csvLoggerContext) << fieldHeader;
    CLOG(INFO, logging::csvLoggerContext) << contextHeader;
    for (const string& o : out)
    {
        CLOG(INFO, logging::csvLoggerContext) << csvSep << o;
    }
}

/*
void logAveragesToCsv(const string& testName,
                      const vector<AdaptiveTestResult>& results,
                      const CsvConfig& csvConfig)
{
    const vector<CsvField>& fields = csvConfig.fields;

    AdaptiveTestResult averageResult;
    averageResult.testName = testName;

    // Supported fields
    float averageDurationMS = 0;

    int numStates = 0;

    for (size_t i = 0; i < results.size(); ++i)
    {
        const AdaptiveTestResult& result = results.at(i);
        if (i == 0)
        {
            numStates = result.numStates;
            averageResult.numStates = numStates;
        }
        else
        {
            // Check if average is requested on a field that makes no sense.
            if (result.numStates != numStates)
            {
                CLOG(FATAL, logging::globalLogger)
                        << "Trying to calculate an average, but number of states differ.";
            }
        }

        if (find(fields.begin(), fields.end(), DURATION_MS) != fields.end())
        {
            averageDurationMS += (result.durationMS - averageDurationMS) / (i + 1);
        }
    }

    averageResult.durationMS = static_cast<long>(averageDurationMS);
    averageResult.avgCount = static_cast<int>(results.size());

    logToCsv(averageResult, csvConfig);
}
*/


void newCsvFile(CsvConfig csvConfig, const string& suffix)
{
    CLOG_IF(VLOG_IS_ON(1), INFO, logging::globalLogger) << "newCsvFile()";
    CLOG_IF(VLOG_IS_ON(1), INFO, logging::globalLogger) << "  suffix: " << suffix;

    if (csvConfig.logEveryIteration)
    {
        logging::setLogfileSuffix(suffix, logging::csvLoggerEveryIteration);
        writeCsvHeader(csvConfig, logging::csvLoggerEveryIteration);
    }
    if (csvConfig.context != TestIteration::END)
    {
        logging::setLogfileSuffix(suffix, logging::csvLoggerContext);
        writeCsvHeader(csvConfig, logging::csvLoggerContext);
    }
}


void printTestBegin(string name)
{
    CLOG(INFO, logging::globalLogger) << "#################### Test " << name << " ####################";
}

void printSummary(const string& testName,
                  const int& executed,
                  const int& passed,
                  const int& notExecuted,
                  const long& durationS,
                  const long& durationM)
{
    CLOG(INFO, logging::globalLogger) << "";
    CLOG(INFO, logging::globalLogger) << "#################### SUMMARY ####################";
    CLOG(INFO, logging::globalLogger) << "# Test name    : " << testName;
    CLOG(INFO, logging::globalLogger) << "# Total tests  : " << executed;
    CLOG(INFO, logging::globalLogger) << "# Passed       : " << passed;
    CLOG(INFO, logging::globalLogger) << "# Failed       : " << executed - passed;
    CLOG(INFO, logging::globalLogger) << "# Not executed : " << notExecuted;
    CLOG(INFO, logging::globalLogger) << "# Duration     : " << durationS << " s (" << durationM << " min).";
    CLOG(INFO, logging::globalLogger) << "#################################################";
}

void printTestConfig(const AdaptiveTestConfig& config)
{
    CLOG(INFO, logging::globalLogger) << "#################### Test Config ####################";
    CLOG(INFO, logging::globalLogger) << "numFsm: " << config.numFsm;
    CLOG(INFO, logging::globalLogger) << "minInput: " << config.minInput + 1;
    CLOG(INFO, logging::globalLogger) << "maxInput: " << config.maxInput + 1;
    CLOG(INFO, logging::globalLogger) << "minOutput: " << config.minOutput + 1;
    CLOG(INFO, logging::globalLogger) << "maxOutput: " << config.maxOutput + 1;
    CLOG(INFO, logging::globalLogger) << "minStates: " << config.minStates + 1;
    CLOG(INFO, logging::globalLogger) << "maxStates: " << config.maxStates + 1;

    CLOG(INFO, logging::globalLogger) << "minOutFaults: " << config.minOutFaults;
    CLOG(INFO, logging::globalLogger) << "maxOutFaults: " << config.maxOutFaults;
    CLOG(INFO, logging::globalLogger) << "minTransFaults: " << config.minTransFaults;
    CLOG(INFO, logging::globalLogger) << "maxTransFaults: " << config.maxTransFaults;

    CLOG(INFO, logging::globalLogger) << "degreeOfCompleteness: " << config.degreeOfCompleteness;
    CLOG(INFO, logging::globalLogger) << "maxDegreeOfNonDeterminism: " << config.maxDegreeOfNonDeterminism;

    CLOG(INFO, logging::globalLogger) << "dontTestReductions: " << std::boolalpha << config.dontTestReductions;
    CLOG(INFO, logging::globalLogger) << "forceTestParameters: " << std::boolalpha << config.forceTestParameters;
    CLOG(INFO, logging::globalLogger) << "seed: " << config.seed;
    CLOG(INFO, logging::globalLogger) << "#####################################################";
}

void printTestResult(AdaptiveTestResult& result, const CsvConfig& csvConfig, const LoggingConfig& loggingConfig)
{
    if (loggingConfig.logTestDetails)
    {
        CLOG(INFO, logging::globalLogger) << testSepLine;
        CLOG(INFO, logging::globalLogger) << "Test                       : " << result.testName;
        CLOG(INFO, logging::globalLogger) << "numStates                  : " << result.numStates;
        CLOG(INFO, logging::globalLogger) << "numInputs                  : " << result.numInputs;
        CLOG(INFO, logging::globalLogger) << "numOutputs                 : " << result.numOutputs;
        CLOG(INFO, logging::globalLogger) << "numDReachableStates        : " << result.numDReachableStates;


        if (loggingConfig.printSetsOfMaximalRDistStates)
        {
            CLOG(INFO, logging::globalLogger) << "setsOfMaximalRDistStates   : ";
            for (auto v : result.setsOfMaximalRDistStates)
            {
                stringstream ss;
                ss << "{";
                for (auto e : v)
                {
                    ss << e->getName() << ", ";
                }
                ss << "}";
                CLOG(INFO, logging::globalLogger) << ss.str();
            }
        }

        CLOG(INFO, logging::globalLogger) << "numSetsOfMaximalRDistStates: " << result.setsOfMaximalRDistStates.size();
        CLOG(INFO, logging::globalLogger) << "numOutFaults               : " << result.numOutFaults;
        CLOG(INFO, logging::globalLogger) << "numTransFaults             : " << result.numTransFaults;
        CLOG(INFO, logging::globalLogger) << "degreeOfCompleteness       : " << result.degreeOfCompleteness;
        CLOG(INFO, logging::globalLogger) << "degreeOfNonDeterminism     : " << result.degreeOfNonDeterminism;
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
        if (loggingConfig.printObservedTraces)
        {
            CLOG(INFO, logging::globalLogger) << "observedTraces             : " << result.observedTraces;
        }
        CLOG(INFO, logging::globalLogger) << "iterations                 : " << result.iterations;
        CLOG(INFO, logging::globalLogger) << "Calculation took " << result.durationMS << " ms ("
                                          << result.durationM << " minutes).";
    }

    if (csvConfig.logEveryIteration)
    {
        logToCsv(result, csvConfig);
    }
}

bool isReduction(Fsm& spec, Fsm& iut, string intersectionName, shared_ptr<Fsm>& intersection)
{
    Fsm inter = spec.intersect(iut, intersectionName);
    intersection = make_shared<Fsm>(inter);
    return !inter.hasFailure();
}

void executeAdaptiveTest(Fsm& spec, Fsm& iut, size_t m, string intersectionName,
                         const bool& toDot, const bool& dontTestReductions, AdaptiveTestResult& result)
{
    Fsm specMin = spec.minimise(false, "", "", false);
    Fsm iutMin = iut.minimise(false, "", "", false);

    result.degreeOfCompleteness = specMin.getDegreeOfCompleteness();
    result.degreeOfNonDeterminism = specMin.getDegreeOfNonDeterminism();

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
    result.numDReachableStates = static_cast<int>(specMin.getDReachableStates().size());
    result.longestObservedTrace = result.observedTraces.getLongestTrace();
}


void createAndExecuteAdaptiveTest(
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
        const bool dontTestReductions,
        const CsvConfig& csvConfig,
        const LoggingConfig loggingConfig,
        AdaptiveTestResult& result)
{

    CLOG_IF(VLOG_IS_ON(2), INFO, logging::globalLogger) << "Creating FSM.";
    shared_ptr<Fsm> spec = Fsm::createRandomFsm(prefix + "-spec",
                                                numInputs,
                                                numOutputs,
                                                numStates,
                                                plSpec,
                                                degreeOfCompleteness,
                                                maxDegreeOfNonDeterminism,
                                                true,
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

    if (loggingConfig.logTestDetails)
    {
        CLOG(INFO, logging::globalLogger) << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
        CLOG(INFO, logging::globalLogger) << "Starting adaptive state counting.";
        CLOG(INFO, logging::globalLogger) << "Name                     : " << result.testName;
        CLOG(INFO, logging::globalLogger) << "numStates                : " << result.numStates;
        CLOG(INFO, logging::globalLogger) << "numInputs                : " << result.numInputs;
        CLOG(INFO, logging::globalLogger) << "numOutputs               : " << result.numOutputs;
        CLOG(INFO, logging::globalLogger) << "numOutFaults             : " << result.numOutFaults;
        CLOG(INFO, logging::globalLogger) << "numTransFaults           : " << result.numTransFaults;
        CLOG(INFO, logging::globalLogger) << "degreeOfCompleteness     : " << degreeOfCompleteness;
        CLOG(INFO, logging::globalLogger) << "maxDegreeOfNonDeterminism: " << maxDegreeOfNonDeterminism;
        CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed      : " << result.createRandomFsmSeed;
        CLOG(INFO, logging::globalLogger) << "createMutantSeed         : " << result.createMutantSeed;
        CLOG(INFO, logging::globalLogger) << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
    }

    executeAdaptiveTest(*spec, *iut, static_cast<size_t>(iut->getMaxNodes()),
                        prefix + "-intersect", loggingConfig.toDot, dontTestReductions, result);

    printTestResult(result, csvConfig, loggingConfig);
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
//TODO Eine CSV-Datei für jede Test-Methode erstellen (ohne den Logger?).
void adaptiveTestRandom(AdaptiveTestConfig& config)
{
    if (config.loggingConfig.logTestDetails)
    {
        printTestBegin(config.testName);
    }
    newCsvFile(config.csvConfig, config.testName);

    std::chrono::steady_clock::time_point totalStart = std::chrono::steady_clock::now();

    if (config.numFsm < 0 ||
            config.maxInput < 0 ||
            config.maxOutput < 0 ||
            config.maxStates < 0 ||
            config.maxOutFaults < 0 ||
            config.maxTransFaults < 0 ||
            config.degreeOfCompleteness < 0 ||
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

    CLOG(INFO, logging::globalLogger) << "Seed: " << config.seed;

    std::mt19937 gen(config.seed);


    const int numberDigits = ((config.numFsm <= 1)? 1 : static_cast<int>(log10(config.numFsm)) + 1);


    const int diffInput = config.maxInput - config.minInput + 1;
    const int diffOutput = config.maxOutput - config.minOutput + 1;
    const int diffStates = config.maxStates - config.minStates + 1;
    const int diffOutFaults = config.maxOutFaults - config.minOutFaults + 1;
    const int diffTransFaults = config.maxTransFaults - config.minTransFaults + 1;

    if (diffInput <= 0 || diffOutput <= 0 || diffStates <= 0 || diffOutFaults <= 0 || diffTransFaults <= 0)
    {
        CLOG(FATAL, logging::globalLogger) << "Please check the test parameters.";
    }

    if (config.loggingConfig.logTestDetails)
    {
        printTestConfig(config);
    }

    float divisor = (diffInput * diffOutput * diffStates * diffOutFaults * diffTransFaults);
    int innerIterations = static_cast<int>(ceil(static_cast<float>(config.numFsm) / divisor));
    int totalIterations = static_cast<int>(innerIterations * divisor);


    CLOG(INFO, logging::globalLogger) << "divisor: " << divisor;
    CLOG(INFO, logging::globalLogger) << "innerIterations: " << innerIterations;
    CLOG(INFO, logging::globalLogger) << "totalIterations: " << totalIterations;
    CLOG(INFO, logging::globalLogger) << "";

    CLOG(INFO, logging::globalLogger) << "diffInput: " << diffInput;
    CLOG(INFO, logging::globalLogger) << "diffOutput: " << diffOutput;
    CLOG(INFO, logging::globalLogger) << "diffStates: " << diffStates;
    CLOG(INFO, logging::globalLogger) << "diffOutFaults: " << diffOutFaults;
    CLOG(INFO, logging::globalLogger) << "diffTransFaults: " << diffTransFaults;

    CLOG(INFO, logging::globalLogger) << "";

    int executed = 0;
    int passed = 0;
    int i = 0;

    const TestIteration& loggingContext = config.csvConfig.context;

    vector<vector<AdaptiveTestResult>> collectedResults;

    vector<AdaptiveTestResult> innerCollectedResults;

    for (int ctInput = config.minInput; ctInput <= config.maxInput; ++ctInput)
    {
        for (int ctOutput = config.minOutput; ctOutput <= config.maxOutput; ++ctOutput)
        {
            for (int ctState = config.minStates; ctState <= config.maxStates ; ++ctState)
            {
                for (int ctOutFault = config.minOutFaults; ctOutFault <= config.maxOutFaults; ++ctOutFault)
                {
                    for (int ctTransFault = config.minTransFaults; ctTransFault <= config.maxTransFaults; ++ctTransFault)
                    {
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
                            logging::setLogfileSuffix(iteration);
                #endif


                            int numInputs = ctInput;
                            int numOutputs = ctOutput;
                            int numStates = ctState;
                            int numOutFaults = ctOutFault;
                            int numTransFaults = ctTransFault;


                            if (numOutputs <= 1 && numOutFaults > 0)
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

                            TIMED_SCOPE(timerBlkObj, "heavy-iter");
                            AdaptiveTestResult result;
                            result.testName = config.testName + "-" + iteration;

                            bool couldCreateMutant = false;
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
                                                iteration,
                                                numStates,
                                                numInputs,
                                                numOutputs,
                                                numOutFaults,
                                                numTransFaults,
                                                config.degreeOfCompleteness,
                                                config.maxDegreeOfNonDeterminism,
                                                createRandomFsmSeed,
                                                createMutantSeed,
                                                plTestSpecCopy,
                                                plTestIutCopy,
                                                config.dontTestReductions,
                                                config.csvConfig,
                                                config.loggingConfig,
                                                result);
                                    couldCreateMutant = true;

                                    if (config.csvConfig.context != TestIteration::END)
                                    {
                                        innerCollectedResults.push_back(result);
                                    }

                                }
                                catch (unexpected_reduction& e)
                                {
                                    CLOG(WARNING, logging::globalLogger) << "IUT is a reduction of the specification.";
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
                            } while (!couldCreateMutant && !abort);
                            if (!couldCreateMutant)
                            {
                                ++i;
                                CLOG(INFO, logging::globalLogger) << "numStates: " << numStates + 1;
                                CLOG(INFO, logging::globalLogger) << "numInput: " << numInputs + 1;
                                CLOG(INFO, logging::globalLogger) << "numOutput: " << numOutputs + 1;
                                CLOG(INFO, logging::globalLogger) << "numOutFaults: " << numOutFaults;
                                CLOG(INFO, logging::globalLogger) << "numTransFaults: " << numTransFaults;
                                CLOG(INFO, logging::globalLogger) << "createRandomFsmSeed: " << createRandomFsmSeed;
                                CLOG(INFO, logging::globalLogger) << "createMutantSeed: " << createMutantSeed;
                                //TODO
                                CLOG(WARNING, logging::globalLogger) << "Could not create requested mutant. Skipping.";
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
                }

                if (loggingContext == TestIteration::STATE)
                {
                    collectedResults.push_back(innerCollectedResults);
                    innerCollectedResults = vector<AdaptiveTestResult>();
                }
            } // End States
        }
    }


    std::chrono::steady_clock::time_point totalEnd = std::chrono::steady_clock::now();
    long durationS = std::chrono::duration_cast<std::chrono::seconds>(totalEnd - totalStart).count();
    long durationM = std::chrono::duration_cast<std::chrono::minutes>(totalEnd - totalStart).count();

    if (config.csvConfig.context != TestIteration::END)
    {
        logToCsv(collectedResults, config.csvConfig);
    }

    if (config.loggingConfig.logTestResultsSummary)
    {
        printSummary(config.testName, executed, passed, i - executed, durationS, durationM);
    }
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

    config.minInput = 3;
    config.maxInput = 3;

    config.minOutput = 3;
    config.maxOutput = 3;

    config.minStates = 10;
    config.maxStates = 10;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 1;
    config.maxOutFaults = 1;

    config.degreeOfCompleteness = 0.6f;
    config.maxDegreeOfNonDeterminism = 1.0f;

    config.dontTestReductions = true;
    config.forceTestParameters = true;

    config.seed = 1337;

    config.csvConfig.logEveryIteration = true;

    config.loggingConfig.logTestDetails = true;
    config.loggingConfig.toDot = true;
    config.loggingConfig.printSetsOfMaximalRDistStates = true;

    writeCsvHeader(config.csvConfig, logging::csvLoggerEveryIteration);

    if (debug)
    {
        CLOG(INFO, logging::globalLogger) << "############## Debugging ##############";

        AdaptiveTestConfigDebug debugConfig;
        debugConfig.numStates = 7;
        debugConfig.numInput = 4;
        debugConfig.numOutput = 4;
        debugConfig.numOutFaults = 1;
        debugConfig.numTransFaults = 0;
        debugConfig.createRandomFsmSeed = 1104848806;
        debugConfig.createMutantSeed = 1837447161;
        debugConfig.degreeOfCompleteness = 0.6f;
        debugConfig.maxDegreeOfNonDeterminism = 1.0f;
        debugConfig.csvConfig.logEveryIteration = false;
        debugConfig.loggingConfig.toDot = true;

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
                    debugConfig.degreeOfCompleteness,
                    debugConfig.maxDegreeOfNonDeterminism,
                    debugConfig.createRandomFsmSeed,
                    debugConfig.createMutantSeed,
                    plTestSpecCopy,
                    plTestIutCopy,
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
    loggingConfig.logTestDetails = true;

    AdaptiveTestResult result;
    result.testName = "00-00";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/00-spec.dot", "00-00-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/00-iut.dot", "00-00-iut");
    result.numOutFaults = 0;
    result.numTransFaults = 0;
    executeAdaptiveTest(spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-00-inter", true, false, result);
    printTestResult(result, csvConfig, loggingConfig);
    assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

void test00_01()
{
    CsvConfig csvConfig;
    csvConfig.logEveryIteration = false;

    LoggingConfig loggingConfig;
    loggingConfig.logTestDetails = true;

    AdaptiveTestResult result;
    result.testName = "00-01";
    printTestBegin(result.testName);
    Fsm spec = Fsm(ascTestDirectory + "00/01-spec.dot", "00-01-spec");
    Fsm iut = Fsm(ascTestDirectory + "00/01-iut.dot", "00-01-iut");

    result.numOutFaults = 1;
    result.numTransFaults = 0;
    executeAdaptiveTest(spec, iut, static_cast<size_t>(iut.getMaxNodes()), "00-01-inter", true, false, result);
    printTestResult(result, csvConfig, loggingConfig);
    assert(result.testName, result.pass);
    CLOG(INFO, logging::globalLogger) << testSepLine;
}

/**
 * Testing FSMs with number of states varying
 * and always one output fault.
 */
void test01_00()
{
    AdaptiveTestConfig config;
    config.testName = "ST-IF";
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

    config.degreeOfCompleteness = 0.6f;
    config.maxDegreeOfNonDeterminism = 1.0f;

    config.dontTestReductions = true;

    config.loggingConfig.toDot = false;
    config.loggingConfig.printSetsOfMaximalRDistStates = false;

    config.csvConfig.logEveryIteration = true;
    config.csvConfig.context = TestIteration::STATE;
    config.csvConfig.fieldsContext.push_back(CsvField::DURATION_MS);
    config.csvConfig.fieldsContext.push_back(CsvField::OBSERVED_TRACES_SIZE);

    config.seed = 1337;

    adaptiveTestRandom(config);
}

void test01_01()
{
    AdaptiveTestConfig config;
    config.testName = "ST-RED";
    config.numFsm = 20000;

    config.minInput = 4;
    config.maxInput = 4;

    config.minOutput = 4;
    config.maxOutput = 4;

    config.minStates = 1;
    config.maxStates = 10;

    config.minTransFaults = 0;
    config.maxTransFaults = 0;

    config.minOutFaults = 0;
    config.maxOutFaults = 0;

    config.degreeOfCompleteness = 0.6f;
    config.maxDegreeOfNonDeterminism = 1.0f;

    config.dontTestReductions = false;

    config.loggingConfig.toDot = false;
    config.loggingConfig.printSetsOfMaximalRDistStates = false;

    config.csvConfig.logEveryIteration = true;
    config.csvConfig.context = TestIteration::STATE;
    config.csvConfig.fieldsContext.push_back(CsvField::DURATION_MS);
    config.csvConfig.fieldsContext.push_back(CsvField::OBSERVED_TRACES_SIZE);

    config.seed = 3537;

    adaptiveTestRandom(config);
}

void runAdaptiveStateCountingTests()
{
    /*
    test00_00();
    test00_01();
    */

    test01_00();
    test01_01();
}

int main(int argc, char* argv[])
{
//    locale::global( locale("") );
    START_EASYLOGGINGPP(argc, argv);
    logging::initLogging();

    initialize();

    CLOG(INFO, logging::globalLogger) << "############## Starting Application ##############";
#ifdef ENABLE_DEBUG_MACRO
    CLOG(INFO, logging::globalLogger) << "This is a debug build!";
#else
    CLOG(INFO, logging::globalLogger) << "This is a release build!";
#endif

    trial(false);

    return 0;

    runAdaptiveStateCountingTests();

    return 0;

    cout << endl << endl;
}



