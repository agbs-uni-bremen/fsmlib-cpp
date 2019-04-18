#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>
#include <random>
#include <iomanip>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"
#include "trees/InputOutputTree.h"
#include "trees/SplittingTree.h"
#include "trees/IOTreeContainer.h"

using namespace std;

unsigned int getRandomSeed() {

    return static_cast<unsigned int>
    (std::chrono::high_resolution_clock::now().time_since_epoch().count());

}

std::mt19937 rand_engine (getRandomSeed());  // replace knuth_b with one of the engines listed below
// std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);

bool random_bool_with_prob( double prob )  // probability between 0.0 and 1.0
{
    std::bernoulli_distribution d(prob);
    return d(rand_engine);
}

shared_ptr<FsmPresentationLayer> createPresentationLayer(const size_t maxInput, const size_t refMin, const size_t maxOutput)
{
    vector<string> in2string = vector<string>();
    for (size_t i = 0; i <= maxInput; i++)
    {
        in2string.push_back(to_string(i));
    }
    vector<string> out2string = vector<string>();
    for (size_t i = 0; i <= maxOutput; i++)
    {
        out2string.push_back(to_string(i));
    }
    vector<string>state2string = vector<string>();
    for (size_t i = 0; i < refMin; i++)
    {
        state2string.push_back(to_string(i));
    }
    shared_ptr<FsmPresentationLayer> pl{
            new FsmPresentationLayer(
                    in2string,
                    out2string,
                    state2string)};

    return pl;
}
bool randomBool() {
    return rand() > (RAND_MAX / 2);
}

void testLeeAds()
{
    cout << "Starting ADS Test for lee94_no_pds.fsm..." << endl;
    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(1,6,1);
    shared_ptr<Dfsm> dfsm = make_shared<Dfsm>("../../../resources/lee94_no_pds.fsm",pl,"lee94_no_pds");

    vector<int> ds = dfsm->createDistinguishingSequence();

    //lee94_no_pds does not possess a pds
    assert(ds.empty());

    shared_ptr<InputOutputTree> ads = dfsm->createAdaptiveDistinguishingSequence();

    //lee94_no_pds does possess an ads
    assert(ads);

    auto nds = dfsm->getNodes();
    auto& nodes = nds;

    auto adsList = make_shared<vector<shared_ptr<InputOutputTree>>>();
    adsList->push_back(ads);
    IOTreeContainer adaptiveTestCases(adsList,dfsm->getPresentationLayer());

    //the ads should distinguish all states from each other
    assert(dfsm->distinguishesAllStates(nodes,nodes,adaptiveTestCases));
    cout << "Finished with SUCCESS!" << endl << endl;
}

void testRandomPdsAndAds(const int numStates,const int numInput,const int numOutput)
{
    cout << "Starting PDS and ADS creation Test for with random DFSM..." << endl;

    int numberOfTests = 100;

    for(int i=0;i<numberOfTests;++i) {
        shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(numInput, numStates, numOutput);
        //create random minimised dfsm with certain number of states, inputs and outputs
        shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl)->minimise());

        vector<int> ds = dfsm->createDistinguishingSequence();
        while(ds.empty()) {
            dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", numStates, numInput, numOutput, pl)->minimise());
            ds = dfsm->createDistinguishingSequence();
        }

        shared_ptr<InputOutputTree> ads = dfsm->createAdaptiveDistinguishingSequence();

        //ads should exist for dfsm with pds
        assert(ads);

        auto nds = dfsm->getNodes();
        auto& nodes = nds;

        InputTrace distinguishingSequence(ds,pl);

        //distinguishing sequence should distinguish all states from each other
        assert(dfsm->distinguishesAllStates(nodes,distinguishingSequence));

        auto adsList = make_shared<vector<shared_ptr<InputOutputTree>>>();
        adsList->push_back(ads);
        IOTreeContainer adaptiveTestCases(adsList,dfsm->getPresentationLayer());

        //adaptive distinguishing sequence should distinguish all states from each other
        assert(dfsm->distinguishesAllStates(nodes,nodes,adaptiveTestCases));
    }
    cout << "Finished with SUCCESS!" << endl << endl;
}

void evaluateDMethodsFaultCoverage() {

    cout << " Start D-Methods Fault Coverage Evaluation..." << endl;

    int numDfsm = 10,
        numMutants = 100;

    vector<vector<int>> dfsm_sizes = {{10,2,2},
                                      {10,3,3},
                                      {20,4,4},
                                      {20,6,6},
                                      {30,6,6},
                                      {30,8,8},
                                      {50,10,10},
                                      {50,13,13},
                                      {75,15,15},
                                      {75,18,18},
                                      {100,20,20},
                                      {100,24,24}
                                      };
    ofstream out("dmethods_fc.csv");
    out << "states  ;"
        << "inputs  ;"
        << "outputs ;"
        << "unequal mutants ;"
        << "not pass. D-Method; "
        << "not pass. D-Method(ADS); "
        << "not pass. Hierons D-Method; "
        << "not pass. Hierons D-Method(ADS); "
        << "D-Method Fault Coverage;"
        << "D-Method(ADS) Fault Coverage;"
        << "Hierons D-Method Fault Coverage;"
        << "Hierons D-Method(ADS) Fault Coverage" << endl;

    for(auto& dfsm_size:dfsm_sizes) {
        cout << endl << "Evaluate Fault Coverage for DFSM with: " << endl
             << "\t -states -> " << to_string(dfsm_size[0]) << endl
             << "\t -inputs -> " << to_string(dfsm_size[1]) << endl
             << "\t -outputs -> " << to_string(dfsm_size[2]) << endl;

        float sidPds_fc = 0,
              sidAds_fc = 0,
              hierPds_fc = 0,
              hierAds_fc = 0;

        int num_nonequal_mut = 0,
                num_not_passing_sidPds = 0,
                num_not_passing_sidAds = 0,
                num_not_passing_hierPds = 0,
                num_not_passing_hierAds = 0;

        for(int i=0;i<numDfsm;++i) {

            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());

            vector<int> ds = dfsm->createDistinguishingSequence();
            while(ds.empty()) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                        dfsm_size[2], pl)->minimise());
                ds = dfsm->createDistinguishingSequence();
            }

            cout << "Created next random DFSM(" << i << ")..." << endl;

            auto a = dfsm->dMethodOnMinimisedDfsm(0,false);
            vector<IOTrace> sidPdsTs;
            for(vector<int> v: *a.getIOLists()) {
                InputTrace i(v,pl);
                IOTrace io = dfsm->applyDet(i);
                sidPdsTs.push_back(io);
            }

            auto b = dfsm->dMethodOnMinimisedDfsm(0,true);
            vector<IOTrace> sidAdsTs;
            for(vector<int> v: *b.getIOLists()) {
                InputTrace i(v,pl);
                IOTrace io = dfsm->applyDet(i);
                sidAdsTs.push_back(io);
            }

            auto c = dfsm->hieronsDMethodOnMinimisedDfsm(false);
            vector<IOTrace> hierPdsTs;
            for(vector<int> v: *c.getIOLists()) {
                InputTrace i(v,pl);
                IOTrace io = dfsm->applyDet(i);
                hierPdsTs.push_back(io);
            }

            auto d = dfsm->hieronsDMethodOnMinimisedDfsm(false);
            vector<IOTrace> hierAdsTs;
            for(vector<int> v: *c.getIOLists()) {
                InputTrace i(v,pl);
                IOTrace io = dfsm->applyDet(i);
                hierAdsTs.push_back(io);
            }
            cout << "Created testsuites for all D-Method variants..." << endl;

            for(int j=0;j<numMutants;++j) {
                unsigned int numOutputFaults = rand() % (dfsm->getMaxNodes()/3),
                        numTransitionFaults = rand() % (dfsm->getMaxNodes()/3);
                if(random_bool_with_prob(0.4)) {
                    numOutputFaults = 0;
                    numTransitionFaults = 0;
                }
                shared_ptr<Dfsm> mutant = dfsm->createMutant("Mutant",numOutputFaults,numTransitionFaults);
                //cout << "Created next random Mutant(" << j << ") with: "
                //     << "\t -num. output faults -> " << numOutputFaults
                //     << "\t -num. transition faults -> " << numTransitionFaults << endl;

                bool result = dfsm->equivalenceCheck(*mutant);
                if (!result) num_nonequal_mut++;

                result = true;
                for (IOTrace io: sidPdsTs) {
                    result &= mutant->pass(io);
                }
                if (!result) num_not_passing_sidPds++;

                result = true;
                for (IOTrace io: sidAdsTs) {
                    result &= mutant->pass(io);
                }
                if (!result) num_not_passing_sidAds++;

                result = true;
                for (IOTrace io: hierPdsTs) {
                    result &= mutant->pass(io);
                }
                if (!result) num_not_passing_hierPds++;

                result = true;
                for (IOTrace io: hierAdsTs) {
                    result &= mutant->pass(io);
                }
                if (!result) num_not_passing_hierAds++;

                //cout << "All testsuites executed on the mutant!" << endl;
            }
        }

        sidPds_fc = (num_not_passing_sidPds>0?(float)num_nonequal_mut/(float)num_not_passing_sidPds*100:
                   ((num_nonequal_mut == 0)?100:0));
        sidAds_fc = (num_not_passing_sidAds>0?(float)num_nonequal_mut/(float)num_not_passing_sidAds*100:
                     ((num_nonequal_mut == 0)?100:0));
        hierPds_fc = (num_not_passing_hierPds>0?(float)num_nonequal_mut/(float)num_not_passing_hierPds*100:
                     ((num_nonequal_mut == 0)?100:0));
        hierAds_fc = (num_not_passing_hierAds>0?(float)num_nonequal_mut/(float)num_not_passing_hierAds*100:
                      ((num_nonequal_mut == 0)?100:0));

        out << dfsm_size[0] << ";"
            << dfsm_size[1] << ";"
            << dfsm_size[2] << ";"
            << num_nonequal_mut << ";"
            << num_not_passing_sidPds << ";"
            << num_not_passing_sidAds << ";"
            << num_not_passing_hierPds << ";"
            << num_not_passing_hierAds << ";"
            << sidPds_fc << ";"
            << sidAds_fc << ";"
            << hierPds_fc << ";"
            << hierAds_fc << endl;
    }
    out.close();
    cout << "Finished!" << endl << endl;
}

void evaluateTestSuiteSizes() {
    cout << " Start Test Suite Size Evaluation..." << endl;

    int numDfsm = 50;

    vector<vector<int>> dfsm_sizes = {{10,1,1}, //20
                                      {10,2,2}, //30
                                      {10,3,3}, //40
                                      {20,2,2}, //60
                                      {20,3,3}, //80
                                      {20,4,4}, //100
                                      {20,5,5}, //120
                                      {30,4,4}, //150
                                      {30,5,5}, //180
                                      {30,6,6}, //210
                                      {50,4,4}, //250
                                      {50,6,6}, //350
                                      {50,7,7}, //400
                                      {50,8,8}, //450
                                      {70,7,7}, //560
                                      {70,8,8}, //630
                                      {70,9,9}, //700
                                      {70,11,11}, //840
                                      {80,13,13}, //1120
                                      {100,9,9}, //1000
                                      {100,11,11}, //1200
                                      {100,12,12}, //1300
    };

    std::cout << std::fixed << std::showpoint;
    std::cout << std::setprecision(2);

    ofstream out("testsuite_sizes.csv");
    out << "states  ;"
        << "inputs  ;"
        << "outputs ;"
        << "transitions ;"
        << "W-Method;"
        << "w Duration;"
        << "Wp-Method;"
        << "wp Duration;"
        << "D-Method;"
        << "D Duration;"
        << "D-Method(ADS);"
        << "D (ADS) Duration;"
        << "Hierons D-Method;"
        << "Hier. D Duration;"
        << "Hierons D-Method(ADS);"
        << "Hier. D (ADS) Duration;"
        << "T-Method;"
        << "T Duration" << endl;

    for(auto& dfsm_size:dfsm_sizes) {
        cout << endl << "Evaluate Test suite size for DFSM with: " << endl
             << "\t -states -> " << to_string(dfsm_size[0]) << endl
             << "\t -inputs -> " << to_string(dfsm_size[1]) << endl
             << "\t -outputs -> " << to_string(dfsm_size[2]) << endl;

        double avgWMethSize = 0,
            avgWMethDuration = 0,
            avgWpMethSize = 0,
            avgWpMethDuration = 0,
            avgDMethSize = 0,
            avgDMethDuration = 0,
            avgADSMethSize = 0,
            avgADSMethDuration = 0,
            avgHDMethSize = 0,
            avgHDMethDuration = 0,
            avgHADSMethSize = 0,
            avgHADSMethDuration = 0,
            avgTMethSize = 0,
            avgTMethDuration = 0;

        cout << "Evaluate Testsuite size for W-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0]) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->wMethodOnMinimisedDfsm(0);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgWMethDuration += elapsed.count();
            avgWMethSize += ts.getFlatSize();
        }
        avgWMethDuration = (double) avgWMethDuration/(double) numDfsm;
        avgWMethSize = (double) avgWMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for Wp-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0]) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->wpMethodOnMinimisedDfsm(0);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgWpMethDuration += elapsed.count();
            avgWpMethSize += ts.getFlatSize();
        }
        avgWpMethDuration = (double) avgWpMethDuration/(double) numDfsm;
        avgWpMethSize = (double) avgWpMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for T-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0]) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->tMethod();

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgTMethDuration += elapsed.count();
            avgTMethSize += ts.getFlatSize();
        }
        avgTMethDuration = (double) avgTMethDuration/(double) numDfsm;
        avgTMethSize = (double) avgTMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for D-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            auto ds = dfsm->createDistinguishingSequence();
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0] || ds.empty()) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
                ds = dfsm->createDistinguishingSequence();
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->dMethodOnMinimisedDfsm(0,false);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgDMethDuration += elapsed.count();
            avgDMethSize += ts.getFlatSize();
        }
        avgDMethDuration = (double) avgDMethDuration/(double) numDfsm;
        avgDMethSize = (double) avgDMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for D-Method(ADS):" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            auto ds = dfsm->createDistinguishingSequence();
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0] || ds.empty()) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
                ds = dfsm->createDistinguishingSequence();
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->dMethodOnMinimisedDfsm(0,true);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgADSMethDuration += elapsed.count();
            avgADSMethSize += ts.getFlatSize();
        }
        avgADSMethDuration = (double) avgADSMethDuration/(double) numDfsm;
        avgADSMethSize = (double) avgADSMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for Hierons D-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            auto ds = dfsm->createDistinguishingSequence();
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0] || ds.empty()) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
                ds = dfsm->createDistinguishingSequence();
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->hieronsDMethodOnMinimisedDfsm(false);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgHDMethDuration += elapsed.count();
            avgHDMethSize += ts.getFlatSize();
        }
        avgHDMethDuration = (double) avgHDMethDuration/(double) numDfsm;
        avgHDMethSize = (double) avgHDMethSize/(double) numDfsm;

        cout << "Evaluate Testsuite size for Hierons D-Method:" << endl;
        for(int i=0;i<numDfsm;++i) {
            shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(dfsm_size[0], dfsm_size[1], dfsm_size[2]);
            //create random minimised dfsm with certain number of states, inputs and outputs
            shared_ptr<Dfsm> dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                                        dfsm_size[2], pl)->minimise());
            auto ds = dfsm->createDistinguishingSequence();
            //generate random dfsm that are already minimised
            while(dfsm->size() < dfsm_size[0] || ds.empty()) {
                dfsm = make_shared<Dfsm>(make_shared<Dfsm>("Dfsm", dfsm_size[0], dfsm_size[1],
                                                           dfsm_size[2], pl)->minimise());
                ds = dfsm->createDistinguishingSequence();
            }
            cout << "Created next random minimised DFSM(" << i << ")..." << endl;

            auto start = chrono::high_resolution_clock::now();

            auto ts = dfsm->hieronsDMethodOnMinimisedDfsm(true);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgHADSMethDuration += elapsed.count();
            avgHADSMethSize += ts.getFlatSize();
        }
        avgHADSMethDuration = (double) avgHADSMethDuration/(double) numDfsm;
        avgHADSMethSize = (double) avgHADSMethSize/(double) numDfsm;

        out << dfsm_size[0] << ";"
            << dfsm_size[1] << ";"
            << dfsm_size[2] << ";"
            << dfsm_size[0] * dfsm_size[1] << ";"
            << avgWMethSize << ";"
            << avgWMethDuration << ";"
            << avgWpMethSize << ";"
            << avgWpMethDuration << ";"
            << avgDMethSize << ";"
            << avgDMethDuration << ";"
            << avgADSMethSize << ";"
            << avgADSMethDuration << ";"
            << avgHDMethSize << ";"
            << avgHDMethDuration << ";"
            << avgHADSMethSize << ";"
            << avgHADSMethDuration << ";"
            << avgTMethSize << ";"
            << avgTMethDuration << endl;
    }
    out.close();
    cout << "Finished!" << endl << endl;
}

int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    //evaluateDMethodsFaultCoverage();
    evaluateTestSuiteSizes();
    //testRandomPdsAndAds(6,2,2);
    //testLeeAds();

}