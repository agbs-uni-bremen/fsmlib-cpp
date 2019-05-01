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

    ofstream out("testsuite_size_hMethod.csv");
    out << "states  ,"
        << "inputs  ,"
        << "outputs ,"
        << "transitions ,"
        << "H-Method ,"
        << "H-Method (duration) " << endl;

    for(auto& dfsm_size:dfsm_sizes) {
        cout << endl << "Evaluate Test suite size for DFSM with: " << endl
             << "\t -states -> " << to_string(dfsm_size[0]) << endl
             << "\t -inputs -> " << to_string(dfsm_size[1]) << endl
             << "\t -outputs -> " << to_string(dfsm_size[2]) << endl;

        double avgHMethSize = 0,
            avgHMethDuration = 0;

        cout << "Evaluate Testsuite size for H-Method:" << endl;
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

            auto ts = dfsm->hMethodOnMinimisedDfsm(0);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            avgHMethDuration += elapsed.count(); //in seconds
            avgHMethSize += ts.getFlatSize();
        }
        avgHMethDuration = (double) avgHMethDuration / (double) numDfsm;
        avgHMethSize = (double) avgHMethSize/(double) numDfsm;


        out << dfsm_size[0] << ","
            << (dfsm_size[1]+1) << ","
            << (dfsm_size[2]+1) << ","
            << (dfsm_size[0]+1) * (dfsm_size[1]+1) << ","
            << avgHMethSize << ","
            << avgHMethDuration << endl;
    }
    out.close();
    cout << "Finished!" << endl << endl;
}

int calcNumAddStatesFCTreshold(IOListContainer& ts,double fcTreshold,int numMutants,
                               shared_ptr<FsmPresentationLayer>& pl,shared_ptr<Dfsm>& dfsm) {

    vector<IOTrace> traces;
    for(vector<int> v: *ts.getIOLists()) {
        InputTrace i(v,pl);
        IOTrace io = dfsm->applyDet(i);
        traces.push_back(io);
    }

    double fc = 100;
    int numAddStates = 0;

    while(fc >= fcTreshold) {
        ++numAddStates;
        //cout << "numAddStates -> " << numAddStates << endl;

        int num_non_equal = 0;
        int num_not_pass = 0;
        fc = 0;
        for(int j=0;j<numMutants;++j) {
            //unsigned int numOutputFaults = rand() % (dfsm->size()/3),
            //       numTransitionFaults = rand() % (dfsm->size()/3);

            shared_ptr<Dfsm> mutant = make_shared<Dfsm>(
                    dfsm->createMutant("Mutant",0,0,numAddStates)->minimise());
            //cout << "mutant number " << j << endl;
            while(mutant->size() < (dfsm->size() + numAddStates)) {
                //cout << "try another mut" << endl;
                mutant = make_shared<Dfsm>(
                        dfsm->createMutant("Mutant",0,0,numAddStates)->minimise());
            }

            bool result = dfsm->equivalenceCheck(*mutant);
            if (!result) num_non_equal++;

            result = true;
            for (IOTrace io: traces) {
                result &= mutant->pass(io);
            }
            if (!result) num_not_pass++;
        }
        //cout << "num not pass -> " << num_not_pass << ", num non equal -> " << num_non_equal << endl;
        fc = (num_non_equal>0?(double)num_not_pass/(double)num_non_equal*100:
              ((num_not_pass == 0)?100:0));
        //cout << "fault coverage -> " << fc << endl;

    }

    return numAddStates;
}

void evaluateFCOutsideFaultDomain() {
    cout << " Start Fault Coverage Evaluation for Additional States(Outside of the Faul Domain)..." << endl;

    int numDfsm = 20;
    int numMuts = 40;
    double faultCoverageTreshhold = 60;

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

    ofstream out("fc_for_add_states_hmeth.csv");
    out << "states  ,"
        << "inputs  ,"
        << "outputs ,"
        << "transitions ,"
        << "H-Method" << endl;

    for(auto& dfsm_size:dfsm_sizes) {
        cout << endl << "Evaluate Test suite size for DFSM with: " << endl
             << "\t -states -> " << to_string(dfsm_size[0]) << endl
             << "\t -inputs -> " << to_string(dfsm_size[1]) << endl
             << "\t -outputs -> " << to_string(dfsm_size[2]) << endl;

        double avgHMethodAddStates = 0;

        cout << "Evaluate Fault Coverage for H-Method:" << endl;
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

            auto ts = dfsm->hMethodOnMinimisedDfsm(0);

            avgHMethodAddStates += calcNumAddStatesFCTreshold(ts,faultCoverageTreshhold,numMuts,pl,dfsm);
        }
        avgHMethodAddStates = (double) avgHMethodAddStates/(double) numDfsm;

        out << dfsm_size[0] << ","
            << (dfsm_size[1]+1) << ","
            << (dfsm_size[2]+1) << ","
            << (dfsm_size[0]) * (dfsm_size[1]+1) << ","
            << avgHMethodAddStates <<  endl;
    }
    out.close();
    cout << "Finished!" << endl << endl;
}

void evaluateTestCaseLength() {
    cout << " Start Test Case Length Evaluation..." << endl;

    int numDfsm = 25;

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

    ofstream out("testcase_average_length_hmethod.csv");
    out << "states  ,"
        << "inputs  ,"
        << "outputs ,"
        << "transitions ,"
        << "H-Method" << endl;

    for(auto& dfsm_size:dfsm_sizes) {
        cout << endl << "Evaluate Testcase length for DFSM with: " << endl
             << "\t -states -> " << to_string(dfsm_size[0]) << endl
             << "\t -inputs -> " << to_string(dfsm_size[1]) << endl
             << "\t -outputs -> " << to_string(dfsm_size[2]) << endl;

        double avgHMethSize = 0;

        cout << "Evaluate Testcase length for H-Method:" << endl;
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

            auto ts = dfsm->hMethodOnMinimisedDfsm(0);

            int avgTCLength = 0;
            for(auto& tc:*ts.getIOLists()) {
                avgTCLength += tc.size();
            }
            avgTCLength = (double) avgTCLength / (double) ts.getIOLists()->size();

            avgHMethSize += avgTCLength;
        }
        avgHMethSize = (double) avgHMethSize/(double) numDfsm;

        out << dfsm_size[0] << ","
            << (dfsm_size[1]+1) << ","
            << (dfsm_size[2]+1) << ","
            << (dfsm_size[0]) * (dfsm_size[1]+1) << ","
            << avgHMethSize <<  endl;
    }
    out.close();
    cout << "Finished!" << endl << endl;
}


int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    evaluateTestCaseLength();
    evaluateFCOutsideFaultDomain();
    evaluateTestSuiteSizes();
}