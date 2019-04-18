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
    out << "states  ;"
        << "inputs  ;"
        << "outputs ;"
        << "transitions ;"
        << "H-Method ;"
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


        out << dfsm_size[0] << ";"
            << dfsm_size[1] << ";"
            << dfsm_size[2] << ";"
            << dfsm_size[0] * dfsm_size[1] << ";"
            << avgHMethSize << ";"
            << avgHMethDuration << endl;
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