#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>

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


void testLeeAds()
{
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
}

void testRandomPdsAndAds(const int numStates,const int numInput,const int numOutput)
{
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
}

int main(int argc, char* argv[])
{

    srand(getRandomSeed());

    testRandomPdsAndAds(6,2,2);
    testLeeAds();
}