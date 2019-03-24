#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>

#include "fsm/Dfsm.h"
#include "fsm/Fsm.h"
#include "fsm/InputTrace.h"
#include "fsm/IOTrace.h"
#include "interface/FsmPresentationLayer.h"
#include "trees/IOListContainer.h"

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

int main(int argc, char* argv[])
{

    //Logging
    //LogCoordinator& logger = LogCoordinator::getStandardLogger();
    //logger.setDefaultStream(cout);

    srand(getRandomSeed());

    shared_ptr<FsmPresentationLayer> pl = createPresentationLayer(6,30,6);
    auto dfsm = make_shared<Dfsm>("Dfsm", 30, 6, 6, pl);
    shared_ptr<Dfsm> dfsmMin = make_shared<Dfsm>(dfsm->minimise());
    //dfsmMin->toDot("dfsm_min_fc");
    //IOListContainer ts = dfsmMin->dMethodOnMinimisedDfsm(0);

}

