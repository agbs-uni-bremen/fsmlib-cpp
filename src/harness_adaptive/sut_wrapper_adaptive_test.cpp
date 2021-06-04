#include "sut_wrapper_adaptive_test.h"


std::shared_ptr<Fsm> sut_fsm;
std::shared_ptr<Fsm> ref_fsm;

std::shared_ptr<FsmNode> currentState;

void sut_init(const std::shared_ptr<Fsm>& spec, bool mutate) { 
    srand((unsigned)time(0)); 

    ref_fsm = spec;

    if (mutate) {
        sut_fsm = spec->createMutant( spec->getName() + "_mutant", 1, 1, true);
    } else {
        int removedTransitions = 0;
        sut_fsm = spec->createReduction( spec->getName() + "_reduction", false, removedTransitions);
    }

    

    sut_reset();
}

void sut_reset() { 
    currentState = sut_fsm->getInitialState();
}


const std::string sut(const std::string input) {
    
    LOG("VERBOSE_SUT_INTERNAL") << "\tSUT called:" << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tCurrent state    : " << currentState->getName();
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tInput            : " << input << std::endl;

    int inputNo = ref_fsm->getPresentationLayer()->in2Num(input);
    
    if (inputNo == -1) {
        return "ERROR: invalid input string: " + input;
    }

    std::vector<std::shared_ptr<FsmTransition>> transitionsForInput;
    for (const auto& t : currentState->getTransitions()) {
        if (t->getLabel()->getInput() == inputNo) {
            transitionsForInput.push_back(t);
        }
    }

    if (transitionsForInput.size() == 0) {
        return "ERROR: no transitions defined for input: " + input;
    }

    std::shared_ptr<FsmTransition> randomChoice = transitionsForInput.at(rand() % transitionsForInput.size());

    currentState = randomChoice->getTarget();

    std::string output = ref_fsm->getPresentationLayer()->getOutId(randomChoice->getLabel()->getOutput());

    LOG("VERBOSE_SUT_INTERNAL") << "\t\tOutput           : " << output << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tPost state       : " << currentState->getName() << std::endl;
    
    return output;
}