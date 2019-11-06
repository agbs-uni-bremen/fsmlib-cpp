#include "sut_wrapper_adaptive.h"
#include <string>

#include "utils/Logger.hpp"


enum SUT_State {
    S1,
    S2,
    S3,
    S4
};

static const std::string sut_state_strings[] = {"S1", "S2", "S3", "S4"};

enum SUT_Input {
    Ia,
    Ib
};

enum SUT_Output {
    O0,
    O1
};

SUT_State currentState;
unsigned int input_length;

void sut_init() { 
    sut_reset();
}

void sut_reset() { 
    currentState = S1;
    input_length = 0;
}



// variables to store the last response to inputs to state that cause nondeterministic behaviour,
// used to cycle through all possible responses
//SUT_Output S1Ia_last_response = O1;
//SUT_Output S3Ib_last_response = O1;

unsigned int s1ia_transitions = 0;
unsigned int s3ib_transitions = 0;

const std::string sut(const std::string input) {
    
    LOG("VERBOSE_SUT_INTERNAL") << "\tSUT called:" << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tCurrent state    : " << sut_state_strings[currentState] << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tInput            : " << input << std::endl;

    SUT_Input sut_input;
    
    if (input.compare("a") == 0) {
        sut_input = Ia;
    } else if (input.compare("b") == 0) {
        sut_input = Ib; 
    } else {
        return "ERROR: invalid input string: " + input;
    }

    SUT_Output output;

    ++input_length;

    switch (currentState) {
        case S1 : 
            switch (sut_input) {
                case Ia: 
                    // switch (S1Ia_last_response) {
                    //     case O0: 
                    //         currentState = S4;
                    //         output = O1;
                    //         break;
                    //     case O1: 
                    //         currentState = S2;
                    //         output = O0;
                    //         break;
                    // };
                    // S1Ia_last_response = output;

                    //if (s1ia_transitions % 3 == 0 || s1ia_transitions % 17 == 0) {
                    if (rand() % 2 == 0) {
                        currentState = S4;
                        output = O1;
                    } else {
                        currentState = S2;
                        output = O0;
                    }

                    ++s1ia_transitions;
                    break;
                case Ib:
                    currentState = S4;
                    output = O1;
                    break;
            };
            break;
        case S2 : 
            switch (sut_input) {
                case Ia: 
                    currentState = S2;
                    output = O0;
                    break;
                case Ib:
                    currentState = S4;
                    output = O1;
                    break;
            };
            break;
        case S3 : 
            switch (sut_input) {
                case Ia: 
                    currentState = S4;
                    output = O1;
                    break;
                case Ib: 
                    // switch (S3Ib_last_response) {
                    //     case O0: 
                    //         currentState = S3;
                    //         output = O1;
                    //         break;
                    //     case O1: 
                    //         currentState = S1;
                    //         output = O0;
                    //         break;
                    // };
                    // S3Ib_last_response = output;

                    // choose different threshold (!= that for s1ia_transitions)
                    //if (s3ib_transitions % 5 == 0 || s3ib_transitions % 7 == 0) {
                    if (rand() % 2 == 0) {    
                        currentState = S1;
                        output = O0;
                    } else {
                        currentState = S3;
                        output = O1;
                    }

                    ++s3ib_transitions;

                    break;
            };
            break;
        case S4 : 
            switch (sut_input) {
                case Ia: 
                    currentState = S3;
                    output = O0;
                    break;
                case Ib:
                    currentState = S1;
                    output = O0;
                    break;
            };
            break;
        default : return "ERROR: invalid state";
    }

    LOG("VERBOSE_SUT_INTERNAL") << "\t\tOutput           : " << output << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\tPost state       : " << sut_state_strings[currentState] << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\ts1ia_transitions : " << s1ia_transitions << std::endl;
    LOG("VERBOSE_SUT_INTERNAL") << "\t\ts3ib_transitions : " << s3ib_transitions << std::endl;

    switch (output) {
        case O0: return "0";
        case O1: return "1";
        default : return "ERROR: invalid output";
    }
}