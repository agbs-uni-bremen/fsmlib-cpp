#include "sut_wrapper_adaptive.h"
#include <string>


enum SUT_State {
    S1,
    S2,
    S3,
    S4
};

enum SUT_Input {
    Ia,
    Ib
};

enum SUT_Output {
    O0,
    O1
};

SUT_State currentState;


void sut_init() { 
    sut_reset();
}

void sut_reset() { 
    currentState = S1;
}



// variables to store the last response to inputs to state that cause nondeterministic behaviour,
// used to cycle through all possible responses
SUT_Output S1Ia_last_response = O1;
SUT_Output S3Ib_last_response = O1;

const std::string sut(const std::string input) {
    
    SUT_Input sut_input;
    
    if (input.compare("a") == 0) {
        sut_input = Ia;
    } else if (input.compare("b") == 0) {
        sut_input = Ib; 
    } else {
        return "ERROR: invalid input string: " + input;
    }

    SUT_Output output;

    switch (currentState) {
        case S1 : 
            switch (sut_input) {
                case Ia: 
                    switch (S1Ia_last_response) {
                        case O0: 
                            currentState = S4;
                            output = O1;
                            break;
                        case O1: 
                            currentState = S2;
                            output = O0;
                            break;
                    };
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
                    switch (S3Ib_last_response) {
                        case O0: 
                            currentState = S3;
                            output = O1;
                            break;
                        case O1: 
                            currentState = S1;
                            output = O0;
                            break;
                    };
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

    switch (output) {
        case O0: return "0";
        case O1: return "1";
        default : return "ERROR: invalid output";
    }

}