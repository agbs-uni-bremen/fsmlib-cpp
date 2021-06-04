/**
 * This "wrapper" simulates the behaviour of an FSM with a 
 * randomly executed transition fault.
 */

#include <string>


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

void sut_reset() { 
    currentState = S1;
    input_length = 0;
}

void sut_init() { 
    srand((unsigned)time(0)); 
    sut_reset();
}



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

    ++input_length;

    switch (currentState) {
        case S1 : 
            switch (sut_input) {
                case Ia: 
                    if (rand() % 2 == 0) {
                        currentState = S4;
                        output = O1;
                    } else {
                        currentState = S2;
                        output = O0;
                    }

                    //randomly injected transition error
                    //if (rand() % 100 == 0) {
                    //    currentState = S2;
                    //    output = O1;
                    //}

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
                    if (rand() % 2 == 0) {    
                        currentState = S1;
                        output = O0;
                    } else {
                        currentState = S3;
                        output = O1;
                    }

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