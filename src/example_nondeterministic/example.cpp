/**
 * Implementation of a simple example SUT representing a state machine with four states.
 */

#include <example_nondeterministic/example.h>

#include <ctime>
#include <stdlib.h>

// the state the SUT currently resides in
SUT_State currentState = S1;

// reset the SUT by setting the current state to the initial state S1
void reset() { 
    currentState = S1;
}

// initialise the SUT by initialising the random number generator and resetting it
void init() { 
    srand((unsigned)time(0)); 
    reset();
}

// Apply an input to the SUT, generating an outputs and possibly modifying the current state.
const SUT_Output apply(const SUT_Input input) {
    
    SUT_Output output;

    switch (currentState) {
        case S1 : 
            switch (input) {
                case Ia: 
                    if (rand() % 2 == 0) {
                        currentState = S4;
                        output = O1;
                    } else {
                        currentState = S2;
                        output = O0;
                    }

                    //randomly injected transition error
                    if (rand() % 100 == 0) {
                        currentState = S2;
                        output = O1;
                    }

                    break;
                case Ib:
                    currentState = S4;
                    output = O1;
                    break;
            };
            break;
        case S2 : 
            switch (input) {
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
            switch (input) {
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
            switch (input) {
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
    }

    return output;
}