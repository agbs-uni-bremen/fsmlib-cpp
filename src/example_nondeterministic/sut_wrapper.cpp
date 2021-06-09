/**
 * This wrapper connects the test harness with an example SUT.
 */

#include <example_nondeterministic/example.h>

#include <string>

// use the reset() operation of the SUT
void sut_reset() { 
    reset();
}

// use the init() operation of the SUT
void sut_init() { 
    init();
}

// translate an input ("a"->Ia; "b"->Ib), apply it to the SUT and translate the output (O0->"0"; O1->"1")
const std::string sut(const std::string input) {
    
    SUT_Input sut_input;
    
    if (input.compare("a") == 0) {
        sut_input = Ia;
    } else if (input.compare("b") == 0) {
        sut_input = Ib; 
    } else {
        return "ERROR: invalid input string: " + input;
    }

    SUT_Output sut_output = apply(sut_input);

    switch (sut_output) {
        case O0: return "0";
        case O1: return "1";
        default: return "ERROR: invalid output observed";
    }
}