/**
 * A simple example SUT representing a state machine with four states.
 * 
 * The SUT almost conforms to the FSM given in resources/harness_adaptive/hierons_fsm.fsm,
 * but exhibits a small chance in applying input a to state S1 that the (S1,a,1,S4) transition 
 * is erroneously modified to reach S2 instead, constituting a transition fault.
 */

#ifndef FSM_EXAMPLE_NONDETERMINISTIC_EXAMPLE_H_
#define FSM_EXAMPLE_NONDETERMINISTIC_EXAMPLE_H_

// Type of states of the SUT
enum SUT_State {
    S1,
    S2,
    S3,
    S4
};

// Type of inputs of the SUT
enum SUT_Input {
    Ia,
    Ib
};

// Type of outputs of the SUT
enum SUT_Output {
    O0,
    O1
};

extern void reset();

extern void init();

extern const SUT_Output apply(const SUT_Input input);

#endif //FSM_EXAMPLE_NONDETERMINISTIC_EXAMPLE_H_
