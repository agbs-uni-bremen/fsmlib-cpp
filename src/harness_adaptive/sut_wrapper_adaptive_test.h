/** 
 * This example wrapper simulates an SUT by using a reduction or a mutant of the specification FSM.
 */


#ifndef SUT_WRAPPER_ADAPTIVE_TEST_H_
#define SUT_WRAPPER_ADAPTIVE_TEST_H_

#include <string>
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmLabel.h"
#include "interface/FsmPresentationLayer.h"
#include "fsm/FsmTransition.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <utility>

#include "utils/Logger.hpp"

/**
 * Create an SUT FSM.
 * 
 * @param spec Specification FSM
 * @param mutate If true, then the SUT is a mutant of the specification FSM.
 *               Otherwise a reduction of the the specification FSM is used as SUT.
 */
void sut_init(const std::shared_ptr<Fsm>& spec, bool mutate);

/**
 * Set the SUT FSM to its initial state.
 */
void sut_reset();

/**
 * Apply an input to the SUT FSM.
 */
const std::string sut(const std::string input);


#endif // SUT_WRAPPER_ADAPTIVE_H_