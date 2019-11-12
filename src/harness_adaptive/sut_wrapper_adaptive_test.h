#ifndef SUT_WRAPPER_ADAPTIVE_TEST_H_
#define SUT_WRAPPER_ADAPTIVE_TEST_H_

#include <string>
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <cstring>
#include <utility>

#include "utils/Logger.hpp"

void sut_init(const std::shared_ptr<Fsm>& spec, bool mutate);

void sut_reset();

//const char* sut(const char* input);
const std::string sut(const std::string input);


#endif // SUT_WRAPPER_ADAPTIVE_H_