#ifndef SUT_WRAPPER_ADAPTIVE_H_
#define SUT_WRAPPER_ADAPTIVE_H_

#include <string>

void sut_init();

void sut_reset();

//const char* sut(const char* input);
const std::string sut(const std::string input);


#endif // SUT_WRAPPER_ADAPTIVE_H_