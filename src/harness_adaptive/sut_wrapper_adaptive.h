/** 
 * This header provides functions required to access the SUT during testing.
 */


#ifndef SUT_WRAPPER_ADAPTIVE_H_
#define SUT_WRAPPER_ADAPTIVE_H_

#include <string>

/**
 * Initialise the SUT.
 */
void sut_init();

/**
 * Reset the SUT to its initial state.
 */
void sut_reset();

/**
 * Apply an input to the SUT and return the observed response.
 */
const std::string sut(const std::string input);


#endif // SUT_WRAPPER_ADAPTIVE_H_