#include "gdclib.h"



/**
 *   Helper data structures and functions for
 *   SUT test wrapper
 */

static const char* outputs[5];
static const char* inputs[4];

static gdc_inputs_t inStr2Enum(const char* input) {
    
    int i;

    for ( i = 0; i < 4; i++ ) {
        if ( strcmp(inputs[i],input) == 0 ) {
            return (gdc_inputs_t)i;
        }
    }
    
    return e1;
}






void sut_init() {
    
    outputs[nop] = strdup("_nop");
    outputs[1] = strdup("a1");
    outputs[2] = strdup("a2");
    outputs[3] = strdup("a3");
    outputs[4] = strdup("a4");
    
    inputs[0] = strdup("e1");
    inputs[1] = strdup("e2");
    inputs[2] = strdup("e3");
    inputs[3] = strdup("e4");
    
    
}

void sut_reset() {
    gdc_reset();
}


const char* sut(const char* input) {
    
    return outputs[gdc(inStr2Enum(input))];
    
}
