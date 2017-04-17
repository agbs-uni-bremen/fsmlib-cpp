#include "gdclib.h"

static gdc_states_t myState = up;

static gdc_outputs handleUp(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = moving_down;
            return a1;
            
        default: return nop;
    }
 
}

static gdc_outputs handleMovingDown(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = suspended_down;
            return a3;
            
        case e2: myState = down;
            return a3;
            
        case e4: myState = moving_up;
            return a4;
            
        default: return nop;
    }
}

static gdc_outputs handleSuspendedDown(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = moving_down;
            return a1;
            
        default: return nop;
    }
}

static gdc_outputs handleDown(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = moving_up;
            return a2;
            
        default: return nop;
    }
}

static gdc_outputs handleMovingUp(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = suspended_up;
            return a3;
            
        case e3: myState = up;
            return a3;
            
        default: return nop;
    }

}

static gdc_outputs handleSuspendedUp(gdc_inputs_t x) {
    
    switch ( x ) {
            
        case e1: myState = moving_up;
            return a2;
            
        default: return nop;
    }

}

void gdc_reset() {
    myState = up;
}

gdc_outputs gdc(gdc_inputs_t x) {
        
    switch ( myState ) {
            
        case up: return handleUp(x);
        case moving_down: return handleMovingDown(x);
        case suspended_down: return handleSuspendedDown(x);
        case down: return handleDown(x);
        case moving_up: return handleMovingUp(x);
        case suspended_up: return handleSuspendedUp(x);
            
        default: return nop;
    }
    
}

