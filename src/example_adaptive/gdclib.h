#include <string.h>


typedef enum {
    up,
    moving_down,
    suspended_down,
    down,
    moving_up,
    suspended_up
} gdc_states_t;

typedef enum {
    e1,
    e2,
    e3,
    e4
} gdc_inputs_t;

typedef enum {
    nop,
    a1,
    a2,
    a3,
    a4
} gdc_outputs;

extern void gdc_reset();
extern gdc_outputs gdc(gdc_inputs_t x);
