set (FSM_EXAMPLE_ADAPTIVE_SOURCES
	gdclib.c
        gdclib.h
        sut_wrapper.c
)

add_library (fsm-example-adaptive OBJECT ${FSM_EXAMPLE_ADAPTIVE_SOURCES})
