set (FSM_EXAMPLE_NONDETERMINISTIC_SOURCES
	example.h
	example.cpp
	sut_wrapper.cpp
)

add_library (fsm-example-nondeterministic OBJECT ${FSM_EXAMPLE_NONDETERMINISTIC_SOURCES})
