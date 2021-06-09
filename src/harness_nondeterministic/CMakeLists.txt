set (FSM_HARNESS_NONDETERMINISTIC_SOURCES
        fsm_test_harness_nondeterministic.cpp
)


add_executable (fsm-harness-nondeterministic ${FSM_HARNESS_NONDETERMINISTIC_SOURCES}
  $<TARGET_OBJECTS:fsm-example-nondeterministic>
)


