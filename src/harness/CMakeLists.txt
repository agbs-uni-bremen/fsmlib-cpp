set (FSM_HARNESS_SOURCES
        fsm-test-harness.c
)


add_executable (fsm-harness ${FSM_HARNESS_SOURCES}
  $<TARGET_OBJECTS:fsm-example>
)


