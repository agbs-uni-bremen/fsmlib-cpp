set (FSM_HARNESS_SOURCES
        fsm-test-harness-adaptive.cpp
)


add_executable (fsm-harness-adaptive ${FSM_HARNESS_SOURCES}
  $<TARGET_OBJECTS:fsm-example-adaptive>
)


