set (FSM_HARNESS_SOURCES
    fsm-test-harness-adaptive.cpp
    sut_wrapper_adaptive.h
    sut_wrapper_adaptive_example.cpp
    #sut_wrapper_adaptive_test.h
    #sut_wrapper_adaptive_test.cpp
    vPrimeEnumerator.cpp
    vPrimeEnumerator.h
)


add_executable (fsm-harness-adaptive ${FSM_HARNESS_SOURCES}
  $<TARGET_OBJECTS:fsm-fsm>
  $<TARGET_OBJECTS:fsm-interface>
  $<TARGET_OBJECTS:fsm-sets>
  $<TARGET_OBJECTS:fsm-trees>
  $<TARGET_OBJECTS:fsm-cloneable>
  $<TARGET_OBJECTS:fsm-utils>
)

target_link_libraries (fsm-harness-adaptive jsoncpp)


