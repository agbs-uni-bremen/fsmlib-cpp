set (METHOD_CHECKER_MAIN_SOURCES
	main.cpp
)
set (EQUIVALENCE_TEST_MAIN_SOURCES
	equivalenceTest.cpp
)
set (TEST_MAIN_SOURCES
	dmethodTest.cpp
)

add_executable (method-checker ${METHOD_CHECKER_MAIN_SOURCES})

add_executable (equivalence-test ${EQUIVALENCE_TEST_MAIN_SOURCES})

add_executable (dmethod-test ${TEST_MAIN_SOURCES})

target_link_libraries (method-checker fsm-fsm fsm-interface fsm-sets fsm-trees jsoncpp)

target_link_libraries (equivalence-test fsm-fsm fsm-interface fsm-sets fsm-trees jsoncpp)

target_link_libraries (dmethod-test fsm-fsm fsm-interface fsm-sets fsm-trees jsoncpp)

