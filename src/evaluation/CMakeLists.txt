set (METHODS_EVALUATION_MAIN_SOURCES
		methods-evaluation.cpp
)
set (HMETHOD_EVALUATION_SOURCES
		hMethod-evaluation.cpp
		)

set (EQUIVALENCE_TEST_MAIN_SOURCES
	equivalenceTest.cpp
)
set (TEST_DMETHOD_SOURCES
	dmethodTest.cpp
)
set (TEST_ADSMETHOD_SOURCES
		adsMethodTest.cpp
		)


add_executable (methods-evaluation ${METHODS_EVALUATION_MAIN_SOURCES}
    	$<TARGET_OBJECTS:fsm-fsm>
    	$<TARGET_OBJECTS:fsm-interface>
    	$<TARGET_OBJECTS:fsm-sets>
    	$<TARGET_OBJECTS:fsm-trees>
		$<TARGET_OBJECTS:fsm-graphs>
		$<TARGET_OBJECTS:fsm-utils>
		$<TARGET_OBJECTS:fsm-cloneable>
		)

add_executable (equivalence-test ${EQUIVALENCE_TEST_MAIN_SOURCES}
		$<TARGET_OBJECTS:fsm-fsm>
    	$<TARGET_OBJECTS:fsm-interface>
    	$<TARGET_OBJECTS:fsm-sets>
    	$<TARGET_OBJECTS:fsm-trees>
		$<TARGET_OBJECTS:fsm-graphs>
		$<TARGET_OBJECTS:fsm-utils>
		$<TARGET_OBJECTS:fsm-cloneable>
		)

add_executable (dmethod-test ${TEST_DMETHOD_SOURCES}
    	$<TARGET_OBJECTS:fsm-fsm>
		$<TARGET_OBJECTS:fsm-interface>
    	$<TARGET_OBJECTS:fsm-sets>
    	$<TARGET_OBJECTS:fsm-trees>
		$<TARGET_OBJECTS:fsm-graphs>
		$<TARGET_OBJECTS:fsm-utils>
		$<TARGET_OBJECTS:fsm-cloneable>
		)

add_executable (adsMethod-test ${TEST_ADSMETHOD_SOURCES}
		$<TARGET_OBJECTS:fsm-fsm>
		$<TARGET_OBJECTS:fsm-interface>
		$<TARGET_OBJECTS:fsm-sets>
		$<TARGET_OBJECTS:fsm-trees>
		$<TARGET_OBJECTS:fsm-graphs>
		$<TARGET_OBJECTS:fsm-utils>
		$<TARGET_OBJECTS:fsm-cloneable>
		)

add_executable (hmethod-evaluation ${HMETHOD_EVALUATION_SOURCES}
		$<TARGET_OBJECTS:fsm-fsm>
		$<TARGET_OBJECTS:fsm-interface>
		$<TARGET_OBJECTS:fsm-sets>
		$<TARGET_OBJECTS:fsm-trees>
		$<TARGET_OBJECTS:fsm-graphs>
		$<TARGET_OBJECTS:fsm-utils>
		$<TARGET_OBJECTS:fsm-cloneable>
		)

target_link_libraries (methods-evaluation jsoncpp)

target_link_libraries (hmethod-evaluation jsoncpp)

target_link_libraries (equivalence-test jsoncpp)

target_link_libraries (dmethod-test jsoncpp)

target_link_libraries (adsMethod-test jsoncpp)



