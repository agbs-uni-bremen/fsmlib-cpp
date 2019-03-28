set (FSM_CHECKER_SOURCES
	fsm-checker.cpp
)

add_executable (fsm-checker ${FSM_CHECKER_SOURCES}
$<TARGET_OBJECTS:fsm-fsm>
$<TARGET_OBJECTS:fsm-interface>
$<TARGET_OBJECTS:fsm-sets>
$<TARGET_OBJECTS:fsm-trees>
$<TARGET_OBJECTS:fsm-graphs>
$<TARGET_OBJECTS:fsm-cloneable>
$<TARGET_OBJECTS:fsm-utils>
)

target_link_libraries (fsm-checker jsoncpp)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()
