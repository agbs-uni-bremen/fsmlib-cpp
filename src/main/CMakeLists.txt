set (FSM_MAIN_SOURCES
	main.cpp
)

add_executable (fsm-main ${FSM_MAIN_SOURCES}
    $<TARGET_OBJECTS:fsm-fsm>
    $<TARGET_OBJECTS:fsm-interface>
    $<TARGET_OBJECTS:fsm-sets>
    $<TARGET_OBJECTS:fsm-trees>
	$<TARGET_OBJECTS:fsm-graphs>
    $<TARGET_OBJECTS:fsm-cloneable>
    $<TARGET_OBJECTS:fsm-utils>
)

target_link_libraries (fsm-main jsoncpp)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()
