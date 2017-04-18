set (FSM_GENERATOR_SOURCES
	fsm-test-generator.cpp
)

add_executable (fsm-generator ${FSM_GENERATOR_SOURCES})

target_link_libraries (fsm-generator fsm-fsm fsm-interface fsm-sets fsm-trees)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()