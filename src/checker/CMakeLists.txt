set (FSM_CHECKER_SOURCES
	fsm-checker.cpp
)

add_executable (fsm-checker ${FSM_CHECKER_SOURCES})

target_link_libraries (fsm-checker fsm-fsm fsm-interface fsm-sets fsm-trees jsoncpp)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()
