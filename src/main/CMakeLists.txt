set (FSM_MAIN_SOURCES
	main.cpp
)

add_executable (fsm-main ${FSM_MAIN_SOURCES})

target_link_libraries (fsm-main fsm-fsm fsm-interface fsm-sets fsm-trees)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()