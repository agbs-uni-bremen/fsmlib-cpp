set (FSM_MAIN_SOURCES
	main.cpp
)

add_executable (fsm-main ${FSM_MAIN_SOURCES})

target_link_libraries (fsm-main fsm-window fsm-fsm fsm-interface fsm-sets fsm-trees Qt5::Widgets ${QT_QTMAIN_LIBRARY})

if(MSVC)
	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
endif()