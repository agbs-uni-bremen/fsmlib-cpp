set (MANAGER_MOC_FILE
	MainWindow.h
	OpenFileWindow.h
)

QT5_WRAP_CPP (MANAGER_MOC ${MANAGER_MOC_FILE})

QT5_WRAP_UI (ui_MainWindow.h MainWindow.ui)

set (FSM_WINDOW_SOURCES
	MainWindow.cpp
	MainWindow.h
	OpenFileWindow.cpp
	OpenFileWindow.h
)

add_library (fsm-window OBJECT ${FSM_WINDOW_SOURCES} ${MANAGER_MOC} MainWindow.ui)

target_link_libraries (fsm-window Qt5::Widgets)


set (FSM_MAIN_SOURCES
	gui.cpp
)

add_executable (fsm-gui OBJECT ${FSM_MAIN_SOURCES})

target_link_libraries (fsm-gui fsm-window fsm-fsm fsm-interface fsm-sets fsm-trees Qt5::Widgets ${QT_QTMAIN_LIBRARY})

if(MSVC)
	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
endif()
