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

add_library (fsm-window ${FSM_WINDOW_SOURCES} ${MANAGER_MOC} MainWindow.ui)

target_link_libraries (fsm-window Qt5::Widgets)