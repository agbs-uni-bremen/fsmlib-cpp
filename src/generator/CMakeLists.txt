set (FSM_GENERATOR_SOURCES
	fsm-test-generator.cpp
  generalized-h-method.hpp
)

add_executable (fsm-generator ${FSM_GENERATOR_SOURCES}
$<TARGET_OBJECTS:fsm-fsm>
$<TARGET_OBJECTS:fsm-interface>
$<TARGET_OBJECTS:fsm-sets>
$<TARGET_OBJECTS:fsm-trees>
$<TARGET_OBJECTS:fsm-cloneable>
$<TARGET_OBJECTS:fsm-utils>
)

target_link_libraries (fsm-generator jsoncpp)

#if(MSVC)
#	set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#endif()
