set (FSM_SETS_SOURCES
	HittingSet.cpp
	HittingSet.h
	HsTreeNode.cpp
	HsTreeNode.h
)

add_library (fsm-sets OBJECT ${FSM_SETS_SOURCES})
