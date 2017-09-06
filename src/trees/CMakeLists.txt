set (FSM_TREES_SOURCES
        AdaptiveTreeNode.cpp
        AdaptiveTreeNode.h
        InputOutputTree.cpp
        InputOutputTree.h
	IOListContainer.cpp
	IOListContainer.h
        IOTreeContainer.cpp
        IOTreeContainer.h
	OutputTree.cpp
	OutputTree.h
	TestSuite.cpp
	TestSuite.h
	Tree.cpp
	Tree.h
	TreeEdge.cpp
	TreeEdge.h
	TreeNode.cpp
	TreeNode.h
)

add_library (fsm-trees ${FSM_TREES_SOURCES})

target_link_libraries(fsm-trees fsm-cloneable)
