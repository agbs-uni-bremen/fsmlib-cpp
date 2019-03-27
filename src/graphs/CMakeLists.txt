set (FSM_GRAPHS_SOURCES
        Node.cpp
        Node.h
        Edge.cpp
        Edge.h
        Graph.cpp
        Graph.h PartitionGraphNode.cpp PartitionGraphNode.h PartitionGraph.cpp PartitionGraph.h)

add_library (fsm-graphs ${FSM_GRAPHS_SOURCES})
