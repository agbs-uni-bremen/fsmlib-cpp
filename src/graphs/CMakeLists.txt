set (FSM_GRAPHS_SOURCES
        Node.cpp
        Node.h
        Edge.cpp
        Edge.h
        Graph.cpp
        Graph.h
        PartitionGraphNode.cpp
        PartitionGraphNode.h
        PartitionGraph.cpp
        PartitionGraph.h
        PartitionGraphEdge.cpp
        PartitionGraphEdge.h)

add_library (fsm-graphs OBJECT ${FSM_GRAPHS_SOURCES})
