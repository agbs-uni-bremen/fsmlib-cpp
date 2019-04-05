/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/Graph.h"
#include "Graph.h"
#include <fstream>
#include <iostream>

Graph::Graph(const vector<shared_ptr<Node>> &nodes)
    :nodes(nodes)
{

}

ostream & operator<<(ostream & out, const Graph & graph) {
    out << "digraph g {" << endl << endl;

    out << "node [shape = ellipse]" << endl;

    for(auto& node:graph.nodes) {
        out << node->getId() << "[label=\"" << node->getId() << "\"];" <<  endl;
    }

    for(auto& node:graph.nodes) {
        for(auto& edge:node->getEdges()) {
            auto targetNode = edge->getTarget();

            string traceString = "";
            for(int input:edge->getTrace()) {
                traceString += to_string(input) + ".";
            }
            traceString = traceString.substr(0,traceString.length()-1);

            out << node->getId() << " -> " << targetNode.lock()->getId() << "[label=\"" << traceString << "\"];"
                << "  //" << node->getId() << " -> " << targetNode.lock()->getId() << endl;
        }
    }

    out << endl << "}" << endl;
    return out;
}

void Graph::toDot(const string &fname) {
    ofstream out(fname + ".dot");
    out << *this;
    out.close();
}
