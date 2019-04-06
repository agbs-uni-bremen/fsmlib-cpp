/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */

#include "graphs/Graph.h"
#include "Graph.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <deque>

Graph::Graph(const vector<shared_ptr<Node>> &nodes)
    :nodes(nodes)
{

}

shared_ptr<vector<shared_ptr<Edge>>>
Graph::shortestPathByBellmanFord(const shared_ptr<Node> &source, const shared_ptr<Node> &target) {
    if(source->getId() >= nodes.size() || target->getId() >= nodes.size()) {
        return shared_ptr<vector<shared_ptr<Edge>>>();
    }

    vector<long double> dist;
    vector<shared_ptr<Edge>> pred;
    for(auto& node:nodes) {
        dist.push_back(numeric_limits<long double>::infinity());
        pred.push_back(nullptr);
    }
    dist[source->getId()] = 0;

    //compute distances and predecessors
    for(int i=0;i<nodes.size()-1;++i) {
        for(auto& node:nodes) {
            for(auto& edge:node->getEdges()) {
                auto targetNode = edge->getTarget().lock();
                if(isinf(dist[node->getId()])) continue;
                if(isinf(dist[targetNode->getId()]) ||
                        dist[targetNode->getId()] > dist[node->getId()] + edge->getCost()) {
                    dist[targetNode->getId()] = dist[node->getId()] + edge->getCost();
                    pred[targetNode->getId()] = edge;
                }
            }
        }
    }
    //check if the target is even reachable
    if(isinf(dist[target->getId()])) {
        return shared_ptr<vector<shared_ptr<Edge>>>();
    }

    //check if there are any negative cost cycles
    for(auto& node:nodes) {
        for(auto& edge:node->getEdges()) {
            auto targetNode = edge->getTarget().lock();
            if(isinf(dist[node->getId()])) continue;
            if(isinf(dist[targetNode->getId()]) ||
               dist[targetNode->getId()] > dist[node->getId()] + edge->getCost()) {
                //apparently the graph contains a negative cost cycle in this case
                return shared_ptr<vector<shared_ptr<Edge>>>();
            }
        }
    }
    deque<shared_ptr<Edge>> path;
    auto currentNode = target;
    path.push_front(pred[target->getId()]);

    while(currentNode != source) {
        //the source node of the edge must be a managed weak_ptr at this point, provided the graph is well formed
        //regarding the pointers
        currentNode = pred[currentNode->getId()]->getSource().lock();
        path.push_front(pred[currentNode->getId()]);
    }

    return make_shared<vector<shared_ptr<Edge>>>(path.begin(),path.end());
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

bool Graph::validateNodeIds() const {
    for(auto& node:nodes) {
        if(node->getId() < nodes.size() && nodes[node->getId()]->getId() == node->getId()) continue;
    }
    return true;
}

