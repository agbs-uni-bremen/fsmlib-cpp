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
#include <unordered_set>
#include <deque>

Graph::Graph(const vector<shared_ptr<Node>> &nodes)
    :nodes(nodes)
{

}

shared_ptr<deque<shared_ptr<Edge>>>
Graph::shortestPathByBellmanFord(const shared_ptr<Node> &source, const shared_ptr<Node> &target) {
    if(source->getId() >= nodes.size() || target->getId() >= nodes.size()) {
        return shared_ptr<deque<shared_ptr<Edge>>>();
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
        return shared_ptr<deque<shared_ptr<Edge>>>();
    }

    //check if there are any negative cost cycles
    for(auto& node:nodes) {
        for(auto& edge:node->getEdges()) {
            auto targetNode = edge->getTarget().lock();
            if(isinf(dist[node->getId()])) continue;
            if(isinf(dist[targetNode->getId()]) ||
               dist[targetNode->getId()] > dist[node->getId()] + edge->getCost()) {
                //apparently the graph contains a negative cost cycle in this case
                return shared_ptr<deque<shared_ptr<Edge>>>();
            }
        }
    }
    auto path = make_shared<deque<shared_ptr<Edge>>>();
    auto currentNode = target;
    // path.push_front(pred[target->getId()]);

    while(currentNode != source) {
        //the source node of the edge must be a managed weak_ptr at this point, provided the graph is well formed
        //regarding the pointers
        path->push_front(pred[currentNode->getId()]);
        currentNode = pred[currentNode->getId()]->getSource().lock();
    }

    return path;
}


shared_ptr<list<shared_ptr<Edge>>> Graph::generateEulerTour() {
    auto eulerTour = make_shared<list<shared_ptr<Edge>>>();
    if(nodes.empty()) return eulerTour;

    unordered_set<shared_ptr<Edge>> unmarkedEdges;

    for(auto& node:nodes) {
        auto& edges = node->getEdges();
        auto& inEdges = node->getInEdges();

        //return nothing if the graph is not symmetric
        if(edges.size() != inEdges.size()) {
            return eulerTour;
        }

        for(auto& edge:edges) {
            unmarkedEdges.insert(edge);
        }
    }

    //create the first tour
    auto currentNode = nodes[0];
    do {
        shared_ptr<Node> nextNode;
        for(auto& edge:currentNode->getEdges()) {
            if(unmarkedEdges.count(edge) > 0) {
                unmarkedEdges.erase(edge);
                nextNode = edge->getTarget().lock();
                eulerTour->push_back(edge);
            }
        }
        if(!nextNode) {
            //should not happen, check must be performed though, otherwise you get stuck in the loop (no valid node ids, graph not well-formed etc.)
            return make_shared<list<shared_ptr<Edge>>>();
        }
        currentNode = nextNode;
    } while(currentNode != nodes[0]);


    while(!unmarkedEdges.empty()) {
        shared_ptr<Node> v;
        auto curr = eulerTour->front()->getSource().lock();
        for(auto& e:curr->getEdges()) {
            if(unmarkedEdges.count(e) > 0) {
                v = curr;
                break;
            }
        }
        if(!v) {
            for (auto &edge:*eulerTour) {
                curr = edge->getTarget().lock();
                for(auto& e:curr->getEdges()) {
                    if(unmarkedEdges.count(e) > 0) {
                        v = curr;
                        break;
                    }
                }
                if(v) break;
            }
            auto firstPart = make_shared<list<shared_ptr<Edge>>>();
            auto secondPart = make_shared<list<shared_ptr<Edge>>>();
            for(auto it = eulerTour->begin(); it == eulerTour->end(); ++it) {

            }
        }
        if(!v) return make_shared<list<shared_ptr<Edge>>>();

        auto nextTour = make_shared<list<shared_ptr<Edge>>>();
        //create the first tour
        auto currentNode = v;
        do {
            shared_ptr<Node> nextNode;
            for(auto& edge:currentNode->getEdges()) {
                if(unmarkedEdges.count(edge) > 0) {
                    unmarkedEdges.erase(edge);
                    nextNode = edge->getTarget().lock();
                    nextTour->push_back(edge);
                }
            }
            if(!nextNode) {
                //should not happen, check must be performed though, otherwise you get stuck in the loop (no valid node ids, graph not well-formed etc.)
                return make_shared<list<shared_ptr<Edge>>>();
            }
            currentNode = nextNode;
        } while(currentNode != nodes[0]);

    }

    return eulerTour;

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


