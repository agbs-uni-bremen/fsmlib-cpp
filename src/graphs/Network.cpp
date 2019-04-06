/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 *
 * Licensed under the EUPL V.1.1
 */
#include "Network.h"
#include "graphs/NetworkEdge.h"

#include <fstream>
#include <iostream>

Network::Network(const vector<shared_ptr<Node>> &nodes,int sourceNodeId, int sinkNodeId)
    : Graph(nodes),sourceNodeId(sourceNodeId),sinkNodeId(sinkNodeId)
{
    residualNetwork = make_shared<Network>(this);
}


Network::Network(Network *network)
    : Graph(network->nodes),sourceNodeId(network->sourceNodeId),sinkNodeId(network->sinkNodeId),residualNetwork(shared_ptr<Network>())
{
    //It is assumed that network itself does not contain any reverse edges and is an actual network and no residual

    for(auto& node:network->nodes) {
        int edgesSize = node->getEdges().size();
        vector<int> delIdx;
        for(int i=0;i<edgesSize;++i) {
            auto& edge = node->getEdges().at(i);
            auto castedEdge = static_pointer_cast<NetworkEdge>(edge);
            auto targetNode = edge->getTarget().lock();
            //mark the maxed out edges of network
            if(castedEdge->getFlow() == castedEdge->getCapacity()) {
                delIdx.push_back(i);
            }
            //add reverse edges to the residual network
            if(castedEdge->getFlow() > 0) {
                auto reverseEdge = make_shared<NetworkEdge>(vector<int>(),nodes[targetNode->getId()],
                        nodes[node->getId()],castedEdge->getFlow(),castedEdge->getCost() * -1);
                nodes[targetNode->getId()]->addEdge(reverseEdge);
                nodes[node->getId()]->addInEdge(reverseEdge);
                reverseEdge->setIsReverse(true);
                //reverse edges are linked from residual to original network and vice versa
                reverseEdge->setReverseEdge(castedEdge);
                castedEdge->setReverseEdge(reverseEdge);
            }
        }
        //delete the edges in the residual network that are maxed out
        for(int idx:delIdx) {
            auto edges = nodes[node->getId()]->getEdges();
            //TODO: vector erase is not that efficient, is there no better way to remove the maxed out edges?
            edges.erase(edges.begin() + idx);
        }
    }
}

void Network::calculateMinimumCostMaximumFlow() {

    auto shortestPath = residualNetwork->shortestPathByBellmanFord(nodes[sourceNodeId],nodes[sinkNodeId]);
    while(!shortestPath->empty()) {
        //find the minimal absolute value by which the flow can be altered along the path
        //for(auto& edge:)
    }
}

ostream & operator<<(ostream & out, const Network & network) {
    out << "digraph g {" << endl << endl;

    out << "node [shape = ellipse]" << endl;

    for(auto& node:network.nodes) {
        out << node->getId() << "[label=\"" << node->getId() << "\"];" <<  endl;
    }

    for(auto& node:network.nodes) {
        for(auto& edge:node->getEdges()) {
            auto targetNode = edge->getTarget();
            auto castedEdge = static_pointer_cast<NetworkEdge>(edge);
            int capacity = castedEdge->getCapacity(),
                flow = castedEdge->getFlow();

            string traceString = "";
            for(int input:edge->getTrace()) {
                traceString += to_string(input) + ".";
            }
            traceString = traceString.substr(0,traceString.length()-1);

            out << node->getId() << " -> " << targetNode.lock()->getId() << "[label=\"" << traceString << ",(" << to_string(flow) << "," << to_string(capacity) << ")\"];"
                << "  //" << node->getId() << " -> " << targetNode.lock()->getId() << endl;
        }
    }

    out << endl << "}" << endl;
    return out;
}

void Network::toDot(const string &fname) {
    ofstream out(fname + ".dot");
    out << *this;
    out.close();
}


