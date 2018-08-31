#include "utils/Logger.hpp"
#include <iostream>

LogCoordinator &LogCoordinator::getStandardLogger() {
    static LogCoordinator standardLogger = LogCoordinator();
    return standardLogger;
}

LogCoordinator::LogCoordinator() : defaultStream(devNull) {
    this->devNull.setstate(std::ios::badbit);
}

LogCoordinator::LogCoordinator(LogCoordinator &&other) : defaultStream(devNull) {
    this->devNull.setstate(std::ios::badbit);
    std::swap(this->streams, other.streams);
}

void LogCoordinator::createLogTargetAndBind(std::string const &name, std::ostream &stream) {
    if(this->streams.count(name) == 0) {
        this->streams.emplace(std::piecewise_construct, std::tuple<std::string const &>(name), std::tuple<std::ostream &>(stream));
    }
}

void LogCoordinator::createLogTarget(std::string const &name) {
    this->createLogTargetAndBind(name, this->defaultStream);
}

void LogCoordinator::setDefaultStream(std::ostream &stream) {
    this->defaultStream = stream;
}

void LogCoordinator::bindToStream(std::string const &name, std::ostream &stream) {
    this->operator[](name) = stream;
}

void LogCoordinator::bindToDevNull(std::string const &name) {
    this->bindToStream(name, this->devNull);
}

void LogCoordinator::bindAllToStream(std::ostream &stream) {
    this->setDefaultStream(stream);
    for(auto &kvp : this->streams) {
        this->bindToStream(kvp.first, stream);
    }
}

void LogCoordinator::bindAllToDevNull() {
    this->bindAllToStream(this->devNull);
}

std::reference_wrapper<std::ostream> &LogCoordinator::operator[](std::string const &name) {
    if(this->streams.count(name) == 0) {
        this->createLogTarget(name);
    }
    return this->streams.at(name);
}

std::ostream &LogCoordinator::operator[](std::string const &name) const {
    return this->streams.at(name).get();
}

