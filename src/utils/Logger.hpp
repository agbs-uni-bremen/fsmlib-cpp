#ifndef __FSMLIB_CPP_UTILS_LOGGER_HPP__
#define __FSMLIB_CPP_UTILS_LOGGER_HPP__

#include <functional>
#include <ostream>
#include <map>
#include <sstream>

class LogCoordinator {
public:
    static LogCoordinator &getStandardLogger();

    void setDefaultStream(std::ostream &);

    void createLogTargetAndBind(std::string const &, std::ostream &);
    void createLogTarget(std::string const &);

    //Single stream redirection
    void bindToStream(std::string const &, std::ostream &);
    void bindToDevNull(std::string const &);

    //Global stream redirection
    void bindAllToStream(std::ostream &);
    void bindAllToDevNull();

    std::reference_wrapper<std::ostream> &operator[](std::string const &);
    std::ostream &operator[](std::string const&) const;
protected:
    LogCoordinator();
    LogCoordinator(LogCoordinator&&);
    LogCoordinator(LogCoordinator&) = default;
    
    std::map<std::string, std::reference_wrapper<std::ostream>> streams;
    std::ostringstream devNull;
    std::reference_wrapper<std::ostream> defaultStream;
};

#define LOG_CONST(x) const_cast<LogCoordinator const&>(LogCoordinator::getStandardLogger())[(x)]
#define LOG(x) (LogCoordinator::getStandardLogger()[(x)].get())

#endif //__FSMLIB_CPP_UTILS_LOGGER_HPP__

