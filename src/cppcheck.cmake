# additional target to perform cppcheck run, requires cppcheck

# get all project files
# HACK this workaround is required to avoid qml files checking ^_^
file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.h)
file(GLOB_RECURSE EXCLUDE
    "${CMAKE_SOURCE_DIR}/externals/jsoncpp-0.10.0/*")
list(REMOVE_ITEM ALL_SOURCE_FILES ${EXCLUDE})
message("${ALL_SOURCE_FILES}")
#foreach (SOURCE_FILE ${ALL_SOURCE_FILES})
#    string(FIND ${SOURCE_FILE} ${PROJECT_TRDPARTY_DIR} PROJECT_TRDPARTY_DIR_FOUND)
#    if (NOT ${PROJECT_TRDPARTY_DIR_FOUND} EQUAL -1)
#        list(REMOVE_ITEM ALL_SOURCE_FILES ${SOURCE_FILE})
#    endif ()
#endforeach ()

add_custom_target(
        cppcheck
        COMMAND /usr/bin/cppcheck
        -I ${CMAKE_SOURCE_DIR}
#        --check-config
        --enable=warning,performance,portability,information,missingInclude
        --std=c++11
        --library=qt.cfg
        --template="[{severity}][{id}] {message} {callstack} \(On {file}:{line}\)"
        --verbose
        --quiet
        ${ALL_SOURCE_FILES}
)
