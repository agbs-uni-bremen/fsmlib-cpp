set(JSONCPP_SOURCES 
	jsoncpp.cpp
	../include/json/json.h
	../include/json/json-forwards.h
)

add_library( jsoncpp ${JSONCPP_SOURCES} )
