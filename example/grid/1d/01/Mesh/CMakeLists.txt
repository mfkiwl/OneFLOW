cmake_minimum_required(VERSION 3.31)

project(Mesh VERSION 1.0)

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

list ( APPEND PRJ_LIBRARIES Qt6::Widgets )
list ( APPEND PRJ_COMPILE_FEATURES cxx_std_23 )

set ( CGNS_INCLUDE_DIRS $ENV{CGNS_INC} )
set ( CGNS_LIBRARIES $ENV{CGNS_LIB_SHARED_NAME} )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    set ( CGNS_ROOT "C:/dev/cgns/4.4.0" )
    set ( CGNS_LIBNAME "cgnsdll.lib" )

    set ( CGNS_INCLUDE_DIRS "${CGNS_ROOT}/include" CACHE PATH "path to CGNS headers" )
    set ( CGNS_LIBRARIES "${CGNS_ROOT}/lib/${CGNS_LIBNAME}" CACHE PATH "path to CGNS library" )
endif()

list ( APPEND PRJ_LIBRARIES ${CGNS_LIBRARIES} )
list ( APPEND PRJ_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} )

if ( WIN32 )
    list ( APPEND PRJ_COMPILE_DEFINITIONS USE_DLL )
endif ()

message( STATUS "PRJ_INCLUDE_DIRS = ${PRJ_INCLUDE_DIRS}")

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

set( PROJECT_SOURCES
    main.cpp
	grid.h grid.cpp
    mainwindow.h mainwindow.cpp mainwindow.ui
	meshdialog.h meshdialog.cpp meshdialog.ui
)

set( PRJ_SOURCES )
foreach(source ${PROJECT_SOURCES})
  #list( APPEND PRJ_SOURCES "codes/${source}" )
  list( APPEND PRJ_SOURCES "${source}" )
endforeach()

foreach(source ${PRJ_SOURCES})
  #message( STATUS "source=${source}" )
endforeach()

add_executable( ${PROJECT_NAME}
    #images/res.qrc
)

target_sources( ${PROJECT_NAME}
  PRIVATE
	${PRJ_SOURCES}
)

target_include_directories ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_INCLUDE_DIRS}
)

target_link_libraries( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_COMPILE_DEFINITIONS} 
)
