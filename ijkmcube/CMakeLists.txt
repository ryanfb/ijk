PROJECT(IJKMCUBE)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMCUBE_DIR "src/ijkmcube")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")
SET(IJK_ISOTABLE_DIR "${IJK_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(expat ITKNrrdIO z)
ADD_DEFINITIONS(-DIJK_ISOTABLE_DIR=\\\"${IJK_ISOTABLE_DIR}\\\")

ADD_EXECUTABLE(ijkmcube ijkmcube_main.cxx ijkmcube.cxx ijktable.cxx ijkoctree.cxx ijkxitIO.cxx)

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS ijkmcube DESTINATION "bin")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKMCUBE_DIR}/ijkmcube.tar ${IJKMCUBE_DIR}/README ${IJKMCUBE_DIR}/INSTALL ${IJKMCUBE_DIR}/CHANGES ${IJKMCUBE_DIR}/*.cxx ${IJKMCUBE_DIR}/*.h ${IJKMCUBE_DIR}/*.txx ${IJKMCUBE_DIR}/CMakeLists.txt ${IJKMCUBE_DIR}/man/*)

