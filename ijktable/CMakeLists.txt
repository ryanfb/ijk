PROJECT(IJKTABLE)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKTABLE_DIR "src/ijktable")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

LINK_LIBRARIES(expat)
ADD_EXECUTABLE(ijkgentable ijkgentable.cxx ijkgenpatch.cxx ijktable.cxx ijkxitIO.cxx clarkson_hull.c)
ADD_EXECUTABLE(ijkdifftable ijkdifftable.cxx ijktable.cxx ijkxitIO.cxx)
ADD_EXECUTABLE(ijktableinfo ijktableinfo.cxx ijktable.cxx ijkxitIO.cxx)
ADD_EXECUTABLE(ijktext2c ijktext2c.cxx)

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS ijkgentable ijkdifftable ijktableinfo DESTINATION "bin")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKTABLE_DIR}/ijktable.tar ${IJKTABLE_DIR}/README ${IJKTABLE_DIR}/INSTALL ${IJKTABLE_DIR}/*.cxx ${IJKTABLE_DIR}/*.h ${IJKTABLE_DIR}/*.txx ${IJKTABLE_DIR}/*.c ${IJKTABLE_DIR}/CMakeLists.txt ${IJKTABLE_DIR}/ijktable.xsd ${IJKTABLE_DIR}/man/* )

