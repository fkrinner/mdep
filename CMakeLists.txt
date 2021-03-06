cmake_minimum_required (VERSION 2.6)

project (Mdep)

set(create_python_links ON)
set(use_auto_diff ON) 



set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )



# Set directories

set( INCLUDE_DIRS
    ${PROJECT_BINARY_DIR} 
    ${PROJECT_SOURCE_DIR} 
) 

if( use_auto_diff )
  add_definitions( -DADOL_ON )
  set( INCLUDE_DIRS
       ${INCLUDE_DIRS}
  )
endif()

if(create_python_links)
  find_package(Boost 1.45.0 COMPONENTS python)
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
  set( INCLUDE_DIRS 
     ${INCLUDE_DIRS} 
     ${Boost_INCLUDE_DIRS} 	
  )
endif()


include_directories(
     ${INCLUDE_DIRS} 
)

# Set libraries

add_library(wave SHARED wave.cxx )
add_library(waveset SHARED waveset.cxx )
add_library(chi2 SHARED chi2.cxx)
add_library(chi2coma SHARED chi2coma.cxx)
add_library(chi2amp SHARED chi2amp.cxx)


if(create_python_links)
  add_library(chi2py SHARED chi2py.cxx)
  add_library(chi2comapy SHARED chi2comapy.cxx)
  add_library(chi2amppy SHARED chi2amppy.cxx)
  add_library(bwpy SHARED bwpy.cxx)
  add_library(pspy SHARED pspy.cxx)
endif()

target_link_libraries(waveset wave)
target_link_libraries(chi2 waveset)
target_link_libraries(chi2amp waveset)
target_link_libraries(chi2coma chi2)

if(create_python_links)
  target_link_libraries(chi2py chi2)
  target_link_libraries(chi2py ${Boost_LIBRARIES})
  target_link_libraries(chi2amppy chi2amp)
  target_link_libraries(chi2amppy ${Boost_LIBRARIES})
  target_link_libraries(chi2comapy chi2coma)
  target_link_libraries(chi2comapy ${Boost_LIBRARIES})
  target_link_libraries(bwpy ${Boost_LIBRARIES})
  target_link_libraries(pspy ${Boost_LIBRARIES})
endif()
if(use_auto_diff) # The links to the auto-diff package are hardcoded here
  target_link_libraries(chi2 /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(chi2amp /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(chi2coma /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(wave /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
endif()



