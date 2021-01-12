cmake_minimum_required(VERSION 3.1)
project(libRandomForest_ZIA)

#确定需要的头文件

include_directories(

    include

)




#确定编译语言

#set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++11)
#设置二进制文件目录

SET(BIN_DESTINATION ${PROJECT_SOURCE_DIR})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DESTINATION})

#如果多个文件相关，则添加源代码（可选）


set(SRC 

    ${PROJECT_SOURCE_DIR}/src/pairnode.cpp

    ${PROJECT_SOURCE_DIR}/src/pairtree.cpp

    ${PROJECT_SOURCE_DIR}/src/pairforest.cpp
    
    ${PROJECT_SOURCE_DIR}/src/forest.cpp

    ${PROJECT_SOURCE_DIR}/src/tree.cpp

    ${PROJECT_SOURCE_DIR}/src/data.cpp

    ${PROJECT_SOURCE_DIR}/src/node.cpp

    ${PROJECT_SOURCE_DIR}/src/utilities.cpp

    ${PROJECT_SOURCE_DIR}/src/nodegini.cpp

    ${PROJECT_SOURCE_DIR}/src/nodeinfogain.cpp

    ${PROJECT_SOURCE_DIR}/src/nodehyperplane.cpp

    #${PROJECT_SOURCE_DIR}/src/RandomForest.cpp


)











#link_directories(${PROJECT_SOURCE_DIR})

#add_library(overload SHARED ${SRC})

#target_link_libraries(overload  -lpthread -lm -lstdc++)

#生成可执行文件
link_directories(${PROJECT_SOURCE_DIR})
#link_libraries(SRC)



#link_directories(${PROJECT_SOURCE_DIR})
#link_libraries(SRC)
AUX_SOURCE_DIRECTORY(./src DIR_SRCS)
#AUX_SOURCE_DIRECTORY(SRC)
#add_executable(RandomForest1 src/RandomForest.cpp ${SRC})

IF(USE_CUDA)
  CUDA_INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/cuda
  )

#  CUDA_ADD_EXECUTABLE(RandomForest
#    RandomForest.cpp
#    cuda/icg_rf_interface.cu
#    ${ICG_RF_SRC}
#  )
  CUDA_ADD_LIBRARY(RandomForestLib
    cuda/icg_rf_interface.cu
    ${SRC}
  )
ELSE(USE_CUDA)
  ADD_EXECUTABLE(RandomForest1
    RandomForest.cpp
    ${SRC}
  )
#  ADD_LIBRARY(RandomForestLib
#    ${ICG_RF_SRC}
#  )
ENDIF(USE_CUDA)



#add_executable(overload src/RandomForest.cpp)
#add_executable(overload ${SRC})
#target_link_libraries(RandomForest1 -libconfig++ -libxml2)
target_link_libraries(RandomForest1 -lpthread -lm -lstdc++)

IF(WIN32)
#TARGET_LINK_LIBRARIES(RandomForest
#  libconfig++ iconv libxml2
#)
TARGET_LINK_LIBRARIES(RandomForestLib
  libconfig++ iconv libxml2
)
ELSE(WIN32)
TARGET_LINK_LIBRARIES(RandomForest1
   xml2 config++ gomp)
ENDIF(WIN32)





