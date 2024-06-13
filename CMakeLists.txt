cmake_minimum_required (VERSION 3.15)

# Setting VCPKG to easy lib installation
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")

# Configuration of project
project (pOptiPharm CXX)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")
set(EXE_NAME pOPShapeSimilarity)

# Configure CUDA compilation
if(ENABLE_CUDA)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        find_package(CUDA REQUIRED)
        include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

        file(GLOB CUDA_FILES "${CMAKE_SOURCE_DIR}/src/functions/*.cu")

        # Disable Eigen vectorization and alignment for CPU SIMD
        add_compile_definitions(EIGEN_DONT_VECTORIZE) # Vectorization
        add_compile_definitions(EIGEN_MAX_ALIGN_BYTES=0) # Alignment
        
        add_compile_definitions(OP_ENABLE_CUDA) # Add CUDA kernels to code
endif()

# Configuring compilers
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-write-strings)

# Set external libs
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# Set header directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Define source files
set(MENSAJES "${CMAKE_SOURCE_DIR}/src/uego/usrintf/msgintf.cc")
file(GLOB FUNCTIONS "${CMAKE_SOURCE_DIR}/src/functions/*.cpp")
file(GLOB MODELS "${CMAKE_SOURCE_DIR}/src/model/*")
file(GLOB READWRITE "${CMAKE_SOURCE_DIR}/src/read_write/*")
file(GLOB UEGOSHAPESIMILARITY "${CMAKE_SOURCE_DIR}/src/uego/usrintf/configur.cc" "${CMAKE_SOURCE_DIR}/src/uego/*.cc"
        "${CMAKE_SOURCE_DIR}/src/uego/project/molvalSS.cc" "${CMAKE_SOURCE_DIR}/src/uego/project/nrealopt.cc" ${MENSAJES})
# file(GLOB UEGOELECTROSTATIC "src/uego/usrintf/configur.cc" "src/uego/*.cc" "src/uego/project/molvalES.cc" "src/uego/project/nrealopt.cc" ${MENSAJES} )

#Executables
if(ENABLE_CUDA)
        add_executable(${EXE_NAME} "${CMAKE_SOURCE_DIR}/src/uego/usrintf/main.cc"
                ${UEGOSHAPESIMILARITY} ${FUNCTIONS} ${MODELS} ${READWRITE} ${CUDA_FILES})        
        set_property(TARGET ${EXE_NAME} PROPERTY CUDA_SEPARABLE_COMPILATION ON)    
else()
        add_executable(${EXE_NAME} "${CMAKE_SOURCE_DIR}/src/uego/usrintf/main.cc"
                ${UEGOSHAPESIMILARITY} ${FUNCTIONS} ${MODELS} ${READWRITE})   
endif()


target_link_libraries(${EXE_NAME} PRIVATE Threads::Threads)
target_link_libraries(${EXE_NAME} PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(${EXE_NAME} PRIVATE Eigen3::Eigen)

if(ENABLE_CUDA)
target_compile_options(${EXE_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        # --default-stream per-thread        
        -dc
        -diag-suppress=815     # Suppress Eigen warnings on device code
        -diag-suppress=20012
        >)        
endif()

