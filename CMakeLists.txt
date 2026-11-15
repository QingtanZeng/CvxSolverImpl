# 1. 声明 CMake 最低版本要求
cmake_minimum_required(VERSION 3.23)

# 指定使用 LLVM 18 (clang-18 / clang++-18) 编译器
set(CMAKE_C_COMPILER clang-18)
set(CMAKE_CXX_COMPILER clang++-18)

# 2. 声明项目名称 and 支持的语言
project(CvxSolverImpl VERSION 1.0 LANGUAGES C CXX)

# 设置默认构建类型为 Debug（启用 -g 和 -O0）
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# ==========================================
# 核心指令：强制 CMake 导出 compile_commands.json
# ==========================================
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 3. 全局设置 C++ 标准为 C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # 禁用编译器特定扩展，保证可移植性

# 4. 包含子目录中的求解器库
add_subdirectory(lib_ecos)
# add_subdirectory(lib_scs/scs-master) # 暂时注释，需要时取消

# 5. 查找 Eigen3 依赖
find_package(Eigen3 REQUIRED)

# 6. 编译测试可执行文件
add_executable(ECOSTest01 test_ecos/ECOSTest01.cpp)
target_link_libraries(ECOSTest01 PRIVATE Eigen3::Eigen ecos)
# 强制 ECOSTest01 编译目标使用 -g 和 -O0
# target_compile_options(ECOSTest01 PRIVATE -g -O0)

# 7. 编译提供给 Julia 的共享库接口
add_library(compute SHARED test_JuliaInterface/lib_compute.cpp)
add_library(mtrctrl SHARED test_JuliaInterface/lib_mtrctrl.cpp)

# ==========================================
# 8. SuiteSparse (AMD, LDL) 测试
# ==========================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib_SuiteSparse/CMakeLists.txt")
    set(SUITESPARSE_ENABLE_PROJECTS "suitesparse_config;amd;ldl" CACHE STRING "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(SUITESPARSE_USE_CUDA OFF CACHE BOOL "" FORCE)
    set(SUITESPARSE_USE_FORTRAN OFF CACHE BOOL "" FORCE)
    add_subdirectory(lib_SuiteSparse EXCLUDE_FROM_ALL)

    add_executable(tst_AMD test_SuiteSparse/tst_AMD.cpp)
    set_target_properties(tst_AMD PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON)
    target_link_libraries(tst_AMD PRIVATE SuiteSparse::AMD SuiteSparse::LDL SuiteSparse::SuiteSparseConfig Eigen3::Eigen)

    # ------------------------------------------
    # Host-Target Two-Stage Build (LDLT_PreCalc -> generated_symbolic.hpp -> LDLT_RunSlv)
    # ------------------------------------------
    set(LDLT_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_SuiteSparse/generated")
    set(LDLT_GEN_HEADER "${LDLT_GEN_DIR}/generated_symbolic.hpp")

    add_executable(LDLT_PreCalc test_SuiteSparse/LDLT_PreCalc.cpp)
    set_target_properties(LDLT_PreCalc PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON)
    target_link_libraries(LDLT_PreCalc PRIVATE SuiteSparse::AMD SuiteSparse::LDL SuiteSparse::SuiteSparseConfig)

    add_custom_command(
        OUTPUT "${LDLT_GEN_HEADER}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${LDLT_GEN_DIR}"
        COMMAND $<TARGET_FILE:LDLT_PreCalc> "${LDLT_GEN_HEADER}"
        DEPENDS LDLT_PreCalc
        COMMENT "[Stage 1] Running LDLT_PreCalc to generate ${LDLT_GEN_HEADER}"
        VERBATIM
    )

    add_executable(LDLT_RunSlv test_SuiteSparse/LDLT_RunSlv.cpp "${LDLT_GEN_HEADER}")
    set_target_properties(LDLT_RunSlv PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON)
    target_include_directories(LDLT_RunSlv PRIVATE "${LDLT_GEN_DIR}")
    target_compile_options(LDLT_RunSlv PRIVATE -fno-exceptions -fno-rtti -Wall -Wextra)
endif()


