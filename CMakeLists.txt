# 1. 声明 CMake 最低版本要求
cmake_minimum_required(VERSION 3.15)

# 指定使用 clang/clang++ 编译器
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

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

# 3. 全局设置 C++ 标准 (这里以 C++20 为例)
set(CMAKE_CXX_STANDARD 20)
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

