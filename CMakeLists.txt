# 1. 声明 CMake 最低版本要求
cmake_minimum_required(VERSION 3.15)

# 2. 声明项目名称和支持的语言
project(OptSolverTest VERSION 1.0 LANGUAGES C CXX)

# ==========================================
# 核心指令：强制 CMake 导出 compile_commands.json
# ==========================================
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 3. 全局设置 C++ 标准 (这里以 C++17 为例)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # 禁用编译器特定扩展，保证可移植性

# 4. 包含子目录中的求解器库
add_subdirectory(lib_ecos)
# add_subdirectory(lib_scs/scs-master) # 暂时注释，需要时取消

# 5. 查找 Eigen3 依赖
find_package(Eigen3 REQUIRED)

# 6. 编译测试可执行文件
add_executable(ECOSTest00 test_ecos/ECOSTest00.cpp)
target_link_libraries(ECOSTest00 PRIVATE Eigen3::Eigen ecos)

# 7. 编译提供给 Julia 的共享库接口
add_library(compute SHARED test_JuliaInterface/lib_compute.cpp)
add_library(mtrctrl SHARED test_JuliaInterface/lib_mtrctrl.cpp)