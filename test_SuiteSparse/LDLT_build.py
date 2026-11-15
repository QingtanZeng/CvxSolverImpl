#!/usr/bin/env python3
"""
Host-Target Two-Stage Build 编排脚本：
  Stage 1 (Host):   使用宿主 Clang C++23 工具链编译并执行 LDLT_PreCalc.cpp (链接 lib_SuiteSparse)，
                    完成 AMD 降填充排序与静态消去树符号分析，产出 generated_symbolic.hpp。
  Stage 2 (Target): 使用目标端 Clang C++23 工具链编译 LDLT_RunSlv.cpp (-fno-exceptions -fno-rtti)，
                    零 SuiteSparse 依赖、零堆内存分配，固化静态符号表并验证数值分解与求解。
"""

import shutil
import subprocess
import sys
from pathlib import Path

# 路径配置
ROOT_DIR = Path(__file__).resolve().parent
WORKSPACE_DIR = ROOT_DIR.parent
SUITESPARSE_SRC_DIR = WORKSPACE_DIR / "lib_SuiteSparse"
CMAKE_BUILD_DIR = WORKSPACE_DIR / "build"
SUITESPARSE_LIB_DIR = CMAKE_BUILD_DIR / "lib_SuiteSparse"

BUILD_HOST_DIR = ROOT_DIR / "build_host"
BUILD_TARGET_DIR = ROOT_DIR / "build_target"
GEN_DIR = ROOT_DIR / "build_artifacts"
GEN_HEADER = GEN_DIR / "generated_symbolic.hpp"

HOST_SRC = ROOT_DIR / "LDLT_PreCalc.cpp"
TARGET_SRC = ROOT_DIR / "LDLT_RunSlv.cpp"

# 编译器与参数配置 (优先使用 clang++-18)
CXX_BIN = shutil.which("clang++-18") or shutil.which("clang++") or "clang++"
HOST_CXX = CXX_BIN
TARGET_CXX = CXX_BIN  # 交叉编译时可替换为目标工具链或追加 --target=...


def run_cmd(cmd: list[str], cwd: Path = ROOT_DIR) -> None:
    print(f"[Orchestrator] Executing: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd, text=True)
    if result.returncode != 0:
        sys.exit(f"[Error] Command failed with exit code {result.returncode}")


def ensure_suitesparse_built() -> None:
    amd_lib = SUITESPARSE_LIB_DIR / "AMD" / "libamd.so"
    ldl_lib = SUITESPARSE_LIB_DIR / "LDL" / "libldl.so"
    cfg_lib = SUITESPARSE_LIB_DIR / "SuiteSparse_config" / "libsuitesparseconfig.so"

    if amd_lib.exists() and ldl_lib.exists() and cfg_lib.exists():
        return

    print("\n=== Prerequisite: Building SuiteSparse (AMD, LDL, SuiteSparseConfig) via CMake ===")
    if not (CMAKE_BUILD_DIR / "CMakeCache.txt").exists():
        run_cmd(["cmake", "-S", str(WORKSPACE_DIR), "-B", str(CMAKE_BUILD_DIR)], cwd=WORKSPACE_DIR)
    run_cmd(
        ["cmake", "--build", str(CMAKE_BUILD_DIR), "--target", "AMD", "LDL", "SuiteSparseConfig", "-j"],
        cwd=WORKSPACE_DIR,
    )


def main() -> None:
    BUILD_HOST_DIR.mkdir(parents=True, exist_ok=True)
    BUILD_TARGET_DIR.mkdir(parents=True, exist_ok=True)
    GEN_DIR.mkdir(parents=True, exist_ok=True)

    ensure_suitesparse_built()

    amd_lib_dir = SUITESPARSE_LIB_DIR / "AMD"
    ldl_lib_dir = SUITESPARSE_LIB_DIR / "LDL"
    cfg_lib_dir = SUITESPARSE_LIB_DIR / "SuiteSparse_config"

    # ---------------------------------------------------------
    # Stage 1: 在宿主机编译并运行 LDLT_PreCalc
    # ---------------------------------------------------------
    print("\n=== Stage 1: Building & Running Host Pre-Calculator (LDLT_PreCalc) ===")
    precalc_bin = BUILD_HOST_DIR / "LDLT_PreCalc"

    host_compile_cmd = [
        HOST_CXX,
        "-std=c++23",
        "-O2",
        "-Wall",
        "-Wextra",
        str(HOST_SRC),
        f"-I{SUITESPARSE_SRC_DIR / 'AMD' / 'Include'}",
        f"-I{SUITESPARSE_SRC_DIR / 'LDL' / 'Include'}",
        f"-I{SUITESPARSE_SRC_DIR / 'SuiteSparse_config'}",
        f"-L{amd_lib_dir}",
        f"-L{ldl_lib_dir}",
        f"-L{cfg_lib_dir}",
        f"-Wl,-rpath,{amd_lib_dir}",
        f"-Wl,-rpath,{ldl_lib_dir}",
        f"-Wl,-rpath,{cfg_lib_dir}",
        "-lamd",
        "-lldl",
        "-lsuitesparseconfig",
        "-o",
        str(precalc_bin),
    ]
    run_cmd(host_compile_cmd)

    # 执行 LDLT_PreCalc，产出静态常量头文件 (generated_symbolic.hpp)
    run_cmd([str(precalc_bin), str(GEN_HEADER)])

    if not GEN_HEADER.exists():
        sys.exit("[Error] Expected generated header was not created!")

    # ---------------------------------------------------------
    # Stage 2: 将预计算结果作为静态常量，编译目标端程序 LDLT_RunSlv
    # ---------------------------------------------------------
    print("\n=== Stage 2: Building Embedded Target Binary (LDLT_RunSlv) ===")
    runslv_bin = BUILD_TARGET_DIR / "LDLT_RunSlv"

    target_compile_cmd = [
        TARGET_CXX,
        "-std=c++23",
        "-O3",
        "-fno-exceptions",
        "-fno-rtti",
        "-Wall",
        "-Wextra",
        f"-I{GEN_DIR}",
        str(TARGET_SRC),
        "-lm",
        "-o",
        str(runslv_bin),
    ]
    run_cmd(target_compile_cmd)

    # ---------------------------------------------------------
    # Verification: 运行目标端二进制验证数值分解与残差
    # ---------------------------------------------------------
    print("\n=== Verification: Running Target Binary (LDLT_RunSlv) ===")
    run_cmd([str(runslv_bin)])


if __name__ == "__main__":
    main()
