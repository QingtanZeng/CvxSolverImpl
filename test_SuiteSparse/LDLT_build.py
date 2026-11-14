#!/usr/bin/env python3
"""
两阶段编译与预计算粘合脚本：
Stage 1: 使用 Host Clang 工具链编译并执行 preCalc01 (链接 SuiteSparse)，产出静态常量头文件与 JSON 快照。
Stage 2: 使用 Target Clang (或交叉编译工具链) 编译 preCalc02，固化 P 与消去树并生成生产环境固件。
"""

import subprocess
import sys
from pathlib import Path

# 路径配置
ROOT_DIR = Path(__file__).resolve().parent
BUILD_HOST_DIR = ROOT_DIR / "build_host"
BUILD_TARGET_DIR = ROOT_DIR / "build_target"
GEN_DIR = ROOT_DIR / "build_artifacts"
GEN_HEADER = GEN_DIR / "generated_symbolic.hpp"

# 编译器与参数配置
HOST_CXX = "clang++"
TARGET_CXX = "clang++"  # 若为交叉编译，可改为目标工具链或追加 --target=arm-none-eabi


def run_cmd(cmd: list[str], cwd: Path = ROOT_DIR):
    print(f"[Orchestrator] Executing: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd, text=True)
    if result.returncode != 0:
        sys.exit(f"[Error] Command failed with exit code {result.returncode}")


def main():
    BUILD_HOST_DIR.mkdir(exist_ok=True)
    BUILD_TARGET_DIR.mkdir(exist_ok=True)
    GEN_DIR.mkdir(exist_ok=True)

    # ---------------------------------------------------------
    # Stage 1: 在主机端编译并运行 preCalc01
    # ---------------------------------------------------------
    print("\n=== Stage 1: Building & Running Host Pre-Calculator (preCalc01) ===")
    precalc01_bin = BUILD_HOST_DIR / "preCalc01"

    host_compile_cmd = [
        HOST_CXX,
        "-std=c++17",
        "-O2",
        str(ROOT_DIR / "src" / "preCalc01.cpp"),
        "-I/usr/include/suitesparse",
        "-I/usr/local/include/suitesparse",
        "-lamd",
        "-lldl",
        "-o",
        str(precalc01_bin),
    ]
    run_cmd(host_compile_cmd)

    # 执行 preCalc01，生成中间预计算结果 (generated_symbolic.hpp)
    run_cmd([str(precalc01_bin), str(GEN_HEADER)])

    if not GEN_HEADER.exists():
        sys.exit("[Error] Expected generated header was not created!")

    # ---------------------------------------------------------
    # Stage 2: 将预计算结果作为静态常量，编译嵌入式目标程序 preCalc02
    # ---------------------------------------------------------
    print("\n=== Stage 2: Building Embedded Target Binary (preCalc02) ===")
    precalc02_bin = BUILD_TARGET_DIR / "preCalc02"

    target_compile_cmd = [
        TARGET_CXX,
        "-std=c++17",
        "-O3",
        "-fno-exceptions",
        "-fno-rtti",  # 典型嵌入式裁剪参数
        f"-I{GEN_DIR}",  # 引入生成的静态头文件目录
        str(ROOT_DIR / "src" / "preCalc02.cpp"),
        "-o",
        str(precalc02_bin),
    ]
    run_cmd(target_compile_cmd)

    # ---------------------------------------------------------
    # 验证运行 (如果在宿主机模拟运行)
    # ---------------------------------------------------------
    print("\n=== Verification: Running Target Binary ===")
    run_cmd([str(precalc02_bin)])


if __name__ == "__main__":
    main()