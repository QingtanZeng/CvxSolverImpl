// src/preCalc02.cpp
#include <cstdio>
// 引入由 preCalc01 预计算生成的静态头文件
#include "generated_symbolic.hpp"

using namespace EmbeddedLDLT;

// 嵌入式端全部使用编译期确定大小的静态/栈内存工作区
static double Lx[L_NNZ > 0 ? L_NNZ : 1]; // L 因子的非零值
static int    Li[L_NNZ > 0 ? L_NNZ : 1]; // L 因子的行索引
static double D[N];                      // 对角阵 D 的数值
static double Y[N];                      // 内部浮点工作区
static int    Pattern[N];                // 消去树非零元可达路径缓冲
static int    Flag[N];                   // 节点访问标记工作区

// 零动态内存的静态消去树数值分解核心函数 (移植自 Timothy Davis 的 LDL 算法)
int embedded_ldl_numeric(const double Ax[]) {
    for (int k = 0; k < N; ++k) {
        Y[k] = 0.0;
        int top = N;
        Flag[k] = k;
        Lnz[k]; // 静态只读引用，实际动态计数用临时数组或原位跟踪
        int lnz_k = 0;

        int kk = P[k];
        for (int p = Ap[kk]; p < Ap[kk + 1]; ++p) {
            int i = Pinv[Ai[p]];
            if (i <= k) {
                Y[i] += Ax[p];
                int len;
                // 沿着静态消去树 Parent 向上回溯至根或已标记节点
                for (len = 0; Flag[i] != k; i = Parent[i]) {
                    Pattern[len++] = i;
                    Flag[i] = k;
                }
                while (len > 0) {
                    Pattern[--top] = Pattern[--len];
                }
            }
        }

        D[k] = Y[k];
        Y[k] = 0.0;
        for (; top < N; ++top) {
            int i = Pattern[top];
            double yi = Y[i];
            Y[i] = 0.0;
            int p2 = Lp[i] + Lnz[i];
            int p;
            for (p = Lp[i]; p < p2; ++p) {
                Y[Li[p]] -= Lx[p] * yi;
            }
            double l_ki = yi / D[i];
            D[k] -= l_ki * yi;
            Li[p] = k;
            Lx[p] = l_ki;
        }
        if (D[k] == 0.0) return -1; // 矩阵奇异
    }
    return 0;
}

int main() {
    // 生产环境中传感器/控制器实时更新的矩阵非零元数值 (与 Ap, Ai 对应)
    const double Ax[A_NNZ] = {
        2.0,        // Col 0: (0,0)
        5.0,        // Col 1: (1,1)
        1.0, 4.0,   // Col 2: (0,2), (2,2)
        2.0, 3.0,   // Col 3: (1,3), (3,3)
        3.0, 1.0, 6.0 // Col 4: (0,4), (2,4), (4,4)
    };

    printf("[preCalc02 - Embedded Target] Loaded static Elimination Tree & Permutation:\n");
    for (int i = 0; i < N; ++i) {
        printf("  Node %d: P[%d]=%d, Parent[%d]=%d, Lp[%d]=%d\n",
               i, i, P[i], i, Parent[i], i, Lp[i]);
    }

    return 0;
}