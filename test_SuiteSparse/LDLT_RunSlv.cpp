/**
 * @file LDLT_RunSlv.cpp
 * @brief Stage 2 (Embedded Target Runtime Solver): 零堆分配、零 SuiteSparse 依赖的
 *        静态消去树 LDL^T 数值分解与线性方程组回代求解器。
 *
 * Host-Target Two-Stage Build 架构特性:
 *   1. 完全脱离 SuiteSparse (无 amd.h / ldl.h 头文件或库依赖)；
 *   2. 仅引入 Stage 1 生成的静态常量头文件 generated_symbolic.hpp；
 *   3. 全部工作区在编译期确定大小 (std::array)，零堆内存分配 (no new / std::vector)；
 *   4. 兼容嵌入式严苛编译选项: -fno-exceptions -fno-rtti；
 *   5. 通过 Stage 1 预计算的 K_From_A_Idx 完成 O(nnz) 数值置换，并在固定拓扑上支持
 *      在线数值更新后的快速重复分解与求解 (典型用于嵌入式 MPC / 内点法 KKT 迭代)。
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>

// 仅依赖 Stage 1 (Host) 预计算生成的静态符号头文件
#include "generated_symbolic.hpp"

namespace cvx::target {

/**
 * @brief 编译期确定容量的静态稀疏 LDL^T 求解内核 (K = P^T * A * P = L * D * L^T)
 *
 * @tparam Case 由 generated_symbolic.hpp 生成的静态算例命名空间包装结构体
 */
template <typename Case>
class StaticLDLTKernel {
public:
    static constexpr int32_t N     = Case::N;
    static constexpr int32_t A_NNZ = Case::A_NNZ;
    static constexpr int32_t L_NNZ = Case::L_NNZ;
    static constexpr size_t  L_CAP = (L_NNZ > 0) ? static_cast<size_t>(L_NNZ) : 1U;

    /**
     * @brief 使用静态映射表将原始上三角矩阵 A 的实时数值 Ax 映射为重排矩阵 K 的数值
     *        时间复杂度 O(A_NNZ)，零排序、零三元组分配
     */
    void update_values(const std::array<double, static_cast<size_t>(A_NNZ)>& a_ax) noexcept {
        a_ax_ = a_ax;
        for (int32_t p = 0; p < A_NNZ; ++p) {
            const auto src = static_cast<size_t>(Case::K_From_A_Idx[static_cast<size_t>(p)]);
            k_ax_[static_cast<size_t>(p)] = a_ax_[src];
        }
    }

    /**
     * @brief 零堆分配的稀疏数值分解 K = L * D * L^T (Up-looking row-by-row LDLT)
     *
     * @return 成功返回 N；若在第 k 列遇到零主元 D[k] == 0，则提前返回 k
     */
    int32_t factorize() noexcept {
        for (int32_t k = 0; k < N; ++k) {
            const auto uk = static_cast<size_t>(k);
            y_[uk]        = 0.0;
            flag_[uk]     = k;
            lnz_curr_[uk] = 0; // 动态跟踪 L 每列当前已写入非零元数 (不可直接修改 constexpr Lnz)

            int32_t top = N;
            const int32_t p_start = Case::K_Ap[uk];
            const int32_t p_end   = Case::K_Ap[uk + 1U];

            for (int32_t p = p_start; p < p_end; ++p) {
                int32_t i = Case::K_Ai[static_cast<size_t>(p)];
                if (i <= k) {
                    y_[static_cast<size_t>(i)] += k_ax_[static_cast<size_t>(p)];
                    int32_t len = 0;
                    while (flag_[static_cast<size_t>(i)] != k) {
                        pattern_[static_cast<size_t>(len++)] = i;
                        flag_[static_cast<size_t>(i)]        = k;
                        i = Case::Parent[static_cast<size_t>(i)];
                    }
                    while (len > 0) {
                        pattern_[static_cast<size_t>(--top)] = pattern_[static_cast<size_t>(--len)];
                    }
                }
            }

            d_[uk] = y_[uk];
            y_[uk] = 0.0;

            for (; top < N; ++top) {
                const int32_t i   = pattern_[static_cast<size_t>(top)];
                const auto    ui  = static_cast<size_t>(i);
                const double  yi  = y_[ui];
                y_[ui]            = 0.0;

                const int32_t p_col_begin = Case::Lp[ui];
                const int32_t p_col_end   = p_col_begin + lnz_curr_[ui];
                int32_t p = p_col_begin;
                for (; p < p_col_end; ++p) {
                    const auto up = static_cast<size_t>(p);
                    y_[static_cast<size_t>(li_[up])] -= lx_[up] * yi;
                }

                const double l_ki = yi / d_[ui];
                d_[uk] -= l_ki * yi;

                const auto up = static_cast<size_t>(p);
                li_[up] = k;
                lx_[up] = l_ki;
                ++lnz_curr_[ui];
            }

            if (d_[uk] == 0.0) {
                return k;
            }
        }
        return N;
    }

    /**
     * @brief 三步回代求解 A * x = b:
     *        1. b_K = P^T * b
     *        2. L * D * L^T * y = b_K
     *        3. x = P * y
     */
    void solve(const std::array<double, static_cast<size_t>(N)>& b,
               std::array<double, static_cast<size_t>(N)>& x_out) const noexcept {
        std::array<double, static_cast<size_t>(N)> work{};

        // 1. b_K[k] = b[P[k]] (对应 ldl_perm)
        for (int32_t k = 0; k < N; ++k) {
            const auto uk = static_cast<size_t>(k);
            work[uk] = b[static_cast<size_t>(Case::P[uk])];
        }

        // 2a. 前代解 L * z = b_K (对应 ldl_lsolve, 单位下三角且对角线未存储)
        for (int32_t j = 0; j < N; ++j) {
            const auto uj    = static_cast<size_t>(j);
            const double xj  = work[uj];
            const int32_t p2 = Case::Lp[uj + 1U];
            for (int32_t p = Case::Lp[uj]; p < p2; ++p) {
                const auto up = static_cast<size_t>(p);
                work[static_cast<size_t>(li_[up])] -= lx_[up] * xj;
            }
        }

        // 2b. 对角求解 D * w = z (对应 ldl_dsolve)
        for (int32_t j = 0; j < N; ++j) {
            const auto uj = static_cast<size_t>(j);
            work[uj] /= d_[uj];
        }

        // 2c. 回代解 L^T * y = w (对应 ldl_ltsolve)
        for (int32_t j = N - 1; j >= 0; --j) {
            const auto uj    = static_cast<size_t>(j);
            const int32_t p2 = Case::Lp[uj + 1U];
            double sum       = work[uj];
            for (int32_t p = Case::Lp[uj]; p < p2; ++p) {
                const auto up = static_cast<size_t>(p);
                sum -= lx_[up] * work[static_cast<size_t>(li_[up])];
            }
            work[uj] = sum;
        }

        // 3. 逆置换还原 x_out[P[j]] = work[j] (对应 ldl_permt)
        for (int32_t j = 0; j < N; ++j) {
            const auto uj = static_cast<size_t>(j);
            x_out[static_cast<size_t>(Case::P[uj])] = work[uj];
        }
    }

    /**
     * @brief 基于上三角 CSC 存储的稀疏对称矩阵向量乘 b = A * x (零稠密矩阵分配)
     */
    void spmv_sym_upper(const std::array<double, static_cast<size_t>(N)>& x,
                        std::array<double, static_cast<size_t>(N)>& b_out) const noexcept {
        b_out.fill(0.0);
        for (int32_t c = 0; c < N; ++c) {
            const auto uc    = static_cast<size_t>(c);
            const int32_t p2 = Case::A_Ap[uc + 1U];
            for (int32_t p = Case::A_Ap[uc]; p < p2; ++p) {
                const auto up   = static_cast<size_t>(p);
                const int32_t r = Case::A_Ai[up];
                const auto ur   = static_cast<size_t>(r);
                const double v  = a_ax_[up];
                b_out[ur] += v * x[uc];
                if (r != c) {
                    b_out[uc] += v * x[ur];
                }
            }
        }
    }

    [[nodiscard]] bool verify_symbolic_counts() const noexcept {
        for (int32_t k = 0; k < N; ++k) {
            const auto uk = static_cast<size_t>(k);
            if (lnz_curr_[uk] != Case::Lnz_Expected[uk]) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] const std::array<double, static_cast<size_t>(N)>& diag() const noexcept {
        return d_;
    }

private:
    std::array<double, static_cast<size_t>(A_NNZ)>  a_ax_{};
    std::array<double, static_cast<size_t>(A_NNZ)>  k_ax_{};
    std::array<int32_t, L_CAP>                      li_{};
    std::array<double, L_CAP>                       lx_{};
    std::array<double, static_cast<size_t>(N)>      d_{};
    std::array<double, static_cast<size_t>(N)>      y_{};
    std::array<int32_t, static_cast<size_t>(N)>     pattern_{};
    std::array<int32_t, static_cast<size_t>(N)>     flag_{};
    std::array<int32_t, static_cast<size_t>(N)>     lnz_curr_{};
};

// 将 generated_symbolic.hpp 生成的命名空间常量适配为模板策略类型
#define CVX_DEFINE_CASE_SPEC(SpecName, NsName)                                             \
    struct SpecName {                                                                      \
        static constexpr int32_t N                = EmbeddedLDLT::NsName::N;                \
        static constexpr int32_t A_NNZ            = EmbeddedLDLT::NsName::A_NNZ;            \
        static constexpr int32_t L_NNZ            = EmbeddedLDLT::NsName::L_NNZ;            \
        static constexpr int32_t UNPERMUTED_L_NNZ = EmbeddedLDLT::NsName::UNPERMUTED_L_NNZ; \
        static constexpr const auto& A_Ap         = EmbeddedLDLT::NsName::A_Ap;             \
        static constexpr const auto& A_Ai         = EmbeddedLDLT::NsName::A_Ai;             \
        static constexpr const auto& A_Ax_Default = EmbeddedLDLT::NsName::A_Ax_Default;     \
        static constexpr const auto& P            = EmbeddedLDLT::NsName::P;                \
        static constexpr const auto& Pinv         = EmbeddedLDLT::NsName::Pinv;             \
        static constexpr const auto& K_Ap         = EmbeddedLDLT::NsName::K_Ap;             \
        static constexpr const auto& K_Ai         = EmbeddedLDLT::NsName::K_Ai;             \
        static constexpr const auto& K_From_A_Idx = EmbeddedLDLT::NsName::K_From_A_Idx;     \
        static constexpr const auto& Parent       = EmbeddedLDLT::NsName::Parent;           \
        static constexpr const auto& Lp           = EmbeddedLDLT::NsName::Lp;               \
        static constexpr const auto& Lnz_Expected = EmbeddedLDLT::NsName::Lnz_Expected;     \
    }

CVX_DEFINE_CASE_SPEC(Matrix5x5Spec, Matrix5x5);
CVX_DEFINE_CASE_SPEC(Arrow6x6Spec, Arrow6x6);
#undef CVX_DEFINE_CASE_SPEC

template <typename Case>
bool run_and_verify_case(const char* title,
                         const std::array<double, static_cast<size_t>(Case::A_NNZ)>& ax_values) noexcept {
    constexpr int32_t N = Case::N;
    StaticLDLTKernel<Case> solver;

    std::printf("\n-----------------------------------------------------------------\n");
    std::printf(" [Stage 2 - Target RunSlv] %s (N=%d, L_NNZ=%d, 未重排填充=%d)\n",
                title, Case::N, Case::L_NNZ, Case::UNPERMUTED_L_NNZ);
    std::printf("-----------------------------------------------------------------\n");

    solver.update_values(ax_values);
    const int32_t rank = solver.factorize();
    if (rank != N) {
        std::printf("[FAIL] 矩阵数值奇异或不定: 在列 %d 处对角元 D[%d] == 0\n", rank, rank);
        return false;
    }

    if (!solver.verify_symbolic_counts()) {
        std::printf("[FAIL] 运行时 L 非零元计数与 Host 端静态消去树预期不一致!\n");
        return false;
    }

    std::printf("因子 D 对角元 = [ ");
    for (int32_t i = 0; i < N; ++i) {
        std::printf("%.4f ", solver.diag()[static_cast<size_t>(i)]);
    }
    std::printf("]\n");

    // 构造精确测试向量 x_exact = [1, 2, ..., N]^T 并用稀疏乘积计算 b = A * x_exact
    std::array<double, static_cast<size_t>(N)> x_exact{};
    std::array<double, static_cast<size_t>(N)> b{};
    std::array<double, static_cast<size_t>(N)> x_sol{};
    std::array<double, static_cast<size_t>(N)> ax_check{};

    for (int32_t i = 0; i < N; ++i) {
        x_exact[static_cast<size_t>(i)] = static_cast<double>(i + 1);
    }

    solver.spmv_sym_upper(x_exact, b);
    solver.solve(b, x_sol);
    solver.spmv_sym_upper(x_sol, ax_check);

    double max_err_x = 0.0;
    double max_res_b = 0.0;
    for (int32_t i = 0; i < N; ++i) {
        const auto ui = static_cast<size_t>(i);
        max_err_x = std::max(max_err_x, std::abs(x_sol[ui] - x_exact[ui]));
        max_res_b = std::max(max_res_b, std::abs(ax_check[ui] - b[ui]));
    }

    std::printf("  真实解 x_exact = [ ");
    for (int32_t i = 0; i < N; ++i) std::printf("%.1f ", x_exact[static_cast<size_t>(i)]);
    std::printf("]\n");

    std::printf("  数值解 x_sol   = [ ");
    for (int32_t i = 0; i < N; ++i) std::printf("%.4f ", x_sol[static_cast<size_t>(i)]);
    std::printf("]\n");

    std::printf("  ||x_sol - x_exact||_inf = %.2e, ||A*x_sol - b||_inf = %.2e\n", max_err_x, max_res_b);

    const bool ok = (max_err_x < 1e-12) && (max_res_b < 1e-12);
    std::printf("  ==> %s\n", ok ? "[SUCCESS] Target 静态消去树数值分解与回代求解高精度通过!"
                                 : "[FAIL] 误差超出容差!");
    return ok;
}

} // namespace cvx::target

int main() {
    std::printf("=================================================================\n");
    std::printf(" [Stage 2 - Target RunSlv] 零堆分配、无 SuiteSparse 依赖求解验证\n");
    std::printf(" 编译契约: C++23 | -fno-exceptions -fno-rtti | 静态工作区\n");
    std::printf("=================================================================\n");

    bool all_ok = true;

    // 1. 基线 5x5 对称稀疏矩阵验证 (对应 tst_AMD 算例 1)
    all_ok &= cvx::target::run_and_verify_case<cvx::target::Matrix5x5Spec>(
        "算例 1: 5x5 对称稀疏矩阵 (默认数值)",
        cvx::target::Matrix5x5Spec::A_Ax_Default
    );

    // 2. 在线数值更新验证: 保持结构不变，仅更新非零元数值 (模拟内点法/MPC KKT 迭代)
    const std::array<double, static_cast<size_t>(cvx::target::Matrix5x5Spec::A_NNZ)> mat5_updated_ax = {
        4.5,         // (0,0)
        8.0,         // (1,1)
        1.2, 6.5,    // (0,2), (2,2)
        2.5, 5.0,    // (1,3), (3,3)
        2.0, 1.5, 9.0 // (0,4), (2,4), (4,4)
    };
    all_ok &= cvx::target::run_and_verify_case<cvx::target::Matrix5x5Spec>(
        "算例 1b: 5x5 矩阵在线数值更新 (复用同一静态消去树)",
        mat5_updated_ax
    );

    // 3. 经典 6x6 箭头矩阵验证 (对应 tst_AMD 算例 2: AMD 将填充从 15 降为 5)
    all_ok &= cvx::target::run_and_verify_case<cvx::target::Arrow6x6Spec>(
        "算例 2: 6x6 箭头矩阵 (Arrowhead, AMD 零额外填充)",
        cvx::target::Arrow6x6Spec::A_Ax_Default
    );

    std::printf("\n=================================================================\n");
    std::printf(" %s\n", all_ok ? "全部 Target 端静态 LDLT 验证顺利通过!" : "Target 端验证存在失败项!");
    std::printf("=================================================================\n");

    return all_ok ? 0 : 1;
}
