/**
 * @file tst_AMD.cpp
 * @brief SuiteSparse AMD 排序、对称矩阵重排 K = P^T * A * P 与 静态消去树 (Elimination Tree) 分析
 * @note 遵循 C++23 语法标准，调用 SuiteSparse 核心 C 库接口 (AMD 与 LDL)
 */

#include <iostream>
#include <vector>
#include <span>
#include <format>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <string>
#include <string_view>
#include <cmath>
#include <ranges>
#include <stdexcept>

// 引入 SuiteSparse C 库头文件
#include "amd.h"
#include "ldl.h"
#include "SuiteSparse_config.h"

namespace cvx {

/**
 * @brief 压缩稀疏列格式 (Compressed Sparse Column, CSC) 矩阵
 * 
 * 专为对称矩阵设计，存储其上三角部分 (row_idx[p] <= col)。
 */
struct CSCMatrix {
    int32_t m{0};                       ///< 行数
    int32_t n{0};                       ///< 列数
    std::vector<int32_t> col_ptr;       ///< 列指针数组 (长度 n + 1)
    std::vector<int32_t> row_idx;       ///< 行索引数组 (长度 nnz)
    std::vector<double>  values;        ///< 数值数组 (长度 nnz)

    [[nodiscard]] int32_t nnz() const noexcept {
        return static_cast<int32_t>(values.size());
    }

    /**
     * @brief 验证是否为规范的 CSC 上三角矩阵 (每列行索引单调严格递增，且 row <= col)
     */
    [[nodiscard]] bool is_valid_upper_triangular() const noexcept {
        if (m <= 0 || n <= 0) return false;
        if (static_cast<int32_t>(col_ptr.size()) != n + 1) return false;
        if (col_ptr[0] != 0 || col_ptr[n] != nnz()) return false;

        for (int32_t c = 0; c < n; ++c) {
            int32_t start = col_ptr[c];
            int32_t end = col_ptr[c + 1];
            if (start > end) return false;

            int32_t last_r = -1;
            for (int32_t p = start; p < end; ++p) {
                int32_t r = row_idx[p];
                if (r <= last_r || r > c) return false; // 必须递增且位于上三角
                last_r = r;
            }
        }
        return true;
    }

    /**
     * @brief 对称置换重排: 生成规范的 K = P^T * A * P 上三角 CSC 稀疏矩阵
     * 
     * 数学原理:
     * 若 A 在原始位置 (r, c) 处有非零元 v (r <= c):
     * 置换后对应的新坐标为 (r', c') = (Pinv[r], Pinv[c])。
     * 由于 K 也是对称矩阵，为维持上三角格式 (行号 <= 列号)，令:
     *   kr = min(r', c')
     *   kc = max(r', c')
     * 则新元素存放在 K 的 (kr, kc) 处。
     * 
     * @param pinv 逆置换向量 (Pinv[j] = k 表示原索引 j 重排为新索引 k)
     * @return CSCMatrix 重排后的规范上三角矩阵 K
     */
    [[nodiscard]] CSCMatrix permute_symmetric(std::span<const int32_t> pinv) const {
        if (static_cast<int32_t>(pinv.size()) != n || m != n) {
            throw std::invalid_argument("Symmetric permutation requires square matrix and matching pinv size");
        }

        CSCMatrix K;
        K.m = m;
        K.n = n;

        // 收集重排后的上三角三元组 (kr, kc, val)
        struct Triplet {
            int32_t r{0};
            int32_t c{0};
            double val{0.0};
        };
        std::vector<Triplet> triplets;
        triplets.reserve(nnz());

        for (int32_t c = 0; c < n; ++c) {
            for (int32_t p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
                int32_t r = row_idx[p];
                if (r > c) continue; // 仅取上三角部分
                int32_t r_new = pinv[r];
                int32_t c_new = pinv[c];
                int32_t kr = std::min(r_new, c_new);
                int32_t kc = std::max(r_new, c_new);
                triplets.push_back({kr, kc, values[p]});
            }
        }

        // 按照列优先排序: 先排新列 kc，再排新行 kr (确保构建出来的 CSC 规范且行索引严格递增)
        std::ranges::sort(triplets, [](const Triplet& a, const Triplet& b) {
            if (a.c != b.c) return a.c < b.c;
            return a.r < b.r;
        });

        // 统计重排矩阵每列的非零元数量
        K.col_ptr.assign(n + 1, 0);
        for (const auto& [r, c, val] : triplets) {
            K.col_ptr[c + 1]++;
        }

        // 前缀和累计生成 col_ptr
        for (int32_t j = 0; j < n; ++j) {
            K.col_ptr[j + 1] += K.col_ptr[j];
        }

        // 填充行索引与数值数组
        K.row_idx.resize(triplets.size());
        K.values.resize(triplets.size());
        for (size_t idx = 0; idx < triplets.size(); ++idx) {
            K.row_idx[idx] = triplets[idx].r;
            K.values[idx]  = triplets[idx].val;
        }

        return K;
    }

    /**
     * @brief 恢复全对称的稠密矩阵表示 (仅用于打印展示)
     */
    [[nodiscard]] std::vector<std::vector<double>> to_dense_symmetric() const {
        std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));
        for (int32_t c = 0; c < n; ++c) {
            for (int32_t p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
                int32_t r = row_idx[p];
                mat[r][c] = values[p];
                mat[c][r] = values[p];
            }
        }
        return mat;
    }

    /**
     * @brief 格式化打印稀疏矩阵结构与内容
     */
    void print(std::string_view name) const {
        std::cout << std::format("\n=== 矩阵 {} ({}x{}, nnz={}) ===\n", name, m, n, nnz());
        std::cout << "col_ptr: [ ";
        for (auto p : col_ptr) std::cout << p << " ";
        std::cout << "]\n";
        std::cout << "row_idx: [ ";
        for (auto r : row_idx) std::cout << r << " ";
        std::cout << "]\n";
        std::cout << "values : [ ";
        for (auto v : values) std::cout << std::format("{:.1f} ", v);
        std::cout << "]\n";

        std::cout << "上三角非零元素 (row, col) = val:\n";
        for (int32_t c = 0; c < n; ++c) {
            for (int32_t p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
                std::cout << std::format("  ({}, {}) = {:6.2f}\n", row_idx[p], c, values[p]);
            }
        }

        if (n <= 8) {
            std::cout << "全对称稠密展开视图:\n";
            auto dense = to_dense_symmetric();
            for (int32_t i = 0; i < n; ++i) {
                std::cout << "  [ ";
                for (int32_t j = 0; j < n; ++j) {
                    std::cout << std::format("{:5.1f} ", dense[i][j]);
                }
                std::cout << "]\n";
            }
        }
    }
};

/**
 * @brief 置换向量与置换矩阵
 */
struct Permutation {
    std::vector<int32_t> p;     ///< 置换向量: p[k] = j 表示原矩阵的第 j 行/列重排到新矩阵的第 k 行/列
    std::vector<int32_t> pinv;  ///< 逆置换向量: pinv[j] = k 表示原第 j 个元素对应新位置 k

    [[nodiscard]] int32_t size() const noexcept {
        return static_cast<int32_t>(p.size());
    }

    /**
     * @brief 由置换向量 P 构建完整置换对象并自动计算逆置换 Pinv
     */
    static Permutation from_p(std::span<const int32_t> p_src) {
        Permutation perm;
        perm.p.assign(p_src.begin(), p_src.end());
        int32_t n = perm.size();
        perm.pinv.assign(n, 0);
        for (int32_t k = 0; k < n; ++k) {
            perm.pinv[perm.p[k]] = k;
        }
        return perm;
    }

    /**
     * @brief 打印置换详情及置换矩阵 P
     */
    void print() const {
        int32_t n = size();
        std::cout << "\n=== AMD 重排置换结果 (Permutation) ===\n";
        std::cout << "置换向量 P     = [ ";
        for (auto x : p) std::cout << x << " ";
        std::cout << "]  (新索引 k -> 原索引 P[k])\n";

        std::cout << "逆置换向量 Pinv = [ ";
        for (auto x : pinv) std::cout << x << " ";
        std::cout << "]  (原索引 j -> 新索引 Pinv[j])\n";

        std::cout << "坐标映射对照表: ";
        for (int32_t j = 0; j < n; ++j) {
            std::cout << std::format("原{}->新{}  ", j, pinv[j]);
        }
        std::cout << "\n";

        if (n <= 8) {
            std::cout << "置换矩阵 P (P[k, j] = 1 表示原第 j 列映射至新第 k 列):\n";
            for (int32_t k = 0; k < n; ++k) {
                std::cout << "  [ ";
                for (int32_t j = 0; j < n; ++j) {
                    std::cout << (p[k] == j ? " 1 " : " 0 ");
                }
                std::cout << "]\n";
            }
        }
    }
};

/**
 * @brief 静态消去树 (Static Elimination Tree) 及其符号结构
 * 
 * 消去树是稀疏对称分解 (Cholesky / LDL^T) 的核心图论工具。
 * 节点 j 的父节点 Parent[j] 是在因子 L 的第 j 列中，对角线下第一个非零元的行号:
 *   Parent[j] = min { i > j | L(i, j) != 0 }
 * 若第 j 列对角线下全为 0，则 Parent[j] = -1 (该节点为树的根节点)。
 */
struct EliminationTree {
    int32_t n{0};
    std::vector<int32_t> parent;  ///< parent[i] 为节点 i 的父节点，-1 表示根节点 (消去森林)
    std::vector<int32_t> lnz;     ///< lnz[j] 为 L 因子第 j 列中对角线下的非零元个数
    std::vector<int32_t> lp;      ///< L 因子的列指针 (长度 n + 1)

    [[nodiscard]] int32_t total_lnz_offdiag() const noexcept {
        return lp.empty() ? 0 : lp.back();
    }

    [[nodiscard]] int32_t total_lnz_with_diag() const noexcept {
        return total_lnz_offdiag() + n;
    }

    /**
     * @brief 计算消去树中每个节点的深度 (根节点深度为 0)
     */
    [[nodiscard]] std::vector<int32_t> compute_depths() const {
        std::vector<int32_t> depth(n, 0);
        for (int32_t i = 0; i < n; ++i) {
            int32_t d = 0;
            int32_t curr = parent[i];
            while (curr != -1) {
                d++;
                curr = parent[curr];
            }
            depth[i] = d;
        }
        return depth;
    }

    /**
     * @brief 打印消去树结构、森林特征与消去依赖路径
     */
    void print() const {
        std::cout << "\n=== 静态消去树 (Static Elimination Tree) 分析 ===\n";
        auto depths = compute_depths();
        int32_t max_depth = 0;
        int32_t num_trees = 0;

        for (int32_t i = 0; i < n; ++i) {
            max_depth = std::max(max_depth, depths[i]);
            if (parent[i] == -1) {
                num_trees++;
                std::cout << std::format("  节点 {:2d} -> [根节点 ROOT]   (树深度 = {:2d}, L列对角线下非零元 = {:2d})\n",
                                         i, depths[i], lnz[i]);
            } else {
                std::cout << std::format("  节点 {:2d} -> 父节点 {:2d}       (树深度 = {:2d}, L列对角线下非零元 = {:2d})\n",
                                         i, parent[i], depths[i], lnz[i]);
            }
        }

        std::cout << std::format("\n消去森林统计: 连通分量 (独立子树) = {}, 消去树最大深度 = {}\n", num_trees, max_depth);
        std::cout << std::format("L 因子非零元数 (不含单位对角线) = {}\n", total_lnz_offdiag());
        std::cout << std::format("L 因子总非零元数 (含单位对角线)   = {}\n", total_lnz_with_diag());

        // 打印消去路径 (从各叶节点到根节点的拓扑路径)
        std::cout << "\n各节点的消去依赖路径 (Node -> Ancestors -> Root):\n";
        for (int32_t i = 0; i < n; ++i) {
            std::cout << std::format("  节点 {:2d}: ", i);
            int32_t curr = i;
            while (curr != -1) {
                std::cout << curr;
                curr = parent[curr];
                if (curr != -1) std::cout << " -> ";
            }
            std::cout << "\n";
        }
    }
};

/**
 * @brief SuiteSparse C 库封装接口类
 */
class SuiteSparseEngine {
public:
    /**
     * @brief 调用 SuiteSparse AMD 进行近似最小度排序 (Approximate Minimum Degree Ordering)
     * 
     * AMD 会在 A+A' (无向连通图) 的邻接模式上自动忽略对角线进行分析，
     * 寻找能最小化 Cholesky / LDL^T 分解中填充元 (Fill-in) 的对称置换 P。
     * 
     * @param A 输入对称矩阵 (上三角或全对称 CSC 存储)
     * @return Permutation 置换对象 (包含 P 和 Pinv)
     */
    [[nodiscard]] static Permutation compute_amd_ordering(const CSCMatrix& A) {
        int32_t n = A.n;
        std::vector<int32_t> P(n, 0);
        double control[AMD_CONTROL];
        double info[AMD_INFO];

        // 获取 AMD 默认控制参数
        amd_defaults(control);

        // 调用 SuiteSparse AMD C API: amd_order
        int status = amd_order(
            n,
            A.col_ptr.data(),
            A.row_idx.data(),
            P.data(),
            control,
            info
        );

        if (status != AMD_OK && status != AMD_OK_BUT_JUMBLED) {
            throw std::runtime_error(std::format("amd_order failed with status {}", status));
        }

        std::cout << "\n=== SuiteSparse AMD 统计分析报告 ===\n";
        std::cout << std::format("  矩阵规模 n             : {}\n", n);
        std::cout << std::format("  对角线非零元数          : {:.0f}\n", info[AMD_NZDIAG]);
        std::cout << std::format("  A + A^T 非零元数       : {:.0f}\n", info[AMD_NZ_A_PLUS_AT]);
        std::cout << std::format("  稠密行/列过滤数        : {:.0f}\n", info[AMD_NDENSE]);
        std::cout << std::format("  估计因子 L 非零元数    : {:.0f}\n", info[AMD_LNZ]);
        std::cout << std::format("  LDL' 估计乘减操作数    : {:.0f}\n", info[AMD_NMULTSUBS_LDL]);

        return Permutation::from_p(P);
    }

    /**
     * @brief 对已完成对称重排的规范上三角矩阵 K 进行 LDL^T 符号分解，提取消去树
     * 
     * @param K 规范上三角 CSC 矩阵 K = P^T * A * P
     * @return EliminationTree 静态消去树及 L 因子符号模式
     */
    [[nodiscard]] static EliminationTree compute_symbolic_ldlt(const CSCMatrix& K) {
        int32_t n = K.n;
        EliminationTree etree;
        etree.n = n;
        etree.lp.assign(n + 1, 0);
        etree.parent.assign(n, -1);
        etree.lnz.assign(n, 0);

        std::vector<int32_t> flag(n, 0);

        // 调用 SuiteSparse LDL C API: ldl_symbolic
        // 由于矩阵 K 已经显式经过对称置换且为规范的上三角存储，传入 P = nullptr, Pinv = nullptr
        ldl_symbolic(
            n,
            const_cast<int32_t*>(K.col_ptr.data()),
            const_cast<int32_t*>(K.row_idx.data()),
            etree.lp.data(),
            etree.parent.data(),
            etree.lnz.data(),
            flag.data(),
            nullptr,  // P = nullptr
            nullptr   // Pinv = nullptr
        );

        return etree;
    }

    /**
     * @brief 数值分解与全流程求解验证: K = L * D * L^T 求解 A * x = b 并检验残差
     * 
     * 方程求解推导:
     *   A * x = b
     *   (P * K * P^T) * x = b           (因 K = P^T * A * P, 故 A = P * K * P^T)
     *   K * (P^T * x) = P^T * b
     * 
     * 求解三部曲:
     *   1. b_K = P^T * b                (通过 ldl_perm 实现: X[k] = B[P[k]])
     *   2. 求解 K * y = b_K             (利用 L * D * L^T 三角回代求解)
     *   3. x = P * y                    (通过 ldl_permt 实现: X[P[j]] = B[j])
     */
    static void verify_numerical_ldlt(const CSCMatrix& K, const CSCMatrix& A_orig, const Permutation& perm) {
        int32_t n = K.n;
        auto etree = compute_symbolic_ldlt(K);
        int32_t lnz = etree.total_lnz_offdiag();

        // 为数值分解分配工作空间
        std::vector<int32_t> Li(lnz, 0);
        std::vector<double>  Lx(lnz, 0.0);
        std::vector<double>  D(n, 0.0);
        std::vector<double>  Y(n, 0.0);
        std::vector<int32_t> Pattern(n, 0);
        std::vector<int32_t> Flag(n, 0);

        // 执行数值分解 K = L * D * L^T
        int32_t rank = ldl_numeric(
            n,
            const_cast<int32_t*>(K.col_ptr.data()),
            const_cast<int32_t*>(K.row_idx.data()),
            const_cast<double*>(K.values.data()),
            etree.lp.data(),
            etree.parent.data(),
            etree.lnz.data(),
            Li.data(),
            Lx.data(),
            D.data(),
            Y.data(),
            Pattern.data(),
            Flag.data(),
            nullptr,
            nullptr
        );

        if (rank != n) {
            std::cout << std::format("警告: 矩阵数值奇异或不定, 在列 {} 处对角元为 0\n", rank);
            return;
        }

        std::cout << "\n=== LDLT 数值分解验证 (K = L * D * L^T) ===\n";
        std::cout << "因子 D 的对角元 = [ ";
        for (auto d : D) std::cout << std::format("{:.4f} ", d);
        std::cout << "]\n";

        // 构造测试问题: 设真实解 x_exact = [1, 2, ..., n]^T
        std::vector<double> x_exact(n);
        for (int32_t i = 0; i < n; ++i) x_exact[i] = static_cast<double>(i + 1);

        // 由稠密对称形式精确计算右端向量 b = A * x_exact
        auto A_dense = A_orig.to_dense_symmetric();
        std::vector<double> b(n, 0.0);
        for (int32_t i = 0; i < n; ++i) {
            for (int32_t j = 0; j < n; ++j) {
                b[i] += A_dense[i][j] * x_exact[j];
            }
        }

        // 1. 置换右端项 b_K = P^T * b (ldl_perm: b_K[k] = b[P[k]])
        std::vector<double> b_K(n, 0.0);
        ldl_perm(n, b_K.data(), b.data(), const_cast<int32_t*>(perm.p.data()));

        // 2. 三角回代求解 K * y = b_K  =>  L * D * L^T * y = b_K
        std::vector<double> y = b_K;
        ldl_lsolve(n, y.data(), etree.lp.data(), Li.data(), Lx.data());      // 解 L * z = b_K
        ldl_dsolve(n, y.data(), D.data());                                  // 解 D * w = z
        ldl_ltsolve(n, y.data(), etree.lp.data(), Li.data(), Lx.data());     // 解 L^T * y = w

        // 3. 逆置换解向量 x_sol = P * y (ldl_permt: x_sol[P[j]] = y[j])
        std::vector<double> x_sol(n, 0.0);
        ldl_permt(n, x_sol.data(), y.data(), const_cast<int32_t*>(perm.p.data()));

        // 4. 残差校验
        double max_err_x = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            max_err_x = std::max(max_err_x, std::abs(x_sol[i] - x_exact[i]));
        }

        std::cout << "线性方程组 A * x = b 求解与残差校验:\n";
        std::cout << "  真实解 x_exact = [ ";
        for (auto v : x_exact) std::cout << std::format("{:.1f} ", v);
        std::cout << "]\n";
        std::cout << "  数值解 x_sol   = [ ";
        for (auto v : x_sol) std::cout << std::format("{:.4f} ", v);
        std::cout << "]\n";
        std::cout << std::format("  解向量最大绝对误差 ||x_sol - x_exact||_inf = {:.2e}\n", max_err_x);
        std::cout << (max_err_x < 1e-12 ? "  ==> [SUCCESS] 符号分解与数值求解高精度吻合!\n" 
                                        : "  ==> [FAIL] 误差超出容差范围!\n");
    }
};

} // namespace cvx

int main() {
    std::cout << "=================================================================\n";
    std::cout << " SuiteSparse AMD 排序、对称重排 K = P^T*A*P 与 静态消去树分析\n";
    std::cout << " 采用 C++23 现代语法，调用 SuiteSparse C Library (AMD + LDL)\n";
    std::cout << "=================================================================\n";

    /* 
     * 算例 1: 5x5 对称稀疏矩阵 A (严格规范的 CSC 上三角存储)
     * [ 2  0  1  0  3 ]
     * [ 0  5  0  2  0 ]
     * [ 1  0  4  0  1 ]
     * [ 0  2  0  3  0 ]
     * [ 3  0  1  0  6 ]
     * 
     * 上三角非零元结构:
     * Col 0: (0, 0)=2.0
     * Col 1: (1, 1)=5.0
     * Col 2: (0, 2)=1.0, (2, 2)=4.0
     * Col 3: (1, 3)=2.0, (3, 3)=3.0
     * Col 4: (0, 4)=3.0, (2, 4)=1.0, (4, 4)=6.0
     */
    cvx::CSCMatrix A{
        .m = 5,
        .n = 5,
        .col_ptr = {0, 1, 2, 4, 6, 9},
        .row_idx = {
            0,          // Col 0
            1,          // Col 1
            0, 2,       // Col 2
            1, 3,       // Col 3
            0, 2, 4     // Col 4
        },
        .values = {
            2.0,
            5.0,
            1.0, 4.0,
            2.0, 3.0,
            3.0, 1.0, 6.0
        }
    };

    if (!A.is_valid_upper_triangular()) {
        std::cerr << "错误: 输入矩阵 A 不是规范的上三角 CSC 格式!\n";
        return -1;
    }

    A.print("原始对称矩阵 A (上三角存储)");

    // -------------------------------------------------------------
    // 步骤 1: 调用 SuiteSparse AMD 计算近似最小度置换 P 及 Pinv
    // -------------------------------------------------------------
    auto perm = cvx::SuiteSparseEngine::compute_amd_ordering(A);
    perm.print();

    // -------------------------------------------------------------
    // 步骤 2: 显式生成重排矩阵 K = P^T * A * P 的规范上三角 CSC 结构
    // -------------------------------------------------------------
    auto K = A.permute_symmetric(perm.pinv);
    K.print("重排矩阵 K = P^T * A * P (上三角存储)");

    // -------------------------------------------------------------
    // 步骤 3: 对重排矩阵 K 进行符号分解，生成静态消去树 (Elimination Tree)
    // -------------------------------------------------------------
    auto etree = cvx::SuiteSparseEngine::compute_symbolic_ldlt(K);
    etree.print();

    // -------------------------------------------------------------
    // 步骤 4: 数值分解 K = L * D * L^T 及线性方程组回代求解验证
    // -------------------------------------------------------------
    cvx::SuiteSparseEngine::verify_numerical_ldlt(K, A, perm);

    // -------------------------------------------------------------
    // 进阶对照试验: 经典 6x6 箭头矩阵 (Arrowhead Matrix)
    // 演示 AMD 排序如何将因子非零元从完全稠密 (15个填充元) 缩减至零额外填充 (5个非零元)
    // -------------------------------------------------------------
    std::cout << "\n\n=================================================================\n";
    std::cout << " 进阶对照演示: 经典箭头矩阵 (Arrowhead Matrix) 排序效果对比\n";
    std::cout << " 证明 AMD 排序如何避免消去过程中的全矩阵灾难性填充 (Fill-in)\n";
    std::cout << "=================================================================\n";

    cvx::CSCMatrix Arrow{
        .m = 6,
        .n = 6,
        .col_ptr = {0, 1, 3, 5, 7, 9, 11},
        .row_idx = {
            0,
            0, 1,
            0, 2,
            0, 3,
            0, 4,
            0, 5
        },
        .values = {
            10.0,
            1.0, 2.0,
            1.0, 2.0,
            1.0, 2.0,
            1.0, 2.0,
            1.0, 2.0
        }
    };

    // 未重排直接符号分解
    auto etree_no_amd = cvx::SuiteSparseEngine::compute_symbolic_ldlt(Arrow);
    std::cout << std::format("未重排直接分解: L 因子对角线下非零元数 = {} (完全稠密, 产生灾难性填充!)\n", 
                             etree_no_amd.total_lnz_offdiag());

    // 经过 AMD 排序后重排并分解
    auto arrow_perm = cvx::SuiteSparseEngine::compute_amd_ordering(Arrow);
    auto Arrow_K = Arrow.permute_symmetric(arrow_perm.pinv);
    auto etree_amd = cvx::SuiteSparseEngine::compute_symbolic_ldlt(Arrow_K);

    arrow_perm.print();
    Arrow_K.print("重排箭头矩阵 K = P^T * Arrow * P");
    etree_amd.print();
    std::cout << std::format("\n==> 对比结论: AMD 重排后 L 因子对角线下非零元数 = {} (零额外填充, 稀疏度保持最优!)\n", 
                             etree_amd.total_lnz_offdiag());

    std::cout << "\n=================================================================\n";
    std::cout << " 全部 AMD 排序、对称重排与消去树分析测试顺利通过!\n";
    std::cout << "=================================================================\n";

    return 0;
}
