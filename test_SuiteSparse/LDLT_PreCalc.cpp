/**
 * @file LDLT_PreCalc.cpp
 * @brief Stage 1 (Host Pre-Calculator): 基于 SuiteSparse (AMD + LDL) 执行图论重排、
 *        对称置换与静态消去树符号分析，并生成 Target 端零依赖静态头文件。
 *
 * Host-Target Two-Stage Build 架构说明:
 *   Stage 1 (本文件, Host 端):
 *     - 运行在宿主机，允许动态内存分配、C++23 标准库及链接 SuiteSparse C 库。
 *     - 负责完成所有“与数值无关的图论与符号计算”：
 *         1. 校验原始稀疏对称矩阵 A 的上三角 CSC 规范性；
 *         2. 调用 amd_order() 计算降填充置换 P 与逆置换 Pinv；
 *         3. 构造规范上三角重排矩阵 K = P^T * A * P 及其 O(nnz) 数值重排映射 k_from_a_idx；
 *         4. 调用 ldl_symbolic() 提取静态消去树 Parent、L 因子列指针 Lp 及列非零计数 Lnz；
 *         5. 输出无任何外部依赖的 C++ constexpr 头文件 (generated_symbolic.hpp)。
 *   Stage 2 (LDLT_RunSlv.cpp, Target 端):
 *     - 不依赖 SuiteSparse、不使用堆分配、兼容 -fno-exceptions -fno-rtti，
 *       直接利用 Stage 1 固化的静态拓扑完成 O(nnz) 数值置换、LDLT 数值分解与回代求解。
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// 引入 SuiteSparse C 库头文件 (仅 Stage 1 Host 端依赖)
#include "SuiteSparse_config.h"
#include "amd.h"
#include "ldl.h"

namespace cvx::host {

/**
 * @brief 压缩稀疏列格式 (CSC) 对称矩阵 (仅存储上三角部分: row_idx[p] <= col)
 */
struct CSCMatrix {
    int32_t m{0};
    int32_t n{0};
    std::vector<int32_t> col_ptr;
    std::vector<int32_t> row_idx;
    std::vector<double>  values;

    [[nodiscard]] int32_t nnz() const noexcept {
        return static_cast<int32_t>(values.size());
    }

    /**
     * @brief 验证是否为规范的 CSC 上三角矩阵 (每列行索引严格递增且 row <= col)
     */
    [[nodiscard]] bool is_valid_upper_triangular() const noexcept {
        if (m <= 0 || n <= 0 || m != n) return false;
        if (static_cast<int32_t>(col_ptr.size()) != n + 1) return false;
        if (static_cast<int32_t>(row_idx.size()) != nnz()) return false;
        if (col_ptr[0] != 0 || col_ptr[n] != nnz()) return false;

        for (int32_t c = 0; c < n; ++c) {
            const int32_t start = col_ptr[c];
            const int32_t end   = col_ptr[c + 1];
            if (start > end) return false;

            int32_t last_r = -1;
            for (int32_t p = start; p < end; ++p) {
                const int32_t r = row_idx[p];
                if (r <= last_r || r > c) return false;
                last_r = r;
            }
        }
        return true;
    }

    /**
     * @brief 对称置换重排: 生成规范的 K = P^T * A * P 上三角 CSC 矩阵与静态 gather 映射
     *
     * 注意: SuiteSparse 的 ldl_symbolic/ldl_numeric 在传入 P 时要求输入矩阵包含完整的
     *       PAP' 上三角，因此若 A 仅存上三角，必须先显式置换为 K，再以 P=nullptr 调用 LDL。
     */
    [[nodiscard]] auto permute_symmetric(std::span<const int32_t> pinv) const;

    [[nodiscard]] std::vector<std::vector<double>> to_dense_symmetric() const {
        std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));
        for (int32_t c = 0; c < n; ++c) {
            for (int32_t p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
                const int32_t r = row_idx[p];
                mat[r][c] = values[p];
                mat[c][r] = values[p];
            }
        }
        return mat;
    }

    void print(std::string_view name) const {
        std::cout << std::format("\n=== 矩阵 {} ({}x{}, nnz={}) ===\n", name, m, n, nnz());
        std::cout << "col_ptr: [ ";
        for (auto p : col_ptr) std::cout << p << ' ';
        std::cout << "]\nrow_idx: [ ";
        for (auto r : row_idx) std::cout << r << ' ';
        std::cout << "]\nvalues : [ ";
        for (auto v : values) std::cout << std::format("{:.1f} ", v);
        std::cout << "]\n";

        if (n <= 8) {
            std::cout << "全对称稠密展开视图:\n";
            const auto dense = to_dense_symmetric();
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
 * @brief 对称置换重排的结果：规范上三角矩阵 K = P^T * A * P 及数值收集索引表
 */
struct PermutedResult {
    CSCMatrix K;
    std::vector<int32_t> k_from_a_idx; ///< K.values[p] == A.values[k_from_a_idx[p]]
};

[[nodiscard]] inline auto CSCMatrix::permute_symmetric(std::span<const int32_t> pinv) const {
    if (static_cast<int32_t>(pinv.size()) != n || m != n) {
        throw std::invalid_argument("Symmetric permutation requires square matrix and matching pinv size");
    }

    struct Triplet {
        int32_t r{0};
        int32_t c{0};
        double  val{0.0};
        int32_t src_idx{0};
    };

    std::vector<Triplet> triplets;
    triplets.reserve(nnz());

    for (int32_t c = 0; c < n; ++c) {
        for (int32_t p = col_ptr[c]; p < col_ptr[c + 1]; ++p) {
            const int32_t r = row_idx[p];
            if (r > c) continue;
            const int32_t r_new = pinv[r];
            const int32_t c_new = pinv[c];
            triplets.push_back({
                std::min(r_new, c_new),
                std::max(r_new, c_new),
                values[p],
                p
            });
        }
    }

    std::ranges::sort(triplets, [](const Triplet& a, const Triplet& b) {
        return (a.c != b.c) ? (a.c < b.c) : (a.r < b.r);
    });

    PermutedResult res;
    res.K.m = m;
    res.K.n = n;
    res.K.col_ptr.assign(n + 1, 0);
    for (const auto& t : triplets) {
        res.K.col_ptr[t.c + 1]++;
    }
    for (int32_t j = 0; j < n; ++j) {
        res.K.col_ptr[j + 1] += res.K.col_ptr[j];
    }

    res.K.row_idx.resize(triplets.size());
    res.K.values.resize(triplets.size());
    res.k_from_a_idx.resize(triplets.size());
    for (size_t idx = 0; idx < triplets.size(); ++idx) {
        res.K.row_idx[idx]      = triplets[idx].r;
        res.K.values[idx]       = triplets[idx].val;
        res.k_from_a_idx[idx]   = triplets[idx].src_idx;
    }
    return res;
}

/**
 * @brief 置换向量与逆置换向量
 */
struct Permutation {
    std::vector<int32_t> p;    ///< p[k] = j: 新位置 k 对应原索引 j
    std::vector<int32_t> pinv; ///< pinv[j] = k: 原索引 j 映射到新位置 k

    [[nodiscard]] static Permutation from_p(std::span<const int32_t> p_src) {
        Permutation perm;
        perm.p.assign(p_src.begin(), p_src.end());
        const auto n = static_cast<int32_t>(perm.p.size());
        perm.pinv.assign(n, 0);
        for (int32_t k = 0; k < n; ++k) {
            perm.pinv[perm.p[k]] = k;
        }
        return perm;
    }

    void print() const {
        const auto n = static_cast<int32_t>(p.size());
        std::cout << "\n=== AMD 重排置换结果 (Permutation) ===\n";
        std::cout << "置换向量 P     = [ ";
        for (auto x : p) std::cout << x << ' ';
        std::cout << "]\n逆置换向量 Pinv = [ ";
        for (auto x : pinv) std::cout << x << ' ';
        std::cout << "]\n坐标映射对照表: ";
        for (int32_t j = 0; j < n; ++j) {
            std::cout << std::format("原{}->新{}  ", j, pinv[j]);
        }
        std::cout << '\n';
    }
};

/**
 * @brief 静态消去树与 L 因子符号结构
 */
struct EliminationTree {
    int32_t n{0};
    std::vector<int32_t> parent;
    std::vector<int32_t> lnz;
    std::vector<int32_t> lp;

    [[nodiscard]] int32_t total_lnz_offdiag() const noexcept {
        return lp.empty() ? 0 : lp.back();
    }

    [[nodiscard]] int32_t total_lnz_with_diag() const noexcept {
        return total_lnz_offdiag() + n;
    }

    [[nodiscard]] std::vector<int32_t> compute_depths() const {
        std::vector<int32_t> depth(n, 0);
        for (int32_t i = 0; i < n; ++i) {
            int32_t d = 0;
            for (int32_t curr = parent[i]; curr != -1; curr = parent[curr]) {
                ++d;
            }
            depth[i] = d;
        }
        return depth;
    }

    void print() const {
        std::cout << "\n=== 静态消去树 (Static Elimination Tree) 分析 ===\n";
        const auto depths = compute_depths();
        int32_t max_depth = 0;
        int32_t num_trees = 0;
        for (int32_t i = 0; i < n; ++i) {
            max_depth = std::max(max_depth, depths[i]);
            if (parent[i] == -1) {
                ++num_trees;
                std::cout << std::format("  节点 {:2d} -> [根节点 ROOT]   (树深度 = {:2d}, L列非零元 = {:2d})\n",
                                         i, depths[i], lnz[i]);
            } else {
                std::cout << std::format("  节点 {:2d} -> 父节点 {:2d}       (树深度 = {:2d}, L列非零元 = {:2d})\n",
                                         i, parent[i], depths[i], lnz[i]);
            }
        }
        std::cout << std::format("消去森林统计: 独立子树 = {}, 最大深度 = {}, L(不含对角) = {}, L(含对角) = {}\n",
                                 num_trees, max_depth, total_lnz_offdiag(), total_lnz_with_diag());
    }
};

/**
 * @brief SuiteSparse Host 端符号预计算封装
 */
class SuiteSparsePreCalculator {
public:
    [[nodiscard]] static Permutation compute_amd_ordering(const CSCMatrix& A) {
        const int32_t n = A.n;
        std::vector<int32_t> P(n, 0);
        double control[AMD_CONTROL]{};
        double info[AMD_INFO]{};

        amd_defaults(control);
        const int status = amd_order(n, A.col_ptr.data(), A.row_idx.data(), P.data(), control, info);
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

    [[nodiscard]] static EliminationTree compute_symbolic_ldlt(const CSCMatrix& K) {
        const int32_t n = K.n;
        EliminationTree etree;
        etree.n = n;
        etree.lp.assign(n + 1, 0);
        etree.parent.assign(n, -1);
        etree.lnz.assign(n, 0);

        std::vector<int32_t> flag(n, 0);
        ldl_symbolic(
            n,
            const_cast<int32_t*>(K.col_ptr.data()),
            const_cast<int32_t*>(K.row_idx.data()),
            etree.lp.data(),
            etree.parent.data(),
            etree.lnz.data(),
            flag.data(),
            nullptr,
            nullptr
        );
        return etree;
    }
};

/**
 * @brief 单个算例的 Host 端预计算完整快照
 */
struct PreCalcArtifact {
    std::string ns_name;
    std::string description;
    CSCMatrix A;
    Permutation perm;
    PermutedResult permuted;
    EliminationTree etree;
    int32_t unpermuted_lnz{0};
};

[[nodiscard]] PreCalcArtifact analyze_case(std::string_view ns_name,
                                           std::string_view description,
                                           const CSCMatrix& A) {
    std::cout << "\n-----------------------------------------------------------------\n";
    std::cout << std::format(" Host 预计算算例 [{}]: {}\n", ns_name, description);
    std::cout << "-----------------------------------------------------------------";

    if (!A.is_valid_upper_triangular()) {
        throw std::runtime_error(std::format("Case '{}' matrix is not valid upper-triangular CSC!", ns_name));
    }

    A.print("原始对称矩阵 A (上三角存储)");

    // 1. 记录未重排直接分解的填充规模 (用于评估 AMD 降填充收益)
    const auto etree_no_amd = SuiteSparsePreCalculator::compute_symbolic_ldlt(A);
    const int32_t unpermuted_lnz = etree_no_amd.total_lnz_offdiag();
    std::cout << std::format("\n[对照] 未重排直接分解 L 非零元数 (不含对角): {}\n", unpermuted_lnz);

    // 2. AMD 排序 -> 显式对称重排 K = P^T * A * P -> 静态消去树符号分析
    auto perm     = SuiteSparsePreCalculator::compute_amd_ordering(A);
    auto permuted = A.permute_symmetric(perm.pinv);
    auto etree    = SuiteSparsePreCalculator::compute_symbolic_ldlt(permuted.K);

    perm.print();
    permuted.K.print("重排矩阵 K = P^T * A * P (上三角存储)");
    etree.print();
    std::cout << std::format("==> AMD 降填充效果: {} -> {}\n", unpermuted_lnz, etree.total_lnz_offdiag());

    return PreCalcArtifact{
        .ns_name        = std::string(ns_name),
        .description    = std::string(description),
        .A              = A,
        .perm           = std::move(perm),
        .permuted       = std::move(permuted),
        .etree          = std::move(etree),
        .unpermuted_lnz = unpermuted_lnz
    };
}

template <typename T>
void emit_std_array(std::ostream& os,
                    std::string_view elem_type,
                    std::string_view name,
                    std::span<const T> data,
                    std::string_view fmt_spec) {
    os << std::format("inline constexpr std::array<{}, {}> {} {{\n    ", elem_type, data.size(), name);
    for (size_t i = 0; i < data.size(); ++i) {
        if constexpr (std::is_same_v<T, double>) {
            os << std::format("{:.17g}", data[i]);
        } else {
            (void)fmt_spec;
            os << std::format("{}", data[i]);
        }
        if (i + 1 < data.size()) {
            os << ((i + 1) % 10 == 0 ? ",\n    " : ", ");
        }
    }
    os << "\n};\n\n";
}

void emit_case_namespace(std::ostream& os, const PreCalcArtifact& art) {
    os << "// -----------------------------------------------------------------------------\n";
    os << std::format("// {}: {}\n", art.ns_name, art.description);
    os << "// -----------------------------------------------------------------------------\n";
    os << std::format("namespace {} {{\n\n", art.ns_name);
    os << std::format("inline constexpr int32_t N                = {};\n", art.A.n);
    os << std::format("inline constexpr int32_t A_NNZ            = {};\n", art.A.nnz());
    os << std::format("inline constexpr int32_t L_NNZ            = {};\n", art.etree.total_lnz_offdiag());
    os << std::format("inline constexpr int32_t UNPERMUTED_L_NNZ = {};\n\n", art.unpermuted_lnz);

    emit_std_array<int32_t>(os, "int32_t", "A_Ap",         art.A.col_ptr,             "{}");
    emit_std_array<int32_t>(os, "int32_t", "A_Ai",         art.A.row_idx,             "{}");
    emit_std_array<double> (os, "double",  "A_Ax_Default", art.A.values,              "{}");
    emit_std_array<int32_t>(os, "int32_t", "P",            art.perm.p,                "{}");
    emit_std_array<int32_t>(os, "int32_t", "Pinv",         art.perm.pinv,             "{}");
    emit_std_array<int32_t>(os, "int32_t", "K_Ap",         art.permuted.K.col_ptr,    "{}");
    emit_std_array<int32_t>(os, "int32_t", "K_Ai",         art.permuted.K.row_idx,    "{}");
    emit_std_array<int32_t>(os, "int32_t", "K_From_A_Idx", art.permuted.k_from_a_idx, "{}");
    emit_std_array<int32_t>(os, "int32_t", "Parent",       art.etree.parent,          "{}");
    emit_std_array<int32_t>(os, "int32_t", "Lp",           art.etree.lp,              "{}");
    emit_std_array<int32_t>(os, "int32_t", "Lnz_Expected", art.etree.lnz,             "{}");

    os << std::format("}} // namespace {}\n\n", art.ns_name);
}

void write_generated_header(const std::filesystem::path& out_path,
                            std::span<const PreCalcArtifact> artifacts) {
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream ofs(out_path);
    if (!ofs) {
        throw std::runtime_error(std::format("Failed to open output header: {}", out_path.string()));
    }

    ofs << "// =============================================================================\n";
    ofs << "// AUTO-GENERATED BY LDLT_PreCalc (Stage 1 Host) - DO NOT EDIT MANUALLY\n";
    ofs << "// Contains static AMD permutations, K = P^T * A * P gather maps, and\n";
    ofs << "// static Elimination Tree symbolic data for zero-allocation Target solving.\n";
    ofs << "// =============================================================================\n";
    ofs << "#pragma once\n\n";
    ofs << "#include <array>\n";
    ofs << "#include <cstdint>\n\n";
    ofs << "namespace EmbeddedLDLT {\n\n";

    for (const auto& art : artifacts) {
        emit_case_namespace(ofs, art);
    }

    ofs << "} // namespace EmbeddedLDLT\n";
    ofs.flush();
    if (!ofs) {
        throw std::runtime_error(std::format("Failed while writing header: {}", out_path.string()));
    }
}

} // namespace cvx::host

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << std::format("Usage: {} <output_generated_header_path>\n", argv[0]);
        return 1;
    }
    const std::filesystem::path out_header = argv[1];

    std::cout << "=================================================================\n";
    std::cout << " [Stage 1 - Host PreCalc] SuiteSparse AMD 排序与静态消去树预计算\n";
    std::cout << "=================================================================\n";

    try {
        // 算例 1: 5x5 对称稀疏正定矩阵 A (规范 CSC 上三角)
        const cvx::host::CSCMatrix mat5{
            .m = 5,
            .n = 5,
            .col_ptr = {0, 1, 2, 4, 6, 9},
            .row_idx = {
                0,
                1,
                0, 2,
                1, 3,
                0, 2, 4
            },
            .values = {
                2.0,
                5.0,
                1.0, 4.0,
                2.0, 3.0,
                3.0, 1.0, 6.0
            }
        };

        // 算例 2: 6x6 经典箭头矩阵 (Arrowhead Matrix)
        const cvx::host::CSCMatrix arrow6{
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

        std::vector<cvx::host::PreCalcArtifact> artifacts;
        artifacts.push_back(cvx::host::analyze_case(
            "Matrix5x5",
            "5x5 Symmetric Sparse Matrix (Baseline AMD + LDLT Verification)",
            mat5
        ));
        artifacts.push_back(cvx::host::analyze_case(
            "Arrow6x6",
            "6x6 Arrowhead Matrix (Zero Extra Fill-in under AMD Ordering)",
            arrow6
        ));

        cvx::host::write_generated_header(out_header, artifacts);
        std::cout << std::format("\n[Stage 1 - Host PreCalc] 成功生成目标端静态符号头文件: {}\n",
                                 out_header.string());
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << std::format("[Stage 1 Error] {}\n", ex.what());
        return 1;
    }
}
