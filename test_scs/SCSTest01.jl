using SCS
using LinearAlgebra, SparseArrays

"""
min x₁² + x₂² 
s.t. [1,1   [x₁    [2
     -1,1] * x₂] =  0]
     
The optimal solution is clearly x₁=1, x₂=1, and cost=2.
"""
# --- Problem Data ---
# 1. Define the problem in SCS standard form
# z = [x₁, x₂]
n = 2 # Number of variables

# Objective: minimize y=(1/2)x^T*P*x
P=sparse(Float64.(I(2))*2)
c=[0.0,0]

# Equality constraints: Ax = b
# x₁ = 1
# x₂ = 2
A = sparse([1.0 1.0; -1.0 1.0])
b = [2.0; 0]
m = length(b)


# Set the problem
# ==========================================
# 2. 定义锥 (Cone Definitions)
# 对应 c_wrapper.jl 中的参数 [cite: 24, 25]
# ==========================================
# z: Zero cone size (对应等式约束的数量)
z = 2 
# l: Linear cone size (对应不等式约束 Ax + s = b, s >= 0)
l = 0
# 其他锥参数初始化为空或 0
bu = Float64[]      # Upper bounds for box cone
bl = Float64[]      # Lower bounds for box cone
q = Int[]           # Second-order cone sizes
s = Int[]           # PSD cone sizes
cs = Int[]          # Complex PSD cone sizes
ep = 0              # Primal exponential cones
ed = 0              # Dual exponential cones
p = Float64[]       # Power cone parameters
d = Int[]           # Log-det cone sizes
nuc_m = Int[]       # Nuclear norm rows
nuc_n = Int[]       # Nuclear norm cols
ell1 = Int[]        # L1-matrix norm sizes
sl_n = Int[]        # Ky-Fan norm cone sizes
sl_k = Int[]        # Ky-Fan norm cone constants

solver_type=SCS.DirectSolver

solution = SCS.scs_solve(
    solver_type,
    m, n,
    A, P, b, c,
    z, l, bu, bl, q, s, cs, ep, ed, p, d, nuc_m, nuc_n, ell1, sl_n, sl_k,
    verbose=true,
    log_csv_filename = "scs_iteration_log.csv" # 指定输出文件名
);

solution.info.status
solution.x