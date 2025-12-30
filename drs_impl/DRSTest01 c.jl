using SCS
using LinearAlgebra, SparseArrays

"""
min x₁² + x₂² 
s.t. [1,1   [x₁    [2
     -1,1] * x₂] =  0]
     x_1>=0, x_2>=0
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
