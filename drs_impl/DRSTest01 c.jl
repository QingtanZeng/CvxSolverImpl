using LinearAlgebra, SparseArrays
using Printf

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

# Objective: minimize y=(1/2)x^T*P*x
P=sparse(Float64.(I(2))*2)
c=[0.0,0]

# Equality constraints: x₁ = 1, x₂ = 2
# Ax + s = b, s=0
# Inequality constraints: x>=0
# -Ix+s=0, s>=0

# Ax + s = b
A = sparse([1.0 1.0; -1.0 1.0; -1 0; 0 -1])
b = [2.0; 0; 0; 0]
m = length(b)
n = length(c) # Number of variables

# Set the problem
# ==========================================
# 2. 定义锥 (Cone Definitions)
# 对应 c_wrapper.jl 中的参数 [cite: 24, 25]
# ==========================================
# z: Zero cone size (对应等式约束的数量)
z = 2 
# l: Linear cone size (对应不等式约束 Ax + s = b, s >= 0)
l = 0

#### Douglas-Rachford Splitting(DRS) ####
M = [P A';-A zeros(m,m)]
q = [c; b]

K=[(I(n)+P) A'; A -I(m)]
F = ldlt(K)
F_LD = sparse(F.LD)

function project_cone(v)
    x_part = v[1:n]
    y_part = v[n+1:end]

    proj_x = x_part;

    y_eq = y_part[1:2]
    y_ineq = max.(0.0, y_part[3:4])

    return [proj_x; y_eq; y_ineq]
end


max_iter = 50
u_nxt = -5*ones(n+m);
u_lst = u_nxt;
w= zeros(n+m)
w[end]=1
println("primal: $(u_nxt[1:n])")
for k in 1:max_iter
     # Step 1
     rhs = w-q;
     rx=rhs[1:n];
     ry=rhs[n+1:end];
     rhs_new = [rx;-ry]
     u_tilde= F \ rhs_new

     # Step 2
     v_input = 2*u_tilde - w
     u_nxt = project_cone(v_input)

     # Step 3
     w = w + u_nxt - u_tilde

     # Stop Critera
     x=u_nxt[1:n]
     y=u_nxt[n+1:end]
     s=b-A*x
     gap=s'*y;
     pres_eq = norm(s[1:2])/norm(b[1:2])
     pres_ineq = norm(max.(0.0, 0 .- s[3:4]))

     pdchange=norm(u_nxt-u_lst)/norm(u_lst)

     if( pres_ineq < 1e-6 &&  pres_eq < 3e-6 && (pdchange < 5E-2 || gap < 1e-2))
          println("Converged at iter $k")
          println("gap: $(@sprintf("%.5E", gap)) ")
          println("pdchange: $(@sprintf("%.5E", pdchange)) ")
          println("primal: [$(join([@sprintf("%.5E", x) for x in u_nxt[1:n]], ", "))]")
          # println("dual: [$(join([@sprintf("%.5E", x) for x in u_nxt[n+1:end]], ", "))]")
          println("primal Equality Residual: $(@sprintf("%.5E",pres_eq)); primal Inequality Residual: $(@sprintf("%.5E",pres_ineq))")
          break;
     end

     u_lst=u_nxt;

     println("Iter $k, ")
     println("gap: $(@sprintf("%.5E", gap)) ")
     println("pdchange: $(@sprintf("%.5E", pdchange)) ")
     println("primal: [$(join([@sprintf("%.5E", x) for x in u_nxt[1:n]], ", "))]")
     # println("dual: [$(join([@sprintf("%.5E", x) for x in u_nxt[n+1:end]], ", "))]")
     println("primal Equality Residual: $(@sprintf("%.5E",pres_eq)); primal Inequality Residual: $(@sprintf("%.5E",pres_ineq))")
     println(" ");
end
