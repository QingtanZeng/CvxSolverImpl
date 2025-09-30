using Plots
using LaTeXStrings

# 使用 GR 绘图后端，并设置一些默认样式
gr(size=(800, 700), legend=false)

# --- 1. 定义数据和函数 ---

# 定义 z₁ 的范围
z1 = -1.2:0.01:1.2

# f₁(z) = z₂ - z₁² = 0  =>  z₂ = z₁²
z2_f1 = z1.^2

# f₂(z) = z₂ + 0.1z₁ = 0.6  =>  z₂ = -0.1z₁ + 0.6
z2_f2 = -0.1 .* z1 .+ 0.6

# 线性化点 z̄ (从图上估计)
z_bar = [0.5, 0.25]

# f₁ 在 z̄ 处的切线: ∇f₁(z̄)ᵀ(z - z̄) = 0
# ∇f₁(z) = [-2z₁, 1]ᵀ, ∇f₁(z̄) = [-1, 1]ᵀ
# -1(z₁ - 0.5) + 1(z₂ - 0.25) = 0  =>  z₂ = z₁ - 0.25
z2_tangent = z1 .- 0.25

# 信任域 (Trust Region) 圆
# ||z - z̄||₂ ≤ η, 从图上估计半径 η ≈ 0.3
η = 0.3
θ = range(0, 2π, length=100)
circle_z1 = z_bar[1] .+ η .* cos.(θ)
circle_z2 = z_bar[2] .+ η .* sin.(θ)

# --- 2. 开始绘图 ---

# 创建基础绘图，设置坐标轴、网格和范围
p = plot(
    xlabel=L"z_1",  # 使用 LaTeXStrings
    ylabel=L"z_2",
    xlims=(-1.2, 1.2),
    ylims=(-1.2, 1.2),
    aspect_ratio=:equal, # 确保圆形看起来是圆的
    framestyle=:origin,  # 坐标轴在原点交叉
    grid=true,
    gridstyle=:dash,
    gridalpha=0.5
)

# --- 3. 逐个添加绘图元素 ---

# 绘制蓝色抛物线 f₁
plot!(p, z1, z2_f1, linewidth=2, color=:blue)

# 绘制绿色直线 f₂
plot!(p, z1, z2_f2, linewidth=2, color=:green)

# 绘制蓝色虚线切线
plot!(p, z1, z2_tangent, linewidth=2, linestyle=:dash, color=:darkblue)

# 绘制 z₁ ≤ 0.85 的紫色虚线边界
vline!(p, [0.85], linewidth=2, linestyle=:dash, color=:magenta)

# 绘制并填充红色的信任域圆
plot!(p, circle_z1, circle_z2, seriestype=[:shape],
    fillcolor=:red, fillalpha=0.3, linecolor=:darkred)

# 标记线性化点 z̄
scatter!(p, [z_bar[1]], [z_bar[2]], markersize=6, markercolor=:black, label=L"\bar{z}")
# 给 z̄ 添加文本标签
annotate!(p, z_bar[1] - 0.1, z_bar[2], text(L"\bar{z}", :black, :right, 12))


# 标记解（黑点）
# 解是切线和绿线的交点：z₁ - 0.25 = -0.1z₁ + 0.6 => 1.1z₁ = 0.85 => z₁ ≈ 0.7727
sol_z1 = 0.85 / 1.1
sol_z2 = sol_z1 - 0.25
scatter!(p, [sol_z1], [sol_z2], markersize=6, markercolor=:black)


# --- 4. 添加文本和箭头标注 ---

# 标注函数 f₁ 和 f₂
annotate!(p, -0.9, 0.9, text(L"f_1(z) \triangleq z_2 - z_1^2 = 0", :blue, :left, 12))
annotate!(p, -0.9, 0.75, text(L"f_2(z) \triangleq z_2 + 0.1z_1 = 0.6", :green, :left, 12))

# 标注切线
annotate!(p, -1.1, -1.0, text(L"f_1(\bar{z}) + \nabla f_1(\bar{z})^\top(z - \bar{z}) = 0", :darkblue, :left, 12))

# 标注信任域
annotate!(p, 0.8, 0.6, text(L"||z - \bar{z}||_2 \le \eta", :darkred, 12))

# 标注 z₁ ≤ 0.85
annotate!(p, 0.9, 0.0, text(L"z_1 \le 0.85", :magenta, 12))
# 添加紫色箭头
annotate!(p, 0.85, 0.2, arrow=Plots.Arrow(:open, :head, 1.5, 1.5))
annotate!(p, 0.85, -0.2, arrow=Plots.Arrow(:open, :head, 1.5, 1.5))
annotate!(p, 0.95, 1.0, arrow=Plots.Arrow(:open, :head, 1.5, 1.5))

# 添加 "Decreasing Cost" 箭头和文本
annotate!(p, 0.5, -0.4, text("Decreasing\nCost", :black, 10))
annotate!(p, 0.5, -0.6, arrow=Plots.Arrow(:closed, :head, 2.0, 2.0))


# --- 5. 显示最终图像 ---
display(p)