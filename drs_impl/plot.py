import matplotlib.pyplot as plt
import numpy as np

# ==========================================
# 1. 原始数据输入
# ==========================================
raw_data = """
primal: [0.0, 0.0]

Iter 1, 
gap: -8.88889E-01 
pdchange: Inf 
primal: [6.66667E-01, 1.00000E+00]
primal Equality Residual: 2.35702E-01; primal Inequality Residual: 0.00000E+00
 
Iter 2, 
gap: -3.88889E-01 
pdchange: 1.66667E-01 
primal: [1.00000E+00, 8.33333E-01]
primal Equality Residual: 1.17851E-01; primal Inequality Residual: 0.00000E+00
 
Iter 3, 
gap: 2.71605E-01 
pdchange: 9.01670E-02 
primal: [1.11111E+00, 1.00000E+00]
primal Equality Residual: 7.85674E-02; primal Inequality Residual: 0.00000E+00
 
Iter 4, 
gap: 3.70370E-01 
pdchange: 6.53720E-02 
primal: [1.11111E+00, 1.05556E+00]
primal Equality Residual: 8.78410E-02; primal Inequality Residual: 0.00000E+00
 
Iter 5, 
gap: 2.70919E-01 
pdchange: 5.04113E-02 
primal: [1.07407E+00, 1.05556E+00]
primal Equality Residual: 6.54729E-02; primal Inequality Residual: 0.00000E+00
 
Iter 6, 
gap: 1.49520E-01 
pdchange: 3.29154E-02 
primal: [1.03704E+00, 1.03704E+00]
primal Equality Residual: 3.70370E-02; primal Inequality Residual: 0.00000E+00
 
Iter 7, 
gap: 6.14236E-02 
pdchange: 1.76676E-02 
primal: [1.01235E+00, 1.01852E+00]
primal Equality Residual: 1.57377E-02; primal Inequality Residual: 0.00000E+00
 
Iter 8, 
gap: 1.22695E-02 
pdchange: 7.96039E-03 
primal: [1.00000E+00, 1.00617E+00]
primal Equality Residual: 4.36486E-03; primal Inequality Residual: 0.00000E+00
 
Iter 9, 
gap: -8.16271E-03 
pdchange: 3.86926E-03 
primal: [9.95885E-01, 1.00000E+00]
primal Equality Residual: 2.90990E-03; primal Inequality Residual: 0.00000E+00
 
Iter 10, 
gap: -1.22949E-02 
pdchange: 2.80138E-03 
primal: [9.95885E-01, 9.97942E-01]
primal Equality Residual: 3.25337E-03; primal Inequality Residual: 0.00000E+00
 
Iter 11, 
gap: -9.58620E-03 
pdchange: 2.06540E-03 
primal: [9.97257E-01, 9.97942E-01]
primal Equality Residual: 2.42492E-03; primal Inequality Residual: 0.00000E+00
 
Iter 12, 
gap: -5.48509E-03 
pdchange: 1.28565E-03 
primal: [9.98628E-01, 9.98628E-01]
primal Equality Residual: 1.37174E-03; primal Inequality Residual: 0.00000E+00
 
Iter 13, 
gap: -2.28666E-03 
pdchange: 6.67004E-04 
primal: [9.99543E-01, 9.99314E-01]
primal Equality Residual: 5.82878E-04; primal Inequality Residual: 0.00000E+00
 
Iter 14, 
gap: -4.57352E-04 
pdchange: 2.95163E-04 
primal: [1.00000E+00, 9.99771E-01]
primal Equality Residual: 1.61661E-04; primal Inequality Residual: 0.00000E+00
 
Iter 15, 
gap: 3.04925E-04 
pdchange: 1.42545E-04 
primal: [1.00015E+00, 1.00000E+00]
primal Equality Residual: 1.07774E-04; primal Inequality Residual: 0.00000E+00
 
Iter 16, 
gap: 4.57317E-04 
pdchange: 1.03165E-04 
primal: [1.00015E+00, 1.00008E+00]
primal Equality Residual: 1.20495E-04; primal Inequality Residual: 0.00000E+00
 
Iter 17, 
gap: 3.55659E-04 
pdchange: 7.61972E-05 
primal: [1.00010E+00, 1.00008E+00]
primal Equality Residual: 8.98119E-05; primal Inequality Residual: 0.00000E+00
 
Iter 18, 
gap: 2.03224E-04 
pdchange: 4.75206E-05 
primal: [1.00005E+00, 1.00005E+00]
primal Equality Residual: 5.08053E-05; primal Inequality Residual: 0.00000E+00
 
Iter 19, 
gap: 8.46749E-05 
pdchange: 2.46863E-05 
primal: [1.00002E+00, 1.00003E+00]
primal Equality Residual: 2.15881E-05; primal Inequality Residual: 0.00000E+00
 
Iter 20, 
gap: 1.69349E-05 
pdchange: 1.09315E-05 
primal: [1.00000E+00, 1.00001E+00]
primal Equality Residual: 5.98746E-06; primal Inequality Residual: 0.00000E+00
 
Iter 21, 
gap: -1.12899E-05 
pdchange: 5.28048E-06 
primal: [9.99994E-01, 1.00000E+00]
primal Equality Residual: 3.99164E-06; primal Inequality Residual: 0.00000E+00
 
Iter 22, 
gap: -1.69350E-05 
pdchange: 3.82173E-06 
primal: [9.99994E-01, 9.99997E-01]
primal Equality Residual: 4.46279E-06; primal Inequality Residual: 0.00000E+00
 
Iter 23, 
gap: -1.31717E-05 
pdchange: 2.82253E-06 
primal: [9.99996E-01, 9.99997E-01]
primal Equality Residual: 3.32637E-06; primal Inequality Residual: 0.00000E+00
 
Converged at iter 24
gap: -7.52670E-06 
pdchange: 1.76015E-06 
primal: [9.99998E-01, 9.99998E-01]
primal Equality Residual: 1.88168E-06; primal Inequality Residual: 0.00000E+00
"""

# ==========================================
# 2. 数据解析
# ==========================================
iters = []
gaps = []
dltper_PDs = []
dltper_pres_eqs = []
pres_ineqs = []

lines = raw_data.strip().split('\n')

for line in lines:
    line = line.strip()
    line = line.replace('\xa0', ' ') 

    if line.startswith("Iter") or line.startswith("Converged at iter"):
        parts = line.split()
        for part in parts:
            clean_part = part.replace(',', '')
            if clean_part.isdigit():
                iters.append(int(clean_part))
                break
                
    elif line.startswith("gap:"):
        val_str = line.split(':')[1].strip()
        gaps.append(float(val_str))
        
    elif line.startswith("pdchange:"):
        val_str = line.split(':')[1].strip()
        if val_str == "Inf":
            dltper_PDs.append(np.nan)
        else:
            dltper_PDs.append(float(val_str))
            
    elif "primal Equality Residual:" in line:
        parts = line.split(';')
        eq_part = parts[0].split(':')[1].strip()
        ineq_part = parts[1].split(':')[1].strip()
        dltper_pres_eqs.append(float(eq_part))
        pres_ineqs.append(float(ineq_part))

# ==========================================
# 3. 绘图 (修改部分)
# ==========================================
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# 修改标题
fig.suptitle('Convergence Analysis of Douglas-Rachford Splitting', fontsize=16)

# --- 子图 1: gap (修改为对数坐标) ---
# 计算绝对值，因为 gap 有负数，不能直接做 log plot
abs_gaps = [abs(g) for g in gaps]

# 使用 semilogy
axs[0, 0].semilogy(iters, abs_gaps, marker='o', markersize=4, linewidth=2, color='blue', label='|gap|')
axs[0, 0].set_title('|Gap| vs. Iteration (Log Scale)')
axs[0, 0].set_xlabel('Iteration')
axs[0, 0].set_ylabel('|Value|')
axs[0, 0].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[0, 0].legend()

# --- 子图 2: dltper_PD ---
axs[0, 1].semilogy(iters, dltper_PDs, marker='s', markersize=4, linewidth=2, color='orange', label='dltper_PD')
axs[0, 1].set_title('dltper_PD vs. Iteration (Log Scale)')
axs[0, 1].set_xlabel('Iteration')
axs[0, 1].set_ylabel('Value')
axs[0, 1].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[0, 1].legend()

# --- 子图 3: dltper_pres_eq ---
axs[1, 0].semilogy(iters, dltper_pres_eqs, marker='^', markersize=4, linewidth=2, color='green', label='dltper_pres_eq')
axs[1, 0].set_title('Primal Eq. Residual vs. Iteration (Log Scale)')
axs[1, 0].set_xlabel('Iteration')
axs[1, 0].set_ylabel('Residual')
axs[1, 0].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1, 0].legend()

# --- 子图 4: pres_ineq ---
axs[1, 1].plot(iters, pres_ineqs, marker='x', markersize=6, linewidth=2, color='red', label='pres_ineq')
axs[1, 1].set_title('Primal Ineq. Residual vs. Iteration')
axs[1, 1].set_xlabel('Iteration')
axs[1, 1].set_ylabel('Residual')
axs[1, 1].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1, 1].set_ylim(-0.1, 1.0) 
axs[1, 1].legend()

# 调整整体布局防止标签重叠
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('convergence_plot.png')