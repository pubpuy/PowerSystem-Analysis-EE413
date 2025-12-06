import numpy as np
import cmath

# ==========================================
# 1. ตั้งค่าตัวแปรและระบบ (System Setup)
# ==========================================
# กำหนดค่าเริ่มต้น
V_mag = np.array([1.05, 1.0, 1.0])  # [V1, V2, V3]
delta = np.array([0.0, 0.0, 0.0])   # [d1, d2, d3] (radians)

# กำหนด Load (ค่าลบ)
P_sch = np.array([0.0, -2.566, -1.386]) 
Q_sch = np.array([0.0, -1.102, -0.452])

# ==========================================
# 2. สร้าง Y-Bus Matrix
# ==========================================
y12 = 1 / (0.02 + 0.04j)
y13 = 1 / (0.01 + 0.03j)
y23 = 1 / (0.0125 + 0.025j)

# สร้าง Matrix 3x3
Y = np.zeros((3, 3), dtype=complex)
Y[0, 1] = Y[1, 0] = -y12
Y[0, 2] = Y[2, 0] = -y13
Y[1, 2] = Y[2, 1] = -y23
Y[0, 0] = y12 + y13
Y[1, 1] = y12 + y23
Y[2, 2] = y13 + y23

G = Y.real
B = Y.imag

# ==========================================
# 3. เริ่มคำนวณ Newton-Raphson (Part a)
# ==========================================
print("Solution Example 1 (Newton-Raphson Method)")
print("=" * 60)
print("a) Iteration Process")
print("-" * 60)

for iteration in range(1, 20):
    # --- 3.1 หาค่า Mismatch ---
    P_calc = np.zeros(3)
    Q_calc = np.zeros(3)
    
    for i in range(3):
        for j in range(3):
            th = delta[i] - delta[j]
            P_calc[i] += V_mag[i] * V_mag[j] * (G[i, j]*np.cos(th) + B[i, j]*np.sin(th))
            Q_calc[i] += V_mag[i] * V_mag[j] * (G[i, j]*np.sin(th) - B[i, j]*np.cos(th))

    # Mismatch (Bus 2, 3)
    mismatch = np.array([
        P_sch[1] - P_calc[1], P_sch[2] - P_calc[2],
        Q_sch[1] - Q_calc[1], Q_sch[2] - Q_calc[2]
    ])
    
    # Check Convergence
    if np.max(np.abs(mismatch)) < 1e-9:
        print(f"Converged in {iteration-1} iterations.")
        break

    # --- 3.2 สร้าง Jacobian Matrix (J) ---
    J = np.zeros((4, 4))
    indices = [1, 2] # Bus 2, 3
    
    for r, i in enumerate(indices):
        for c, j in enumerate(indices):
            # J1: dP/dd
            if i == j: J[r, c] = -Q_calc[i] - (B[i,i] * V_mag[i]**2)
            else:      J[r, c] = V_mag[i] * V_mag[j] * (G[i,j]*np.sin(delta[i]-delta[j]) - B[i,j]*np.cos(delta[i]-delta[j]))
            
            # J2: dP/dV
            if i == j: J[r, c+2] = (P_calc[i]/V_mag[i]) + (G[i,i] * V_mag[i])
            else:      J[r, c+2] = V_mag[i] * (G[i,j]*np.cos(delta[i]-delta[j]) + B[i,j]*np.sin(delta[i]-delta[j]))
            
            # J3: dQ/dd
            if i == j: J[r+2, c] = P_calc[i] - (G[i,i] * V_mag[i]**2)
            else:      J[r+2, c] = -V_mag[i] * V_mag[j] * (G[i,j]*np.cos(delta[i]-delta[j]) + B[i,j]*np.sin(delta[i]-delta[j]))
            
            # J4: dQ/dV
            if i == j: J[r+2, c+2] = (Q_calc[i]/V_mag[i]) - (B[i,i] * V_mag[i])
            else:      J[r+2, c+2] = V_mag[i] * (G[i,j]*np.sin(delta[i]-delta[j]) - B[i,j]*np.cos(delta[i]-delta[j]))

    # --- 3.3 Solve Correction ---
    dx = np.linalg.solve(J, mismatch)
    delta[1] += dx[0]; delta[2] += dx[1]
    V_mag[1] += dx[2]; V_mag[2] += dx[3]

# --- สรุปผลข้อ a ---
print("-" * 60)
print("The final solution (Part a)\n")

# แปลงเป็น Rectangular -> ปัดเศษ -> แปลงกลับเป็น Polar (เพื่อให้ตรงเฉลย)
V1_rect = cmath.rect(V_mag[0], delta[0])
V2_rect = cmath.rect(V_mag[1], delta[1])
V3_rect = cmath.rect(V_mag[2], delta[2])

V2_ans = complex(round(V2_rect.real, 4), round(V2_rect.imag, 4))
V3_ans = complex(round(V3_rect.real, 4), round(V3_rect.imag, 4))

mag2, ang2 = cmath.polar(V2_ans)
mag3, ang3 = cmath.polar(V3_ans)

print(f"V2 = {V2_ans.real:.4f} - j{abs(V2_ans.imag):.4f} = {mag2:.5f}∠{np.degrees(ang2):.4f}° pu")
print(f"V3 = {V3_ans.real:.4f} - j{abs(V3_ans.imag):.4f} = {mag3:.5f}∠{np.degrees(ang3):.4f}° pu")

# ==========================================
# 4. คำนวณ Slack Bus (Part b)
# ==========================================
print("\n" + "=" * 60)
print("b) The slack bus calculation")
print("-" * 60)

# สูตร: S1* = V1* [ I1_injected ]
# I1_injected = V1(y12 + y13) - (y12*V2 + y13*V3)

# ใช้ค่า V1, V2_ans, V3_ans (ที่ปัดเศษแล้ว)
V1_val = V1_rect # V1 ไม่ต้องปัดเศษเพราะเป็นค่าคงที่ (1.05)

term_I = V1_val*(y12 + y13) - (y12*V2_ans + y13*V3_ans)
S1_conj = np.conj(V1_val) * term_I

P1 = S1_conj.real
Q1 = -S1_conj.imag # S* = P - jQ -> Q = -Imag

print(f"P1 - jQ1 = {P1:.3f} - j{abs(Q1):.3f}")
print(f"P1 = {P1:.3f} pu = {P1 * 100:.1f} MW")
print(f"Q1 = {Q1:.3f} pu = {Q1 * 100:.1f} Mvar")
print("=" * 60)