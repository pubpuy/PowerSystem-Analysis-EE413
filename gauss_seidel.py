import numpy as np
import cmath

# -------------------------------------------------
# 1. กำหนดค่าตัวแปร (Input Data)
# -------------------------------------------------
V1 = 1.05 + 0j      # Slack Bus
V2 = 1.0 + 0j       # เริ่มต้น Flat Start
V3 = 1.0 + 0j       # เริ่มต้น Flat Start

S2 = -2.566 - 1.102j 
S3 = -1.386 - 0.452j 

y12 = 1 / (0.02 + 0.04j)
y13 = 1 / (0.01 + 0.03j)
y23 = 1 / (0.0125 + 0.025j)

Y22 = y12 + y23
Y33 = y13 + y23

# -------------------------------------------------
# 2. เริ่มคำนวณ (Part a)
# -------------------------------------------------
print("Solution (E1)")
print("=" * 60)

for k in range(1, 8):
    # หา V2
    term_S2 = np.conj(S2) / np.conj(V2)
    V2 = (term_S2 + y12*V1 + y23*V3) / Y22

    # หา V3
    term_S3 = np.conj(S3) / np.conj(V3)
    V3 = (term_S3 + y13*V1 + y23*V2) / Y33

    # --- ส่วนแสดงผล (Printing Logic) แก้ไขใหม่ ---
    if k == 1:
        # รอบที่ 1: โชว์แค่ V3 ตามรูป
        print(f"V3^(1) = {V3.real:.4f} - j{abs(V3.imag):.4f}")
    
    elif k == 2:
        # รอบที่ 2: โชว์ข้อความหัวข้อ และค่า V2, V3
        print("\nFor the 2nd iteration,")
        print(f"V2^(2) = {V2.real:.4f} - j{abs(V2.imag):.4f}")
        print("and")
        print(f"V3^(2) = {V3.real:.4f} - j{abs(V3.imag):.4f}")
        print("-" * 60) # ขีดเส้นคั่นเพื่อเข้าสู่ตารางสรุป
        
    else:
        # รอบที่ 3-7: โชว์เป็นตารางคู่กัน
        print(f"V2^({k}) = {V2.real:.4f} - j{abs(V2.imag):.4f}      "
              f"V3^({k}) = {V3.real:.4f} - j{abs(V3.imag):.4f}")

print("-" * 60)
print("The final solution is\n")

# --- ปัดเศษทศนิยมให้เหลือ 4 ตำแหน่ง ---
V2_ans = complex(round(V2.real, 4), round(V2.imag, 4))
V3_ans = complex(round(V3.real, 4), round(V3.imag, 4))

mag2, ang2 = cmath.polar(V2_ans)
mag3, ang3 = cmath.polar(V3_ans)

print(f"V2 = {V2_ans.real:.4f} - j{abs(V2_ans.imag):.4f} = {mag2:.5f}∠{np.degrees(ang2):.4f}° pu")
print(f"V3 = {V3_ans.real:.4f} - j{abs(V3_ans.imag):.4f} = {mag3:.5f}∠{np.degrees(ang3):.4f}° pu")

# -------------------------------------------------
# 3. คำนวณ Slack Bus (Part b)
# -------------------------------------------------
print("\n" + "=" * 60)
print("b) The slack bus calculation")
print("-" * 60)

term_I = V1*(y12 + y13) - (y12*V2_ans + y13*V3_ans)
S1_conj = np.conj(V1) * term_I

P1 = S1_conj.real
Q1 = -S1_conj.imag 

print(f"P1 - jQ1 = {P1:.3f} - j{abs(Q1):.3f}")
print(f"P1 = {P1:.3f} pu = {P1 * 100:.1f} MW")
print(f"Q1 = {Q1:.3f} pu = {Q1 * 100:.1f} Mvar")
print("=" * 60)