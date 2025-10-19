import numpy as np
from math import cos, sqrt, atan, e, pi
import matplotlib.pyplot as plt

class ISA:
    def __init__(self) -> None:
        self.h_trans = (0, 11_000, 20_000, 32_000, 47_000, 51_000, 71_000, 86_000, 90_000, 105_000)
        self.Tlapse = (-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0, 0.004)

    def get(self, h, T_sea=288.15, P_sea=101325):
        i = 0
        T = T_sea
        P = P_sea
        while True:
            h_cap = min(h, self.h_trans[i+1])
            T_old = T
            T = T + (h_cap - self.h_trans[i])*self.Tlapse[i]

            if abs(self.Tlapse[i]) < 1e-6:
                P = P * e**(-g/R/T*(h_cap-self.h_trans[i]))
            else:
                P = P * (T/T_old)**(-g/self.Tlapse[i]/R)

            rho = P/R/T

            if h <= self.h_trans[i+1]:
                break
            i += 1
        
        return (P, rho, T)

def GetWingWeight(ZFM, n_ult):
    kw = 6.67*10**-3
    bref = 1.905 # m
    bs = b / cos(lambda_halfc) # m
    return (kw*bs**0.75 * (1+sqrt(bref/bs)) * n_ult**0.55 * ((bs/tr)/(ZFM/S))**0.3) * ZFM

atm = ISA()
g = 9.80665 # m/s²
R = 287.05287
# P, rho, T = atm.get(19_000)
rho, T = 1.225, 288.15


lambda_halfc = -3.746203854 *pi/180 # rad
b = 31.4 # m
S = 92.9 # m²
thickness = 0.10 # fraction
MTOM = 18144 # kg
MTOW = MTOM*g # N
ZFM = 7610 # kg
ZFW = ZFM*g # N
EOM = 7250 # kg
EOW = EOM*g # N
M_des = 13481 # kg
taper = 0.18 # fraction
Cr = 5.014574112 # m
tr = Cr * thickness # m
cl_alpha = 5.01545857
cl_max_clean = 1.1925 # -
V_b = 89.1
V_cr = 0.68 * sqrt(1.4*R*T) # m/s (=M*a)
V_d = 1.5 * V_cr
V_s = sqrt(MTOW/S*2/1.225*1/cl_max_clean) # m/s
MGC = (1+taper)/2*Cr

# statistical estimates
# ground
u_hat_d = 25 *0.3048 # m/s
u_hat_cr = 50 *0.3048 # m/s
u_hat_alpha = 66 *0.3048 # m/s

# cruise alt
# u_hat_d = 12.5 *0.3048 # m/s
# u_hat_cr = 25 *0.3048 # m/s
# u_hat_alpha = 38 *0.3048 # m/s

mu_g = 2 * MTOW/S / (rho * g * cl_alpha * MGC)
K = 0.88 * mu_g / (5.3 + mu_g)

u_d = K * u_hat_d
u_cr = K * u_hat_cr
u_alpha = K * u_hat_alpha

delta_n_d = rho * V_d * cl_alpha * u_d / ( 2 * MTOW/S) # check W
delta_n_cr = rho * V_cr * cl_alpha * u_cr / ( 2 * MTOW/S) # check W
delta_n_alpha = 3.04-1

n_max = max(delta_n_d+1, delta_n_cr+1, delta_n_alpha+1)

v_vals = np.linspace(0, 400, 250)
n_s = (v_vals/V_s)**2
delta_n_d_line = rho * v_vals * cl_alpha * u_d / ( 2 * MTOW/S) # check W
delta_n_cr_line = rho * v_vals * cl_alpha * u_cr / ( 2 * MTOW/S) # check W
delta_n_alpha_line = rho * v_vals * cl_alpha * u_alpha / ( 2 * MTOW/S) # check W


plt.figure(figsize=(9, 7))
plt.plot(v_vals, n_s, lw=0.5, color='gray')
plt.plot(v_vals, 1+delta_n_d_line, '--k', lw=0.5)
plt.plot(v_vals, 1+delta_n_cr_line, '--k', lw=0.5)
plt.plot(v_vals, 1+delta_n_alpha_line, '--k', lw=0.5)
plt.plot(v_vals, 1-delta_n_d_line, '--k', lw=0.5)
plt.plot(v_vals, 1-delta_n_cr_line, '--k', lw=0.5)
plt.plot(v_vals, 1-delta_n_alpha_line, '--k', lw=0.5)
plt.plot([0, V_b, V_cr, V_d, V_d, V_cr, V_b , 0], [1, 1+delta_n_alpha, 1+delta_n_cr, 1+delta_n_d, 1-delta_n_d, 1-delta_n_cr, 1-delta_n_alpha, 1], color='black')
plt.scatter([0, V_b, V_cr, V_d, V_b, V_cr, V_d], [1, 1+delta_n_alpha, 1+delta_n_cr, 1+delta_n_d, 1-delta_n_alpha, 1-delta_n_cr, 1-delta_n_d], color='black')
# label the key velocity points
plt.annotate('$V_B$', xy=(V_b, 1+delta_n_alpha), xytext=(V_b, 1+delta_n_alpha+0.25),fontsize=9, ha='center', va='bottom')
plt.annotate('$V_{CR}$', xy=(V_cr, 1+delta_n_cr), xytext=(V_cr, 1+delta_n_cr+0.25),fontsize=9, ha='center', va='bottom')
plt.annotate('$V_D$', xy=(V_d, 1+delta_n_d), xytext=(V_d, 1+delta_n_d+0.25), fontsize=9, ha='center', va='bottom')
plt.ylim(bottom=-4, top=6)
plt.xlim(left=0, right=375)
plt.xlabel('Velocity (m/s)')
plt.ylabel('n')
plt.title('Load factor vs Velocity')
plt.grid(True)
plt.tight_layout()
plt.show(block=False)

WM0 = GetWingWeight(ZFM, 1.5*n_max)
const_term = ZFM - WM0
fuel = MTOM - ZFM

n_ult_new = 1.5*n_max * 0.7

i = 0
iter = []
mtom_list = []
err = 1
while err>1e-6:
    old = ZFM
    WM = GetWingWeight(ZFM, n_ult_new)
    ZFM = WM + const_term
    err = abs(ZFM-old)
    mtom_list.append(ZFM+fuel)
    iter.append(i)
    i += 1

print(mtom_list[-1])
plt.figure(figsize=(9, 7))
plt.plot(iter, mtom_list, 'k', lw=1)
plt.xlim(left=0)
plt.xlabel('Iteration')
plt.ylabel('MTOM [kg]')
plt.title('Algorithm iteration vs Maximum take-off mass')
plt.grid(True)
plt.tight_layout()
plt.show(block=True)