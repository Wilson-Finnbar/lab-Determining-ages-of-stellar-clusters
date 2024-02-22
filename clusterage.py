#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uncertainties as uc
from uncertainties import unumpy

A_V = uc.ufloat(14.952,0.007)
B_V = uc.ufloat(18.304,0.008)
C_V = uc.ufloat(17.510,0.007)
D_V = uc.ufloat(17.914,0.007)
E_V = uc.ufloat(15.875,0.006)
F_V = uc.ufloat(16.925,0.007)
G_V = uc.ufloat(18.494,0.008)

zp = pd.read_fwf('zeropointstars.txt', header=None)
print(zp)
AV_zp = uc.ufloat(zp[1][0],zp[2][0])
BV_zp = uc.ufloat(zp[1][1],zp[2][1])
CV_zp = uc.ufloat(zp[1][2],zp[2][2])
DV_zp = uc.ufloat(zp[1][3],zp[2][3])
EV_zp = uc.ufloat(zp[1][4],zp[2][4])
FV_zp = uc.ufloat(zp[1][5],zp[2][5])
GV_zp = uc.ufloat(zp[1][6],zp[2][6])

m_instr = np.array([AV_zp,BV_zp,CV_zp,DV_zp,EV_zp,FV_zp,GV_zp])
m_star = np.array([A_V,B_V,C_V,D_V,E_V,F_V,G_V])

m_zp = []

for i in range(len(m_star)):
    m_zp.append( m_star[i] - m_instr[i])
print(sum(m_zp)/len(m_zp))

for i in range(len(m_star)):
    print(f"& {m_star[i].n : .3f} & {m_star[i].s : .3f} & {m_instr[i].n : .3f} & {m_instr[i].s : .3f} \u005C\u005C")

# %%
vstar = pd.read_fwf('vstars.txt', header=None)
bstar = pd.read_fwf('bstars.txt', header=None)
print(bstar.head())
vmag = vstar[3] - 5*np.log10(6800/10)
vmag_err = vstar[4]
bmag = (bstar[3] - 5*np.log10(6800/10)) - vmag
bmag_err = bstar[4]

plt.errorbar(bmag, vmag,
             yerr = vmag_err,
             xerr = bmag_err, 
             fmt ='r.',
             ecolor = 'k')
plt.xlim(-0.7, 1.7)
plt.ylim(max(vmag), min(vmag))
plt.xlabel('B-V')
plt.ylabel('V')
plt.savefig('B-V.pdf')

# %%
zpbdf = pd.read_fwf('zeropointb.txt', header=None)
print(zpbdf)
bm_instr = unumpy.uarray(zpbdf[3],zpbdf[4])
bm_star = unumpy.uarray([14.952+0.793,18.304+0.440,17.510+0.647,17.914+0.538,15.875-0.058,16.925+0.677,18.494+0.434],[0.010,0.010,0.011,0.011,0.010,0.010,0.010])
bm_zp = bm_star - bm_instr
print(np.mean(bm_zp))

for i in range(len(bm_star)):
    print(f"& {bm_star[i].n : .3f} & {bm_star[i].s : .3f} & {bm_instr[i].n : .3f} & {bm_instr[i].s : .3f} \u005C\u005C")
# %%
plt.errorbar(bmag, vmag,
             yerr = vmag_err,
             xerr = bmag_err, 
             fmt ='r.',
             ecolor = 'k')
plt.xlim(-0.7, 1.7)
plt.ylim(max(vmag), min(vmag))
plt.xlabel('B-V')
plt.ylabel('V')

ubv_iso_logAge, ubv_iso_Bmag, ubv_iso_Vmag = np.loadtxt('UBV.dat',usecols=(2,27,28),unpack=True)

uniq = np.unique(ubv_iso_logAge)

for i in range(0,len(uniq),2):
    age = np.where(ubv_iso_logAge == uniq[i])
    Bmags = ubv_iso_Bmag[age]
    Vmags = ubv_iso_Vmag[age]
    Bmags = Bmags - Vmags
    plt.plot(Bmags[Bmags<1.1],Vmags[Bmags<1.1],'k')
plt.savefig('isofitwhole.pdf')
# %%
print(10**uniq[17]/1000000000)
age = np.where(ubv_iso_logAge == uniq[17])
Bmags = ubv_iso_Bmag[age]
Vmags = ubv_iso_Vmag[age]
Bmags = Bmags - Vmags
plt.plot(Bmags[Bmags<1.1],Vmags[Bmags<1.1],'k')
plt.errorbar(bmag, vmag,
             yerr = vmag_err,
             xerr = bmag_err, 
             fmt ='r.',
             ecolor = 'k')
plt.xlim(-0.7, 1.7)
plt.ylim(max(vmag), min(vmag))
plt.xlabel('B-V')
plt.ylabel('V')
plt.savefig('isofit17.pdf')
# %%
rey = pd.read_fwf('rey_data.txt')
print(rey.columns)

reyv = rey['Vmag'] - 5*np.log10(6800/10)
reyb = ((rey['B-V'] + rey['Vmag']) - 5*np.log10(6800/10)) - reyv

plt.errorbar(reyb, reyv,
             yerr = rey['e_Vm'],
             xerr = rey['e_B-V'], 
             fmt ='r.',
             ecolor = 'k',
             zorder=1)
plt.xlim(-0.7, 1.7)
plt.ylim(max(reyv), min(reyv))
plt.xlabel('B-V')
plt.ylabel('V')

age = np.where(ubv_iso_logAge == uniq[15])
Bmags = ubv_iso_Bmag[age]
Vmags = ubv_iso_Vmag[age]
Bmags = Bmags - Vmags
plt.plot(Bmags[Bmags<1.1],Vmags[Bmags<1.1],'k',zorder=2)
plt.savefig('rey.pdf')
# %%
print(10**uniq[17]/1000000000)
print(10**uniq[15]/1000000000)
print(Bmags[Bmags<1.1])
# %%
vabs = vstar[3] - 5*np.log10(6800/10)
vapp = vstar[3]
v_err = vstar[4]
babs = (bstar[3] - 5*np.log10(6800/10))
bapp = bstar[3]
b_err = bstar[4]

for i in range(len(bmag)):
    print(f"{i+1} & \cellcolor\u007B gray!10 \u007D {vapp[i] : .3f} & {v_err[i] : .3f} & \cellcolor\u007B gray!10 \u007D {bapp[i] : .3f} & {b_err[i] : .3f} & \cellcolor\u007B gray!10 \u007D {vabs[i] : .3f} & {v_err[i] : .3f} & \cellcolor\u007B gray!10 \u007D {babs[i] : .3f} & {b_err[i] : .3f} \u005C\u005C")
# %%
