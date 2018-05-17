import json
from high_throughput.vasp.plot import plt_set

ms = json.load(open('./seebeck_eff_mass_vs_Ef_tensor.json'))
comp = json.load(open('./Fermi_surface_complexicity_factor_vs_Ef_tensor.json'))
mc = json.load(open('./eff_mass_vs_Ef_tensor.json'))

vbm = ms['vbm']

import matplotlib.pyplot as plt
plt_set(plt,width=30,height=40,scale=23)
fig, (ax0,ax2,ax1) = plt.subplots(3,1)
xlim=(-0.4,0.0)


carr_Efs = ms['carr_Ef_pairs']

n_mc = len(mc['mu_steps'])
ax0.plot(mc['mu_steps'],[mc['600'][i][0][0] for i in range(n_mc)],label = '$xx, yy$')
#ax0.plot(mc['mu_steps'],[mc['600'][i][1][1] for i in range(n_mc)],label = '$yy $',linestyle='dashed')
ax0.plot(mc['mu_steps'],[mc['600'][i][2][2] for i in range(n_mc)],label = '$zz $')
ax0.set_xlim(xlim)
ax0.set_ylim(0,9)
ax0.set_ylabel('$m_c (m_e)$')
ax0.legend(loc=1)

#ax0r = ax0.twinx()
n_ms = len(mc['mu_steps'])
ax2.plot(ms['mu_steps'],[ms['600'][i][0][0] for i in range(n_ms)],label = '$xx, yy$')
#ax0.plot(ms['mu_steps'],[ms['600'][i][1][1] for i in range(n_ms)],label = '$yy $',linestyle='dashed')
ax2.plot(ms['mu_steps'],[ms['600'][i][2][2] for i in range(n_ms)],label = '$zz$')
ax2.set_xlim(xlim)
ax2.set_ylim(0,0.6)
ax2.set_ylabel('$m_s (m_e)$')
ax2.legend(loc=2)


n_com = len(comp['mu_steps'])
ax1.plot(comp['mu_steps'],[comp['600'][i][0][0] for i in range(n_com)],label = '$xx, yy$')
#ax1.plot(comp['mu_steps'],[comp['600'][i][1][1] for i in range(n_com)],label = '$yy $',linestyle='dashed')
ax1.plot(comp['mu_steps'],[comp['600'][i][2][2] for i in range(n_com)],label = '$zz $')
ax1.set_xlabel('$E - E_{VBM}(eV)$')
ax1.set_ylabel('$N_v^*K^*$')
ax1.set_xlim(xlim)
ax1.set_ylim(0,1.5)
ax1.legend(loc=2)


size = 90
ax0.axvline(carr_Efs['1e+20'],linewidth=4,linestyle='dashed',color='black')
ax0.text(carr_Efs['1e+20'],1.8,'$10^{20} cm^{-3}$', fontsize=size)
ax0.axvline(carr_Efs['1e+21'],linewidth=4,linestyle='dashed',color='black')
ax0.text(carr_Efs['1e+21'],6.0,'$10^{21} cm^{-3}$', fontsize=size)

ax2.axvline(carr_Efs['1e+20'],linewidth=4,linestyle='dashed',color='black')
ax2.text(carr_Efs['1e+20'],0.1,'$10^{20} cm^{-3}$', fontsize=size)
ax2.axvline(carr_Efs['1e+21'],linewidth=4,linestyle='dashed',color='black')
ax2.text(carr_Efs['1e+21'],0.1,'$10^{21} cm^{-3}$', fontsize=size)

ax1.axvline(carr_Efs['1e+20'],linewidth=4,linestyle='dashed',color='black')
ax1.text(carr_Efs['1e+20'],0.35,'$10^{20} cm^{-3}$', fontsize=size)
ax1.axvline(carr_Efs['1e+21'],linewidth=4,linestyle='dashed',color='black')
ax1.text(carr_Efs['1e+21'],0.35,'$10^{21} cm^{-3}$', fontsize=size)

fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
fig.savefig('Fermi_surface_complexity_power_factor.png')
