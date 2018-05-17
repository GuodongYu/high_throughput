import os
import json
from high_throughput.vasp.plot import plt_set

root_dir = '/home/gyu/workplace/intrinsic_defects/Zintl_Antimonides/Yb14AlSb11/Yb_2/bulk/sta/8x8x8/boltztrap'
root_dir = '/home/gyu/workplace/intrinsic_defects/Zintl_Antimonides/Yb14MnSb11/Yb_2/AF/magnetic_phase1/U_Mn_d/U_3eV/8x8x8/transport/boltztrap/semi/boltztrap'
seeb_eff_mass_file = 'seebeck_eff_mass_vs_Ef_doping_level.json'
cond_eff_mass_file = 'eff_mass_vs_Ef.json'
fermi_surf_comp_fact_file = 'Fermi_surface_complexicity_factor_vs_Ef.json'

def plot_eff_mass_and_Fermi_surf_complexity_factor():
    import matplotlib.pyplot as plt
    plt_set(plt)
    fig, (ax0, ax1) = plt.subplots(2,1)
    seeb_eff_mass = json.load(open(os.path.join(root_dir,seeb_eff_mass_file)))
    cond_eff_mass = json.load(open(os.path.join(root_dir,cond_eff_mass_file)))
    fermi_surf_fact = json.load(open(os.path.join(root_dir,fermi_surf_comp_fact_file)))
    Efs = seeb_eff_mass['carr_Ef_pairs']
    legends=[]
    legend_mc,=ax0.plot(cond_eff_mass['mu_steps'],cond_eff_mass['600'],label='$m_c (m_0)$')
    legend_ms,=ax0.plot(seeb_eff_mass['mu_steps'],seeb_eff_mass['600'],label='$m_s (m_0)$')
    legends.append(legend_mc);legends.append(legend_ms)
    ax0.legend(handles=legends,loc=2)
    legend_fac,=ax1.plot(fermi_surf_fact['mu_steps'],fermi_surf_fact['600'],label='$N_v^*K^*$')
    ax1.legend(handles=[legend_fac],loc=2)
    ax0.axvline(Efs['1e+20'],color='black',linestyle='dashed',linewidth=3);ax0.text(Efs['1e+20'],1,'$10^{20}cm^{-3}$')
    ax0.axvline(Efs['1e+21'],color='black',linestyle='dashed',linewidth=3);ax0.text(Efs['1e+21'],1,'$10^{21}cm^{-3}$')
    ax1.axvline(Efs['1e+20'],color='black',linestyle='dashed',linewidth=3);ax1.text(Efs['1e+20'],0.2,'$10^{20}cm^{-3}$')
    ax1.axvline(Efs['1e+21'],color='black',linestyle='dashed',linewidth=3);ax1.text(Efs['1e+21'],0.2,'$10^{21}cm^{-3}$')
    ax0.set_xlim([-0.5,0]);ax1.set_xlim([-0.5,0])
    ax0.set_ylim([0,5]);ax1.set_ylim([0,0.8])
    #ax0.set_xlabel('Fermi energy (eV)');ax1.set_xlabel('Fermi energy (eV)')
    ax1.set_xlabel('$E-VBM (eV)$')
    fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})
    plt.tight_layout()
    plt.savefig('Fermi_surface_complexity_factor.svg')
    plt.clf()
    
plot_eff_mass_and_Fermi_surf_complexity_factor()
