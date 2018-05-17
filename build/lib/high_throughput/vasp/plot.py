from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Rectangle
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.electronic_structure.plotter import BSPlotter, \
                    BSPlotterProjected, plot_brillouin_zone, HighSymmKpath,DosPlotter
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.units import *
from pymatgen.core.periodic_table import Element
from high_throughput.vasp.transport import seebeck_eff_mass_from_seebeck_carr
from high_throughput.config import *
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import OrbitalType,Spin
import json
import glob
import os
from monty.json import jsanitize
import numpy as np
from operator import itemgetter
from os.path import isfile, join, isdir
from scipy.constants import e,m_e


def plt_set(plt,nrow=1,ncolumn=1,width=20,height=20,scale=10):
    fig_width = width*ncolumn
    fig_height = height*nrow
    fig_size = [fig_width, fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 4*scale,
              'font.size': 3*scale,
              'font.weight': 'normal',
              'legend.fontsize': 4*scale,
              'xtick.labelsize': 4*scale,
              'ytick.labelsize': 4*scale,
              'text.usetex': False,
              'figure.figsize': fig_size,
              'lines.markersize': 2*scale,
              'lines.linewidth':  0.7*scale
              }
    plt.rcParams.update(params)
    #plt.tight_layout()

def show_charge_curves(chargefile):
    """
    Args:
        chargefile: charge formated file, such as CHGCAR and PARCHG 
    """
    import matplotlib.pyplot as plt
    data=Chgcar.from_file(chargefile)
    fig,axes=plt.subplots(3,1)
    plt_set(plt,3,1)
    axes[0].plot(data.get_axis_grid(0),data.get_average_along_axis(0),label='a_direction')
    axes[1].plot(data.get_axis_grid(1),data.get_average_along_axis(1),label='b_direction')
    axes[2].plot(data.get_axis_grid(2),data.get_average_along_axis(2),label='c_direction')
    plt.show()
    plt.clf()

def show_brillouin_zone(stru_source):
    import matplotlib
    matplotlib.use('Tkagg')
    try:
        stru = Structure.from_file(stru_source)
    except:
        stru=m.get_structure_by_material_id(stru_source)
    kpath = HighSymmKpath(stru)
    labels = kpath.kpath['kpoints']
    plot_brillouin_zone(stru.lattice.reciprocal_lattice,labels=kpath.kpath['kpoints']).show()


def show_BZ_from_BS_run(bs_vasprun_file):
    import matplotlib
    matplotlib.use('Tkagg')
    xml=Vasprun(bs_vasprun_file)
    bs=xml.get_band_structure(line_mode=True)
    stru=xml.final_structure
    kpts={k.label:k.frac_coords for k in bs.kpoints if k.label}
    plot_brillouin_zone(stru.lattice.reciprocal_lattice,labels=kpts).show()


class Dosplot(object):
    def __init__(self,dir):
        self.dir = dir
        self.run = Vasprun(join(self.dir,'vasprun.xml'))
        pdos = dict(zip(self.run.final_structure,self.run.pdos))
        self.complete_dos = CompleteDos(self.run.final_structure,self.run.complete_dos,pdos)

    def get_OrbitalType(self,orb):
        if orb == 's':
            return OrbitalType.s
        elif orb == 'p':
            return OrbitalType.p
        elif orb == 'd':
            return OrbitalType.d
        elif orb == 'f':
            return OrbitalType.f


    def plot_elt_spd(self,orbits,with_tot=True,output_dir='./',xlim=None,fmt='png'):
        """
        orbits: ['Yb_p','Yb_d','Sb_p'] like
        """
        dosp = DosPlotter()
        if with_tot:
            dosp.add_dos('Total',self.run.tdos)
        for i in orbits:
            ele,orb=(i.split('_')[0],i.split('_')[1])
            orbitaltype = self.get_OrbitalType(orb)
            dosp.add_dos(i,self.complete_dos.get_element_spd_dos(ele)[orbitaltype])
        if xlim is None:
            xlim=[-4,4]
        dosp.save_plot(os.path.join(output_dir,'pdos')+'.%s' % fmt,img_format=fmt,xlim=xlim)

    def plot_elt(self,with_tot=True,output_dir='./',xlim=None,fmt='png'):
        dos = DosPlotter()
        if with_tot:
            dos.add_dos('Total',self.run.tdos)
        dos_elt = self.complete_dos.get_element_dos()
        
        for i in [Element(j) for j in ['Sb','Al','Ca']]:
            elt = i.symbol
            dos.add_dos(elt,dos_elt[i])
        if xlim is None:
            xlim=[-4,4]
        plt = dos.get_plot(xlim)
        plt_set(plt)
        plt.savefig(os.path.join(output_dir,'dos_on_elt')+'.%s' % fmt, format=fmt)
        plt.clf()


class BandPlot(object):

    def __init__(self, dir):
        """
        Class to many kinds of band structure and density of states plots. 
        Args:
            dir: directory for Vasp band structure calculation
            project: whether plot projected band structure on elt or orbital
            ylim: energy range for plot
            fmt: picture format to save the band structure
        """
        self.dir = dir
        self.run=Vasprun(join(self.dir,'vasprun.xml'))
        self.project=False
        self.bs = self.run.get_band_structure(line_mode=True)
        self.structure = self.run.final_structure
        self.formula = self.structure.composition.reduced_formula

    def _open_project(self):
        self.run=Vasprun(join(self.dir,'vasprun.xml'),parse_projected_eigen=True)
        self.bs = self.run.get_band_structure(line_mode=True)
        self.project=True

    def plot_picked_bands(self,nth=[],ylim=(-0.5,1.5),spin='up',fmt='png',line_val=[]):
        """
        Args:
            nth: a list for choosed bands counted from 1 
            spin: 'up' or 'down' for spin_up and spin_down
            fmt: picture format
        """
        import matplotlib.pyplot as plt
        plt_set(plt)
        labels=[i.label for i in self.bs.kpoints]
        nk=len(labels)
        inds_del=[] ### the inds_del will be removed
        for i in range(len(labels)-1):
            if labels[i] and labels[i+1] and labels[i]==labels[i+1]:
                inds_del.append(i+1)
            if labels[i] and labels[i+1] and labels[i]!=labels[i+1]:
                labels[i]=labels[i]+'|'+labels[i+1]
                labels[i+1]=labels[i]
        labels_new=[labels[i] for i in range(nk) if i not in inds_del]

        nk_new=len(labels_new)
        i,j=(0,0)
        xlim=[]
        while j<=nk_new-1:
            if labels_new[i] and labels_new[i+1]:
                xlim.append(i);xlim.append(i)
                labels_new[i+1]=None
                j+=2
            else:
                xlim.append(i)
                j+=1
            i+=1

        inds_break=[] ### at some points band breaks, inds_break save the inds just after the break
        inds_break.append(0)
        for i in range(len(xlim)-1):
            if xlim[i]==xlim[i+1]:
                inds_break.append(i+1)
        inds_break.append(None)

        labels_final=[]
        for label in labels_new:
            if label:
                labels_final.append('$'+label+'$')
            else:
                labels_final.append('')

        if spin == 'up':
            bs = self.bs.bands[Spin.up]
        elif spin == 'down':
            bs = self.bs.bands[Spin.down]
        else:
            raise KeyError("Only 'up' or 'donw' can be accepted! ")
        bs_pick = []
        for i in nth:
            bsi = bs[i-1]
            bsi_pick = [bsi[i] for i in range(len(bsi)) if i not in inds_del]
            bs_pick.append(bsi_pick)
        legends=[]
        lines=['-','--','-.',':']
        colors=['red','blue','orange','green','black','yellow']
        j=0;k=0
        for i in range(len(nth)):
            nbs=np.array(bs_pick[i])-self.bs.efermi
            for nk in range(len(inds_break)-1):
                start = inds_break[nk]
                end = inds_break[nk+1]
                legend,= plt.plot(xlim[start:end],nbs[start:end],color=colors[k],linestyle=lines[j],label=str(nth[i])+'$^{th}$')
            legends.append(legend)
            if not j%len(lines) and j:
                j=0
            else:
                j+=1
            if not k%len(colors) and k:
                k=0
            else:
                k+=1
        for i in range(nk_new):
            if labels_final[i]:
                plt.axvline(xlim[i],color='black',linestyle="-",linewidth=2)
        for i in line_val:
            plt.axhline(i,color='red',linestyle="--",linewidth=2)
        plt.xticks(xlim, labels_final)
        plt.xlim(0,len(set(xlim))-1)
        plt.ylim(ylim)
        plt.legend(handles=legends, loc=1)
        plt.savefig(join(self.dir,'%s_picked_up_bs.%s' % (str(self.formula),fmt)),format=fmt)
        plt.clf()
        

    def plot_band_structure(self,fmt='png',ylim=(-1,1),smooth=True,line_val=None,vbm_cbm_marker=False):
        """
        fmt: figure format for save band structure
        ylim: energy range for plot
        smooth: whether smooth band lines
        line_val: plot an addition line at this energy if not None
        """
        from matplotlib import pyplot as plt
        plt_set(plt)
        fig,(ax) = plt.subplots()
        ymajorLocator   = MultipleLocator(0.5)
        pl = BSPlotter(self.bs).get_plot(vbm_cbm_marker=vbm_cbm_marker, smooth=smooth)
        ax = pl.axes()
        if line_val:
            plt.plot([-100,100],[line_val,line_val])
        ax.set_ylim(ylim)
        ax.set_xlabel('Wave Vector')
        ax.set_ylabel('$E-E_f (eV)$')
        ax.yaxis.set_major_locator(ymajorLocator)
        plt.tight_layout()
        plt.savefig(join(self.dir,'%s_bs.%s' % (str(self.formula),fmt)),format=fmt)
        plt.clf()
    
    def plot_projected_bs_on_elt(self, gap_scissor=False, gap_after_scissor=None, color_mode=False, fmt='pdf', ylim=(-1,1)):
        """
        Args:
            color_mode: True, the different colors represent the different projections
                        False, one subfigure santds for one projection
        """
        from matplotlib import pyplot as plt
        if not self.project:
            self._open_project()
        if gap_scissor:
            bs_ = self.bs.apply_scissor(gap_after_scissor)
        else:
            bs_ = self.bs
        if color_mode:
            mark='color'
            proj_bands = BSPlotterProjected(bs_).get_elt_projected_plots_color(elt_ordered=[Element('Bi'),Element('Ba'),Element('Te')])
        else:
            mark='subplots'
            proj_bands = BSPlotterProjected(bs_).get_elt_projected_plots()
        plt.ylim(ylim)
        plt.savefig(join(self.dir,'%s_projected_bs_on_elt_%s.%s' % (str(self.formula), mark, fmt)), format=fmt)
 
    def plot_projected_bs_on_orbit(self, dictio, fmt='jpg', ylim=(-1,1),smooth=True):
        """
        Args:
            dictio: orbital dict for projection. for example {'Cu':['d','s'],'O':['p']}
        """
        if self.project==False:
            self._open_project()
        proj_bands = BSPlotterProjected(self.bs).get_projected_plots_dots(dictio,ylim=ylim)
        proj_bands.savefig(join(self.dir,'%s_projected_bs_on_orbit.%s' % (str(self.formula), fmt)), format=fmt,smooth=smooth)

class BoltztrapDataPlot(object):
    def __init__(self,dir):
        self.dir = dir
        self.an = BoltztrapAnalyzer.from_files(self.dir)
        temps=self.an.mu_doping['p'].keys()
        temps.sort()
        self.temps=temps
        outfile = glob.glob(join(self.dir,'*outputtrans'))[0]
        with open(outfile) as f:
            data = f.read().split('\n')
            vbm_ry = float([i for i in data if 'VBM:' in i][0].split()[1])
            self.vbm = vbm_ry*Ry_to_eV
        self.cbm = self.vbm+self.an.gap
        self.gap = self.an.gap
        self.nstep = len(self.an.mu_steps)
        self.colors=['red','blue','green','orange','purple','cyan','pink']*3
        self.marks=['D','o','^','H','s','p']*3

    def get_fermi_energy_from_doping_level(self,n,T):
        doping_steps = self.an.get_carrier_concentration() 
        abs_diff=np.absolute(np.array(doping_steps[T])-n)
        ind=min(enumerate(abs_diff), key=itemgetter(1))[0]
        mu = self.an.mu_steps[ind]
        #print 'Doping level is %s at Fermi level %s' % (str(doping_steps[T][ind]),str(mu))
        return mu, ind

    def get_Fermi_energy_from_Hall_carr(self,nH,T):
        hall_carr = self.an.get_hall_carrier_concentration()
        ind = self._get_index_of_hall_carr_at_T(nH, T,self.vbm,self.cbm)
        mu = self.an.mu_steps[ind]
        print 'Hall carr is %s at Fermi level %s' % (str(-hall_carr[T][ind]),str(mu))
        return mu

    def _get_quantity_and_unit(self, prop, output, doping_levels,relaxation_time=1.e-14,kl=1.0):
        """
        For effective mass, please know:
        (1) Only constant time relaxation time calculation can be applied.
        (2) Effective mass plot of only semiconductor is right because of 
            the disaggrement between carrier concentration and doping level.
        """
        if prop == 'seebeck':
            quantity=self.an.get_seebeck(output, doping_levels)
            unit= '($\mu VK^{-1}$)'
        elif prop == 'conductivity':
            quantity = self.an.get_conductivity(output, doping_levels, relaxation_time)
            unit = '($1/\Omega m$)'
        elif prop == 'eff_mass':
            quantity = self.an.get_average_eff_mass(output,doping_levels)
            unit = '$m_e$'
        elif prop == 'resistivity':
            cond = self.an.get_conductivity(output, doping_levels,relaxation_time)
            if not doping_levels:
                quantity = {T:1.0e+5/np.array(cond[T]) for T in self.temps}
            else:
                quantity = {}
                for type_dope in ['p','n']:
                    quantity[type_dope] = {T:1.0e+5/np.array(cond[type_dope][T]) for T in self.temps}
            unit = '($m\Omega cm$)'
        elif prop == 'thermal_conductivity':
            quantity = self.an.get_thermal_conductivity(output, doping_levels,relaxation_time)
            unit = '(W/mK)'
        elif prop == 'power_factor':
            quantity = self.an.get_power_factor(output,doping_levels,relaxation_time)
            unit ='($\mu W/mK^2$)'
        elif prop == 'zT':
            quantity = self.an.get_zt(output, doping_levels,relaxation_time,kl)
            unit = ''
        else:
            raise ValueError('Not allowed prop!')

        ### convert to np.array format ###
        if not doping_levels:
            quantity={T:np.array(quantity[T]) for T in self.temps}
        else:
            for doping_type in ['p', 'n']:
                quantity[doping_type]={T:np.array(quantity[doping_type][T]) for T in self.temps}
        return (quantity,unit)

    def plot_property_vs_temperature_fixdoping(self, prop, dope_type, exp_data_file=None, fmt='jpg',xlim=None, ylim=None):
        """
        Plot property average over three directions against temperature for the fixed doping mode.
        Fixed doping mode means the doping levels are just the ones given in the Boltztrap inputfile
        boltztrap.intrans. And physical quantities are read from the output files ended by _fixdoping. 
        Args:
            prop: physical quantity to plot. Followings are allowed: 
                          seebeck, conductivity, resistivity, thermal_conductivity, power_factor, zt, eff_mass
            dope_type: n or p
            exp_data_file: file storing experimental data with the formate like this:
                            100, 50.1
                            200, 80.1
                            .....    
                            first (temperateur) and second (quantity) columns are separated by ','.  
                                    
        """
        import matplotlib.pyplot as plt
        doping=self.an.doping[dope_type]
        if dope_type == 'n':
            legend_dope = [-1*i for i in doping]
        else:
            legend_dope =doping
        if doping==[]:
            raise ValueError('No %s type doping was calculated!' % dope_type)
        num_doping=len(doping)
        quantity, unit = self._get_quantity_and_unit(prop,doping_levels=True)
        quantity=quantity[dope_type]
        quantity_avg={T:np.array(quantity[T]).mean(axis=1) for T in self.temps}
        legends=[]
        for i in range(num_doping):
            legend,=plt.plot(self.temps,[quantity_avg[T][i] for T in self.temps],
                             color=self.colors[i],marker=self.marks[i],label=legend_dope[i])
            legends.append(legend)
        if exp_data_file:
            with open(exp_data_file) as f:
                contents = f.read().split('\n')
                temps_exp = [float(i.split(',')[0]) for i in contents if i]
                quantity_exp = [float(i.split(',')[1]) for i in contents if i]
            legend,=plt.plot(temps_exp, quantity_exp,
                                    color='black',marker='o',label='experiment')
            legends.append(legend)
        plt.legend(handles=legends, loc=1)
        plt.xlabel('Temperature (K)')
        plt.ylabel('%s %s' % (prop, unit))
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.savefig('%s_vs_T_fixdoping.%s' % (prop,fmt),format=fmt)
        plt.clf()

    def plot_property_vs_temperature_input_doping(self, prop, output, dopings, exp_data_files=[], fmt='jpg',relaxation_time=1e-14,kl=1.0,output_json=True,xlim=None,ylim=None,width=20,height=20):
        """
        Args:

            dopings: list giving the doping levels you want. For example: [1e+19, 1e+20, -1e+20]
                    plus and negative mean the holes and electrons doping respectively.
            exp_data_file: file saving experimental data with the formate like this:
                            100, 50.1
                            200, 80.1
                            .....    
                            first (temperateur) and second (seebeck) columns are separated by ','.  
        """
        import matplotlib.pyplot as plt
        plt_set(plt,width=width, height=height)
        fig,ax = plt.subplots(1,1)
        avg = False
        eig = False
        if output in ['avg','avg+eig','eig+avg']:
            avg = True
            if prop == 'zT':
                quantity_steps_avg,unit=self._get_quantity_and_unit(prop, 'average', doping_levels=False,relaxation_time=relaxation_time,kl=kl)
            else:
                quantity_steps_avg,unit=self._get_quantity_and_unit(prop, 'average', doping_levels=False,relaxation_time=relaxation_time)
        if output in ['eig','avg+eig','eigs','eig+avg']:
            eig = True
            if prop == 'zT':
                quantity_steps_eig,unit=self._get_quantity_and_unit(prop, 'eig', doping_levels=False,relaxation_time=relaxation_time,kl=kl)
            else:
                quantity_steps_eig,unit=self._get_quantity_and_unit(prop, 'eig', doping_levels=False,relaxation_time=relaxation_time)
        doping_steps=self.an.get_carrier_concentration()
        if avg:
            quantity_avg = {}
        if eig:
            quantity_eig = {}
        Ef_comms = {}
        for T in self.temps:
            if avg:
                quantity_avg[T]=[]
            if eig:
                quantity_eig[T]=[]
            Ef_comms[T]=[]
            for i in range(len(dopings)):
                abs_diff=np.absolute(np.array(doping_steps[T])-dopings[i])
                ind=min(enumerate(abs_diff), key=itemgetter(1))[0]
                mu = self.an.mu_steps[ind]
                if dopings[i] > 0:
                    Ef_comm = 'v+%s' % str(round(mu - self.vbm,2)) if mu - self.vbm >0 else 'v%s' % str(round(mu - self.vbm,2))
                else:
                    Ef_comm = 'c+%s' % str(round(mu - self.cbm,2)) if mu - self.cbm >0 else 'c%s' % str(round(mu - self.cbm,2))
                print 'T = %sK, input doping = %s cm-3, real doping = %s cm-3' % (T, str(dopings[i]), str(doping_steps[T][ind]))
                diff_min = abs(dopings[i]-doping_steps[T][ind])/abs(dopings[i])
                if diff_min > 0.1:
                    print 'Warning: doping difference is large, please consider to rerun boltztrap after decreasing energygrid!!'
                if eig:
                    s_eig=quantity_steps_eig[T][ind]
                    quantity_eig[T].append(s_eig)
                if avg:
                    s_avg=quantity_steps_avg[T][ind]
                    quantity_avg[T].append(s_avg)
                Ef_comms[T].append(Ef_comm)
        for i in range(len(dopings)):
            if output in ['avg','avg+eig','average']:
                ax.plot(self.temps, [quantity_avg[T][i] for T in self.temps], color=self.colors[i],
                        marker=self.marks[i],label='avg %.4g $cm^{-3}$' % dopings[i])
                for ll in self.temps:
                    ax.text(ll,quantity_avg[ll][i],Ef_comms[ll][i],fontsize=15,horizontalalignment='center',verticalalignment='center')
            if output in ['eig','eigs','avg+eig','eig+avg']:
                if avg:
                    w = 3*i+1
                else:
                    w = 3*i
                if not avg:
                    ax.plot([0],[0],color='white',label = "doping: %.4g $cm^{-3}$" % dopings[i])
                ax.plot(self.temps, [quantity_eig[T][i][0] for T in self.temps], color=self.colors[w],
                            marker=self.marks[w],label='0')
                ax.plot(self.temps, [quantity_eig[T][i][1] for T in self.temps], color=self.colors[w+1],
                            marker=self.marks[w+1],label='1')
                ax.plot(self.temps, [quantity_eig[T][i][2] for T in self.temps], color=self.colors[w+2],
                            marker=self.marks[w+2],label='2')
        if exp_data_files:
            shapes = ['o','^']
            ii = 1
            for exp_data_file in exp_data_files:
                with open(exp_data_file) as f:
                    contents = f.read().split('\n')
                    temps_exp = [float(i.split(',')[0]) for i in contents if i]
                    quantity_exp = [float(i.split(',')[1]) for i in contents if i]
                ax.plot(temps_exp, quantity_exp,
                                    color='black',marker=shapes[ii-1],linewidth=0,label='exp %i' % ii)
                ii+=1
        if output_json:
            json_data={}
            json_data['Temps']=self.temps
            json_data['dopings']=dopings
            json_data['relaxation_time']=relaxation_time
            json_data[prop]={}
            for i in range(len(dopings)):
                if avg:
                    json_data[prop]['average'] = []
                    json_data[prop]['average'].append([quantity_avg[T][i] for T in self.temps])
                if eig:
                    json_data[prop]['eigs'] = []
                    json_data[prop]['eigs'].append([quantity_eig[T][i] for T in self.temps])
            with open(prop+'.json','w') as f:
                f.write(json.dumps(jsanitize(json_data)))

        ax.legend()
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('%s %s' % (prop.replace('eff','effective'),unit))
        ax.set_title('$\\tau$ %e' % relaxation_time)
        plt.tight_layout()
        plt.savefig('%s_vs_T_doping_level.%s' % (prop,fmt),format=fmt)
        plt.clf()


    def _get_index_of_hall_carr_at_T(self, HC0, T0,vbm_fake,cbm_fake):
        """
        Get the index to make BoltztrapAnalyzer.get_hall_carrier_concentration()[T0][index]=HC0
        """
        hall_carr_T0 = -1.0*np.array(self.an.get_hall_carrier_concentration()[T0])
        mu_steps_vbm_ref = np.array(self.an.mu_steps)-vbm_fake
        mu_steps_cbm_ref = np.array(self.an.mu_steps)-cbm_fake
        ind_vbm = min(enumerate(np.absolute(mu_steps_vbm_ref)), key=itemgetter(1))[0]
        print 'vbm index: %i' % ind_vbm
        ind_cbm = min(enumerate(np.absolute(mu_steps_cbm_ref)), key=itemgetter(1))[0]
        print 'cbm index: %i' % ind_cbm
        ######### choose ind range for p type #########
        if HC0>0:
            ind_tmp = ind_vbm
            while hall_carr_T0[ind_tmp] > 0.0:
                ind_begin_p = ind_tmp
                ind_tmp = ind_tmp-1
                if hall_carr_T0[ind_begin_p] > HC0:
                    break
            #while hall_carr_T0[ind_tmp] < 0.0:
            #    ind_begin_p = ind_tmp
            #    ind_tmp = ind_tmp+1
            ind_tmp = ind_begin_p +1
            while hall_carr_T0[ind_tmp] > hall_carr_T0[ind_tmp+1]:
                ind_end_p = ind_tmp +1
                ind_tmp = ind_tmp+1
                if hall_carr_T0[ind_end_p] < HC0:
                    break
            ind_range_p = [ind_begin_p, ind_end_p]
        print ind_range_p
        ################################################
        ######### choose ind range for n type ##########
        if HC0 <0:
            ind_tmp = ind_cbm
            while hall_carr_T0[ind_tmp] < 0.0:
                ind_end_n = ind_tmp
                ind_tmp = ind_tmp+1
                if hall_carr_T0[ind_end_n]<HC0:
                    break
            ind_tmp = ind_cbm
            while hall_carr_T0[ind_tmp] < hall_carr_T0[ind_tmp-1]:
                ind_begin_n = ind_tmp -1
                ind_tmp = ind_tmp-1
                if hall_carr_T0[ind_begin_n]>HC0:
                    break
            ind_range_n = [ind_begin_n, ind_end_n]
        ################################################
        if HC0>0:
            hall_carr_T0_cut=hall_carr_T0[ind_range_p[0]:ind_range_p[1]+1]
            ind=min(enumerate(np.absolute(hall_carr_T0_cut-HC0)),key=itemgetter(1))[0]
            ind=ind+ind_range_p[0]
        if HC0<0:
            hall_carr_T0_cut=hall_carr_T0[ind_range_n[0]:ind_range_n[1]+1]
            ind=min(enumerate(np.absolute(hall_carr_T0_cut-HC0)),key=itemgetter(1))[0]
            ind=ind+ind_range_n[0]
        print 'T=%sK, input and real hall carrs are %s and %s cm-3' % (str(T0), str(HC0), str(hall_carr_T0[ind]))
        return ind

    def plot_property_vs_temperature_input_hall_carrier_shot(self, prop, output, HC0, T0=300, exp_data_file=None, fmt='jpg',relaxation_time=1e-14,output_json=True,xlim=None,ylim=None,**kws):
        """
        from
                 f(T,mu,doping)=0 and g(T,mu,hall_carr)=0, 
        it gives
                 w(T,doping,hall_carr)=0
                 
        After fixing hall_carr and T to be HC0 and T0, doping level can be determined.
        Then the function will plot the property at this doping level over all temperature.
        Args:
            HC0: choosed hall carrier concentration
            T0: choosed temperature
        """
        try:
            vbm=kws['vbm']
            cbm=kws['cbm']
        except:
            vbm=self.vbm
            cbm=self.cbm
        if T0 not in self.temps:
            raise ValueError('T0 = %s K was not calculated!!' % T0)
        ind = self._get_index_of_hall_carr_at_T( HC0, T0, vbm, cbm)
        carr_T0=self.an.get_carrier_concentration()[T0]
        doping=carr_T0[ind]
        print 'T=%sK, doping=%s \n' % (str(T0), str(doping))
        self.plot_property_vs_temperature_input_doping(prop,output,[doping],exp_data_file=exp_data_file, fmt=fmt,relaxation_time=relaxation_time,output_json=output_json,xlim=xlim,ylim=ylim)

    def plot_property_vs_temperature_input_hall_carrier(self, prop, hall_carrs, exp_data_file=None, fmt='jpg',relaxation_time=1e-14,**kws):
        """
        plot the property at Hall carrier concentration over the whole temperature.
        Args:
            hall_carrs: Hall carrier concentration list.
        """
        try:
            vbm=kws[vbm];cbm=kws['cbm']
        except:
            vbm=self.vbm;cbm=self.cbm
        import matplotlib.pyplot as plt
        plt_set(plt)
        #quantity_steps,unit=self._get_quantity_and_unit(prop,doping_levels=False)
        quantity_avg_steps,unit=self._get_quantity_and_unit(prop, output='average', doping_levels=False,relaxation_time=relaxation_time)
        #quantity_avg_steps={T:quantity_steps[T].mean(axis=1) for T in self.temps}
        quantity_avg={}
        for HC in hall_carrs:
            for T in self.temps:
                quantity_avg[T]=[]
                try:
                    ind=self._get_index_of_hall_carr_at_T(HC,T,vbm,cbm)
                    quantity_avg[T].append(quantity_avg_steps[T][ind])
                except:
                    quantity_avg[T].append(None)
        legends=[]
        for i in range(len(hall_carrs)):
            legend,=plt.plot(self.temps, [quantity_avg[T][i] for T in self.temps], color=self.colors[i],
                            marker=self.marks[i],label=hall_carrs[i])
            legends.append(legend)

        if exp_data_file:
            with open(exp_data_file) as f:
                contents = f.read().split('\n')
                temps_exp = [float(i.split(',')[0]) for i in contents if i]
                quantity_exp = [float(i.split(',')[1]) for i in contents if i]
            legend,=plt.plot(temps_exp, quantity_exp,
                                    color='black',marker='o',linewidth=0, label='experiment')
            legends.append(legend)
        plt.legend(handles=legends, loc=1)
        plt.xlabel('Temperature (K)')
        plt.ylabel('%s %s' % (prop,unit))
        plt.savefig('%s_vs_T_Hall_carr.%s' % (prop,fmt),format=fmt)
        plt.clf()

    def plot_prop_vs_doping_level(self,prop='zT',output='avg',inds=[0,1,2],Temps=[300],relaxation_time = 1.0e-14,kl=1.0,xlim=[1e+19,1e+21],ylim=None,write_json=True,width=20,height=20,scale=10):
        """
        Args:
            prop: can be 'seebeck', 'resistivity', 'conductivity', 'power_factor', 'zT'
            output: average avg or eig
            index: the directory of the matrix. It will be ignored when output is average
            Temps: a list with all temperatures for plot
            relaxation_time: relaxation_time scale timing on the boltztrap calculation
            kl: lattice thermal conductivity for zT, in W/(m*K)
        """
        from matplotlib import pyplot as plt
        if output in ['average','avg']:
            otuput = 'average'
        elif output in ['eig','eigs']:
            output = 'eigs'
        plt_set(plt,width=width,height=height,scale=scale)
        fig,ax = plt.subplots(1,1)
        if prop == 'zT':
            data,unit=self._get_quantity_and_unit(prop, output, doping_levels=False,relaxation_time=relaxation_time,kl=kl)
        else:
            data,unit=self._get_quantity_and_unit(prop, output, doping_levels=False,relaxation_time=relaxation_time)
        carrs = self.an.get_carrier_concentration()
        for T in Temps:
            if output == 'average':
                ax.plot(carrs[T],data[T],label = '%i K' % T)
            elif output == ['eigs','eig']:
                for ii in inds:
                    ax.plot(carrs[T],[i[ii] for i in data[T]],label = '%i, %iK' % (ii,T))
            elif output == 'tensor':
                ii_labels = ['xx','yy','zz']
                for ii in inds:
                    ax.plot(carrs[T],[i[ii][ii] for i in data[T]],label = '%s, %iK' % (ii_labels[ii], T))
            else:
                raise KeyError('output can not be accepted')

        ax.legend()
        ax.set_xlabel('Doping level ($cm^{-3}$)')
        ax.set_ylabel(prop+' %s' % unit)
        ax.set_xlim(xlim)
        if write_json:
            json_out = {}
            json_out['relaxation_time'] = relaxation_time
            json_out['units']={'Temperature':'K', 'carrier_concentration':'cm-3', 'power_factor':'microW m-1 K-2','relaxation_time':'s'}
            for T in Temps:    
                json_out[T] = {}
                json_out[T]['carriers']=carrs[T]
                json_out[T][prop] = data[T]
            with open(prop+'_carr_%s.json' % output,'w') as f:
                f.write(json.dumps(jsanitize(json_out)))
                                            
        if ylim is not None:
            ax.set_ylim(ylim)
        plt_set(plt)
        plt.savefig(prop+'_vs_carr_%s.svg' % output)
        plt.clf()

    def plot_Propmax_vs_doping_level(self,prop,output='avg', index=0,doping_level_range=[1.0e+18,1.0e+22], T_range=[300,1300],relaxation_time=1.0e-14,xlim=None,ylim=None,write_json=True,width=20,height=20):
        """
        prop: zT, power_factor ...
        output: avg or eigs
        index: 0, 1, 2 when output is eigs. This tag will be ignored when output is 'avg' or 'average'
        doping_level_range: the xlim range
        """
        from matplotlib import pyplot as plt
        plt_set(plt,width=width,height=height)
        fig,ax = plt.subplots(1,1)
        if output in ['avg','average']:
            output = 'average'
        if output in ['eig','eigs']:
            output = 'eigs'
        data,unit=self._get_quantity_and_unit(prop, output, doping_levels=False,relaxation_time=relaxation_time)
        carrs = self.an.get_carrier_concentration()
        data_max = []
        carr_max = []
        temps = [t for t in self.temps if t <= max(T_range) and t >= min(T_range)]
        for T in temps:
            if output == 'average':
                data_T = data[T]
            elif output == 'eigs':
                data_T = [i[index] for i in data[T]]
            carrs_T = carrs[T]
            pair = [i for i in zip(data_T,carrs_T) if i[1] < max(doping_level_range) and i[1]>min(doping_level_range)]
            pair_sort = sorted(pair,key=lambda x:x[0],reverse=True)
            data_max_T, carr_max_T = pair_sort[0]
            data_max.append(data_max_T)
            carr_max.append(carr_max_T)

        ax.scatter(carr_max,data_max)
        for i in range(len(carr_max)):
            ax.text(carr_max[i],data_max[i],'%i K' % temps[i])
        ax.set_xlabel('Doping level ($cm^{-3}$)')
        ax.set_ylabel(prop+'_max %s' % unit)
        ax.set_xlim(xlim)
        if write_json:
            json_out={}
            json_out['units']={'Temperature':'K', 'carrier_concentration':'cm-3', prop:unit,'relaxation_time':'s'}
            json_out['units']={'Temperature':'K', 'carrier_concentration':'cm-3', prop:unit,'relaxation_time':'s'}
            json_out['temperatures'] = temps
            json_out['carriers']= carr_max
            json_out[prop+"_max"] = data_max
            json_out['relaxation_time'] = relaxation_time
            with open(prop+'_max_vs_carr.json','w') as f:
                f.write(json.dumps(jsanitize(json_out)))
                                            
        if ylim is not None:
            ax.set_ylim(ylim)
        plt_set(plt)
        plt.savefig(prop+'_max_vs_carr.jpg')
        plt.clf()
        

    def get_seebeck_eff_mass_vs_Ef(self,temps,Lambda,output):
        carr = self.an.get_carrier_concentration()
        seebeck = self.an.get_seebeck(output, doping_levels=False)
        quantity_steps = {}
        for T in temps:
            quantity_steps[T]=[]
            seebeck_T = seebeck[T]
            carr_T = carr[T]
            for i in range(self.nstep):
                if output == 'average':
                    seeb_mass = seebeck_eff_mass_from_seebeck_carr(abs(seebeck_T[i]), carr_T[i], T, Lambda)
                elif output in ['eigs','eig']:
                    seeb_mass = []
                    for j in [0,1,2]:
                        seeb_mass_j = seebeck_eff_mass_from_seebeck_carr(abs(seebeck_T[i][j]), carr_T[i], T, Lambda)
                        seeb_mass.append(seeb_mass_j)
                elif output == 'tensor':
                    seeb_mass = [[0,0,0] for __ in range(3)]
                    for j in [0,1,2]:
                        for k in [0,1,2]:
                            seeb_mass_jk = seebeck_eff_mass_from_seebeck_carr(abs(seebeck_T[i][j][k]), carr_T[i], T, Lambda)
                            seeb_mass[j][k] = seeb_mass_jk
                quantity_steps[T].append(seeb_mass)
        return quantity_steps

    def plot_seebeck_eff_mass_vs_Ef(self,output='average',inds=[0,1,2],temps=[300],Lambda=0.5,carr_markers=[1e+20,-1e+20],carr_marker_T=300,fmt='png',xlim=None,ylim=None,ref_vbm=True,json_write=True):
        from matplotlib import pyplot as plt
        plt_set(plt)
        fig,ax = plt.subplots(1,1)
        if xlim is None:
            xlim=(self.vbm-0.5,self.cbm+0.5)
        unit = '($m_e$)'
        mu_steps=np.array(self.an.mu_steps)-self.vbm if ref_vbm else np.array(self.an.mu_steps)
        carr_marker_coll = {}
        for n in carr_markers:
            mu = self.get_fermi_energy_from_doping_level(n,carr_marker_T)[0]
            carr_marker_coll[n] = mu - self.vbm if ref_vbm else mu
        quantity_steps = self.get_seebeck_eff_mass_vs_Ef(temps, Lambda, output)
        for T in temps:
            if output == 'average':
                ax.plot([mu_steps[i] for i in range(self.nstep)],quantity_steps[T],label='avg, %sK' % T)
            elif output in ['eig','eigs']:
                output = 'eigs'
                for ii in inds:
                    ax.plot([mu_steps[i] for i in range(self.nstep)],[j[ii] for j in quantity_steps[T]],label='%i, %sK' % (ii,T))
            elif output == 'tensor':
                ii_labels = ['xx','yy','zz']
                for ii in inds:
                    ax.plot([mu_steps[i] for i in range(self.nstep)],[j[ii][ii] for j in quantity_steps[T]],label='%s, %sK' % (ii_labels[ii],T))
            else:
                raise KeyError('%s is unknown' % output)
                    
        if json_write:
            json_data={}
            json_data['mu_steps']=mu_steps
            json_data['vbm'] = 0.0 if ref_vbm else self.vbm
            json_data['carr_Ef_pairs'] = carr_marker_coll
            for T in temps:
                json_data[T]=quantity_steps[T]
            with open('seebeck_eff_mass_vs_Ef_%s.json' % output,'w') as f:
                f.write(json.dumps(jsanitize(json_data)))

        for i in carr_marker_coll:
            carr  = i
            mu = carr_marker_coll[i]
            ax.axvline(mu,color='black',linewidth=1,linestyle='dashed')
            if not ylim:
                ylim = plt.axes().get_ylim()
            ax.text(mu,(ylim[0]+ylim[1])/2.0,str(carr))
            
        ax.legend(loc=1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if ref_vbm:
            ax.axvline(0,color='black',linewidth=0.5)
            ax.axvline(self.gap,color='black',linewidth=0.5)
        else:
            ax.axvline(self.vbm,color='black',linewidth=0.5)
            ax.axvline(self.cbm,color='black',linewidth=0.5)
        ax.axhline(0.0,color='black',linewidth=0.5)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('Seebeck effective mass')
        fig.savefig('Seebeck_effective_mass_vs_Ef_%s.%s' % (output,fmt),format=fmt)
        plt.clf()

    def plot_Fermi_surface_complexity_factor_vs_Ef(self,output='average',inds=[0,1,2],temps=[300],Lambda=0.5,fmt='png',xlim=None,ylim=None,ref_vbm=True,json_write=True):
        from matplotlib import pyplot as plt
        plt_set(plt)
        fig,ax=plt.subplots(1,1)
        if xlim is None:
            xlim=(self.vbm-0.5,self.cbm+0.5)
        if ylim is None:
            ylim = (0,10)
        mu_steps=np.array(self.an.mu_steps)-self.vbm if ref_vbm else np.array(self.an.mu_steps)
        seebeck_eff_mass = self.get_seebeck_eff_mass_vs_Ef(temps,Lambda,output)
        cond_eff_mass,unit = self._get_quantity_and_unit('eff_mass', output=output, doping_levels=False)
        comp_fac={}
        for T in temps:
            comp_fac[T]=(np.array(seebeck_eff_mass[T])/np.abs(np.array(cond_eff_mass[T])))**1.5
            if output =='average':
                ax.plot(mu_steps,comp_fac[T],label='avg, T=%sK' % T)
            elif output in ['eig','eigs']:
                output = 'eigs'
                for i in inds:
                    ax.plot(mu_steps,[j[i] for j in comp_fac[T]],label='%i, T=%sK' % (i,T))
            elif output == 'tensor':
                ii_labels = ['xx','yy','zz']
                for i in inds:
                    ax.plot(mu_steps,[j[i][i] for j in comp_fac[T]],label='%s, T=%sK' % (ii_labels[i],T))
            else:
                raise KeyError('%s is unknown' % output)

        if json_write:
            with open('Fermi_surface_complexicity_factor_vs_Ef_%s.json' % output,'w') as f:
                json_data={}
                for T in temps:
                    json_data[T]= comp_fac[T]
                json_data['mu_steps']=mu_steps
                f.write(json.dumps(jsanitize(json_data)))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.legend(loc=1)
        if ref_vbm:
            ax.axvline(0,color='black',linewidth=0.5)
            ax.axvline(self.gap,color='black',linewidth=0.5)
        else:
            ax.axvline(self.vbm,color='black',linewidth=0.5)
            ax.axvline(self.cbm,color='black',linewidth=0.5)
        ax.axhline(0.0,color='black',linewidth=0.5)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('complexity factor')
        plt.savefig('Fermi_surface_complexcity_factor_vs_Ef_%s.%s' % (output,fmt),format=fmt)
        plt.clf()



    def plot_property_vs_Ef(self, prop, output='average', inds=[0,1,2],temps=[300], fmt='png', xlim=None, ylim=None, relaxation_time=1.0e-14,ref_vbm=True,json_write=True):
        """
        Args:
            prop: the property plotted, 'zT', 'eff_mass', 'conductivity', 'resistivity', 'seebeck'
            output: 'average', 'eigs', 'tensor'
            inds: [0,1,2], [0], [0,1] the directions you want to plot when output is 'eigs' or 'tensor'.
            temps: a list with temperatures
            fmt: the format of the plot
        """
        from matplotlib import pyplot as plt
        fig,ax = plt.subplots(1,1)
        plt_set(plt)
        if xlim is None:
            xlim=(self.vbm-0.5,self.cbm+0.5)
        if ylim is None:
            ylim=(0,2)
        if prop == 'seebeck_eff_mass':
            print 'Plese use function plot_seebeck_eff_mass_vs_Ef !'
            return
        quantity_steps,unit=self._get_quantity_and_unit(prop,output,doping_levels=False,relaxation_time=relaxation_time)
        mu_steps=np.array(self.an.mu_steps)-self.vbm if ref_vbm else np.array(self.an.mu_steps)
        for T in temps:
            quantity=quantity_steps[T]
            if output == 'average':
                ax.plot(mu_steps,quantity,label='T=%sK' % T)
            elif output in ['eig','eigs']:
                output = 'eigs'
                for ii in inds:
                    ax.plot(mu_steps,[quantity[i][ii] for i in range(self.nstep)],label='%d T=%sK' % (ii,T))
            elif output =='tensor':
                ii_labels = ['xx','yy','zz']
                for ii in inds:
                    ax.plot(mu_steps,[quantity[i][ii][ii] for i in range(self.nstep)],label='%s T=%sK' % (ii_labels[ii],T))
        if json_write:
            json_data={}
            json_data['vbm'] = 0.0 if ref_vbm else self.vbm
            json_data['gap'] = self.gap
            json_data['mu_steps']=mu_steps
            for T in temps:
                json_data[T] = quantity_steps[T]
            f=open(prop+'_vs_Ef_%s.json' % output,'w')
            f.write(json.dumps(jsanitize(json_data)))
        ax.legend(loc=1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if ref_vbm:
            ax.axvline(0,color='black',linewidth=0.5)
            ax.axvline(self.gap,color='black',linewidth=0.5)
        else:
            ax.axvline(self.vbm,color='black',linewidth=0.5)
            ax.axvline(self.cbm,color='black',linewidth=0.5)
        ax.axhline(0.0,color='black',linewidth=0.5)
        ax.set_xlabel('Fermi level (eV)')
        ax.set_ylabel('%s %s' % (prop.replace('eff_','effective '),unit))
        fig.savefig('%s_vs_Ef_%s.%s' % (prop,output,fmt),format=fmt)
        plt.clf()

    def plot_LorenzNumber_vs_Ef(self, temps=[300], fmt='png', xlim=None, ylim=None):
        from matplotlib import pyplot as plt
        plt_set(plt)
        if not xlim:
            xlim=(self.vbm-0.5,self.cbm+0.5)
        thermal_cond=self._get_quantity_and_unit('thermal_conductivity',doping_levels=False)[0]
        thermal_cond_avg={T:thermal_cond[T].mean(axis=1) for T in self.temps}
        cond=self._get_quantity_and_unit('conductivity',doping_levels=False)[0]
        cond_avg={T:cond[T].mean(axis=1) for T in self.temps}
        legends=[]
        if temps == 'all':
            temps=self.temps
        for T in temps:
            ln=np.array(thermal_cond_avg[T])/(float(T)*np.array(cond_avg[T]))
            legend,=plt.semilogy(self.an.mu_steps,ln,label='T=%sK' % T)
            legends.append(legend)
        plt.legend(handles=legends,loc=1)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.axvline(self.vbm,color='black',linewidth=0.5)
        plt.axvline(self.cbm,color='black',linewidth=0.5)
        plt.plot([xlim[0],xlim[1]],[0,0],'k',linewidth=0.5)
        plt.xlabel('Fermi level (eV)')
        plt.ylabel('Lorenz number (%s)' % '$W\Omega K^{-2}$')
        plt.tight_layout()
        plt.savefig('Lorenz_number_vs_Ef.%s' % fmt,format=fmt)
        plt.clf()

    def plot_doping_hall_carr_comparison(self,temp,fmt='png',xlim=None,ylim=None,vlines=[],hlines=[],width=20, height=20, scale=10):
        from matplotlib import pyplot as plt
        fig,ax = plt.subplots(1,1)
        plt_set(plt,width=width, height=height, scale=scale)
        if not xlim:
            xlim=(self.vbm-0.5,self.cbm+0.5)
        if not ylim:
            ylim=(0,1e+22)
        dopings = np.abs(self.an.get_carrier_concentration()[temp])
        hall_carr = np.abs(self.an.get_hall_carrier_concentration()[temp])
        ax.semilogy(self.an.mu_steps,dopings,label='$n$')
        ax.semilogy(self.an.mu_steps,hall_carr,label='$n_H$')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.axvline(self.vbm,color='black',linewidth=2.0)
        ax.axvline(self.cbm,color='black',linewidth=2.0)
        for ver in vlines:
            if ver<0:
                xv = self.vbm +ver
            elif ver > 0:
                xv = self.cbm+ver
            ax.axvline(xv,color='red',linewidth=1.0)
        for hl in hlines:
            ax.axhline(hl,color='red',linewidth=1.0)
        ax.legend()
        ax.set_xlabel('Fermi level (eV)')
        fig.tight_layout()
        fig.savefig('doping_hall_carr_comparison.%s' % fmt,format=fmt)
        plt.clf()

def relaxation_time(scattering_type,T,tau0,energy_shift):
    KbT0 = 8.6173303e-5*300. # T0 is 300 K
    Kb = 8.6173303e-5
    if scattering_type == 'ac':
        tau = tau0 * (energy_shift/(Kb*T))**-0.5 * (T/300.)**-1.5
    elif scattering_type == 'ionic':
        tau = tau0* (energy_shift/(KbT0))**1.5
    return tau


def plot_average_relaxation_time_vs_temperature(dirs,dir0='',method='tau_ef',dopings=[1e+20],taus=[],xlim=None,ylim=None,marks=[],write_json=True):
    """
    To get the average of relaxation time over energy vs temperature.
    Args:
        dirs: a list for directories of the boltztrap calculations with nonCRTA
        dir0: the directory of the boltztrap calculation with the 1s CRTA, if method is 'tau_ef', this parameter can be ignored
        method: the method to get average relaxation time, 'tau_ef' and 'cond'
            'tau_ef', average of relaxation is chosen to be the relaxation at Ef. Good for heavy-doped semoconductor. 
            'cond', average of relaxation is choosen to be the cond_nonCRTA/cond0_1s_CRTA
        dopings: a list of doping levels for getting the average of relaxation time
        taus: the times affecting on the boltztrap input. For example, in boltztrap calculation
            tau = tau0 (e/kT0)^r (T/T0)^s, the final tau = relaxatin_time*tau
        marks: the label for each relaxation time model in dirs, 'ac': acoustic phonon scattering; 'ionic': ionic impurity scattering
    """
    if len(dirs) != len(marks) or len(dirs) != len(taus):
        raise KeyError("Number of marks != Number of boltztrap calculations or the numbers of taus!!")
    from matplotlib import pyplot as plt
    fig,ax = plt.subplots(1,1)
    plt_set(plt)
    if method == 'cond':
        bzp0 = BoltztrapDataPlot(dir0)
        temps0 = bzp0.temps
    tau_avg = {}
    i = 0
    json_out = {}
    json_out['unit']='fs'
    for dir1 in dirs:
        bzp1 = BoltztrapDataPlot(dir1)
        temps1 = bzp1.temps
        if method == 'cond':
            temps = [ll for ll in temps0 if ll in temps1]
        else:
            temps = temps1
        json_out[marks[i]]={}
        json_out[marks[i]]['temperatures']=temps
        for dope in dopings:
            json_out[marks[i]][dope]=[]
            C = str(dope).split('e')[0][:3]
            e = str(dope).split('e')[-1].split('+')[-1]
            dope_str = '$'+C+'\\times 10^{%s} cm^{-3}$' % e
            tau_avg[dope]=[]
            for T in temps:
                mu1, ind1 = bzp1.get_fermi_energy_from_doping_level(dope,T)
                real_dope1 = bzp1.an.get_carrier_concentration()[T][ind1]
                print "%s: %i K,  %e" % (marks[i], T, real_dope1)
                if method == 'cond':
                    mu0, ind0 = bzp0.get_fermi_energy_from_doping_level(dope,T)
                    real_dope0 = bzp0.an.get_carrier_concentration()[T][ind0]
                    print "CRTA: %i K,  %e" % (T, real_dope0)
                    cond1 = bzp1.an.get_conductivity(output='average',doping_levels=False,relaxation_time=taus[i]*1.e+15)[T][ind1] 
                    cond0 = bzp0.an.get_conductivity(output='average',doping_levels=False,relaxation_time=1)[T][ind0] 
                    tau_avg_T = cond1/cond0
                if method == 'tau_ef':
                    if dope > 0:
                        ene_shift = abs(bzp1.vbm - mu1)
                    else:
                        ene_shift = abs(bzp1.cbm - mu1)
                    tau_avg_T = relaxation_time(marks[i],T,taus[i],ene_shift)
                tau_avg[dope].append(tau_avg_T)
                json_out[marks[i]][dope].append(tau_avg_T)
            ax.plot(temps,tau_avg[dope],label = marks[i] + ' %s' % dope_str)
        i += 1
            
    if write_json:
        with open('average_tau_vs_temperature.json','w') as f:
            f.write(json.dumps(jsanitize(json_out)))
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.legend()
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('relaxation time average (fs)')
    plt.savefig('Average_of_relaxation_time_vs_temperature_%s.jpg' % method)



class LobsterPlot(object):
    def __init__(self):
        pass

    def ReadCOHPData(self,cohpfile):
        with open(cohpfile) as f:
            data = f.read().split('\n')
        label_start_OK = False; label_end_OK = False
        for line_num, line in enumerate(data):
            if 'No.' in line and 'No.' not in data[line_num-1]:
                label_start = line_num
                label_start_OK = True
            if 'No.' in line and 'No.' not in data[line_num+1]:
                label_end = line_num
                label_end_OK = True
            if label_start_OK and label_end_OK:
                break
        #labels = ['average']+[i.split('(')[0][3:].replace('>','') for i in data[label_start:label_end+1]]
        labels = ['average']+[i.split(':')[0][3:].replace('>','') for i in data[label_start:label_end+1]]
        npairs = len(labels)
        cohps_orig = [[float(j) for j in i.split()] for i in data[label_end+1:]]
        energies = [i[0] for i in cohps_orig]
        pcohps = [];ipcohps=[]
        for i in range(npairs):
            pcohps.append([j[2*i+1] for j in cohps_orig])
            ipcohps.append([j[2*(i+1)] for j in cohps_orig])
        out = {}
        out['labels'] = labels
        out['pcohp'] = pcohps
        out['icohp'] = ipcohps
        out['energies'] = energies
        return out

    def plot_coph(self,cohpfile,Nos,out_to='./',fmt='png',xlim=None,ylim=None):
        """
        cohpfile: COHPCAR.lobster file
        Nos: list including the bond pairs, like [3,4,5,8] or ['1-5',6,9], 0 means the average 
        """
        nums = []
        for No in Nos:
            if '-' in No:
                start = int(No.split('-')[0])
                end = int(No.split('-')[1])
                nums += range(start,end+1)
            else:
                nums += int(No)
        import matplotlib.pyplot as plt
        plt_set(plt)
        data = self.ReadCOHPData(cohpfile)
        legends = []
        cols = ['red','blue','green','orange','purple','gray','yellow']*10
        i=0
        for n in nums:
            legend,=plt.plot(-np.array(data['pcohp'][n]),data['energies'],label = data['labels'][n],color=cols[i])
            legends.append(legend)
            i+=1
        plt.xlabel('-pCOHP')
        plt.ylabel('$E-E_f(eV)$')
        plt.legend(handles=legends,loc=1)
        plt.axvline(0,color='black',linewidth=2)
        plt.axhline(0,color='black',linewidth=2)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.tight_layout()
        plt.savefig(os.path.join(out_to,'pcohp.'+fmt),format=fmt)
        plt.clf()




