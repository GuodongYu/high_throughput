#!/usr/bin/env python

__author__ = "Geoffroy Hautier, Bharat Medasani, Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

import numpy as np
import matplotlib
matplotlib.use('agg')
from pymatgen.util.plotting_utils import get_publication_quality_plot
import json

class DefectPlotter(object):
    """
    Class performing all the typical plots from a defects study
    """

    def __init__(self, analyzer):
        """
        Args:
            analyzer: DefectsAnalyzer object 
        """

        self._analyzer = analyzer

    def get_plot_form_energy(self, xlim=None, ylim=(-1,3), ncolumn=3):
        """
        Formation energy vs Fermi energy plot
        Args:
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
        Returns:
            a matplotlib object

        """
        if xlim is None:
            xlim = (-0.2, self._analyzer._band_gap+0.5)
        max_lim = xlim[1]
        min_lim = xlim[0]
        nb_steps = 10000
        x = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nb_steps)
        y = {}
        trans_level_pt = {}
        trans_level_with_ends_pt = {}
        text_q_labels_pt = {}
        text_q_xy_pt = {}
        for t in self._analyzer._get_all_defect_types():
            y_tmp = []
            trans_level = []
	    text_q_labels = []
	    text_q_y_ends = []
            prev_min_q, cur_min_q = None, None
            for x_step in x:
                min = 10000
                for i, dfct in enumerate(self._analyzer._defects):
                    if dfct.name == t:
                        val = self._analyzer._formation_energies[i] + \
                                dfct.charge*x_step
                        if val < min:
                            min = val
                            cur_min_q = dfct.charge
		if x_step==x[0]:
		    for i, dfct in enumerate(self._analyzer._defects):
			if dfct.name == t and dfct.charge==cur_min_q:
			    text_q_y_begin=self._analyzer._formation_energies[i] + \
                                dfct.charge*x_step
		if x_step==x[-1]:
		    for i, dfct in enumerate(self._analyzer._defects):
                       if dfct.name == t and dfct.charge==cur_min_q:
                            text_q_y_end=self._analyzer._formation_energies[i] + \
                                dfct.charge*x_step
		if cur_min_q != prev_min_q:
		    text_q_labels.append(cur_min_q)
                if prev_min_q is not None:
                    if cur_min_q != prev_min_q:
                        trans_level.append((x_step, min))
                prev_min_q = cur_min_q
                y_tmp.append(min)
		y[t]=y_tmp
		
            trans_level_pt[t] =  trans_level
	    trans_level_with_ends = [(x[0],text_q_y_begin)] + trans_level + [(x[-1],text_q_y_end)]
	    trans_level_with_ends_pt[t] = trans_level_with_ends
	    text_q_labels_pt[t] = text_q_labels
	    text_q_xy_pt[t]=[( (trans_level_with_ends[i][0]+trans_level_with_ends[i-1][0])/2., 
                            (trans_level_with_ends[i][1]+trans_level_with_ends[i-1][1])/2.)  
			     for i in range(1,len(trans_level_with_ends))]
	    for i in range(1,len(trans_level_with_ends)):
		check1 = trans_level_with_ends[i][1] < ylim[0] and trans_level_with_ends[i-1][1] < ylim[0]
		check2 = trans_level_with_ends[i][1] > ylim[1] and trans_level_with_ends[i-1][1] > ylim[1]
		if check1 or check2:
		    text_q_labels_pt[t][i-1]=''
		check3 = ylim[0] < trans_level_with_ends[i][1] < ylim[1] and trans_level_with_ends[i-1][1] < ylim[0]
		check4 = ylim[0] < trans_level_with_ends[i][1] < ylim[1] and trans_level_with_ends[i-1][1] > ylim[1]
		check5 = ylim[0] < trans_level_with_ends[i-1][1] < ylim[1] and trans_level_with_ends[i][1] < ylim[0]
		check6 = ylim[0] < trans_level_with_ends[i-1][1] < ylim[1] and trans_level_with_ends[i][1] > ylim[1]
		checkallout1 = trans_level_with_ends[i-1][1] < ylim[0] and trans_level_with_ends[i][1] > ylim[1]
		checkallout2 = trans_level_with_ends[i][1] < ylim[0] and trans_level_with_ends[i-1][1] > ylim[1]
		if check3 or check5:
		    for idfc, dfct in enumerate(self._analyzer._defects):
                        if dfct.name == t and dfct.charge==text_q_labels[i-1]:
                            xedge=(ylim[0]-self._analyzer._formation_energies[idfc])/dfct.charge
			    if check3:
			    	x0=(xedge+trans_level_with_ends[i][0])/2.0
			    	y0=(ylim[0]+trans_level_with_ends[i][1])/2.0
			    if check5:
				x0=(xedge+trans_level_with_ends[i-1][0])/2.0
				y0=(ylim[0]+trans_level_with_ends[i-1][1])/2.0
                            text_q_xy_pt[t][i-1]=(x0,y0)
		if check4 or check6:
		    for idfc, dfct in enumerate(self._analyzer._defects):
                        if dfct.name == t and dfct.charge==text_q_labels[i-1]:
                            xedge=(ylim[1]-self._analyzer._formation_energies[idfc])/dfct.charge
			    if check4:
			    	x1=(xedge+trans_level_with_ends[i][0])/2.0
                            	y1=(ylim[1]+trans_level_with_ends[i][1])/2.0
			    if check6:
				x1=(xedge+trans_level_with_ends[i-1][0])/2.0
                                y1=(ylim[1]+trans_level_with_ends[i-1][1])/2.0
                            text_q_xy_pt[t][i-1]=(x1,y1)
		if checkallout1 or checkallout2:
		    for idfc, dfct in enumerate(self._analyzer._defects):
                        if dfct.name == t and dfct.charge==text_q_labels[i-1]:
                            xedge1=(ylim[0]-self._analyzer._formation_energies[idfc])/dfct.charge
			    xedge2=(ylim[1]-self._analyzer._formation_energies[idfc])/dfct.charge
			    x0=(xedge1+xedge2)/2.0
			    y0=(ylim[0]+ylim[1])/2.0
			    text_q_xy_pt[t][i-1]=(x0,y0)
        def del_num(name):
            new=''
            for i in name:
                if not i.isdigit():
                    new+=i
            return new
	def get_min_lines(keys):
	    keys_min=[i for i in keys]
	    for i in keys:
		for j in keys:
		    if j==i:
			continue
		    yi=np.array(y[i]);yj=np.array(y[j])
	            if set(yi>=yj)==set([True]):
	                keys_min.remove(i)
			break 
	    keys_min_absolute=[i for i in keys_min]

	    all=set(keys_min_absolute)
	    done=set([])
	    left= all-done
	    for i in all:
		done.add(i)
		for j in all-done:
		    yi=np.array(y[i])+0.1;yj=np.array(y[j])
		    if set(yi>=yj)==set([True]):
			keys_min.remove(i)
			break
		    else:
			done=done-set([i])
	    return keys_min
        main_keys=set([del_num(i) for i in y.keys()])
        key_group={}
        for i in main_keys:
            key_group[i]=[]
        for key in y.keys():
            main_key=del_num(key)
            key_group[main_key].append(key)
        for main_key in key_group:
	    left=get_min_lines(key_group[main_key])
	    for key in key_group[main_key]:
		if key not in left:
		    del y[key];del trans_level_pt[key];del text_q_labels_pt[key]
##################
        from matplotlib import pyplot as plt
        from math import ceil
        nline_in_subplot=20
        adict=locals()
        nsubplot=int(ceil(len(y)/float(nline_in_subplot)))  #### subplot numbers
        nrow=int(ceil(nsubplot/float(ncolumn)))
        #plt = get_publication_quality_plot(12, 3*nrow)
        #fig_width = 18
        #fig_height = nrow*3
        fig_width = 8*ncolumn
        fig_height = 6*nrow
        if nsubplot <= ncolumn:
            fig_width = 8*nsubplot
            fig_height = 6
        fig_size = [fig_width, fig_height]
        params = {'backend': 'ps',
                'axes.labelsize': 25,
                'font.size': 25,
                'font.weight': 'light',
                'legend.fontsize': 15,
                'xtick.labelsize': 25,
                'ytick.labelsize': 25,
                'text.usetex': False,
                'figure.figsize': fig_size,
                'lines.markersize': 5}
        plt.rcParams.update(params)
        if nsubplot > 1:
            fig,axes=plt.subplots(nrow,ncolumn,sharey=True)
            n=0
            for i in range(nrow):
                for j in range(ncolumn):
                    if nrow>1:
                        adict['ax%s' % n]=axes[i][j]
                    elif nrow==1:
                        adict['ax%s' % n]=axes[j]
                    adict['legend%s' % n]=[]
                    n=n+1
            if nrow*ncolumn>nsubplot:
                for i in range(nsubplot,nrow*ncolumn):
                    fig.delaxes(adict['ax%s' % str(i)])
        elif nsubplot==1:
            fig,adict['ax0']=plt.subplots(1,1)
            adict['legend0']=[]
        import matplotlib.cm as cm

        def get_legends(types):
            legends = []
            for name in types:
                legends.append(name)
            return legends


	ykeys=y.keys()

	#ykeys.sort(key=lambda x:(x.split('_')[-1][:-1],x.split('_')[0]))
	ykeys.sort()
	colors=['blue','green','red','cyan','pink','gray']
	linetypes=['-','--','-.',':']
	#colors=cm.Dark2(np.linspace(0, 1, 6))
	i=0
        for c,cnt in zip(ykeys,range(len(y))):
	    j=i%len(colors)
	    k=i/len(colors)
            x_trans = [pt[0] for pt in trans_level_pt[c]]
            y_trans = [pt[1] for pt in trans_level_pt[c]]
            ind=int(ceil((cnt+1)/float(nline_in_subplot)))-1
            adict['legend%s' % str(ind)].append(c)
	    #adict['ax%s' % str(ind)].scatter(x_trans, y_trans,c=colors[j],marker='*',s=60)
            adict['ax%s' % str(ind)].plot(x, y[c], linewidth=2,color=colors[j],linestyle=linetypes[k])#, color=colors[cnt])
            #adict['ax%s' % str(ind)].scatter(x_trans, y_trans,c=colors[j],marker='*',markersize=8, linestyle='none',label=None)
 	    for ll in range(len(text_q_labels_pt[c])):
	    	adict['ax%s' % str(ind)].text(text_q_xy_pt[c][ll][0],text_q_xy_pt[c][ll][1],str(text_q_labels_pt[c][ll]),color=colors[j],horizontalalignment='center')
	    if i==nline_in_subplot-1:
		i=0
		continue
	    i=i+1
            #plt.plot(x, y[c], linewidth=3, color=colors[cnt])

        #i=0
        #for c,cnt in zip(ykeys,range(len(y))):
        #    j=i%len(colors)
        #    ind=int(ceil((cnt+1)/float(nline_in_subplot)))-1
        #    x_trans = [pt[0] for pt in trans_level_pt[c]]
        #    y_trans = [pt[1] for pt in trans_level_pt[c]]
        #    adict['ax%s' % str(ind)].plot(x_trans, y_trans,c=colors[j],marker='*',markersize=8, linestyle='none')
 	#    for ll in range(len(text_q_labels_pt[c])):
	#    	adict['ax%s' % str(ind)].text(text_q_xy_pt[c][ll][0],text_q_xy_pt[c][ll][1],str(text_q_labels_pt[c][ll]),color=colors[j],horizontalalignment='center')
        #    i=i+1
	#
        for i in range(nsubplot):
            adict['ax%s' % i].plot([min_lim, max_lim], [0, 0], 'k-')
            adict['ax%s' % i].axvline(x=0.0, linestyle='--', color='k', linewidth=2)
            adict['ax%s' % i].axvline(x=self._analyzer._band_gap, linestyle='--', color='k',
                        linewidth=2)
            if ylim is not None:
                adict['ax%s' % i].set_ylim(ylim)
                adict['ax%s' % i].set_xlim(xlim)
                adict['ax%s' % i].legend(get_legends(adict['legend%s' % i]),loc='upper right')
            #adict['ax%s' % i].set_xlabel("Fermi energy (eV)"s
            #adict['ax%s' % i].set_ylabel("Formation Energy (eV)")
        adict['ax%s' % 0].set_ylabel("Formation energy (eV)")
        adict['ax%s' % str(nsubplot-1)].set_xlabel('Fermi energy (eV)')
        fig.set_tight_layout({'pad':0.2, 'h_pad':0.1, 'w_pad':0.0})

        #for cnt, c in enumerate(y):
        ##for c in y:
        #   # plt.plot(x, y[c], next(linecycler), linewidth=6, color=colors[cnt])
        #    x_trans = [pt[0] for pt in trans_level_pt[c]]
        #    y_trans = [pt[1] for pt in trans_level_pt[c]]
        #    plt.plot(x_trans, y_trans,  marker='*', color=colors[cnt], markersize=12, fillstyle='full')
        return plt

    def plot_conc_temp(self, me=[1.0, 1.0, 1.0], mh=[1.0, 1.0, 1.0]):
        """
        plot the concentration of carriers vs temperature both in eq and non-eq after quenching at 300K
        Args:
            me:
                the effective mass for the electrons as a list of 3 eigenvalues
            mh:
                the effective mass for the holes as a list of 3 eigenvalues
        Returns;
            a matplotlib object

        """
        temps = [i*100 for i in range(3,20)]
        qi = []
        qi_non_eq = []
        for t in temps:
            qi.append(self._analyzer.get_eq_Ef(t, me, mh)['Qi']*1e-6)
            qi_non_eq.append(
                    self._analyzer.get_non_eq_Ef(t, 300, me, mh)['Qi']*1e-6)

        plt = get_publication_quality_plot(12, 8)
        plt.xlabel("temperature (K)")
        plt.ylabel("carrier concentration (cm$^{-3}$)")
        plt.semilogy(temps, qi, linewidth=3.0)
        plt.semilogy(temps, qi_non_eq, linewidth=3)
        plt.legend(['eq','non-eq'])
        return plt

    def plot_carriers_ef(self, temp=300, me=[1.0, 1.0, 1.0], mh=[1.0, 1.0, 1.0]):
        """
        plot carrier concentration in function of the fermi energy
        Args:
            temp:
                temperature
            me:
                the effective mass for the electrons as a list of 3 eigenvalues
            mh:
                the effective mass for the holes as a list of 3 eigenvalues
        Returns:
            a matplotlib object
        """
        plt = get_publication_quality_plot(12, 8)
        qi = []
        efs = []
        for ef in [x * 0.01 for x in range(0, 100)]:
            efs.append(ef)
            qi.append(self._analyzer.get_Qi(ef, temp, me, mh)*1e-6)
        plt.ylim([1e14, 1e22])
        return plt.semilogy(efs, qi)
