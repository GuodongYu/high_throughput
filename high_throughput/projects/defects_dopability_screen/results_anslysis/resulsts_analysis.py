
three_five_file = '/home/gyu/workplace/My_Work/defects_thermoelectrons_HT/3-5.dat'
TE_file = '/home/gyu/workplace/My_Work/defects_thermoelectrons_HT/thermoelectrics_list.dat'


e_nega = {'B':2.04,'Al':1.61,'Ga':1.81,'In':1.78,'N':3.04,'P':2.19,'As':2.18,'Sb':2.05,'Bi':2.02}


class analyzer(object):
    def __init__(self,system='3-5'):
        self.colors = ['red','blue','green','black','orange','cyan','yellow','magenta','purple','pink','brown','gray']*20
        self.system = system
        if system == '3-5':
            file_ = three_five_file
        elif system == 'TE':
            file_ = TE_file
        results = {}
        with open(file_) as f:
            file_cont = f.read().split('\n')
        for i in file_cont:
            line = i.split('\t')
            if line[0] == 'mpid' or line == ['']:
                continue
            result = {}
            mpid = line[0]
            result['formula']=line[1]
            result['space_group']=line[2]
            try:
                result['hse_gap']=float(line[3])
                result['pbe_gap']=float(line[4])
                result['vbm_shift']=float(line[5])
                result['cbm_shift']=float(line[7])
                result['BBD_ready']=line[8]
                if result['BBD_ready'] == 'TTT':
                    result['epsilon'] = float(line[9])
                    result['p_dopability'] = line[10]
                    result['plus_lost'] = line[11]
                    result['plus_potalign_bad'] = line[12]
                    result['n_dopability'] = line[13]
                    result['minus_lost'] = line[14]
                    result['minus_potalign_bad'] = line[15]
            except:
                pass
            results[mpid]=result
        self.results = results

    def get_elts_from_formula(self,formula):
        elts = []
        for i in range(len(formula)-1):
            if formula[i].isupper():
                if formula[i+1].islower():
                    elt = formula[i]+formula[i+1]
                else:
                    elt = formula[i]
                elts.append(elt)
        if formula[-1].isupper():
            elts.append(formula[-1])
        return elts

    def vbm_shift_vs_enegativity(self,elts=['N','P','As','Sb','Bi'],space_group='all'):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)

    def gap_vs_electronegativity_diff(self,space_group='all'):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        for i in self.results:
            result = self.results[i]
            spg = result['space_group']
            if space_group not in ['spg','all']:
                continue
            formula = result['formula']
            elts = self.get_elts_from_formula(formula)
            for elt in elts:
                if elt in ['N','P','As','Sb']:
                    elt_minus = elt
                elif elt in ['B','Al','Ga','In']:
                    elt_plus = elt
            e_nega_diff = e_nega[elt_minus] - e_nega[elt_plus]
            if 'hse_gap' in result:
                ax.scatter(e_nega_diff,result['hse_gap'],color='red')
                ax.scatter(e_nega_diff,result['pbe_gap'],color='blue')
            ax.set_xlabel('electronegativity difference')
            ax.set_ylabel('gap (eV)')
            fig.savefig('gap_vs_electronegativity.png')

    def band_edge_shift_plot(self,elts=['N','P','As','Sb'],space_group=['all'],xlim=None,ylim=None):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        colors = ['red','blue','green','black','orange','cyan','yellow','magenta','purple','pink','brown','gray']*20

        for i in self.results:
            result = self.results[i]
            if 'all' not in space_group:
                if result['space_group'] not in space_group:
                    continue

            formula = result['formula']
            if 'vbm_shift' not in result:
                continue
            vbm_shift = result['vbm_shift']
            cbm_shift = result['cbm_shift']
            for i in range(len(elts)):
                all_elts = self.get_elts_from_formula(formula)
                if elts[i] in all_elts and len(all_elts) == 2:
                    c = colors[i]
                    ax.scatter(vbm_shift,cbm_shift,s=20, color=c)
                    break
        if xlim is None:
            xlim = ax.get_xlim()
        if ylim is None:
            ylim = ax.get_ylim()

        for ii in range(len(elts)):
            ax.scatter(xlim[0]-1,ylim[0]-1,s=20,color=colors[ii],label=elts[ii])
        xs=[x*0.02 for x in range(200)]
        ax.plot(xs,[0.5*i for i in xs])
        ax.set_xlabel('vbm_shift (eV)')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel('cbm_shift (eV)')
        ax.legend()
        fig.savefig(self.system+'-band_edge_shift.png')

    def gap_comparison_plot(self):
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        for i in self.results:
            result = self.results[i]
            try:
                ax.scatter(result['hse_gap'],result['pbe_gap'],s=20,color='black')
            except:
                pass
        xs = [x*0.02 for x in range(200)]
        ax.plot(xs,[i*0.5 for i in xs])
        ax.set_xlabel('hse_gap (eV)')
        ax.set_ylabel('pbe_gap (eV)')
        fig.savefig(self.system+'_gap_comparison.png')

def main():
    system = 'TE'
    system = '3-5'
    #system = 'import'
    analy = analyzer(system)
    if system == '3-5':
        #analy.band_edge_shift_plot(space_group = 'all')
        analy.gap_vs_electronegativity_diff()
    elif system == 'TE':
        #analy.band_edge_shift_plot(['F','O','Cl','N','Br','I','S','Se','P','As','Te','Sb'])
        analy.gap_comparison_plot()


if __name__ == '__main__':
    main()
            
        

