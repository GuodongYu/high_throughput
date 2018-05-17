import os
import json
import csv

def plot_eff_mass_resis(compound,scale):
    root = os.path.join('.',compound)
    m = json.load(open(os.path.join(root,'eff_mass.json')))
    effm = m['eff_mass']['average'][0]
    Ts_m = m['Temps']

    Ts_r = []
    res =[]
    for  i in csv.reader(open(os.path.join(root,'resis_exp.csv'))):
        Ts_r.append(float(i[0]))
        res.append(float(i[1]))

    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(1,1)
    ax.plot(Ts_m,[i*scale for i in effm],label='eff_mass * %f' % scale)
    ax.scatter(Ts_r,res,label='experiment')
    ax.legend()
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('resistivity')
    fig.savefig('eff_mass_vs_resis.svg')

def main():
    plot_eff_mass_resis('Yb14MnSb11',2.4)


if __name__ == '__main__':
    main()

