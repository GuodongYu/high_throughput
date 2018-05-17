import numpy as np


def fit_linear_seebeck_vs_T(exp_data_file):
    import matplotlib.pyplot as plt
    with open(exp_data_file) as f:
        contents = f.read().split('\n')
        temps_exp = [float(i.split(',')[0]) for i in contents if i]
        quantity_exp = [float(i.split(',')[1]) for i in contents if i]
    plt.scatter(temps_exp,quantity_exp)
    temps=[i for i in temps_exp if i<=700]
    z1 = np.polyfit(temps, quantity_exp[0:len(temps)], 1)
    p1 = np.poly1d(z1)
    plt.plot(temps_exp, [p1(i) for i in temps_exp])
    plt.savefig('slope.png')
    print z1, p1

