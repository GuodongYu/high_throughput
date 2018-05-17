from pymatgen.io.vasp.outputs import Outcar
import numpy as np

class ElasticAnalyzer(object):
    def __init__(self,OUTCAR_file):
        self.out = Outcar(OUTCAR_file)
        self.OUTCAR = OUTCAR_file

    def read_elastic_modulus(self):
        C = np.mat(np.zeros((6,6)))
        string = 'TOTAL ELASTIC MODULI'
        with open(self.OUTCAR) as f:
            out = f.read().split('\n')
        for line_num, line in enumerate(out):
            if string in line:
                start_line = line_num
                break
        block = [i.split()[1:] for i in out[start_line+3:start_line+9]]
        C = [[float(j) for j in i] for i in block]
        return np.mat(C)
    
    def get_elastic_compliances(self):
        C = self.read_elastic_modulus()
        S = np.linalg.inv(C)
        return S
        

            




