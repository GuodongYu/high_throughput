from pymatgen import SETTINGS, SETTINGS_FILE
from pymatgen import MPRester
try:
    m = MPRester(SETTINGS['PMG_MAPI_KEY'])
except:
    m = MPRester(SETTINGS['MAPI_KEY'])


