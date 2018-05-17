"""
This module computes finite size supercell charge corrections for
defects in anistropic systems using extended Freysoldt (or Kumagai) method 
developed by Kumagai and Oba.
Kumagai method includes
   a) anisotropic PC energy
   b) potential alignment by atomic site averaging at Wigner Seitz cell
      edge
If you use the corrections implemented in this module, cite
 a) Kumagai and Oba, Phys. Rev. B. 89, 195205 (2014) and
 b) Freysoldt, Neugebauer, and Van de Walle,
    Phys. Status Solidi B. 248, 1067-1076 (2011)  and
in addition to the pycdt paper
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import math
import logging

import numpy as np

from pymatgen.io.vasp.outputs import Locpot, Outcar
from pymatgen.core.lattice import Lattice

from pycdt.corrections.utils import *
from pycdt.utils.units import hart_to_ev

norm = np.linalg.norm


def kumagai_init(structure, dieltens):
    angset = structure.lattice.get_cartesian_coords(1)

    dieltens = np.array(dieltens)
    if not len(dieltens.shape):
        dieltens = dieltens*np.identity(3)
    elif len(dieltens.shape) == 1:
        dieltens = np.diagflat(dieltens)

    logging.getLogger(__name__).debug('Lattice constants (in Angs): '
                                      + str(cleanlat(angset)))
    [a1, a2, a3] = ang_to_bohr * angset  # convert to bohr
    bohrset = [a1, a2, a3]
    vol = np.dot(a1, np.cross(a2, a3))

    logging.getLogger(__name__).debug('Lattice constants (in Bohr): '
                                      + str(cleanlat([a1, a2, a3])))
    determ = np.linalg.det(dieltens)
    invdiel = np.linalg.inv(dieltens)
    logging.getLogger(__name__).debug('inv dielectric tensor: ' + str(invdiel))
    return angset, bohrset, vol, determ, invdiel


def real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance):
    invdiel = np.linalg.inv(dieltens)
    determ = np.linalg.det(dieltens)
    realpre = q / np.sqrt(determ)
    tolerance /= hart_to_ev

    #Real space sum by converging with respect to real space vectors
    #create list of real space vectors that satisfy |i*a1+j*a2+k*a3|<=N
    Nmaxlength = 40  #tolerance for stopping real space sum convergence
    N = 2
    r_sums = []
    while N < Nmaxlength:  
        r_sum = 0.0
        if norm(r):
            for i in range(-N, N+1):
                for j in range(-N, N+1):
                    for k in range(-N, N+1):
                        r_vec = i*a1 + j*a2 + k*a3 - r
                        loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                        nmr = math.erfc(gamma * np.sqrt(loc_res))
                        dmr = np.sqrt(determ * loc_res)
                        r_sum += nmr / dmr  
        else:
            for i in range(-N, N+1):
                for j in range(-N, N+1):
                    for k in range(-N, N+1):
                        if i == j == k == 0:
                            continue
                        else:
                            r_vec = i*a1 + j*a2 + k*a3 
                            loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                            nmr = math.erfc(gamma * np.sqrt(loc_res))
                            dmr = np.sqrt(determ * loc_res)
                            r_sum += nmr / dmr  
        r_sums.append([N, realpre * r_sum])

        if N == Nmaxlength-1:
            logging.getLogger(__name__).warning(
                'Direct part could not converge with real space translation '
                'tolerance of {} for gamma {}'.format(Nmaxlength-1, gamma))
            return
        elif len(r_sums) > 3:
            if abs(abs(r_sums[-1][1]) - abs(r_sums[-2][1])) < tolerance:
                r_sum = r_sums[-1][1]
                logging.debug("gamma is {}".format(gamma))
                logging.getLogger(__name__).debug(
                    "convergence for real summatin term occurs at step {} "
                    "where real sum is {}".format(N,  r_sum * hart_to_ev))
                break

        N += 1
    return r_sum


def get_g_sum_at_r(g_sum, structure, dim, r):
    """
    Args:
        g_sum: Reciprocal summation calculated from reciprocal_sum method
        structure: Bulk structure pymatgen object
        dim : ngxf dimension
        r: Position relative to defect (in cartesian coords)
    Returns:
        reciprocal summ value at g_sum[i_rx,j_ry,k_rz]
    """

    fraccoord = structure.lattice.get_fractional_coords(r)
    i, j, k = getgridind(structure, dim, fraccoord)

    return g_sum[i, j, k]


def anisotropic_madelung_potential(structure, dim, g_sum, r, dieltens, q,  
                                   gamma, tolerance):
    """
    Compute the anisotropic Madelung potential at r not equal to 0.
    For r=(0,0,0) use anisotropic_pc_energy function
    Args:
        structure: Bulk pymatgen structure type
        dim : ngxf dimension
        g_sum: Precomputed reciprocal sum for all r_vectors
        r: r vector (in cartesian coordinates) relative to defect position. 
           Non zero r is expected
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        tolerance: Tolerance parameter for numerical convergence
        gamma (float): Convergence parameter 
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens)

    recippartreal = q * get_g_sum_at_r(g_sum, structure, dim, r)
    directpart = real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q * np.pi / (vol * (gamma ** 2))
    logging.getLogger(__name__).debug('self interaction piece is {}'.format(
            selfint * hart_to_ev))

    pot = hart_to_ev * (directpart + recippartreal - selfint)
    return pot


def anisotropic_pc_energy(structure, g_sum, dieltens, q, gamma, tolerance):
    """
    Compute the anistropic periodic point charge interaction energy.
    Args:
        structure: Bulk pymatgen structure type
        g_sum : comes from KumagaiBulkInit class
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        gamma : convergence parameter optimized in KumagaiBulkInit class
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens)

    g_part = q*g_sum[0,0,0]
    r_part = real_sum(a1, a2, a3, [0,0,0], q, dieltens, gamma, tolerance)
    selfint = q*np.pi / (vol * (gamma**2)) #self interaction term
    #surface term (only for r not at origin)
    surfterm = 2*gamma*q / np.sqrt(np.pi*determ)

    logger = logging.getLogger(__name__)
    logger.debug('reciprocal part: {}'.format(g_part * hart_to_ev))
    logger.debug('real part: {}'.format(r_part * hart_to_ev))
    logger.debug('self interaction part: {}'.format(selfint * hart_to_ev))
    logger.debug('surface term: {}'.format(surfterm * hart_to_ev))

    pc_energy = -(q*0.5*hart_to_ev) * (r_part + g_part - selfint - surfterm)
    logging.debug('Final PC Energy term: {} eV'.format(pc_energy))

    return pc_energy


def getgridind(structure, dim, r, gridavg=0.0):
    """
    Computes the index of a point, r, in the locpot grid
    Args:
        structure:
            Pymatgen structure object
        dim:
            dimension of FFT grid (NGXF dimension list in VASP)
        r:
            Relative co-ordinates with respect to abc lattice vectors
        gridavg:
            If you want to do atomic site averaging, set gridavg to
            the radius of the atom at r
    Returns:
        [i,j,k]: Indices as list
    TODO: Once final, remove the getgridind inside disttrans function
    """
    abc = structure.lattice.abc
    grdind = []

    if gridavg:
        radvals = [] #radius in terms of indices
        dxvals = []

    for i in range(3):
        if r[i] < 0:
            while r[i] < 0:
                r[i] += 1
        elif r[i] >= 1:
            while r[i] >= 1:
                r[i] -= 1
        r[i] *= abc[i]
        num_pts = dim[i]
        x = [now_num / float(num_pts) * abc[i] for now_num in range(num_pts)]
        dx = x[1] - x[0]
        x_rprojection_delta_abs = np.absolute(x - r[i])
        ind = np.argmin(x_rprojection_delta_abs)
        if x_rprojection_delta_abs[ind] > dx*1.1: #to avoid numerical errors
            logger = logging.getLogger(__name__)
            logger.error("Input position not within the locpot grid")
            logger.error("%d, %d, %f", i, ind, r)
            logger.error("%f", x_rprojection_delta_abs)
            raise ValueError("Input position is not within the locpot grid")
        grdind.append(ind)
        if gridavg:
            radvals.append(int(np.ceil(gridavg/dx)))
            dxvals.append(dx)

    if gridavg:
        grdindfull = []
        for i in range(-radvals[0], radvals[0]+1):
            for j in range(-radvals[1], radvals[1]+1):
                for k in range(-radvals[2], radvals[2]+1):
                    dtoc = [i*dxvals[0], j*dxvals[1], k*dxvals[2]]
                    if norm(dtoc) < gridavg:
                        ival = (i+grdind[0]) % dim[0]
                        jval = (j+grdind[1]) % dim[1]
                        kval = (k+grdind[2]) % dim[2]
                        grdindfull.append((ival, jval, kval))
        grdind = grdindfull

    return grdind


def disttrans(struct, defstruct, dim):
    """
    To calculate distance from defect to each atom and finding NGX grid
    pts at each atom.
    Args:
        struct: Bulk structure object
        defstruct: Defect structure object
        dim: dimensions of FFT grid
    """

    #Find defect location in bulk and defect cells
    blksite, defsite = find_defect_pos(struct, defstruct)
    logger = logging.getLogger(__name__)
    if blksite is None and defsite is None:
        logger.error('Not able to determine defect site')
        return
    if blksite is None:
        logger.debug('Found defect to be Interstitial type at %s',
                      repr(defsite))
    elif defsite is None:
        logger.debug('Found defect to be Vacancy type at %s', repr(blksite))
    else:
        logger.debug('Found defect to be antisite/subsitution type at %s ' \
                      ' in bulk, and %s in defect cell', 
                      repr(blksite), repr(defsite))

    if blksite is None:
        blksite = defsite
    elif defsite is None:
        defsite = blksite

    def_ccoord = blksite[:]
    defcell_def_ccoord = defsite[:]

    sitelist = struct.sites[:]

    #better image getter since pymatgen wasnt working well for this
    def returnclosestr(vec):
        from operator import itemgetter
        listvals = []
        abclats = defstruct.lattice.matrix
        trylist = [-1, 0, 1]
        for i in trylist:
            for j in trylist:
                for k in trylist:
                    transvec = i*abclats[0] + j*abclats[1] + k*abclats[2]
                    #rnew = vec - (defcell_def_ccoord + transvec)
                    rnew = vec + transvec - defcell_def_ccoord
                    listvals.append([norm(rnew), rnew, transvec])
        listvals.sort(key=itemgetter(0))
        return listvals[0] #will return [dist,r to defect, and transvec for defect]

    grid_sites = {}  # dictionary with indices keys in order of structure list
    for i in sitelist:
        if np.array_equal(i.coords, def_ccoord):
            logging.debug('Site {} is defect! Skipping '.format(i))
            continue

        blksite, defsite = closestsites(struct, defstruct, i.coords)

        # blkindex=struct.index(blksite[0])
        # defindex=defstruct.index(defsite[0])
        blkindex = blksite[-1]
        defindex = defsite[-1]

        dcart_coord = defsite[0].coords ## orig
        #dcart_coord = blksite[0].coords ## yugd
        closeimage = returnclosestr(dcart_coord)
        cart_reldef = closeimage[1]
        defdist = closeimage[0]

        if abs(norm(cart_reldef) - defdist) > 0.1:
            logger.warning('Image locater issue encountered for site = %d',
                            blkindex)
            logger.warning('In defect supercell')
            logger.warning('Distance should be %f', defdist)
            logger.warning('But, calculated distance is %f', norm(cart_reldef))
            #dont want to break the code here, but want flag to exist...what to do?
            #return 

        if blkindex in grid_sites:
            logger.warning('Index %d already exists in potinddict!', blkindex)
            logger.warning('Overwriting information.')

        grid_sites[blkindex] = {
                'dist': defdist, ### distance of site from defect
                'cart': dcart_coord, ### absulute vector of site 
                'cart_reldef': cart_reldef,  ### relative vector of site to defect
                'siteobj': [i.coords, i.frac_coords, i.species_string], 
                'bulk_site_index': blkindex, #### index of site in bulk structure
                'def_site_index': defindex}  ### index of site in defect structure

    return grid_sites

def get_rcut(site_pots,shell_thickness=0.1,npts_in_sample=5):
    sites = [i['dist'] for i in site_pots]
    sites = list(set([round(i,3) for i in sites]))
    sites_remove = []
    scan =[]
    for i in sites:
        scan.append(i)
        noscan=[k for k in sites if k not in scan]    
        for j in noscan:
            if abs(i-j) <= shell_thickness:
                sites_remove.append(j)
                scan.append(j)
    sites_short= [i for i in sites if i not in sites_remove]
    sites_short.sort(reverse=True)
    return (sites_short[4]+sites_short[5])/2.0
    

def check_converge(site_pots, rcut):
    pots=[(i['dist'],i['potdiff']) for i in site_pots]
    rmax = max([i[0] for i in pots])
    steps = int(round((rmax - rcut)/0.1))
    
    
            

def wigner_seitz_radius(structure):
    """
    Calculate the Wigner Seitz radius for the given structure.
    Args:
        structure: pymatgen Structure object
    """
    wz = structure.lattice.get_wigner_seitz_cell()

    dist = []
    for facet in wz:
        midpt = np.mean(np.array(facet), axis=0)
        dist.append(norm(midpt))
    wsrad = min(dist)
    wsrad_max = max(dist)

    return (wsrad + wsrad_max)/2.0


def read_ES_avg(location_outcar):
    """
    Reads NGXF information and Electrostatic potential at each atomic 
    site from VASP OUTCAR file.
    """
    with open(location_outcar,'r') as file:
        tmp_out_dat = file.read()

        out_dat = tmp_out_dat.split('\n')
        start_line = 0
        end_line = 0
        ngxf_line = 0
        for line_num, line in enumerate(out_dat):
            if "dimension x,y,z NGXF" in line:
                ngxf_line = line_num
            elif "average (electrostatic) potential at core" in line:
                start_line = line_num
                end_line = 0 # For multiple electrostatic read outs
            elif start_line and not end_line and not len(line.split()):
                end_line = line_num

        ngxlineout = out_dat[ngxf_line].split()
        ngxf_dims = list(map(int, ngxlineout[3:8:2]))

        rad_line = out_dat[start_line+1].split()
        radii = list(map(float, rad_line[5:]))
        #Above would be better to do as dictionary but no structure available

        ES_data = {'sampling_radii': radii, 'ngxf_dims': ngxf_dims}
        pot = []
        for line_num in range(start_line+3, end_line):
            line = out_dat[line_num].split()
            #TO BHARAT: Outcar can sometimes produces lines that get read in here like:
            #['16-105.4492', '17-105.4453', '18-105.4643', '19-105.4629', '20-105.4482']
            #so I added this stupid hack to get around when the numbers were together in Outcar...
            badflag=0
            for entry in line:
                if ('-' in entry) and (entry[0]!='-'):
                    badflag=1
            if badflag:
                line = []
                for inval in out_dat[line_num].split():
                    if (len(inval.split('-'))!=1) and inval[0]!='-':
                        # print 'rewriting',inval
                        [x,y] = inval.split('-')
                        line.append(x)
                        line.append('-'+y)
                    else:
                        line.append(inval)
                # print 'made ',line
            ############
            avg_es = map(float, line[1::2])
            pot += avg_es
        ES_data.update({'potential': pot})

        return ES_data

    return None


def read_ES_avg_fromlocpot(locpot):
    """
    TODO: smarter radii based on ENAUG to reproduce the core

    Reads Electrostatic potential at each atomic
    site from Locpot Pymatgen object
    """
    structure = locpot.structure
    radii = {specie: 1.0 for specie in set(structure.species)}
    # The above needs to be smarter (related to ENAUG?)

    ES_data = {'sampling_radii': radii, 'ngxf_dims': locpot.dim}
    pot = []
    for site in structure.sites:
        indexlist = getgridind(structure, locpot.dim,  site.frac_coords,
                               gridavg=radii[site.specie])
        samplevals = []
        for u,v,w in indexlist:
            samplevals.append(locpot.data["total"][u][v][w])
        pot.append(np.mean(samplevals))

    ES_data.update({'potential': pot})

    return ES_data


class KumagaiBulkInit(object):
    """
    Compute the anisotropic madelung potential array from the bulk 
    locpot. This helps in evaluating the bulk supercell related part 
    once to speed up the calculations.
    """
    def __init__(self, structure, dim, epsilon, encut=520, tolerance=0.0001,
                 optgamma=False):
        """
        Args
            structure: 
                Pymatgen structure object of bulk cell
            dim: 
                Fine FFT grid dimensions as a list 
                For vasp this is NGXF grid dimensions
            epsilon: 
                Dielectric tensor
            encut (float): 
                Energy cutoff for optimal gamma
            tolerance (float): 
                Accuracy parameter
            optgamma: 
                if you know optimized gamma, give its value. 
                Otherwise it will be computed.
        """
        self.structure = structure
        self.dim = dim
        self.epsilon = epsilon
        self.encut = encut
        self.tolerance = tolerance
        #self.silence = silence
        if not optgamma:
            self.gamma = self.find_optimal_gamma()
        else:
            self.gamma = optgamma
        print 'optimized gamma %s' % str(self.gamma)
        if self.gamma is None:
            self.gamma = 20.
            print "Warnning: set gamma to bt 20.0 by hand, please check for sure!"
        self.g_sum = self.reciprocal_sum()
        logging.getLogger(__name__).info('optimized gamma: %f', self.gamma)

    def find_optimal_gamma(self):
        """
        Find optimal gamma by evaluating the brute force reciprocal
        summation and seeing when the values are on the order of 1,
        This calculation is the anisotropic Madelung potential at r = (0,0,0).
        Note this only requires the STRUCTURE not the LOCPOT object.
        """
        print 'Optimizing gamma until reci_summation > 1 eV ......'
        print '%-14s %-17s' % ('gamma', 'reci_summation:')
        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                self.structure, self.epsilon)
        optgam = None

        #do brute force recip summation
        def get_recippart(encut, gamma):
            recippart = 0.0
            for rec in genrecip(a1, a2, a3, encut):
                Gdotdiel = np.dot(rec, np.dot(self.epsilon, rec))
                summand = math.exp(-Gdotdiel / (4 * (gamma ** 2))) / Gdotdiel
                recippart += summand
            recippart *= 4*np.pi/vol
            return recippart, 0.0

        def do_summation(gamma):
            # Do recip sum until it is bigger than 1eV
            # First do Recip space sum convergence with respect to encut for 
            # this gamma
            encut = 20  #start with small encut for expediency
            recippartreal1, recippartimag1 = get_recippart(encut, gamma)
            encut += 10
            recippartreal, recippartimag = get_recippart(encut, gamma)
            converge = [recippartreal1, recippartreal]
            logger = logging.getLogger(__name__)
            encut_step=10.
            while abs(abs(converge[0]) - abs(converge[1])) * hart_to_ev > \
                    self.tolerance:
                if encut >=self.encut-10:
                    encut_step = 4.
                encut += encut_step
                recippartreal, recippartimag = get_recippart(encut, gamma)
                converge.reverse()
                converge[1] = recippartreal
                if encut > self.encut:
                    msg = 'Optimal gamma not found at {} eV cutoff'.format(
                                self.encut)
                    logger.error(msg)
                    raise ValueError(msg)
            #print 'gamma',gamma, 'sum', abs(recippartreal)*hart_to_ev
            if abs(recippartimag) * hart_to_ev > self.tolerance:
                logger.error("Imaginary part of reciprocal sum not converged.")
                logger.error("Imaginary sum value is {} (eV)".format(
                    recippartimag * hart_to_ev))
                return None, None, None
            logger.debug('Reciprocal sum converged to %f eV',
                         recippartreal * hart_to_ev)
            logger.debug('Convergin encut = %d eV', encut)

            if (abs(converge[1]) * hart_to_ev < 1 and not optgam):
                #logger.warning('Reciprocal summation value is less than 1 eV.')
                #logger.warning('Might lead to errors')
                #logger.warning('Change gamma.')
                #return None, 'Try Again'
                return recippartreal*hart_to_ev, gamma ,'Try Again'
            
            return recippartreal*hart_to_ev, gamma, 'OK'

        logger = logging.getLogger(__name__)
        #start with gamma s.t. gamma*L=5 (some paper said this is optimal)
        #optimizing gamma for the reciprocal sum to improve convergence 
        gamma = 5.0/(vol ** (1/3.0))
        optimal_gamma_found = False

        opti_hist=[]
        predict = False
        step = 10.
        while not optimal_gamma_found:
            recippartreal, optgamma, state = do_summation(gamma)
            print '%-14s %-17s' % (str(gamma), str(recippartreal))
            if state == 'OK':
                logger.debug('optimized gamma found to be %f', optgamma)
                optimal_gamma_found = True
            elif 'Try Again' in state:
                if gamma < 10.:
                    gamma *= 1.5
                else:
                    opti_hist.append((gamma,recippartreal))
                    gamma += step 
                    if len(opti_hist)>1:
                        x1,y1=opti_hist[0]
                        x2,y2=opti_hist[1]
                        a = (y2-y1)/(x2-x1)
                        b = y1 - a*x1
                        gamma_should= (1.-b)/a
                        print '(predicted gamma: %s)' % str(gamma_should + 0.5)
                        step = 2.
                        predict = True
                        gamma = gamma_should +0.5
                        if gamma > 60:
                            return None
            else:
                logger.error('Had problem in gamma optimization process.')
                return None

            #if gamma > 80:
            #    logger.error('Could not optimize gamma before gamma = %d', 80)
            #    return None

        return optgamma 

    def reciprocal_sum(self):
        """
        Compute the reciprocal summation in the anisotropic Madelung 
        potential.

        TODO: Get the input to fft cut by half by using rfft instead of fft
        """
        logger = logging.getLogger(__name__)
        logger.debug('Reciprocal summation in Madeling potential')
        over_atob = 1.0 / ang_to_bohr
        atob3 = ang_to_bohr ** 3

        latt = self.structure.lattice
        vol = latt.volume * atob3 # in Bohr^3

        reci_latt = latt.reciprocal_lattice
        [b1, b2, b3] = reci_latt.get_cartesian_coords(1)
        b1 = np.array(b1) * over_atob # In 1/Bohr
        b2 = np.array(b2) * over_atob
        b3 = np.array(b3) * over_atob

        nx, ny, nz = self.dim
        logging.debug('nx: %d, ny: %d, nz: %d', nx, ny, nz)
        ind1 = np.arange(nx)
        for i in range(int(nx/2), nx):
            ind1[i] = i - nx
        ind2 = np.arange(ny)
        for i in range(int(ny/2), ny):
            ind2[i] = i - ny
        ind3 = np.arange(nz)
        for i in range(int(nz/2), nz):
            ind3[i] = i - nz

        g_array = np.zeros(self.dim, np.dtype('c16'))
        gamm2 = 4*(self.gamma**2)
        for i in ind1:
            for j in ind2:
                for k in ind3:
                    g = i*b1 + j*b2 + k*b3
                    g_eps_g = np.dot(g, np.dot(self.epsilon, g))
                    if i == j == k == 0:
                        continue
                    else:
                        g_array[i,j,k] = math.exp(-g_eps_g/gamm2) / g_eps_g

        r_array = np.fft.fftn(g_array)
        over_vol = 4*np.pi/vol # Multiply with q later
        r_array *= over_vol
        r_arr_real = np.real(r_array)
        r_arr_imag = np.imag(r_array)

        max_imag = r_arr_imag.max()
        logger.debug('Max imaginary part found to be %f', max_imag)

        return r_arr_real


class KumagaiCorrection(object):
    """
    Extended freysoldt correction developed by Kumagai and Oba.
    """
    def __init__(self, dim_dfc, dielectric_tensor, q, gamma, g_sum, bulk_structure,
                 defect_structure, puredat, defdat, energy_cutoff=520, madetol=0.0001, 
                 lengths=None, **kw):
        """
        Args:
            dim_dfc: dim for defect calculation
            dielectric_tensor: 
                Macroscopic dielectric tensor
                Include ionic also if defect is relaxed, othewise ion clamped.
                Can be a matrix array or scalar.
            q: 
                Charge associated with the defect. Typically integer
            gamma:  
                Convergence parameter. Obtained from KumagaiBulkPart
            g_sum: 
                value that is dependent on the Bulk only. 
                Obtained from KumagaiBulkPart
            bulk_structure: 
                bulk Pymatgen structure object. Need to specify this if 
                using Outcar method for atomic site avg.
                (If you specify outcar files for bulk_file_path but dont 
                specify structure then code will break)
                (TO DO: resolve this dumb dependency by being smarter 
                about where structure comes from?)
            defect_structure: 
                defect structure. Needed if using Outcar method
            puredat: output of read_ES_avg from bulk OUTCAR
            defdat: output of read_ES_avg from defect OUTCAR
            energy_cutoff: 
                Energy for plane wave cutoff (in eV).
                If not given, Materials Project default 520 eV is used.
            madetol: 
                Tolerance for convergence of energy terms in eV 
            lengths:
                Lengths of axes, for speeding up plotting slightly
            keywords:
                1) bulk_locpot: Bulk Locpot file path OR Bulk Locpot 
                   defect_locpot: Defect Locpot file path or defect Locpot 
                2) (Or) bulk_outcar:   Bulk Outcar file path 
                   defect_outcar: Defect outcar file path 
        """
        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self.dieltens = np.identity(3) * dielectric_tensor
        else:
            self.dieltens = np.array(dielectric_tensor)

        self.puredat = puredat
        self.defdat = defdat
        self.dim = dim_dfc
        self.madetol = madetol
        self.q = q
        self.encut = energy_cutoff
        self.structure = bulk_structure
        self.defstructure = defect_structure
        self.gamma = gamma
        self.g_sum = g_sum

        self.lengths=lengths

    def correction(self, title=None, partflag='All'):
        """
        Computes the extended Freysoldt correction for anistropic systems 
        developed by Y. Kumagai and F. Oba (Ref: PRB 89, 195205 (2014)
        Args:
            title:
                If plot of potential averaging process is wanted set title 
            partflag:
                Specifies the part of correction computed
                'pc': periodic interaction of defect charges (point charge) only
                'potalign': potential alignmnet correction only, 
                'All' (default): pc and potalign combined into one value, 
                'AllSplit' for correction in form [PC, potterm, full]
        """
        logger = logging.getLogger(__name__)
        logger.info('This is Kumagai Correction.')

        if not self.q:
            if partflag == 'AllSplit':
                return [0., 0., 0., 0.]
            else:
                return 0.0

        if partflag != 'potalign':
            energy_pc = self.pc()

        if partflag != 'pc':
            potalign = self.potalign(title)

        #logger.info('Kumagai Correction details:')
        #if partflag != 'potalign':
        #    logger.info('PCenergy (E_lat) = %f', round(energy_pc, 5))
        #if partflag != 'pc':
        #    logger.info('potential alignment (-q*delta V) = %f',
        #                 round(potalign, 5))
        if partflag in ['All','AllSplit']:
            logger.info('Total Kumagai correction = %f',
                         round(energy_pc+potalign[0], 5))

        if partflag == 'pc':
            return round(energy_pc, 5)
        elif partflag == 'potalign':
            return round(potalign[0], 5)
        elif partflag == 'All':
            return round(energy_pc+potalign[0], 5)
        else:
            corr = map(lambda x: round(x, 5),
                        [energy_pc, potalign[0], energy_pc+potalign[0], potalign[1]])
            return corr 
    def pc(self):

        energy_pc = anisotropic_pc_energy(
                self.structure, self.g_sum, self.dieltens, self.q,
                self.gamma, self.madetol)

        logger = logging.getLogger(__name__)
        logger.info('PC energy determined to be %f eV (%f Hartree)',
                     energy_pc, energy_pc/hart_to_ev)

        return energy_pc

    def potalign(self, title=None):
        """
        Potential alignment for Kumagai method
        Args:
            title: Title for the plot. None will not generate the plot
        """
        logger = logging.getLogger(__name__)
        logger.info('\nRunning potential alignment (atomic site averaging)')

        
        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                self.structure, self.dieltens)

        potinddict = disttrans(self.structure, self.defstructure, self.dim)

        minlat = min(norm(a1), norm(a2), norm(a3))
        lat_perc_diffs = [100 * abs(norm(a1) - norm(lat)) / minlat for lat \
                          in [a2, a3]]
        lat_perc_diffs.append(100 * abs(norm(a2) - norm(a3)) / minlat)
        if not all(i < 45 for i in lat_perc_diffs):
            logger.warning('Detected that cell was not very cubic.')
            logger.warning('Sampling atoms outside wigner-seitz cell may '\
                            'not be optimal')
        wsrad = wigner_seitz_radius(self.structure)
        logger.debug('wsrad %f', wsrad)

        # get ES potential from either Outcar (VASP) or Locpot pymatgen object

        site_pots =[]
        for i in potinddict.keys():
            j = potinddict[i]['def_site_index'] #assuming zero defined
            k = potinddict[i]['bulk_site_index']
            v_qb = self.defdat['potential'][j] - self.puredat['potential'][k]


            cart_reldef = potinddict[i]['cart_reldef']
            v_pc = anisotropic_madelung_potential(
                    self.structure, self.dim, self.g_sum, cart_reldef, 
                    self.dieltens, self.q, self.gamma, self.madetol)
            v_qb *= -1 #change charge sign convention

            potinddict[i]['Vpc'] = v_pc
            potinddict[i]['Vqb'] = v_qb
            
            site_pots.append({'elt':potinddict[i]['siteobj'][2],'dist':potinddict[i]['dist'],'potdiff':v_qb-v_pc})

            logger.debug('Atom: %d, anisotropic madelung potential: %f',
                          i, v_pc)
            logger.debug('Atom: %d, bulk/defect difference = %f', i, v_qb)
        
        wsrad = get_rcut(site_pots)
        for i in potinddict.keys():
            logger.debug("Atom %d, distance: %f", i, potinddict[i]['dist'])
            if potinddict[i]['dist'] > wsrad:
                potinddict[i]['OutsideWS'] = True
            else:
                potinddict[i]['OutsideWS'] = False


        if title:
            fullspecset = self.structure.species
            specset = list(set(fullspecset))
            shade, forplot = {}, {}
            for i in specset:
                shade[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': []}
                forplot[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': [],'sites':[]}

        forcorrection = []
        for i in potinddict.keys():
            if (not title and not potinddict[i]['OutsideWS']):
                continue
            if potinddict[i]['OutsideWS']:
                forcorrection.append(potinddict[i]['Vqb']-potinddict[i]['Vpc'])
                if title:
                    elt = fullspecset[i].symbol
                    shade[elt]['r'].append(potinddict[i]['dist'])
                    shade[elt]['Vpc'].append(potinddict[i]['Vpc'])
                    shade[elt]['Vqb'].append(potinddict[i]['Vqb'])
            if title:
                elt = fullspecset[i].symbol
                forplot[elt]['r'].append(potinddict[i]['dist'])
                forplot[elt]['Vpc'].append(potinddict[i]['Vpc'])
                forplot[elt]['Vqb'].append(potinddict[i]['Vqb'])
                forplot[elt]['sites'].append(potinddict[i]['siteobj'])

        potalign = np.mean(forcorrection)
        std_devia = np.std(forcorrection)
        if title:
            forplot['EXTRA'] = {'wsrad': wsrad, 'potalign': potalign}
            try:
                forplot['EXTRA']['lengths']=self.structure.lattice.abc
            except:
                forplot['EXTRA']['lengths']=self.lengths

            if title != 'written':
                KumagaiCorrection.plot(forplot, title, std_devia)
            else:
                #TODO: use a more descriptive fname that describes the defect
                from monty.serialization import dumpfn
                from monty.json import MontyEncoder
                fname = 'KumagaiData.json'
                dumpfn(forplot, fname, cls=MontyEncoder)

        logger.info('potential alignment (site averaging): %f',
                     np.mean(forcorrection))
        logger.info('Potential correction energy: %f eV',
                     -self.q * np.mean(forcorrection))
        return (-self.q*potalign, std_devia)

    @classmethod
    def plot(cls, forplot, title, std_devia):
        """
        Plotting of locpot data
        TODO: Rename forplot to a more descriptive name
        """
        import matplotlib.pyplot as plt
        text='std_devia:%s' % str(round(std_devia,3))
        plt.figure()
        plt.clf()
        collis = ['b', 'g', 'c', 'm', 'y', 'w', 'k']
        ylis = []
        rlis = []
        for i in range(len(forplot.keys())):
            inkey = list(forplot.keys())[i]
            if inkey == 'EXTRA':
                continue
            for k in forplot[inkey]['r']:
                rlis.append(k)
            for k in ['Vqb', 'Vpc']:
                for u in forplot[inkey][k]:
                    ylis.append(u)
            plt.plot(forplot[inkey]['r'], forplot[inkey]['Vqb'], 
                     color=collis[i], marker='^', alpha=0.5, linestyle='None',
                     label=str(inkey) + ': $V_{q/b}$')
            plt.plot(forplot[inkey]['r'], forplot[inkey]['Vpc'], 
                     color=collis[i], marker='o', alpha=0.5, linestyle='None',
                     label=str(inkey) + ': $V_{pc}$')
        full = []
        for i in forplot.keys():
            if i == 'EXTRA':
                continue
            for k in range(len(forplot[i]['Vpc'])):
                full.append([
                    forplot[i]['r'][k], 
                    forplot[i]['Vqb'][k] - forplot[i]['Vpc'][k]
                    ])
        realfull = sorted(full, key=lambda x: x[0])
        r, y = [], []
        for i in realfull:
            r.append(i[0])
            y.append(i[1])
        wsrad = forplot['EXTRA']['wsrad']
        potalign = forplot['EXTRA']['potalign']
        plt.plot(r, y, color=collis[-1], marker='x', linestyle='None',
                 label='$V_{q/b}$ - $V_{pc}$')
        plt.xlabel('Distance from defect ($\AA$)',fontsize=20)
        plt.ylabel('Potential (V)',fontsize=20)

        x = np.arange(wsrad, max(forplot['EXTRA']['lengths']), 0.01)
        plt.fill_between(x, min(ylis) - 1, max(ylis) + 1, facecolor='red', 
                         alpha=0.15, label='sample')
        plt.axhline(y=potalign, linewidth=0.5, color='red',
                    label='potalign')
        plt.legend(loc='lower left')
        ymin = min(ylis) - 0.5
        ymax = max(ylis) + 0.5
        ymin = -1; ymax =1
        y0 = 0.125*ymin+0.875*ymax
        x0 = (max(rlis) + 3)/2.0
        plt.text(x0,y0,text,fontsize=20)
        plt.axhline(y=0, linewidth=0.2, color='black')
        plt.ylim([ymin,ymax])
        plt.xlim([0, max(rlis) + 3])

        plt.title('%s atomic site potential plot' % title)
        plt.savefig('%s_kumagaisiteavgPlot.pdf' % title)
        plt.clf()

    @classmethod
    def plot_from_datfile(cls, name='KumagaiData.json', title='default'):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()

        """
        from monty.serialization import loadfn
        from monty.json import MontyDecoder

        forplot = loadfn(name, cls=MontyDecoder)
        cls.plot(forplot, title=title)
