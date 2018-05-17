from high_throughput.defects.database import DefectOperater, TasksOperater, \
                                        HseGapOperater , closest_defect_distance
from high_throughput.config import *
from high_throughput.defects.utils.util import *
from high_throughput.defects.defectsmaker import DefectRemaker
from high_throughput.defects.utils.vasp import make_vasp_defect_files
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pycdt.corrections.kumagai_correction import find_defect_pos,closestsites
import numpy as np
import json, os
from numpy import mat

def get_mpids(mpidfile):
    if os.path.split(mpidfile)[-1].split('.')[-1] == 'json':
        data = json.load(open(mpidfile))
        return data['mpids']
    with open(mpidfile,'r') as f:
        data = f.read()
    mpids = [i for i in data.split('\n') if i]
    return mpids


class BatchTools(object):
    def __init__(self):
        self.do = DefectOperater()

    def _check_gap_vbm(self, mpid):
        try:
	    pbe_vbm, pbe_gap = self.do._get_vbm_gap(mpid, 'dft')
            pbe_bs = {'vbm':pbe_vbm,'gap':pbe_gap}
            gap_dft = True
        except:
            print "%s dft gap: not ready!!" % mpid
            gap_dft = False
            pbe_bs = {'vbm':None,'gap':None}
        try:
	    hse_vbm, hse_gap = self.do._get_vbm_gap(mpid, 'hse')
            hse_bs = {'vbm':hse_vbm,'gap':hse_gap}
            gap_hse = True
        except:
            print "%s hse gap: not ready!!" % mpid
            gap_hse = False
            hse_bs = {'vbm':None,'gap':None}
        return {'dft':pbe_bs,'hse':hse_bs}

    def _check_rawdata(self, mpid):
        keys=['bulk', 'dielectric', 'defect']
        rawdata={i:True for i in keys}
        if mpid not in self.do.to.groups['bulk']:
            print "%s bulk: not ready!!" % mpid
            rawdata['bulk'] = False
        if mpid not in self.do.to.groups['dielectric']:
            print "%s dielectric: not ready!!" % mpid
            rawdata['dielectric'] = False
        if mpid not in self.do.to.groups['defect']:
            print "%s defects: not ready!!" % mpid
            rawdata['defect'] = False
        return rawdata    

    def refresh_database(self):
        self.do.to._classify()
        self.do._gather_all_mpids()
        self.do.hgo._gather_all_mpids()

    def plot_potalign(self,mpid,path='./'):
        self.do.to.calculate_defect_correction(mpid,path_save_potalign=path,db_update=False)

    def get_formation_energies_without_corr(self,mpid,vbmgap_from='hse'):
        rec = self.do.collection.find_one({'mpid':mpid})
        Ecorr={}
        forms_nocorr={}
        for gp in rec['defects']:
            for dfc in gp['computed_defects']:
                Ecorr[dfc['full_name']] = dfc['charge_correction']
        results = rec['results'][vbmgap_from]['regions']
        gap = rec['results'][vbmgap_from]['gap']
        for region in results:
            forms_nocorr[region]={}
            for dfc in results[region]['formation_energies']:
                name = dfc['name']
                q = dfc['charge']
                full_name=name+'_'+str(q)
                energy=dfc['energy']
                energy_nocorr = energy-Ecorr[full_name]
                Evbm = energy_nocorr
                Ecbm = q*gap+Evbm
                forms_nocorr[region][full_name]={'VBM':Evbm,'CBM':Ecbm}
        return forms_nocorr

    def get_epsilon(self,mpid):
        idd = self.do.to.groups['dielectric'][mpid][0]
        rec = self.do.to.collection.find_one({'_id':idd},['calculations.output.epsilon_ionic','calculations.output.epsilon_static'])
        epsilon_ionic= rec['calculations'][0]['output']['epsilon_ionic']
        eps_ion_avg = (epsilon_ionic[0][0]+epsilon_ionic[1][1]+epsilon_ionic[2][2])/3.0
        epsilon_static = rec['calculations'][0]['output']['epsilon_static']
        eps_ele_avg = (epsilon_static[0][0]+epsilon_static[1][1]+epsilon_static[2][2])/3.0
        return eps_ion_avg+eps_ele_avg, mat(epsilon_ionic)+mat(epsilon_static)

    def vbm_cbm_shift(self,mpid,**kws):
        try:
            pbe_bs = kws['pbe_bs']
            hse_bs = kws['hse_bs']
        except:
            check_gap = self._check_gap_vbm(mpid)
            pbe_bs = check_gap['pbe_bs']
            hse_bs = check_gap['hse_bs']
        if pbe_bs['gap'] is None or hse_bs['gap'] is None:
            return '', ''
        pbe_vbm = pbe_bs['vbm'];pbe_cbm = pbe_vbm + pbe_bs['gap']
        hse_vbm = hse_bs['vbm'];hse_cbm = hse_vbm + hse_bs['gap']
        vbm_shift = pbe_vbm - hse_vbm
        cbm_shift = hse_cbm - pbe_cbm
        return vbm_shift, cbm_shift

    def screen_materials(self, mpidfile='/home/gyu/scripts/mpid.dat',output_filename=None,vbmgap_from='hse', potalign_tol=0.06):
        mpids=get_mpids(mpidfile)
        if output_filename is None:
            f=open('screen_output.dat','w')
        else:
            f=open(output_filename,'w')
        keywords=['mpid','formula','spaceGroup','hse_gap','pbe_gap','vbm_shift','percent','cbm_shift','BDD?','epsilon','p_dopability','q+_lost',
                'q+_potalign_bad','n_dopability','q-_lost','q-_potalign_bad']
        f.write('\t'.join(keywords)+'\n')
        screen = {}
        for mpid in mpids:
            print mpid
            stru = self.do.to.get_structure(mpid)
            spa=SpacegroupAnalyzer(stru)
            sg = spa.get_space_group_symbol()
            formula = stru.composition.reduced_formula.encode()
            gap_check = self._check_gap_vbm(mpid)
            pbe_bs = gap_check['dft']; pbe_gap = pbe_bs['gap']
            hse_bs = gap_check['hse']; hse_gap = hse_bs['gap']
            kws = {'pbe_bs':pbe_bs,'hse_bs':hse_bs}
            vbm_shift,cbm_shift = self.vbm_cbm_shift(mpid,**kws)
            vbm_shift_perc = str(round(vbm_shift/pbe_gap*100,1))+'%' if vbm_shift else ''
            vbm_shift = str(round(vbm_shift,2)) if vbm_shift else vbm_shift
            cbm_shift = str(round(cbm_shift,2)) if cbm_shift else cbm_shift
            pbegap = '' if pbe_bs['gap'] is None else str(round(pbe_bs['gap'],2))
            hsegap = '' if hse_bs['gap'] is None else str(round(hse_bs['gap'],2))
            rawdata = self._check_rawdata(mpid)
            defect = str(rawdata['defect'])[0];bulk=str(rawdata['bulk'])[0];dielectric=str(rawdata['dielectric'])[0]
            if not hsegap or set([defect,bulk,dielectric])!=set(['T']):
                f.write('\t'.join([mpid,formula,sg,hsegap,pbegap,vbm_shift,vbm_shift_perc,cbm_shift,bulk+dielectric+defect])+'\n')
                continue
            try:
                report=self.report_dopability(mpid,vbmgap_from,potalign_tol)
            except:
                try:
                    self.update_database( mpid, from_rawdata=False)
                    report=self.report_dopability(mpid,vbmgap_from,potalign_tol)
                except:
                    comment = 'Kumagai_corr: Failed'
                    f.write('\t'.join([mpid,formula,sg,hsegap,pbegap,vbm_shift,vbm_shift_perc,cbm_shift,bulk+dielectric+defect])+'\t'+comment+'\n')
                    continue
            epsilon_avg = self.get_epsilon(mpid)[0]
            pdopability = report['p_type_dopability']
            ndopability = report['n_type_dopability']
            defects_need = self.defects_need_rerun(mpid,vbmgap_from,potalign_tol)
            lost = defects_need['lost']
            potalign_bad =defects_need['potalign_bad']
            lost_plus = {name:[i for i in lost[name] if i>0] for name in lost}
            lost_plus = {name:lost_plus[name] for name in lost_plus if lost_plus[name]}
            if not lost_plus:
                lost_plus = ''

            lost_minus = {name:[i for i in lost[name] if i<0] for name in lost}
            lost_minus = {name:lost_minus[name] for name in lost_minus if lost_minus[name]}
            if not lost_minus:
                lost_minus = ''

            potalign_bad_plus = {name:[i for i in potalign_bad[name] if i>0] for name in potalign_bad}
            potalign_bad_plus = {name:potalign_bad_plus[name] for name in potalign_bad_plus if potalign_bad_plus[name]}
            if not potalign_bad_plus:
                potalign_bad_plus = ''

            potalign_bad_minus = {name:[i for i in potalign_bad[name] if i<0] for name in potalign_bad}
            potalign_bad_minus = {name:potalign_bad_minus[name] for name in potalign_bad_minus if potalign_bad_minus[name]}
            if not potalign_bad_minus:
                potalign_bad_minus = ''
            
            record_dict =[mpid,formula,sg,hsegap,pbegap,vbm_shift,vbm_shift_perc,cbm_shift,bulk+dielectric+defect,str(epsilon_avg),pdopability,
                    str(lost_plus).encode(),str(potalign_bad_plus).encode(),ndopability,
                    str(lost_minus).encode(),str(potalign_bad_minus).encode()] 
            f.write('\t'.join(record_dict)+'\n')
        f.close()

    def defects_need_rerun(self, mpid, vbmgap, potalign_tol=0.06):
        """
        vbmgap: 'hse' or 'dft', the source of vbm and gap
        label: hse or dft. If bad p-dopability, positive lost and potalign unconverged  charged defects are not necessary to rerun
        potalign_tol: defects with std_devia more than potalign_tol are thought of being potalign unconverged
        """
        dfcs_lost = self.do.get_lost_defects(mpid)
        dfcs_unconverged = self.do.get_potalign_unconverged_defects(mpid,potalign_tol)
        dopability = self.report_dopability(mpid, vbmgap, True, potalign_tol)
        pdopability = dopability['p_type_dopability']
        ndopability = dopability['n_type_dopability']
        if pdopability == 'Bad':
            dfcs_lost = {i:[j for j in dfcs_lost[i] if j <= 0] for i in dfcs_lost}
            dfcs_unconverged = {i:[j for j in dfcs_unconverged[i] if j <= 0] for i in dfcs_unconverged}
        if ndopability == 'Bad':
            dfcs_lost = {i:[j for j in dfcs_lost[i] if j >= 0] for i in dfcs_lost}
            dfcs_unconverged = {i:[j for j in dfcs_unconverged[i] if j >= 0] for i in dfcs_unconverged}
        qrange_shrinked = self.do.shrink_defect_charge_range(mpid, vbmgap, corr=True, potalign_tol=potalign_tol)
        dfcs_lost_end = {}
        for i in dfcs_lost:
            if dfcs_lost[i]:
                if i not in qrange_shrinked:
                    dfcs_lost_end[i]=dfcs_lost[i]
                else:
                    qrange=qrange_shrinked[i]
                    qmin = qrange[0]; qmax=qrange[1]
                    new_qs=[q for q in dfcs_lost[i] if qmin<= q <=qmax]
                    if new_qs:
                        dfcs_lost_end[i]=[q for q in dfcs_lost[i] if qmin<= q <=qmax]
        ##### shrink #####
        dfcs_unconverged_end={}
        for i in dfcs_unconverged:
            if dfcs_unconverged[i]:
                if i not in qrange_shrinked:
                    dfcs_unconverged_end[i]=dfcs_unconverged[i]
                else:
                    qrange=qrange_shrinked[i]
                    qmin = qrange[0]; qmax=qrange[1]
                    new_qs=[q for q in dfcs_unconverged[i] if qmin<=q<=qmax]
                    if new_qs:
                        dfcs_unconverged_end[i]=new_qs
        #### just leave the ones whose formation energy less than 0 befort correction ######
        forms_nocorr = self.get_formation_energies_without_corr(mpid,vbmgap_from=vbmgap)
        dfcs_unconverged_lt0_nocorr={}
        for dfc in dfcs_unconverged_end:
            for q in dfcs_unconverged[dfc]:
                where = 'VBM' if q >=0 else 'CBM'
                full_name=dfc+'_'+str(q)
                forms=[]
                for reg in forms_nocorr:
                    forms.append(forms_nocorr[reg][full_name][where])
                if min(forms)<0:
                    if dfc not in dfcs_unconverged_lt0_nocorr:
                        dfcs_unconverged_lt0_nocorr[dfc]=[]
                    dfcs_unconverged_lt0_nocorr[dfc].append(q)
        return {'lost':dfcs_lost_end,'potalign_bad':dfcs_unconverged_lt0_nocorr} 

    def update_database(self, mpid, from_rawdata=True, path_save_potalign=None):
        """
        Args:
            mpid: Materials Project id
            path_potalign_save: the directory for saving the potential alignment plots
            potalign_tol: for all points of POT[qb]-POT[model], when standard deviation is less than potalign_tol
                          potential alignment is thought of as converged, otherwise unconverged.
        """
        self.refresh_database()
        rawdata = self._check_rawdata(mpid)
        check_gap = self._check_gap_vbm(mpid)
        if set(rawdata.values()) != set([True]):
            return
        if [check_gap[i]['gap'] for i in ['hse','dft']] == [None,None]:
            return 
        if from_rawdata:
            self.do.to.calculate_defect_correction(mpid, db_update=True, path_save_potalign=path_save_potalign)
        hsegap = False if check_gap['hse']['gap'] is None else True
        pbegap = False if check_gap['dft']['gap'] is None else True
        self.do._update_defects(mpid, hsegap, pbegap)
        self.refresh_database()
        print "Update successfully!"

    def update_database_all(self, mpidsfile='/home/gyu/scripts/mpid.dat', from_rawdata=True, path_potalign_save=None):
        mpids = get_mpids(mpidsfile)
        for mpid in mpids:
            try:
                self.update_database(mpid,from_rawdata,path_potalign_save)
            except:
                pass

    def remove_from_database(self,mpid):
        """
        Will remove the raw data in tasks collection and the record in defect collection
        """
        self.do.to.del_rawdata(mpid)
        r=self.do.collection.find_one({'mpid':mpid},[''])
        try:
            self.do.collection.delete_one({'_id':r['_id']})
        except:
            pass
        

    def report_potalign_convergence(self,mpid,potalign_tol=0.06):
        ids_dfc = self.do.to.groups['defect'][mpid]
        state = {}
        state['vac']={}
        state['subst']={}
        for idd in ids_dfc:
            record = self.do.to.collection.find_one({'_id':idd},['transformations','kumagai_correction'])
            charge = record['transformations']['charge']
            defect = record['transformations']['defect_type']
            std_devia = record['kumagai_correction']['std_devia']
            if std_devai<=potalign_tol:
                potalign_conv=True
            else:
                potalign_conv=False
            if 'vac' in defect:
                if defect not in state['vac']:
                    state['vac'][defect]={}
                state['vac'][defect][charge]=potalign_conv
            if 'subst' in defect:
                if defect not in state['subst']:
                    state['subst'][defect]={}
                state['subst'][defect][charge]=potalign_conv
        return state

    def report_potalign_convergence_summary(self,mpidfile):
        mpids = get_mpids(mpidfile)
        report={}
        report['vac']={};report['subst']={}
        for tp in ['vac','subst']:
            report[tp]['num_total'] = 0
            report[tp]['num_converged'] = 0
        report['materials']={}
        report['materials']['num_converged']={}
        for tp in ['vac','subst','vac_subst']:
            report['materials']['num_converged'][tp]=0
        for i in ['num_inlist','num_with_defects']:
            report['materials'][i]=0

        def get_state_indetail(state):
            ste={}
            for tp in ['vac','subst']:
                state_coll=[]
                for dfc in state[tp]:
                    for chrg in state[tp][dfc]:
                        state_coll.append(state[tp][dfc][chrg])
                ste[tp]=True if set(state_coll)==set([True]) else False
            ste['vac_subst']=True if ste['vac'] and ste['subst'] else False
            return ste

        for mpid in mpids:
            report['materials']['num_inlist'] += 1
            try:
                state = self.report_potalign_convergence(mpid) 
                state_indetail = get_state_indetail(state)
                mats={}
                for tp in state_indetail:
                    if state_indetail[tp]:
                        report['materials']['num_converged'][tp]+=1
                report['materials']['num_with_defects'] += 1

            except:
                continue
            for tp in ['vac','subst']:
                for dfc in state[tp]:
                    for chrg in state[tp][dfc]:
                        if chrg:
                            report[tp]['num_total']+=1
                            if state[tp][dfc][chrg]:
                                report[tp]['num_converged']+=1
        
        return report


    def report_dopability(self,mpid, label, corr=True, potalign_tol=0.06, in_detail=False):
        """
        Args:
            mpid: the material id in Materials Project
            label: can be hse or dft
        """
        report = self.do.report_dopability(mpid, label, corr, potalign_tol, in_detail)
        return report
    
    def plot_formation_energies(self, mpid, label, corr=True, ext_elts=[], potalign_tol=0.06, xlim=None, ylim=(-1,3), path='.', fmt='jpg'):
        self.do.plot_formation_energies(mpid, label, corr, ext_elts, potalign_tol, xlim, ylim, path, fmt)
        #try:
        #    self.do.plot_formation_energies(mpid, label, corr, potalign_tol, xlim, ylim, path, fmt)
        #except:
        #    print 'Please update database for %s first!' % mpid

    def _plot_formation_energies_old_method(self,mpid, label, xlim=None, ylim=(-1,3), path='', fmt='png'): 
        """
        Args:
            mpid: the material id in Materials Project
            label: only can be hse or dft
            xlim: electronic fermi energy range for defect formation energies plot
            ylim: formation energy range for show in defects plot
            path: the path for saving plots
            vbm:  the value you want, not the dft vbm from Materials Project
                    or the hse vbm from our own database
            gap: as vbm
        """
        if label not in ['hse','dft']:
            raise KeyError("label must be 'hse' or 'dft'!!")
        else:
            try:
                rawdata = self._check_rawdata(mpid, label)
                if set(rawdata.values()) != set([True]):
                    return
                self.do._plot_formation_energies_backup(mpid, label, xlim=xlim, ylim=ylim, path=path, fmt=fmt, vbm=vbm, gap=gap)
            except:
                self.update_database(mpid, label, from_rawdata=False, report_update = True)
                self.do._plot_formation_energies_backup(mpid, label, xlim=xlim, ylim=ylim, path=path, fmt=fmt, vbm=vbm, gap=gap)

    def generate_input_files_further(self,mpid, vbmgap='hse', dist_step=1.0, hse=False, potalign_tol=0.06,path='./'):
        """
        Generate defect input files further for the defects lost \
        and potalign unconverged with larger defect distance
        """
        self._insert_structures_tag_into_db(mpid)
        ### just get the composition of bulk with any size ### 
        id_bulk = self.do.to.groups['bulk'][mpid][0]
        record_blk = self.do.to.collection.find_one({'_id':id_bulk},['pretty_formula'])
        from pymatgen.core.composition import Composition
        comp = Composition(record_blk['pretty_formula'])
        #######################################################
        def defects_bad():
            structs={}
            ids_dfc = self.do.to.groups['defect'][mpid]
            for idd in ids_dfc:
                record =self.do.to.collection.find_one({'_id':idd},['structures','transformations'])
                name = record['transformations']['defect_type']
                charge = record['transformations']['charge']
                structures = record['structures']
                if name not in structs:
                    structs[name]={}
                structs[name][charge]=structures
            dfcs_rerun = self.defects_need_rerun(mpid, vbmgap, potalign_tol) 
            potalign_bad = dfcs_rerun['potalign_bad']
            lost = dfcs_rerun['lost']
            dfcs = self.do.get_defects(mpid)
            dfc_ok={}
            for name in dfcs:
                if name in potalign_bad:
                    if set(potalign_bad[name])!=set(dfcs[name]):
                        dfc_ok[name] = [i for i in dfcs[name] if i not in potalign_bad[name]]
                else:
                    dfc_ok[name]=dfcs[name]
            out=[]
            for name in potalign_bad:
                if name in dfc_ok and 0 in dfc_ok[name]:
                    if len(dfc_ok[name])==1:
                        out.append({'name':name,'charge':0, 'structures':structs[name][0],'bigger_dist':True})
                for q in potalign_bad[name]:
                    out.append({'name':name,'charge':q, 'structures':structs[name][q],'bigger_dist':True})
            for name in lost:
                for q in lost[name]:
                    if not q:
                        if name in dfc_ok:
                            q0=dfc_ok[name][0]
                            structures = structs[name][q0]
                            out.append({'name':name,'charge':0, 'structures':structures,'bigger_dist':False})
                        else:
                            q0=dfcs[name][0]
                            structures = structs[name][q0]
                            out.append({'name':name,'charge':0, 'structures':structures,'bigger_dist':True})
                    else:
                        sign = np.sign(q)
                        qs = [i for i in dfcs[name] if np.sign(i)==sign or not i]
                        qs.append(q)
                        if sign == 1:
                            qs.sort()
                        elif sign == -1:
                            qs.sort(reverse=True)
                        if qs[0]!=q:
                            q_closest = qs[qs.index(q)-1]
                        else:
                            q_closest = qs[1]
                        structures=structs[name][q_closest]
                        if name in dfc_ok and q_closest in dfc_ok[name]:
                            out.append({'name':name,'charge':q, 'structures':structures,'bigger_dist':False})
                        else:
                            out.append({'name':name,'charge':q, 'structures':structures,'bigger_dist':True})
            return out
        defects_rerun = defects_bad()
        defects = []
        print 'Generate input files for:'
        bulks=[]
        ids_blk = self.do.to.groups['bulk'][mpid]
        for idd in ids_blk:
            rec = self.do.to.collection.find_one({'_id':idd},['input.crystal'])
            stru_blk = Structure.from_dict(rec['input']['crystal'])
            bulks.append(stru_blk)
        diels=[]
        ids_diel = self.do.to.groups['dielectric'][mpid]
        for idd in ids_diel:
            rec = self.do.to.collection.find_one({'_id':idd},['input.crystal'])
            stru_diel = Structure.from_dict(rec['input']['crystal'])
            diels.append(stru_diel)
        for dfc in defects_rerun:
            print '\t %s %s' % (dfc['name'],str(dfc['charge']))
            bulk_sc = Structure.from_dict(dfc['structures']['bulk_supercell'])
            stru_dfc = Structure.from_dict(dfc['structures']['defect_no_relaxation'])
            if  dfc['bigger_dist']:
                defects_one = DefectRemaker(bulk_sc,stru_dfc,dfc['name'],dfc['charge'],dfc_dist_add=dist_step).defects
            else:
                bulk_uc = Structure.from_dict(dfc['structures']['bulk_unitcell'])
                kws={'stru_bulk_uc':bulk_uc}
                defects_one = DefectRemaker(bulk_sc,stru_dfc,dfc['name'],dfc['charge'],dfc_dist_add=0.0,**kws).defect_backup
            for dfc in defects_one:
                if dfc['short_name'] == 'bulk':
                    bulk = dfc['supercells'][0]['structure']
                    existed = False
                    for blk in bulks:
                        if are_structures_equal(blk, bulk):
                            existed = True
                            break
                    if not existed:        
                        bulks.append(bulk)
                        defects.append(dfc)
                elif dfc['short_name'] == 'dielectric':
                    diel = dfc['structure']
                    existed = False
                    for diel_ in diels:
                        if are_structures_equal(diel, diel_):
                            existed = True
                            break
                    if not existed:        
                        diels.append(diel)
                        defects.append(dfc)
                else:
                    defects.append(dfc)
        make_vasp_defect_files(defects,path,mpid,comp,hse)

    def _insert_structures_tag_into_db(self,mpid):
        grps=self.do.to.get_group_on_mpid(mpid)
        for grp in grps:
            id_bulk_uc = grp['dielectric']
            rec_blk_uc = self.do.to.collection.find_one({'_id':id_bulk_uc},['input.crystal'])
            blk_uc_dict = rec_blk_uc['input']['crystal']
            id_bulk_sc = grp['bulk']
            rec_blk_sc = self.do.to.collection.find_one({'_id':id_bulk_sc},['input.crystal'])
            blk_sc_dict = rec_blk_sc['input']['crystal']
            blk_sc = Structure.from_dict(blk_sc_dict)
            for dfc_id in grp['defect']:
                structures={}
                structures['bulk_unitcell'] = blk_uc_dict
                structures['bulk_supercell'] = blk_sc_dict
                rec_dfc = self.do.to.collection.find_one({'_id':dfc_id},['input.crystal',
                    'transformations.defect_type','structures'])
                if 'structures' in rec_dfc:
                    continue
                name_ref = rec_dfc['transformations']['defect_type']
                def del_num(string):
                    out = ''
                    for i in string:
                        if i.isdigit():
                            continue
                        out+=i
                    return out
                name_ref = unicode(del_num(name_ref))
                dfc_stru_dict = rec_dfc['input']['crystal']
                stru_dfc = Structure.from_dict(dfc_stru_dict)
                pos_blk, pos_dfc=find_defect_pos(blk_sc, stru_dfc)
                if pos_blk is None and pos_dfc is None:
                    raise KeyError("Bulk_sc doesn't match defect structure!!")
                elif pos_blk is None:
                    dfc_type ='inter'
                elif pos_dfc is None:
                    dfc_type = 'vac'
                else:
                    dfc_type = 'subst'
                blk_sc_tmp = blk_sc.copy()
                if dfc_type == 'vac':
                    site = closestsites(blk_sc,stru_dfc,pos_blk)[0][0]
                    index = blk_sc.sites.index(site)
                    name = site.specie.symbol+'_vac'
                    if name != name_ref:
                        raise KeyError("Defect structure donsn't match defect name!!")
                    blk_sc_tmp.remove_sites([index])
                elif dfc_type == 'subst':
                    site_blk = closestsites(blk_sc,stru_dfc,pos_blk)[0][0]
                    index = blk_sc.sites.index(site_blk)
                    site_blk_symbol = site_blk.specie.symbol
                    site_dfc = closestsites(blk_sc,stru_dfc,pos_blk)[1][0]
                    site_dfc_specie = site_dfc.specie
                    site_dfc_symbol = site_dfc.specie.symbol
                    name = site_blk_symbol+'_subst_'+site_dfc_symbol
                    if name != name_ref:
                        raise KeyError("Defect structure donsn't match defect name!!")
                    blk_sc_tmp.replace(index,site_dfc_specie)
                elif dfc_type == 'inter':
                    raise KeyError("Can not deal with interstial now!!")
                structures['defect_no_relaxation'] = blk_sc_tmp.as_dict()
                self.do.to.collection.update({'_id':dfc_id},{'$set':{'structures':structures}})

