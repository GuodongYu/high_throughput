from high_throughput.defects.database import HseGapOperater
import glob
hg=HseGapOperater()
dirs=glob.glob('/home/acad/ucl-naps/yugd/workplace/hse_gap_pbe_stru/3-5_semiconductors/new/*')

for dir in dirs:
    try:
    	hg.insert_new_record(dir)
    except:
	pass

