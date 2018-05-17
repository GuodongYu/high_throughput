from high_throughput.utils.query import common_binary_semiconductors
from high_throughput.defects.utils.vasp import generate_input_files
have = ['mp-1639','mp-2172','mp-2624','mp-1700','mp-2490','mp-4824','mp-830','mp-20351','mp-2524','mp-8883']

data = common_binary_semiconductors('3-5')
mpids = [i['task_id'] for i in data]

for mpid in mpids:
    if mpid in have:
        continue
    generate_input_files(mpid,path='/home/gyu/3_5_semiconductor',min_dis=11.0)


