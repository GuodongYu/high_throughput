from pymongo.mongo_client import MongoClient
from pymongo.collection import Collection,ObjectId
from pymongo.database import Database

try:
    tasks_collection_name='tasks'
    defect_collection_name='defect'
    hsegap_collection_name='HSE_gaps_pbe_structure'
    db_info={"aliases": {},  
        "database": "results_GY", 
        "host": "marilyn.pcpm.ucl.ac.be", 
        "port": 27017,
        "admin_user": "gyu", 
        "admin_password": "pulco",
        "readonly_user": "gyu", 
        "readonly_password": "pulco"}
    
    clt = MongoClient(host=db_info['host'],port=db_info['port'])
    db = Database(clt,db_info['database'])
    db.authenticate(db_info['admin_user'],db_info['admin_password'])
except:
    print 'Can not connect database'
