
from tqdm.auto import tqdm
from dgfi_if import SWORD_IF, SWORD_Reach, MVA_IF
from pathlib import Path
import warnings
from tqdm.contrib.concurrent import process_map, thread_map
tqdm.pandas()


SAVE_ROOT = Path('./reaches/')
MAX_WORKERS = 100
SWORD_VERSION = '2.0'
along_angle_th = 65
along_conf_th = 300

sword = SWORD_IF(sword_version=SWORD_VERSION)

def get_reach_data_helper(reach_id : int):
    reach_id = int(reach_id)
    path = SWORD_Reach.get_pickle_path_from_reach_id(reach_id, SAVE_ROOT)

    r = SWORD_Reach(reach_id, sword_version = SWORD_VERSION)
    if r.n_nodes < 2:
        del r
        return None

    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore')
        r.get_icesat2_slope(silent=True,along_angle_th=along_angle_th,along_conf_th=along_conf_th)

    r.save(path)
    del r

def run():
    reach_id_dic = sword.list_reaches()
    for nc_filename, reach_ids in tqdm(reach_id_dic.items(), desc='Iterating Continents'):
        reach_ids = [x for x in reach_ids if x % 10 < 4]
        process_map(get_reach_data_helper,reach_ids,desc='Retrieving Reach Data', max_workers=MAX_WORKERS, chunksize=10)

if __name__ == "__main__":
    run()