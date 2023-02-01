import logging, os
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)
# import tum_design as tum

mission_groups = {
    'jason': {
        'passlocator': 'jason2',
        'mva': ['jason2_hf','jason3_hf', 'jason3_igdr_f_hf'],
        'color': '#ff0000',
        'output_name': 'Jason-2/3'
    },
    'envisat': {
        'passlocator': 'envisat',
        'mva': ['envisat_v3_hf'],
        'color': '#ff0000',
        'output_name': 'Envisat'
    },
    'saral': {
        'passlocator': 'saral',
        'mva': ['saral_hf'],
        'color': '#ff0000',
        'output_name': 'Saral'
    },
    'sentinel3a': {
        'passlocator': 'sentinel3a',
        'mva': ['sentinel3al_hf','sentinel3al_stc_hf'],
        'color': '#ff0000',
        'output_name': 'Sentinel-3A'
    },
    'sentinel3b': {
        'passlocator': 'sentinel3b_nominal',
        'mva': ['sentinel3bl_hf','sentinel3bl_stc_hf'],
        'color': '#ff0000',
        'output_name': 'Sentinel-3B'
    },
    'cryosat2': {
        'passlocator': 'cryosat2',
        'mva': ['cryosatD_sar_hf', 'cryosatD_lrm_hf', 'cryosatD_sin_hf'],
        'color': '#ff0000',
        'output_name': 'CryoSat-2'
    },
    'icesat': {
        'passlocator': 'icesat',
        'mva': ['icesat_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat'
    },
    'icesat2_gt1l': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt1l_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    },
    'icesat2_gt1r': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt1r_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    },
    'icesat2_gt2l': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt2l_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    },
    'icesat2_gt2r': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt2r_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    },
    'icesat2_gt3l': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt3l_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    },
    'icesat2_gt3r': {
        'passlocator': 'icesat2_nominal',
        'mva': ['icesat2_gt3r_atl13v5_hf'],
        'color': '#ff0000',
        'output_name': 'ICESat-2'
    }
}

waterfall_override = {
    # Sword Version : Reaches that are actually no waterfalls or dams
    '11': [74210000264]
}
