import subprocess
import numpy as np
import json
import sys
import datetime
import random
import os
import openslide
from bson import json_util
from pymongo import MongoClient

is_shifted = False

# something need to be changed
# wsi, perdfile, output
# studyid,  executionid
# cell_idx
# heat_list

# columns info: 0-x, 1-y, 2-epithelium, 3-stroma, 4-tumor, 5-necrosis, 6-dysplasia, 7-tumor/non-tumor

wsi = '/data10/shared/tcga_all/coad/'
predfile = '/data05/shared/yuwei/coad_CM-5348/cancer_txt/coad_model/'
output = '/data05/shared/yuwei/coad_CM-5348/cancer_json/coad_model/tumor'

def write_json(svs, pred, output):
    fullname = pred.split(predfile)[1].split('prediction-')[1].split('.')[0]
    subjectid = '-'.join(fullname.split('-')[0:3])
    """if subjectid == caseid"""
    #caseid = subjectid
    """if subjectid != caseid"""
    caseid = '-'.join(fullname.split('-')[3:6])

    print("subject id:",subjectid, "caseid:",caseid)
    # initial cell type: epithelium, stroma, tumor...
    cell_idx = 2
    # the number of cell types
    n_heat=2
    heat_list = ['tumor','lymph']
    weight_list = [0.5,0.5]
    
    heatmapfile = output + '/heatmap_' + pred.split('prediction-')[1] + '.json'
    #metafile = output + '/meta_' + pred.split('prediction-')[1] + '.json'

    oslide = openslide.OpenSlide(svs)
    slide_width_openslide = oslide.dimensions[0]
    slide_height_openslide = oslide.dimensions[1]

    #print("Loaded wsi and predictionid: ", svs, pred)

    # initialize param
    x_arr = np.zeros(10000000, dtype = np.float64)
    y_arr = np.zeros(10000000, dtype = np.float64)
    score_arr = np.zeros(10000000, dtype = np.float64)
    score_set_arr = np.zeros(shape=(10000000, n_heat), dtype = np.float64)
    idx = 0

    #read info in prediction file
    with open(pred) as infile:
        for line in infile:
            parts = line.split(' ')
            x = int(parts[0])
            y = int(parts[1])
            score = float(parts[cell_idx])
            score_set = score

            x_arr[idx] = x
            y_arr[idx] = y
            score_arr[idx] = score
            score_set_arr[idx,0]=score_set
            idx += 1

    x_arr = x_arr[:idx]
    y_arr = y_arr[:idx]
    score_arr = score_arr[:idx]   # related to metric_value
    score_set_arr = score_set_arr[:idx]    #related to metric_array

    patch_width = max(x_arr[1] - x_arr[0], y_arr[1] - y_arr[0])
    patch_height = patch_width

    slide_width = int(slide_width_openslide)
    slide_height = int(slide_height_openslide)

    x_arr = x_arr / slide_width
    y_arr = y_arr / slide_height
    patch_width = patch_width / slide_width
    patch_height = patch_height / slide_height
    patch_area = patch_width * slide_width * patch_height * slide_height

    dict_img = {}
    dict_img['case_id'] = caseid
    dict_img['subject_id'] =  subjectid

    analysis_execution_date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

    dict_analysis = {}
    dict_analysis['cancer_type'] = 'quip'
    dict_analysis['study_id'] = 'TCGA:coad'
    dict_analysis['execution_id'] = 'tumor_coad_model'
    dict_analysis['source'] = 'computer'
    dict_analysis['computation'] = 'heatmap'
    dict_analysis['exexcution_date'] = analysis_execution_date


    if (is_shifted == True):
        shifted_x = -3*patch_width / 4.0
        shifted_y = -3*patch_height / 4.0
    else:
        shifted_x = 0
        shifted_y = 0

    with open(heatmapfile, 'w') as f:
        for i_patch in range(idx):
            dict_patch = {}
            dict_patch['type'] = 'Feature'
            dict_patch['parent_id'] = 'self'
            dict_patch['normalized'] = 'true'
            dict_patch['object_type'] = 'heatmap_multiple'
            dict_patch['footprint'] = patch_area
            dict_patch['x'] = x_arr[i_patch] + shifted_x
            dict_patch['y'] = y_arr[i_patch] + shifted_y
            
            x1 = dict_patch['x'] - patch_width/2
            x2 = dict_patch['x'] + patch_width/2
            y1 = dict_patch['y'] - patch_height/2
            y2 = dict_patch['y'] + patch_height/2
            dict_patch['bbox'] = [x1,y1,x2,y2]

            dict_geo = {}
            dict_geo['type'] = 'Polygon'
            dict_geo['coordinates'] = [[[x1,y1],[x2,y1],[x2,y2],[x1,y2],[x1,y1]]]
            dict_patch['geometry'] = dict_geo

            dict_prop = {}
            dict_prop['metric_value'] = score_arr[i_patch]
            dict_prop['metric_type'] = 'tile_dice'
            dict_prop['human_mark'] = -1
            
            dict_multiheat = {}
            dict_multiheat['human_weight'] = -1
            dict_multiheat['weight_array'] = weight_list
            dict_multiheat['heatname_array'] = heat_list
            dict_multiheat['metric_array'] = score_set_arr[i_patch].tolist()
            
            dict_prop['multiheat_param'] = dict_multiheat
            dict_patch['properties'] = dict_prop
            
            dict_provenance = {}
            dict_provenance['image'] = dict_img
            dict_provenance['analysis'] = dict_analysis
            dict_patch['provenance'] = dict_provenance

            dict_patch['date'] = datetime.datetime.now()
            
            json.dump(dict_patch, f, default = json_util.default)
            f.write('\n')

if __name__ == "__main__":
    wsi_cases = os.listdir(wsi)
    pred_txt_files = os.listdir(predfile)
    json_files = os.listdir(output)
    for p in pred_txt_files:
        if p.startswith('prediction-'):
            pred_path = predfile + p
            svs_name = p.split('prediction-')[1]+'.svs'
            json_name = 'heatmap_'+p.split('prediction-')[1]+'.json'
            if json_name in json_files:
                print(json_name + "exist, pass!")
            else:
                if svs_name in wsi_cases:
                    svs_path = wsi + svs_name
                    write_json(svs_path,pred_path,output)
                else:
                    print(svs_name+' doesn\'t exist.')
        else:
            pass
