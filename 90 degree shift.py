# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:48:44 2019

@author: Yuan Ji

input: raw signal

output: 90 degree shifted and write into txt files
"""
import os
import numpy as np
#import qt



   
def _shift(signal):
    #shift raw cosine-like data to sine-like curve for a better fitting
    '''
    sometimes the maxima does not appear at the first data point
    '''
    return signal[28:]+signal[:28]

def _find_txt(root):
    pure_txt_name=[]
    txt_full_dir=[]
    g = os.walk(root) 
    for par_dir, _, files in g: 
        for file in files: 
            if file.endswith(".txt"):
                name=file.strip('.txt')
                #print(name)
                pure_txt_name.append(name)
                filepath = os.path.join(par_dir, file) 
                #print(filepath)
                txt_full_dir.append(filepath)

    #    
        #print(pure_txt_name)
    return txt_full_dir, pure_txt_name  
        

def read_shift_write(txt_full_path, name, output_path):
    #get 1w and 2w signal
    angle=[]
    R1_raw = []
    V1_phase = []
    V1_raw = []
    R2_raw = []
    V2_phase = []
    V2_raw = []
    V1=[]
    R1=[]
    V2=[]
    R2=[]
    
    file_to_read=open(txt_full_path , 'r') 
    
    lines = file_to_read.readlines() # 整行读取数据
    
    '''
    definition:
        R1_raw, R2_raw: read from SR830
        gain: set in SR 560, belongs to [10, 100, 1000]
        R1_raw = V1_raw/I
        R2_raw = \sqrt(2)*V2_raw/(I^2)
        V1_actual = V1_raw * gain
        R1_actual = R1_raw * gain
        V2_actual = V2_raw * gain
        R2_actual = V1_raw * gain
    '''
    
    for line in lines:
        angle_temp,R1_raw_, V1_phase_, V1_raw_, R2_raw_, V2_phase_, V2_raw_, R1_actual, V1_actual, R2_actual, V2_actual  = [float(i) for i in line.split()] 
        angle.append(angle_temp)
        
        R1_raw.append(R1_raw_)
        V1_phase.append(V1_phase_)
        V1_raw.append(V1_raw_)
        
        R2_raw.append(R2_raw_)
        V2_phase.append(V2_phase_)
        V2_raw.append(V2_raw_)
        
        V1.append(V1_actual)
        R1.append(R1_actual)
        V2.append(V2_actual)
        R2.append(R2_actual)
    
    #caculation
    current1 = np.true_divide(V1_raw, R1_raw)
    currentSquare = np.true_divide(V2_raw, R2_raw)*np.sqrt(2)
    #print(current1)
    #print(currentSquare)
    current = np.average(current1)
    currentsquare = np.average(currentSquare)
    print('Current: '+str(int(current*1000000))+'uA', '\nCurrent^2(A^2): '+str(currentsquare))
    
    gain_1 = np.true_divide(V1_raw, V1)
    gain_2 = np.true_divide(V2_raw, V2)
    gain1 = np.average(gain_1)
    gain2 = np.average(gain_2)
    print('gain: '+str(gain1)+', '+str(gain2))
    
    shifted_V1_raw = _shift(V1_raw)
    shifted_V1_phase = _shift(V1_phase)
    
    shifted_V2_raw = _shift(V2_raw)
    shifted_V2_phase = _shift(V2_phase)
    
    shifted_R1 = _shift(R1)
    shifted_V1 = _shift(V1)
    shifted_R2 = _shift(R2)
    shifted_V2 = _shift(V2)
    
    
    #write into txt
    f = open(os.path.join(output_path, name+'.txt'), 'w')
    for _angle,sV1_raw, sV1_phase, sV2_raw, sV2_phase, sR1, sV1, sR2, sV2 in zip(angle, shifted_V1_raw, shifted_V1_phase, shifted_V2_raw, shifted_V2_phase, shifted_R1, shifted_V1, shifted_R2, shifted_V2):
        f.write(str(_angle)+'\t'+str(sV1_raw)+'\t'+str(sV1_phase)+'\t'+str(sV2_raw)+'\t'+str(sV2_phase)+'\t'+str(sV1)+'\t'+str(sV2)+'\n')
    f.close()
    
def _process(input_root, output_root):
    txt_full_dir, pure_txt_name = _find_txt(input_root)
    for dire, name in zip(txt_full_dir, pure_txt_name):
        read_shift_write(dire, name, output_root)



input_root = r"D:\BaiduNetdiskDownload\Cr2O3 magnon\20190625 B264\Cr2O3 11-20 8w 18kV 21nm 20190418\10 pin 9-10-7-8-11-12 2K 500uA\data"
output_root = r"D:\BaiduNetdiskDownload\Cr2O3 magnon\20190625 B264\Cr2O3 11-20 8w 18kV 21nm 20190418\10 pin 9-10-7-8-11-12 2K 500uA\shifted_added"
_process(input_root, output_root)