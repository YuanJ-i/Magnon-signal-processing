# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:48:44 2019

@author: Yuan Ji

input: raw signal

output: 90 degree shifted and write into txt files
"""
import os
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
    V_omega=[]
    R_omega=[]
    V_2omega=[]
    R_2omega=[]
    
    file_to_read=open(txt_full_path , 'r') 
    
    lines = file_to_read.readlines() # 整行读取数据
    
    for line in lines:
        angle_temp,a, b, c,d,e,f ,R_omega_temp, V_omega_temp, R_2omega_temp, V_2omega_temp  = [float(i) for i in line.split()] 
        angle.append(angle_temp)
        V_omega.append(V_omega_temp)
        R_omega.append(R_omega_temp)
        V_2omega.append(V_2omega_temp)
        R_2omega.append(R_2omega_temp)
    
    
    shifted_R1 = _shift(R_omega)
    shifted_V1 = _shift(V_omega)
    shifted_R2 = _shift(R_2omega)
    shifted_V2 = _shift(V_2omega)
    
    #write into txt
    f = open(output_path+'//'+name+'.txt', 'w')
    for _angle, sR1, sV1, sR2, sV2 in zip(angle, shifted_R1, shifted_V1, shifted_R2, shifted_V2):
        f.write(str(_angle)+'\t'+str(sR1)+'\t'+str(sV1)+'\t'+str(sR2)+'\t'+str(sV2)+'\n')
    f.close()
    
def _process(input_root, output_root):
    txt_full_dir, pure_txt_name = _find_txt(input_root)
    for dire, name in zip(txt_full_dir, pure_txt_name):
        read_shift_write(dire, name, output_root)



input_root = "F:\\BaiduNetdiskDownload\\Measurement\\20190618 B264\\Cr2O3 11-20 8w 18kV 21nm 20190418\\4-pin 11-12-9-10-13-14 2K 500uA\\data"
output_root = "F:\\BaiduNetdiskDownload\\Measurement\\20190618 B264\\Cr2O3 11-20 8w 18kV 21nm 20190418\\4-pin 11-12-9-10-13-14 2K 500uA\shifted_data"
_process(input_root, output_root)