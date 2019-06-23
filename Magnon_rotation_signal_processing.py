# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 19:01:43 2019

@author: Yuan Ji

NEED optimization for sine fitting
and abnormal signal exclusion
"""
import os
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t
import itertools
import numpy as np
#import qt
import re

class signal_processing(object):
    '''
    class for single data file in rotation measurement, eg. 2K 9T 500uA gain=1000.txt
    '''
    def __init__(self, folder, filename):
        self.root=folder
        self.file_full_path=self.root+'//'+filename+'.txt'
        self.fig_folder = self.root+'\\figs'
        self.temeprature, self.current, self.field, self.gain = self._parse(filename)
        self.angle, self.V_omega, self.R_omega, self.V_2omega, self.R_2omega = self._get_data()
        self.IsGoodFitting=True
        self.amplitude = 0
        self.amplitude_std = 0
        self.cleaned_shifted_signal = []
        self.fitted_signal = []
    
    def _get_info(self):
        return self.IsGoodFitting, self.amplitude, self.amplitude_std
    
    def _parse(self, filename):
        #parse the filename, extract [Temp, Curr, Field, gain]
        
        number_list=re.findall(r"\d+\.?\d*",filename)
        temperature_sting=re.findall(r"\d+\.?\d*K",filename)[0]
        temperature=float(re.findall(r"\d+\.?\d*",temperature_sting)[0])
        
    #    print(type(temperature))
        
        current_sting=re.findall(r"\d+\.?\d*uA",filename)[0]
        current=int(re.findall(r"\d+\.?\d*",current_sting)[0])
        field_sting=re.findall(r"\d+\.?\d*Oe",filename)[0]
        field=float(re.findall(r"\d+\.?\d*",field_sting)[0])
        gain=int(number_list[-1])
        
        return temperature, current, field, gain
    

    def _get_data(self):
        #get 1w and 2w signal
        angle=[]
        V_omega=[]
        R_omega=[]
        V_2omega=[]
        R_2omega=[]
        
        file_to_read=open(self.file_full_path , 'r') 
        
        lines = file_to_read.readlines() # 整行读取数据
        
        for line in lines:
            angle_temp,a, b, c,d,e,f ,R_omega_temp, V_omega_temp, R_2omega_temp, V_2omega_temp  = [float(i) for i in line.split()] 
            angle.append(angle_temp)
            V_omega.append(V_omega_temp)
            R_omega.append(R_omega_temp)
            V_2omega.append(V_2omega_temp)
            R_2omega.append(R_2omega_temp)
        
        return angle, V_omega, R_omega, V_2omega, R_2omega
    
    
    def _abnormal_data_detect_and_deletion(self, signal):
        #delete data out of range, which deviates more than 3 sigma away from average
        #return clean data
        
        #sudden jump detection
        signal_bar=np.mean(signal)
        signal_std=np.std(signal)
        cleaned_angle=[]
        cleaned_signal=[]
        for angle, signal in zip(self.angle, signal):
            if abs(signal-signal_bar)<=2*signal_std:
                cleaned_angle.append(angle)
                cleaned_signal.append(signal)
    #    print(cleaned_angle)
    #    print("raw data number: "+str(len(angle_list)))        
    #    print("cleaned data number: "+str(len(cleaned_angle)))
        return cleaned_angle, cleaned_signal

    def _sine_curve_fit(self, angle_list, signal_list):
    #    n=len(angle_list)
        #sine curve fit
        def f_fit(x,a, b, c):
            return a*np.sin(math.pi*x/180+b)+c
        def f_predict(x,a,b,c):
            return [a*np.sin(kx*math.pi/180+b)+c for kx in x]
        
        
        a_guess=max(signal_list)
        p_fit,pcov=curve_fit(f_fit,angle_list,signal_list, p0=(a_guess,0, 0))#曲线拟合
        a,b,c=p_fit.tolist()
        a_std=math.sqrt(np.diag(pcov)[0])
        #    print(p_fit)#optmized a,b,c
    #    print(pcov)#最优参数的协方差估计矩阵
        return a, a_std, f_predict(angle_list,a, b, c)    
    
    def _get_cod(self, raw_signal, predicted_signal):
        # compute adjusted coeffienct of determination
        aver=np.mean(raw_signal)
        sse=0
        sst=0
        ssr=0
        for actual, predict in zip(raw_signal, predicted_signal):
            sse=sse+(actual-predict)**2
            sst=sst+(actual-aver)**2
            ssr=ssr+(predict-aver)**2
        cod=ssr/sst
    #    cod2=1-sse/sst
    #    print("Curr: "+str(current)+" cod: "+str(cod))
    #    print("Curr: "+str(current)+" cod2: "+str(cod2))
    #    print()
        return cod
    
    
    def _noise(self):
        sequence = np.array(self.cleaned_shifted_signal, dtype = float)
        sine = np.array(self.fitted_signal, dtype = float)
        #find noise value of cleaned signal
        s=np.fft.rfft(sequence)
#        print(s)
        
        
        #naive method
        residual = sequence-sine
        
        difference = []
        for i in range(len(residual)):
            if i==0:
                pass
            difference.append(abs(residual[i]-residual[i-1]))
            
        noise_max = max(difference)
        noise_aver = np.average(np.array(difference))
        noise_median = np.median(np.array(difference))

#        print(str(self.current)+'uA: max noise'+str(noise_max))
#        print(str(self.current)+'uA: aver noise'+str(noise_aver))
#        print(str(self.current)+'uA: median noise'+str(noise_median))
        
        return noise_median
        
    
    def _shift(self, signal):
        #shift raw cosine-like data to sine-like curve for a better fitting
        '''
        sometimes the maxima does not appear at the first data point
        '''
        return signal[28:]+signal[:28]
    
    def _process(self):
        '''core method'''
        #90 degree shift of V2w
        V_2omega_shifted=self._shift(self.V_2omega)
        #out of range data deletion
        angle_cleaned, V_2omega_shifted_cleaned=self._abnormal_data_detect_and_deletion(V_2omega_shifted)
        raw_data=V_2omega_shifted_cleaned
        self.cleaned_shifted_signal = V_2omega_shifted_cleaned
        #sine curve fit
        self.amplitude, self.amplitude_std , predicted_value = self._sine_curve_fit(angle_cleaned, raw_data) 
        
        self.fitted_signal = predicted_value
        cod = self._get_cod(raw_data, predicted_value)
        #set amplitude to 0 if cod<0.5
        if cod<0.5:
            self.IsGoodFitting=False
            self.amplitude=0
            
        return angle_cleaned, V_2omega_shifted_cleaned, predicted_value
    
    def _plot_raw(self, fig_folder):
        #plot raw data
        plt.figure()
        plt.plot(self.angle , self.V_2omega, 'r',label='raw',marker='.')
        plt.title(str(self.temeprature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-RAW')
        plt.xlabel('Angle')
        plt.ylabel('V2w')
        plt.legend()
        plt.savefig(fname=fig_folder+'\\'+str(self.temperature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-RAW.jpg')
        plt.close()
        
    def _plot_shifted(self):
        plt.figure()
        plt.plot(self.angle , self._shift(self.V_2omega), 'r',label='raw',marker='.')
        plt.title(str(self.temeprature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-Shifted')
        plt.xlabel('Angle')
        plt.ylabel('V2w')
        plt.legend()
        plt.show()
        
    def _plot_fitted(self, fig_folder):
        angle_cleaned, V_2omega_shifted_cleaned, predicted_value = self._process()
        #data washed
        plt.figure()
        plt.plot(angle_cleaned,V_2omega_shifted_cleaned,'r',label='cleaned',marker='.')
        plt.plot(angle_cleaned,predicted_value,'b--',label='fitting',marker='.')
        plt.title(str(self.temeprature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-Fitted with cod='+str(self.cod))
        plt.xlabel('Angle')
        plt.ylabel('V2w')
        plt.legend()
        plt.savefig(fname=fig_folder+'\\'+str(self.temperature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-Fit.jpg')
    #    print(len(angle_cleaned))
    #    plt.savefig()
        plt.close()


class signal_merge(object):
    #merge different signals, eg. current dependence txt data at certain field
    def __init__(self, root):
        #判定属于哪种数据： curr dep , field dep or temp dep
        self.root =  root
       
        self.dependence = 'current'
        self.dependence_unit = 'uA'
        self.temperatrue = 300.0
        self.field = 0.00
        self.current = 500.0
#        self.temperature, self.field, self.current, self.gain
        self.folder_name, self.dictionary = self._parse_path(root) 
        self.variable=[]
        self.amplitude=[]
        self.amplitude_std=[]
        
        aquired_para = 0
        self.para_constant =''
        for k,v in self.dictionary.items():
            if(len(v)==0):
                self.dependence = k
#                print(self.dependence)
                switch = {'temperature': 'K', 'current':'uA', 'field':'T'}
                self.dependence_unit = switch[k]
            else:
                if(k=='temperature'):
                    self.temperatrue = float(v[0].strip('K'))
                    self.para_constant = self.para_constant + v[0] +' '
                    aquired_para = aquired_para+1
                if(k=='field'):
                    self.field = float(v[0].strip('T'))
                    self.para_constant = self.para_constant + v[0] +' '
                    aquired_para = aquired_para+1
                if(k=='current'):
                    self.current = float(v[0].strip('uA'))
                    self.para_constant = self.para_constant + v[0]+' ' 
                    aquired_para = aquired_para+1
        
        if(aquired_para == 3):
            print('This is a single run file, not dependence series!')
        if(aquired_para == 2):
            print(self.dependence+'-dep '+self.para_constant+'data processing begins.')
        if(aquired_para == 1):
            print('Is this a 2-variable data folder?Only one variable dependence series can be processed!')
        if(aquired_para == 0):
            print('Not enough info in folder name')
                

    def _parse_path(self, path):
        #parse the filename, extract [Temp, Curr, Field, gain]
        
        #extract the last string of the root path
        
        project_string = self.root.split('\\')[-1]
        temperature_pattern = re.findall(r'\d+K', project_string)
        field_pattern = re.findall(r'\d+T', project_string)
        current_pattern = re.findall(r'\d+uA', project_string)
        
        return project_string, {'temperature':temperature_pattern, 'field':field_pattern, 'current':current_pattern}
    
    def _find_txt(self):
        #find all txt data file under certain path
        
        pure_txt_name=[]
        txt_full_dir=[]
        g = os.walk(self.root) 
        for par_dir, _, files in g: 
            for file in files: 
                if file.endswith(".txt"):
                    name=file.strip('.txt')
        #            print(name)
                    pure_txt_name.append(name)
                    filepath = os.path.join(par_dir, file) 
        #            print(filepath)
                    txt_full_dir.append(filepath)
    
    #    
    #    print(pure_log_name)
    #    print(txt_full_dir)  
        return pure_txt_name, txt_full_dir 
    

    def _normalize(self, signal):
        #normalize by the maximum and minimum
        normalizer=max(signal)-min(signal)    
        return [(s-min(signal))/normalizer for s in signal] 
    
    
    def _process(self):
        '''
        #core method
        '''
        txt_name_list, txt_dirs = self._find_txt()
        for name, dire in zip(txt_name_list, txt_dirs):
            data_model = signal_processing(self.root, name)
            temperature, current, field, gain = data_model._parse(name)
            if len(txt_name_list)<5:
                print("Not enough data point in ["+self.folder_name+"]")
            data_model._process()
            flag, amplitude, a_std = data_model._get_info()
            data_model._noise()
            #change var according to dependence
            switch = {'temperature': temperature, 'current':current, 'field':field}
            var = switch[self.dependence]
            
            if flag:              #delete the noisy point
                self.variable.append(var)
                self.amplitude.append(amplitude)
                self.amplitude_std.append(a_std)
                
#            self.variable.append(var)
#            self.amplitude.append(amplitude)
#            self.amplitude_std.append(a_std)
#        
        #sort by current
        
        zipped = sorted(zip(self.variable,  self.amplitude, self.amplitude_std))
        self.variable=[k[0] for k in zipped]
        self.amplitude=[k[1] for k in zipped]
        self.amplitude_std=[k[2] for k in zipped]
        print(self.dependence+'-dep '+self.para_constant+'data processing ends.')

    
    def _plot(self):
        #create destination folder for figs
        res_folder = self.root+'\\result'
        if not os.path.exists(res_folder):
            os.mkdir(res_folder)
        
        plt.figure()
        
        plt.errorbar(self.variable, self.amplitude, yerr=self.amplitude_std, fmt='-o', barsabove=True)
        plt.title(self.dependence+" dependence @"+self.para_constant)
        plt.xlabel(self.dependence+'('+self.dependence_unit+')')
        plt.ylabel('V2w(V)')
        plt.savefig(fname=res_folder+'\\'+self.dependence+' dependence at '+self.para_constant+'.jpg')
        plt.show() 
        plt.close()

class folder_merge(object):
    def __init__(self, root, meta_var = 'field', meta_constant = 'temperature', assigned_constant = '2K', IsReverse = False):
        self.root = root
        self.metaVar = meta_var
        self.meta_constant = meta_constant
        self.assigned_constant = assigned_constant
        self.IsReverse = IsReverse
        self.folder_name_list = []
        self.folder_path_list = []
        self._find_folder()
        self.selected_folder = self._select_folder()
        self.model_list = []
        
        
        
    def _select_folder(self):
        #return data folder for certain dependence; field, current or temperature
        switch = {'temperature': 'K', 'current':'uA', 'field':'T'}
        
        selected_folders = []
        meta_var_list = []
        
        for folder_name in self.folder_name_list:
            dictionary = self._parse_name(folder_name)
            main_var = dictionary[self.metaVar]
            if main_var=='':
                print('FOLDER: ['+folder_name+'] does not have '+self.metaVar+' variable')
                
            main_constant = dictionary[self.meta_constant]
            if main_constant=='':
                print('FOLDER: ['+folder_name+'] does not have '+self.meta_constant+' variable')
            if((not main_constant=='') and (not main_var=='')):
                if main_constant == self.assigned_constant:
                    
                    meta_var_list.append(int(main_var.strip(switch[self.metaVar])))
                    selected_folders.append(folder_name)
        #sort according to main var
        a=self.IsReverse
        zipped = sorted(zip(meta_var_list, selected_folders), reverse = a)
        selected_folders=[k[1] for k in zipped]
        
        
        print("Selected Folder: ")
        for folder in selected_folders:
            print(folder)
            
        
        return selected_folders
    
    def _parse_name(self, folder_name):
        #parse the filename, extract [Temp, Curr, Field, gain]
        
        
        
        temperature_pattern = re.findall(r'\d+K', folder_name)
        temperature_string = temperature_pattern[0] if len(temperature_pattern)>0 else ''
        field_pattern = re.findall(r'\d+T', folder_name)
        field_string = field_pattern[0] if len(field_pattern)>0 else ''
        current_pattern = re.findall(r'\d+uA', folder_name)
        current_string = current_pattern[0] if len(current_pattern)>0 else ''
        return {'temperature':temperature_string, 'field':field_string, 'current':current_string}
    
    def _find_folder(self):
        #find all folder name under root: depth==1
            
        
        g = os.walk(self.root) 
        for par_dir, dirnames, files in g: 
            for folder in dirnames:
                # only 1st-order folder is needed under root
                folder_full_path = os.path.join(self.root,folder)
                if os.path.exists(folder_full_path):
#                    print(folder)
                    self.folder_name_list.append(folder)
                    self.folder_path_list.append(folder_full_path)
    
    def _process(self, normalize = False):
        
        for folder in self.selected_folder:
            full_path = os.path.join(self.root,folder)
            model = signal_merge(full_path)
            model._process()
            self.model_list.append(model)
            model._plot()
       
    
    def _plot(self):
        
        res_folder = self.root+'\\result'
        if not os.path.exists(res_folder):
            os.mkdir(res_folder)
        
        plt.figure()
        for model in self.model_list:
            label = model.dictionary[self.metaVar][0]
            #Current square
            CurrentSquarre = [(i/1000000)**2 for i in model.variable]
            plt.errorbar(CurrentSquarre, model.amplitude, yerr=model.amplitude_std, fmt='-o', barsabove=True, label = label)
        
        plt.title('Current dependence at '+self.assigned_constant+' for different '+self.metaVar)
        plt.xlabel('I^2(A)')
        plt.ylabel('V2w(V)')
        plt.legend()
        plt.savefig(fname=res_folder+'\\Current dependence at '+self.assigned_constant+' for different '+self.metaVar+'.jpg')
        plt.show() 
        plt.close()
        
                

 
if __name__ == "__main__":
      
      root="F:\\BaiduNetdiskDownload\\Measurement\\from XWY\\20190415-20190319Cr2O3(11-20) 20nm-S23 10um"
      folder_root="F:\\BaiduNetdiskDownload\\Measurement\\from XWY\\20190415-20190319Cr2O3(11-20) 20nm-S23 10um\\4-Pin3-4-5-6 9T 2K Curr"
      filename="1 S23 2.00K90000Oe (7Hz) 0_360deg  700uA1000"
      '''
      signal = signal_processing(folder_root, filename)
      signal._process()
      signal._noise()
      
      signal_pack=signal_merge(root)
      signal_pack._process()
      
      '''
      
      #test for folder_marge class
      
      model = folder_merge(root, meta_var = 'temperature', meta_constant = 'field', assigned_constant = '9T')
      model._process()
      model._plot()
     
            


'''
nineTesla="F:\\BaiduNetdiskDownload\\Measurement\\from XWY\\20190415-20190319Cr2O3(11-20) 20nm-S23 10um\\4-Pin3-4-5-6 9T 2K Curr"
model = signal_merge(nineTesla)
model._process()
model._plot()
'''