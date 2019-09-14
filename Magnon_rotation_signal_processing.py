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
        self.temperature, self.current, self.field, self.gain = self._parse(filename)
        self.angle, self.V_2omega_W, self.R_2omega_W, self.V_2omega_w, self.R_2omega_w = self._get_data()
        
        self.IsGoodFitting_w=True
        self.amplitude_w = 0
        self.amplitude_std_w = 0
        self.cod_w = 0
        self.angle_cleaned_w = []
        self.cleaned_shifted_signal_w = []
        self.fitted_signal_w = []
        
        self.IsGoodFitting_W=True
        self.amplitude_W = 0
        self.amplitude_std_W = 0
        self.cod_W = 0
        self.angle_cleaned_W = []
        self.cleaned_shifted_signal_W = []
        self.fitted_signal_W = []
    
    
    def _get_info(self):
        return self.IsGoodFitting_W, self.amplitude_W, self.amplitude_std_W, \
        self.IsGoodFitting_w, self.amplitude_w, self.amplitude_std_w
    
    def _parse(self, filename):
        #parse the filename, extract [Temp, Curr, Field, gain]
        
        number_list=re.findall(r"\d+\.?\d*",filename)
        temperature_sting=re.findall(r"\d+\.?\d*K",filename)[0]
        temperature=float(re.findall(r"\d+\.?\d*",temperature_sting)[0])
#        print(temperature)
    #    print(type(temperature))
        
        current_sting=re.findall(r"\d+\.?\d*uA",filename)[0]
        current=int(re.findall(r"\d+\.?\d*",current_sting)[0])
        field_sting=re.findall(r"\d+\.?\d*Oe",filename)[0]
        field=float(re.findall(r"\d+\.?\d*",field_sting)[0])
        gain=int(number_list[-1])
        
        return temperature, current, field, gain
    

    def _get_data(self):
        #get 2w signal from 2 SR830 simoutaneously
        angle=[]
        V_2omega_W=[]
        R_2omega_W=[]
        V_2omega_w=[]
        R_2omega_w=[]
        
        file_to_read=open(self.file_full_path , 'r') 
        
        lines = file_to_read.readlines() # 整行读取数据
        
        for line in lines:
            angle_temp,a, b, c,d,e,f ,R_omega_temp, V_omega_temp, R_2omega_temp, V_2omega_temp  = [float(i) for i in line.split()] 
            angle.append(angle_temp)
            V_2omega_W.append(V_omega_temp)
            R_2omega_W.append(R_omega_temp)
            V_2omega_w.append(V_2omega_temp)
            R_2omega_w.append(R_2omega_temp)
        
        return angle, V_2omega_W, R_2omega_W, V_2omega_w, R_2omega_w
    
    
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
                cleaned_signal.append(signal-signal_bar)    #normalize to 0
            else:
                print('angle, signal: '+str(angle)+' , '+str(signal)+' is omitted!')
#        print(cleaned_angle)
#        print("raw data number: "+str(len(angle_list)))        
#        print("cleaned data number: "+str(len(cleaned_angle)))
        return cleaned_angle, cleaned_signal
    
    
    def _normalize(self, signal):
        #minus the mean
        return [s-np.mean(signal) for s in signal]
        
        
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
    
    
    def _process_w(self, threshold_w = 0.7):
        '''SR 830, measure shorter Pt'''
        #90 degree shift of V2w
        V_2omega_w_shifted=self._shift(self.V_2omega_w)
        #out of range data deletion
        self.angle_cleaned_w, self.cleaned_shifted_signal_w = self._abnormal_data_detect_and_deletion(V_2omega_w_shifted)

        
        #sine curve fit
        self.amplitude_w, self.amplitude_std_w , self.fitted_signal_w = self._sine_curve_fit(self.angle_cleaned_w, self.cleaned_shifted_signal_w) 
        
        
        cod_w = self._get_cod(self.cleaned_shifted_signal_w, self.fitted_signal_w)
        self.cod_w = cod_w
        #set amplitude to 0 if cod<0.5
        if cod_w<threshold_w:
            self.IsGoodFitting_w=False
#            self.amplitude_w=0
        
        print("Is it good fitting with cod_w>"+str(threshold_w)+"?\n"+str(self.IsGoodFitting_w)+"\ncod_w = "+str(self.cod_w))
            
        return self.angle_cleaned_w, self.cleaned_shifted_signal_w, self.fitted_signal_w
        
    
    def _process_W(self, threshold_W = 0.7):
        '''SR 830, measure longer Pt'''
        #90 degree shift of V2w
        V_2omega_W_shifted=self._shift(self.V_2omega_W)
        #out of range data deletion
        self.angle_cleaned_W, self.cleaned_shifted_signal_W = self._abnormal_data_detect_and_deletion(V_2omega_W_shifted)

        
        #sine curve fit
        self.amplitude_W, self.amplitude_std_W , self.fitted_signal_W = self._sine_curve_fit(self.angle_cleaned_W, self.cleaned_shifted_signal_W) 
        
        
        cod_W = self._get_cod(self.cleaned_shifted_signal_W, self.fitted_signal_W)
        self.cod_W = cod_W
        #set amplitude to 0 if cod<0.5
        if cod_W<threshold_W:
            self.IsGoodFitting_W=False
#            self.amplitude_W=0
        
        print("Is it good fitting with cod_W>"+str(threshold_W)+"?\n"+str(self.IsGoodFitting_W)+"\ncod_W = "+str(self.cod_W))
            
        return self.angle_cleaned_W, self.cleaned_shifted_signal_W, self.fitted_signal_W
    
    def _process(self):
        '''define threshold'''
        print('****************************************************************************')
        print('plot raw data')
        self._plot_raw(fig_folder = self.root)
        print('Anaylization begins.')
        print('Long Pt')
        self._process_w()
        print('----------------------------------------------------------------------')
        print('short Pt')
        self._process_W()
        print('----------------------------------------------------------------------')
        self._plot_shifted(fig_folder = self.root)
        print('----------------------------------------------------------------------')
        print('plot fitted data')
        self._plot_fitted(fig_folder = self.root)
        
        print('Amplitude for long Pt: ', self.amplitude_W, 'std: ',self.amplitude_std_W)
        print('Amplitude for short Pt: ', self.amplitude_w, 'std: ', self.amplitude_std_w)
        if self.amplitude_W>0 and self.amplitude_w:
            ratio = self.amplitude_W/self.amplitude_w
            print('V2w(W)/V2w(w) =', ratio)
        
    def _plot_raw(self, fig_folder,savefig = False):
        #plot raw data
        plot_name = str(self.temperature)+'K '+str(self.field/10000)+'T '+str(self.current)+'uA gain'+str(self.gain)+'-RAW'
        plt.figure()
        plt.plot(self.angle , self.V_2omega_W, 'r',label='long Pt-raw',marker='.')
        plt.plot(self.angle , self.V_2omega_w, 'b',label='short Pt-raw',marker='.')
        plt.title(plot_name)
        plt.xlabel('Angle')
        plt.ylabel('V2w/V')
        plt.legend()
        if savefig:
            plt.savefig(fname=fig_folder+'\\'+plot_name+'.jpg')
            plt.close()
        else:
            plt.show()
        
    def _plot_shifted(self, fig_folder, savefig = False):
        plot_name = str(self.temperature)+'K '+str(self.field/10000)+'T '+str(self.current)+'uA gain'+str(self.gain)+'-Shifted'
        plt.figure()
        plt.plot(self.angle , self._normalize(self._shift(self.V_2omega_W)), 'r',label='long Pt-shifted',marker='.')
        plt.plot(self.angle , self._normalize(self._shift(self.V_2omega_w)), 'b',label='short Pt-shifted',marker='.')
        plt.title(plot_name)
        plt.xlabel('Angle')
        plt.ylabel('V2w/V')
        plt.legend()
        
        if savefig:
            plt.savefig(fname=fig_folder+'\\'+plot_name+'.jpg')
            plt.close()
        else:
            plt.show()
        
    def _plot_fitted(self, fig_folder, savefig = False):
        #angle_cleaned, V_2omega_shifted_cleaned, predicted_value = self._process()
        #data washed
        plot_name = str(self.temperature)+'K '+str(self.field/10000)+'T '+str(self.current)+'uA gain'+str(self.gain)
        plt.figure()
        plt.plot(self.angle_cleaned_W,self.cleaned_shifted_signal_W,'r',label='long Pt-cleaned',marker='.')
        plt.plot(self.angle_cleaned_W,self.fitted_signal_W,'cyan',label='long Pt-fitting',marker='.')
        
        plt.plot(self.angle_cleaned_w,self.cleaned_shifted_signal_w,'b',label='short Pt-cleaned',marker='.')
        plt.plot(self.angle_cleaned_w,self.fitted_signal_w,'cyan',label='short Pt-fitting',marker='.')
        
        plt.title(plot_name+'\nFitted with cod_W='+str(format(self.cod_W, '.2f')) +' and cod_w=' + str(format(self.cod_w, '.2f')) )
        plt.xlabel('Angle')
        plt.ylabel('V2w/V')
        plt.legend()
        if savefig:
            plt.savefig(fname=fig_folder+'\\'+str(self.temperature)+'K '+str(self.field)+'Oe '+str(self.current)+'uA-Fit.jpg')
    #    print(len(angle_cleaned))
    #    plt.savefig()
            plt.close()
        else:
            plt.show()


class signal_merge(object):
    #merge different signals, eg. current dependence txt data at certain field
    def __init__(self, root):
        #判定属于哪种数据： curr dep , field dep or temp dep
        self.root =  root
        
        #create destination folder for figs and output txts
        res_folder = self.root+'\\result'
        if not os.path.exists(res_folder):
            os.mkdir(res_folder)
        self.res_folder = res_folder
       
        self.dependence = 'current'
        self.dependence_unit = 'uA'
        self.temperatrue = 300.0
        self.field = 0.00
        self.current = 500.0
#        self.temperature, self.field, self.current, self.gain
        self.folder_name, self.dictionary = self._parse_path(root) 
        self.variable_W=[]    # can be list of T, H or I: long Pt
        self.variable_w = []     #short Pt
        self.amplitude_W=[]
        self.amplitude_std_W=[]
        self.amplitude_w=[]
        self.amplitude_std_w=[]
        
        aquired_para = 0
        self.para_constant =''
        for k,v in self.dictionary.items():
            if(len(v)==0):
                self.dependence = k
#                print(self.dependence)
                switch = {'temperature': 'K', 'current':'uA', 'field':'Oe'}
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
        for parent_dir, _, files in g: 
            for file in files: 
                if os.path.exists(os.path.join(self.root, file)) and file.endswith(".txt"):   #ensure ther is a txt under root path
                    name=file.strip('.txt')
#                    print(name)
                    pure_txt_name.append(name)
                    filepath = os.path.join(self.root, file) 
#                    print(filepath)
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
        matrix = []      # matrix for txt write, each row is a slice of data
        
        txt_name_list, txt_dirs = self._find_txt()
        for name, dire in zip(txt_name_list, txt_dirs):
            data_model = signal_processing(self.root, name)
            temperature, current, field, gain = data_model._parse(name)
            if len(txt_name_list)<5:
                print("Not enough data point in ["+self.folder_name+"]")
            data_model._process()
            flag_W, amplitude_W, amplitude_W_std, flag_w, amplitude_w, amplitude_w_std = data_model._get_info()
            
            #data_model._noise()
            #change var according to dependence
            switch = {'temperature': temperature, 'current':current, 'field':field}
            var = switch[self.dependence]
            
            matrix.append([var, amplitude_W, amplitude_W_std, amplitude_w, amplitude_w_std, amplitude_W/amplitude_w])
            
            
            if flag_W:
                self.variable_W.append(var)
                self.amplitude_W.append(amplitude_W)
                self.amplitude_std_W.append(amplitude_W_std)
            else:
                print(var,self.dependence_unit,amplitude_W, ' is dropped.')
            if flag_w:
                self.variable_w.append(var)
                self.amplitude_w.append(amplitude_w)
                self.amplitude_std_w.append(amplitude_w_std)
            else:
                print(var,self.dependence_unit, amplitude_w, ' is dropped.')
          
            
        #sort data by var and write to txt
        #print(matrix)
        name = self.dependence+' dependence at '+self.para_constant
        self.write2txt(matrix, name)
#            self.variable.append(var)
#            self.amplitude.append(amplitude)
#            self.amplitude_std.append(a_std)
#         
        #sort by current, or field, or temp
        
        zipped_W = sorted(zip(self.variable_W,  self.amplitude_W, self.amplitude_std_W))
        self.variable_W=[k[0] for k in zipped_W]
        self.amplitude_W=[k[1] for k in zipped_W]
        self.amplitude_std_W=[k[2] for k in zipped_W]
        
        zipped_w = sorted(zip(self.variable_w,  self.amplitude_w, self.amplitude_std_w))
        self.variable_w=[k[0] for k in zipped_w]
        self.amplitude_w=[k[1] for k in zipped_w]
        self.amplitude_std_w=[k[2] for k in zipped_w]
        
        print('\n\n*********************************************************************************\n')
        print(self.dependence+'-dep '+self.para_constant+'data processing ends.')
        
    def write2txt(self, matrix, name):
        #np.savetxt(os.path.join(self.res_folder, name+'.txt'), matrix);
        matrix = np.array(matrix)
        i = np.argsort(matrix[:,0])
        data = matrix[i]
        data = data.tolist()
        
        
        #write into txt
        f = open(os.path.join(self.res_folder, name+'.txt'), 'w')
        #write to txt
        for u in data:
            for v in u:
                f.write(str(v)+' ')
            f.write('\n')
        f.close()
        
        
    
    def _plot(self):
        
        #write to txt
        
        
        
        plt.figure()
        
        plt.errorbar(self.variable_W, self.amplitude_W, yerr=self.amplitude_std_W, fmt='-o', barsabove=True)
        plt.errorbar(self.variable_w, self.amplitude_w, yerr=self.amplitude_std_w, fmt='-o', barsabove=True)
        plt.title(self.dependence+" dependence @"+self.para_constant)
        plt.xlabel(self.dependence+'('+self.dependence_unit+')')
        plt.ylabel('V2w(V)')
        plt.savefig(fname=self.res_folder+'\\'+self.dependence+' dependence at '+self.para_constant+'.jpg')
        plt.show() 
        plt.close()
        
    def _plot_ratio(self):
        #compute the ratio of signal from long Pt and short Pt
        
        dict_W = dict(zip(self.variable_W, self.amplitude_W))
        dict_w = dict(zip(self.variable_w, self.amplitude_w))
        Var_commom = list(set(dict_W.keys()).intersection(set(dict_w.keys())))
        ratio = []
        for var in Var_commom:
            ratio.append(dict_W[var]/dict_w[var])
        zip_ratio = sorted(zip(Var_commom, ratio))
        Var_commom=[k[0] for k in zip_ratio]
        ratio=[k[1] for k in zip_ratio]
        
        plt.figure()
        plt.plot(Var_commom, ratio, 'b',label='V2w(long Pt)/v2w(short Pt)',marker='.')
        plt.title(self.dependence+" dependence @"+self.para_constant)
        plt.xlabel(self.dependence+'('+self.dependence_unit+')')
        plt.ylabel('ratio')
        res_folder = self.root+'\\result'
        plt.savefig(fname=self.res_folder+'\\'+self.dependence+' dependence at '+self.para_constant+'-ratio.jpg')
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
      
      root=r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Rotator Thermometer"
      folder_root_1 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Rotator Thermometer\1-5K 9T 800uA"
      folder_root_2 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Rotator Thermometer\2-5K 800uA gain 100-field dep"
      folder_root_3 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\8-2K 800uA gain 100-field dep"
      folder_root_4 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\7-2K 9T gain 100-curr dep"
      folder_root_5 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\4-2K 800uA gain 1000 field dep"
      folder_root_6 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\9-9T 800uA gain 100-temp dep"
      folder_root_7 = r"D:\Data\20190909to0911 B264\Cr2O3 11-20 20190608 18nm S12-100um 5um\2-2K 500uA gain 100-field dep"
      folder_root_8 = r"D:\Data\20190909to0911 B264\Cr2O3 11-20 20190608 18nm S12-100um 5um\3-2K 9T gain 100-curr dep\data"
      folder_root_9 = r"D:\Data\20190909to0911 B264\Cr2O3 11-20 20190608 18nm S12-100um 5um\4-2K 800uA gain 100-field dep"
      folder_root_10 = r"D:\Data\20190909to0911 B264\Cr2O3 11-20 20190608 18nm S12-100um 5um\5-9T 500uA gain 100-temp dep"
      folder_root_11  = r"D:\Data\20190909to0911 B264\Cr2O3 11-20 20190608 18nm S12-100um 5um\6-9T 800uA gain 100-temp dep"
      folder_root_12 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\10-2K 500uA gain 100-field dep"
      folder_root_13 = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Sys Thermometer\11-9T 500uA gain 100-temp dep"
      
      fig_folder = r"D:\Data\20190903to0905 B264\Cr2O3 11-20 20190507 18nm, Rotator Thermometer\1-5K 9T 800uA\test"
      filename_0 = "3 20190507 Cr2O3 11-20 18nm- S1 200um 4um 4.98K90000Oe (AC) 0_360deg  800uA10"
      filename_1 = "1 20190507 Cr2O3 11-20 18nm- S1 200um 4um 4.99K90000Oe (AC) 0_360deg  800uA100"
      filename_2 = "2 20190507 Cr2O3 11-20 18nm- S1 200um 4um 4.98K90000Oe (AC) 0_360deg  800uA1000"
      '''
      signal = signal_processing(folder_root, filename_0)
      
      signal._process(fig_folder)
      '''
      signal_m = signal_merge(folder_root_13)
      signal_m._process()
      signal_m._plot()
      signal_m._plot_ratio()


      
      
    
            

