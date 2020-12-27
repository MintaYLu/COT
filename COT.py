import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

class COT:
    def __init__(self, df_raw=None, df_mean=None, silent=False):
        self.df_raw = df_raw
        self.df_mean = df_mean
        self.df_cos = None
        self.silent = silent
        
        if self.df_raw is None and self.df_mean is None:
            raise ValueError("df_raw and df_mean cannot be None at the same time.")
        if not self.silent:
            print(f"COT: package initiated.")
    
    # def load_data(self, filename, logarithmic_data=False):
    #   self.df_raw = pd.read_csv(filename, index_col=0)
    #    if logarithmic_data:
    #        self.df_raw = self.df_raw.apply(lambda x: 2 ** x)
    #    
    #    if not self.silent:
    #        print(f"COT: data imported from \'{filename}\'.")
    
    def generate_subtype_means(self):
        map_subtypes = {}
        for col in self.df_raw.columns:
            map_subtypes.setdefault(col.split('.')[0], []).append(col)
        
        self.df_mean = pd.DataFrame(index=self.df_raw.index)
        for subtype, cols in map_subtypes.items():
            self.df_mean[subtype] = self.df_raw[cols].mean(axis=1)
        
        if not self.silent:
            print(f"COT: subtype means generated.")

    def generate_cos_values(self):
        df_mean_cos = self.df_mean.apply(lambda x: x / np.linalg.norm(x), axis=1)
        
        self.df_cos = pd.DataFrame(index=df_mean_cos.index)
        self.df_cos['cos'] = df_mean_cos.apply(lambda x: np.max(x), axis=1)
        self.df_cos['subtype'] = df_mean_cos.apply(lambda x: np.argmax(x), axis=1)
        
        if not self.silent:
            print(f"COT: cos values generated.")
    
    def save_cos_values(self, filename, threshold=None, sorted=True):
        df_output = self.df_cos.copy()
        if threshold:
            df_output = df_output[df_output['cos'] >= threshold]
            filename = '.'.join(filename.split('.')[:-1]) + f"_threshold={threshold}.{filename.split('.')[-1]}"
        if sorted:
            df_output = df_output.sort_values(by='cos', ascending=False)
        
        df_output.to_csv(filename, index=True)
        
        if not self.silent:
            print(f"COT: cos values exported to \'{filename}\'.")
    
    def cos_pipeline(self, input_file, output_file, logarithmic_input=False, sorted_output=True, output_threshold=None):
        self.load_data(filename=input_file, logarithmic_data=logarithmic_input)
        self.generate_subtype_means()
        self.generate_cos_values()
        
        self.save_cos_values(filename=output_file, threshold=None, sorted=sorted_output)
        if output_threshold:
            self.save_cos_values(filename=output_file, threshold=output_threshold, sorted=sorted_output)
