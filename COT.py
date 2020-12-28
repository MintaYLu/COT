import collections
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
# import warnings
# warnings.filterwarnings("ignore")

class COT:
    def __init__(self, df_raw=None, df_mean=None, logarithmic_data=False, normalization=True, silent=False):
        self.df_raw = df_raw
        self.df_mean = df_mean
        self.subtypes = collections.defaultdict(list)
        self.df_cos = None
        self.markers = {}
        
        self.silent = silent
        
        
        if self.df_raw is None and self.df_mean is None:
            raise ValueError("df_raw and df_mean cannot be None at the same time.")
            
        if logarithmic_data:
            if self.df_raw is not None:
                self.df_raw = self.df_raw.apply(lambda x: 2 ** x)
            if self.df_mean is not None:
                self.df_mean = self.df_mean.apply(lambda x: 2 ** x)
                
        if normalization:
            if self.df_raw is not None:
                nor_mean = self.df_raw.sum(axis=0).mean()
                nor_digit = math.floor(math.log10(nor_mean))
                nor_const = math.floor(nor_mean / (10**nor_digit)) * (10**nor_digit)
                for sample in self.df_raw:
                    self.df_raw[sample] = self.df_raw[sample] / self.df_raw[sample].sum() * nor_const
                
            if self.df_mean is not None:
                nor_mean = self.df_mean.sum(axis=0).mean()
                nor_digit = math.floor(math.log10(nor_mean))
                nor_const = math.floor(nor_mean / (10**nor_digit)) * (10**nor_digit)       
                for sample in self.df_mean:
                    self.df_mean[sample] = self.df_mean[sample] / self.df_mean[sample].mean() * nor_const
        
        if not self.silent:
            print(f"COT: package initiated.")
    
    # def load_data(self, filename, logarithmic_data=False):
    #   self.df_raw = pd.read_csv(filename, index_col=0)
    #    if logarithmic_data:
    #        self.df_raw = self.df_raw.apply(lambda x: 2 ** x)
    #    
    #    if not self.silent:
    #        print(f"COT: data imported from \'{filename}\'.")
    
    def generate_subtype_means(self, subtype_label):
        # map_subtypes = {}
        # for col in self.df_raw.columns:
        #   map_subtypes.setdefault(col.split('.')[0], []).append(col)
        for sample in range(len(subtype_label)):
            self.subtypes[subtype_label[sample]].append(self.df_raw.columns[sample])
        
        self.df_mean = pd.DataFrame(columns=self.subtypes.keys())
        for subtype in self.subtypes:
            self.df_mean[subtype] = self.df_raw[self.subtypes[subtype]].mean(axis=1)
        
        self.markers = {i: [] for i in self.subtypes.keys()}
        
        if not self.silent:
            print(f"COT: subtype means generated.")

            
    def generate_cos_values(self):
        df_mean_cos = self.df_mean.apply(lambda x: x / np.linalg.norm(x), axis=1)
        
        self.df_cos = pd.DataFrame(index=df_mean_cos.index)
        self.df_cos['cos'] = df_mean_cos.apply(lambda x: np.max(x), axis=1)
        self.df_cos['subtype'] = df_mean_cos.apply(lambda x: df_mean_cos.columns[np.argmax(x)], axis=1)
        
        if not self.silent:
            print(f"COT: cos values generated.")
            
            
    def obtain_subtype_markers(self, pThre=None, top=None, per=None):
        if top is not None:
            cos_sorted = self.df_cos.sort_values(by='cos', ascending=False)
            for i in range(top):
                self.markers[cos_sorted.iloc[i, 1]].append(cos_sorted.index[i])
   

    def plot_simplex(self):
        X = self.df_mean
        Xproj = X.divide(X.sum(axis=1), axis=0)
        Xproj = Xproj.to_numpy().transpose()

        K = len(X.columns)
        A = np.identity(K)
        PS = np.vstack((np.cos(np.arange(K) * 2 * np.pi / K), np.sin(np.arange(K) * 2 * np.pi / K)))
        tmp = np.matmul(PS, np.linalg.pinv(A))
        Xp = pd.DataFrame(np.matmul(tmp, Xproj), columns=X.index)
        
        
        plt.axes().set_aspect('equal', 'datalim')
        plt.scatter(Xp.iloc[0, ], Xp.iloc[1, ], marker = 'o', s=10,
                    color='#%02x%02x%02x' % (200, 200, 200), facecolors='none', alpha=0.3)
        mg_col = ['red','blue','orange','green','purple']

        for i, cell in enumerate(self.markers):
            plt.scatter(Xp.loc[0, self.markers[cell]], Xp.loc[1, self.markers[cell]], marker = 'o', s=10,
                        color=mg_col[i], facecolors='none')

        plt.scatter(PS[0], PS[1], marker='^', color='k', s=15)
    
    # def save_cos_values(self, filename, threshold=None, sorted=True):
    #    df_output = self.df_cos.copy()
    #    if threshold:
    #        df_output = df_output[df_output['cos'] >= threshold]
    #        filename = '.'.join(filename.split('.')[:-1]) + f"_threshold={threshold}.{filename.split('.')[-1]}"
    #    if sorted:
    #        df_output = df_output.sort_values(by='cos', ascending=False)
    #    
    #    df_output.to_csv(filename, index=True)
    #    
    #    if not self.silent:
    #        print(f"COT: cos values exported to \'{filename}\'.")
    
    def cos_pipeline(self, input_file, output_file, logarithmic_input=False, sorted_output=True, output_threshold=None):
        self.load_data(filename=input_file, logarithmic_data=logarithmic_input)
        self.generate_subtype_means()
        self.generate_cos_values()
        
        self.save_cos_values(filename=output_file, threshold=None, sorted=sorted_output)
        if output_threshold:
            self.save_cos_values(filename=output_file, threshold=output_threshold, sorted=sorted_output)
