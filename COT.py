import collections
import math
import time
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

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
        
        if not self.silent:
            print(f"COT: subtype means generated.")
    
    def generate_cos_values(self):
        df_mean_cos = self.df_mean.apply(lambda x: x / np.linalg.norm(x), axis=1)
        
        self.df_cos = pd.DataFrame(index=df_mean_cos.index)
        self.df_cos['cos'] = df_mean_cos.apply(lambda x: np.max(x), axis=1)
        self.df_cos['subtype'] = df_mean_cos.apply(lambda x: np.argmax(x), axis=1)
        
        if not self.silent:
            print(f"COT: cos values generated.")
    
    def estimate_p_values(self):
        print("COT: estimating p-values ...")
        count = 0
        dataFit = self.df_cos.sort_values(by='cos', ascending=False)
        dataFit = np.array(dataFit['cos'])
        pThreL = [0, 0.001, 0.005, 0.01, 0.05]
        numExp = np.zeros([len(pThreL) - 1])
        numObs = np.ones([len(pThreL) - 1])
        dataNum = len(dataFit)
        pvalFit = np.zeros(dataNum)
        subNum = len(self.subtypes)
        
        def mix_norm_cdf(x, weights, means, covars, lower, upper):
            mcdf = np.zeros([1, len(x)])
            mnor = np.zeros([1, len(x)])
            for i in range(len(weights)):
                mcdf += weights[i] * (norm.cdf(x, loc=means[i], scale=covars[i]**(0.5)) - 
                                      norm.cdf(lower, loc=means[i], scale=covars[i]**(0.5))) 
                mnor += weights[i] * (norm.cdf(upper, loc=means[i], scale=covars[i]**(0.5)) - 
                                      norm.cdf(lower, loc=means[i], scale=covars[i]**(0.5)))
            mcdf /= mnor
            return mcdf.reshape(-1)
        
        while ((sum(pvalFit < pThreL[-1]) - dataNum * pThreL[-1]) / sum(pvalFit < pThreL[-1])) > 0.001:
            tic = time.perf_counter()
            count += 1
            cluster = AgglomerativeClustering(n_clusters=subNum, linkage='ward')
            cluster.fit_predict(dataFit.reshape(-1, 1))
            clusDic = collections.defaultdict(list)
            for i in range(len(cluster.labels_)):
                clusDic[cluster.labels_[i]].append(dataFit[i])
            
            gmW = np.zeros(subNum)
            gmM = np.zeros([subNum, 1])
            gmP = np.zeros([subNum, 1, 1])
            for i, clus in enumerate(clusDic):
                gmW[i] = len(clusDic[clus]) / len(dataFit)
                gmM[i] = sum(clusDic[clus]) / len(clusDic[clus])
                gmP[i] = 1 / np.array(clusDic[clus]).var()
            
            gm = GaussianMixture(n_components=subNum, tol=1e-5, max_iter=10000,
                                 weights_init=gmW, means_init=gmM, precisions_init=gmP).fit(dataFit.reshape(-1, 1))
            pvalFit = 1 - mix_norm_cdf(dataFit, gm.weights_, gm.means_, gm.covariances_, 1/(subNum**0.5), 1)
            
            dataFitNew = []   
            for i in range(1, len(pThreL)):       
                numExp[i - 1] = int(dataNum * (pThreL[i] - pThreL[i-1]))
                numObs[i - 1] = sum(pvalFit < pThreL[i]) - sum(pvalFit < pThreL[i-1])       
                if int(dataNum * (pThreL[i] - pThreL[i-1])) < (sum(pvalFit < pThreL[i]) - sum(pvalFit < pThreL[i-1])):          
                    picked = (np.linspace(sum(pvalFit < pThreL[i-1]), sum(pvalFit < pThreL[i]), int(dataNum * (pThreL[i] - pThreL[i-1])), 
                                      endpoint=False)).astype(int)
                    dataFitNew = np.hstack([dataFitNew, dataFit[picked]])
                else:
                    dataFitNew = np.hstack([dataFitNew, dataFit[sum(pvalFit < pThreL[i-1]):sum(pvalFit < pThreL[i])]])
            dataFit = np.hstack([dataFitNew, dataFit[sum(pvalFit < pThreL[-1]):]])    
            dataNum = len(dataFit)
            toc = time.perf_counter()
            print(f"Iteration {count}: {toc - tic:0.4f} seconds")
        
        self.df_cos['p.value'] = 1 - mix_norm_cdf(self.df_cos['cos'], gm.weights_, gm.means_, gm.covariances_, 1/(subNum**0.5), 1)
        self.df_cos['q.value'] = multipletests(self.df_cos['p.value'], method='fdr_bh')[1]
        print("COT: p-values estimated.")
    
    def obtain_subtype_markers(self, pThre=None, qThre=0.05, top=None, per=None):
        self.markers = {i: [] for i in self.subtypes.keys()}
        cos_sorted = self.df_cos.sort_values(by='cos', ascending=False)  
        
        for i in range(top if top else len(cos_sorted)):
            self.markers[cos_sorted.iloc[i, 1]].append(cos_sorted.index[i])
        
        if per is not None:
            for subtype in self.markers:
                if len(self.markers[subtype]) > per:
                    self.markers[subtype] = self.markers[subtype][:per]
        
        if pThre is not None:
            for subtype in self.markers:
                self.markers[subtype] =\
                cos_sorted.loc[self.markers[subtype]].index[cos_sorted.loc[self.markers[subtype], 'p.value'] <= pThre]
        
        if qThre is not None:
            for subtype in self.markers:
                self.markers[subtype] =\
                cos_sorted.loc[self.markers[subtype]].index[cos_sorted.loc[self.markers[subtype], 'q.value'] <= qThre]
    
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
        mg_col = ['red', 'blue', 'orange', 'green', 'purple', 'pink', 'cyan', 'brown', 'lime', 'magenta']
        
        handles = [None] * K
        for i, cell in enumerate(self.markers):
            handles[i] = plt.scatter(Xp.loc[0, self.markers[cell]], Xp.loc[1, self.markers[cell]], marker = 'o', s=10, 
                                     color=mg_col[i], facecolors='none')
        plt.scatter(PS[0], PS[1], marker='^', color='k', s=15)
        
        legText = [''] * K
        for i, cell in enumerate(self.subtypes):
            legText[i] = f"{cell} ({len(self.markers[cell])})"
        
        leg = plt.legend(handles, legText, prop={'size': 8})
        for i, text in enumerate(leg.get_texts()):
            text.set_color(mg_col[i])
    
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
