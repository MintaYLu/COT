# COT Readme
**COT**

**License**

Copyright <2020> 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

**Code**

    import numpy as np
    import pandas as pd
    import warnings
    warnings.filterwarnings("ignore")
    
    class COT:
        def __init__(self, silent=False):
            self.df_raw = None
            self.df_mean = None
            self.df_cos = None
            self.silent = silent
            
            if not self.silent:
                print(f"COT: package initiated.")
        
        def load_data(self, filename, logarithmic_data=False):
            self.df_raw = pd.read_csv(filename, index_col=0)
            if logarithmic_data:
                self.df_raw = self.df_raw.apply(lambda x: 2 ** x)
            
            if not self.silent:
                print(f"COT: data imported from \'{filename}\'.")
        
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

README



**COT** is a Python module for machine learning built on top of NumPy and Pandas, and is distributed under the MIT license. The project was developed and maintained by Virginia Tech CBIL Group.


## **Installation**

**Dependencies**
COT requires:

- Python (>= 3.6)
- NumPy (>= 1.13.3)
- Pandas (>= 0.25.0)

**Installation**
If you already have a working installation of numpy and pandas, the easiest way to install COT is downloading the COT.py file at: <Add link here>.

**Usage**
Once the COT.py is downloaded, the user need to either put the COT.py file into the working directory or add it to the Python module search path by:

    import sys
    sys.path.append("the-path-of-the-COT.py-file")

Then, the user should be able to use the COT module. See the following code for an example:

    from COT import COT
    
    cot = COT()
    cot.cos_pipeline(input_file, output_file, logarithmic_input=False, sorted_output=True, output_threshold=None)

The module will read the original data from the csv file (input_file), then generate the output csv file (output_file), which contains three columns: gene_id, cos_value, subtype_label.
If the output_threshold is set to some finite value between 0 and 1, then another csv file, which only stores the genes with a cos_value >= output_threshold, will be generated.

**Description**
Cosine based One-sample Test (COT) Python package is designed to detect marker genes among many subtypes using subtype-enriched expression profiles, which is pplicable to multi-omics data. 


## Functions
- COT(silent): COT class constructor.
    - Parameters:
        - silent (True/False): If “True”, turn off the standard output. Default value is “False”.
- load_data(filename, logarithmic_data): Import the raw gene data from the csv file. The csv file has the following format: each row represents a gene, and each column represents the expression level of a sample; and the first column represents the gene names, and the first row represents the subtypes of the corresponding samples.
    - Parameters:
        - filename (str): The filename of the gene data.
        - logarithmic_data (True/False): If “True”, then the input data in the csv file is logarithmic, and we perform the exponential transformation to restore the raw gene features during the data loading process. Default value is “False”.
    - Outputs:
        - Store the raw gene data at the class attribute, df_raw. In the dataframe, df_raw, the index is the gene name, and the column is the sample name, which is “subtype_name”, “subtype_name.1”, “subtype_name.2”, etc.
- generate_subtype_means(): Generate a dataframe to restore the subtype means, which are the mean values of the gene data within the same subtypes.
    - Parameters:
        - N/A
    - Ouputs:
        - Store the subtype mean data at the class attribute, df_mean, where the index is the gene name, and the column represents the subtype name.
- generate_cos_values(): Compute the cos value for each subtype, then find the maximum cos value and the corresponding subtype. 
    - Parameters:
        - N/A
    - Outputs:
        - Store the maximum cos value and the corresponding subtype for each gene at the class attribute, df_cos, where the index is the gene name, and the column “cos”, “subtype” represent the maximum cos value and the corresponding subtype, repectively.
- save_cos_values(filename, threshold=None, sorted=True): Save the gene names, maximum cos values and the corresponding subtypes into a csv file.
    - Parameters:
        - filename (str): The filename of the output data.
        - threshold (float, [0, 1]): If threshold is set to some finite value, then another csv file will be generated, which only stores the genes with a maximum cos value >= threshold. Default value is “None”.
        - sorted (True/False): If “true“, the gene will be sorted by the cos value in the descending order. Default value is “True”.
    - Outputs:
        - A csv file (filename) stores all the genes with their cos values and the corresponding subtypes.
        - Another csv file (filename + threshold information) when threshold is set to some finite value, which stores the data when cos value >= threshold.
- cos_pipeline(input_file, output_file, logarithmic_input, sorted_output, output_threshold): A pipeline to process the COT algorithm.
    - Parameters:
        - input_file (str): Same as filename in the class method load_data.
        - output_file (str): Same as filename in the class method save_cos_values.
        - logarithmic_input (True/False): Same as logarithmic_data in the class method load_data.
        - sorted_output (True/False): Same as sorted in the class method save_cos_values.
        - output_threshold (float, [0, 1]): Same as threshold in the class method save_cos_values.
    - Outputs:
        - csv file(s) which stores the input genes with their cos values and the corresponding subtypes. Same as the outputs of the class method save_cos_values.
## Example:


    from COT import COT
    
    cot = COT()
    cot.cos_pipeline("input_gene.csv", "output_cos.csv", output_threshold=0.9)

input_gene.csv:

| gene | A   | A   | B   | B   |
| ---- | --- | --- | --- | --- |
| 0    | 0.5 | 0.7 | 0.7 | 0.9 |
| 1    | 1.0 | 1.0 | 0.0 | 0.0 |

output_cos.csv:

| gene | cos | subtype |
| ---- | --- | ------- |
| 1    | 1.0 | A       |
| 0    | 0.8 | B       |

  output_cos_threshold=0.9.csv:

| gene | cos | subtype |
| ---- | --- | ------- |
| 1    | 1.0 | A       |


 
 

