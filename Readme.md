# Readme
**COT**


## **Description**

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
- COT(df_raw, df_mean, logarithmic_data, normalization, silent): COT class constructor. During this step either the raw gene data or the subtype mean data will be imported and preprocessed, i.e., performing exponential transformation or sample-wise normalization.
    - Parameters:
        - df_raw (pandas.DataFrame): The raw gene data stored as a pandas Dataframe. df_raw has the following format: each row represents a gene, and each column represents a sample. Default value is None.
        - df_mean (pandas.DataFrame): The subtype mean data stored as a pandas Dataframe. df_mean has the following format: each row represents a gene, and each column represents a subtype. Default value is None.
        - logarithmic_data (True/False): If “True”, then the input data, df_raw and df_mean, is logarithmic, and we perform the exponential transformation to restore the raw gene features during the data loading process. Default value is “False”.
        - normalization (True/False): If “True”, then the feature-wise normalizayion will be applied to the input data, df_raw and df_mean. Default value is “True”.
        - silent (True/False): If “True”, turn off the standard output. Default value is “False”.
    - Outputs:
        - Store the raw gene data at the class attribute, df_raw.
        - Store the subtype mean data at the class attribute, df_mean.
- generate_subtype_means(subtype_label): Generate a dataframe to restore the subtype means, which are the mean values of the gene data within the same subtypes.
    - Parameters:
        - subtype_label (list(str)): The subtype label for each sample stored as a python list. Default value is None.
    - Outputs:
        - Store the subtype mean data at the class attribute, df_mean, where the index is the gene name, and the column represents the subtype name.
- generate_cos_values(): Compute the cos value for each subtype, then find the maximum cos value and the corresponding subtype.
    - Parameters:
        - N/A
    - Outputs:
        - Store the maximum cos value and the corresponding subtype for each gene at the class attribute, df_cos, where the index is the gene name, and the column “cos”, “subtype” represent the maximum cos value and the corresponding subtype, repectively.
- estimate_p_values(): Estimate the p-value for the maximum cos value of each gene. During this step, the distribution of the maximum cos value is fitted to $k$ Gaussian functions, $$f_\text{dist}(x)$$, where $k$ is the number of subtypes. Then the p-value can be estimated as $p = \int_{x_\text{min}}_{x_\text{max}}{f_\text{dist}(x) dx}$
    - Parameters:
        - N/A
    - Outputs:
        - Store the p-value for the maximum cos value of each gene at the column 'p.value' of the class attribute, df_cos.
        - Store the multiple test q-value for the maximum cos value of each gene at the column 'q.value' of the class attribute, df_cos.
- obtain_subtype_markers(pThre, qThre, top, per): 
    - Parameters:
        - pThre (): Default value is None.
        - qThre (): Default value is 0.05.
        - top (): Default value is None.
        - per (): Default value is None.
    - Outputs:
        - N/A
- plot_simplex(): 
    - Parameters:
        - N/A
    - Outputs:
        - N/A



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
1. Perform COT calculation step by step:
    1. Import the COT package
    from COT import COT
    cot = COT()
    2. Load the csv file (“input_gene.csv”)
    cot.load_data(filename="input_gene.csv", logarithmic_data=False)

input: “input_gene.csv”

| gene | A   | A   | B   | B   |
| ---- | --- | --- | --- | --- |
| 0    | 0.5 | 0.7 | 0.7 | 0.9 |
| 1    | 1.0 | 1.0 | 0.0 | 0.0 |

output: df_raw

| gene | A   | A.1 | B   | B.1 |
| ---- | --- | --- | --- | --- |
| 0    | 0.5 | 0.7 | 0.7 | 0.9 |
| 1    | 1.0 | 1.0 | 0.0 | 0.0 |

    3. Generate the subtype mean values
    cot.generate_subtype_means()

output: df_mean

| gene | A   | B   |
| ---- | --- | --- |
| 0    | 0.6 | 0.8 |
| 1    | 1.0 | 0.0 |

    4. Generate the cosine values
    cot.generate_cos_values()

output: df_cos

| gene | cos | subtype |
| ---- | --- | ------- |
| 0    | 0.8 | B       |
| 1    | 1.0 | A       |

    5. Export the cosine values (without a threshold)
    cot.save_cos_values(filename="output_cos.csv", threshold=None, sorted=True)

output: “output_cos.csv”

| gene | cos | subtype |
| ---- | --- | ------- |
| 1    | 1.0 | A       |
| 0    | 0.8 | B       |

    6. Export the cosine values (with a threshold)
    cot.save_cos_values(filename="output_cos.csv", threshold=0.9, sorted=True)

output: “output_cos_threshold=0.9.csv”

| gene | cos | subtype |
| ---- | --- | ------- |
| 1    | 1.0 | A       |

2. Or use the pipeline:
    from COT import COT
    
    cot = COT()
    cot.cos_pipeline(input_file="input_gene.csv",output_file="output_cos.csv", logarithmic_input=False, sorted_output=True, output_threshold=0.9)

Then we will obtain the same output with a single step.
