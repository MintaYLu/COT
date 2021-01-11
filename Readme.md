# Readme
**COT**


## **Description**

**COT** is a Python module for machine learning built on top of NumPy and Pandas, and is distributed under the MIT license. Cosine based One-sample Test (COT) is an accurate and efficient method -to detect SMG among many subtypes using subtype-enriched expression profiles. Basically, COT uses the cosine similarity between a molecule’s cross-subtype expression pattern and the exact mathematical definition of SMG as the test statistic, and formulates the detection problem as a one-sample test. Under the assumption that a significant majority of genes are associated with the null hypothesis, COT approximates the empirical null distribution for calculating p values. The project was developed and maintained by Virginia Tech CBIL Group.


## **Installation**

**Dependencies**
COT requires:

- python (>= 3.7.4)
- numpy (>= 1.16.5)
- scipy (>= 1.3.1)
- pandas (>= 0.25.1)
- statsmodels (>= 0.10.1)
- scikit-learn (>= 0.21.3)
- seaborn (>= 0.9.0)
- matplotlib (>= 3.1.1)

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
        - logarithmic_data (True/False): If True, then the input data, df_raw and df_mean, is logarithmic, and we perform the exponential transformation to restore the raw gene features during the data loading process. Default value is False.
        - normalization (True/False): If True, then the feature-wise TPM normalization will be applied to the input data, df_raw and df_mean. Default value is True.
        - silent (True/False): If True, turn off the standard output. Default value is False.
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
        - Store the maximum cos value and the corresponding subtype for each gene at the class attribute, df_cos, where the index is the gene name, and the column "cos", "subtype" represent the maximum cos value and the corresponding subtype, respectively.
- estimate_p_values(): Estimate the p-value for the maximum cos value of each gene. During this step, the distribution function of the maximum cos value, <img src="https://render.githubusercontent.com/render/math?math=f_\text{dist}(x)">, is fitted to <img src="https://render.githubusercontent.com/render/math?math=k"> Gaussian functions, where <img src="https://render.githubusercontent.com/render/math?math=k"> is the number of subtypes. Then the p-value can be estimated as <img src="https://render.githubusercontent.com/render/math?math=p(x) = \frac{\int_{x}^{x_\text{max}}{f_\text{dist}(t) dt}}{\int_{x_\text{min}}^{x_\text{max}}{f_\text{dist}(t) dt}}">, where <img src="https://render.githubusercontent.com/render/math?math=x"> represents the maximum cos value for each gene, and <img src="https://render.githubusercontent.com/render/math?math=x_\text{min} = \frac{1}{\sqrt{k}}">, <img src="https://render.githubusercontent.com/render/math?math=x_\text{max} = 1"> are its minimum and maximum values, respectively. Furthermore, the multiple test q-values are calculated using the FDR control method. Note that the true positive instances, which appear close to <img src="https://render.githubusercontent.com/render/math?math=x = 1"> and above the fitted distribution function <img src="https://render.githubusercontent.com/render/math?math=f_\text{dist}(x)">, are removed using an iterative process.
    - Parameters:
        - N/A
    - Outputs:
        - Store the p-value for the maximum cos value of each gene at the column "p.value" of the class attribute, df_cos.
        - Store the multiple test q-value for the maximum cos value of each gene at the column "q.value" of the class attribute, df_cos.
- obtain_subtype_markers(pThre, qThre, top, per): Obtain a dictionary showing the marker genes for each subtype, and filtered by p-value, q-value and maximum cos value thresholds.
    - Parameters:
        - pThre (float, [0, 1]): If set to be a finite value, then only the instances with "p.value" <= pThre will be selected. Default value is None.
        - qThre (float, [0, 1]): If set to be a finite value, then only the instances with "q.value" <= qThre will be selected. Default value is 0.05.
        - top (int): If set to be a finite value, then only the instances with the top n (n = parameter "top") maximum cos values will be selected. Default value is None.
        - per (int): If set to be a finite value, then only the instances with the top n (n = parameter "per") maximum cos values per subtype will be selected. Default value is None.
    - Outputs:
        - Store the marker genes for each subtype as a python dictionary, where the subtypes and corresponding marker genes are the keys and values, respectively, at the class attribute, markers.
- plot_simplex(): Simplex plot of the cos values for each gene, where the radius and angle represent the cos values and subtypes, respectively. The marker genes for seperate subtypes, which are obtained from the class method, obtain_subtype_markers(), are labeled with different colors.
    - Parameters:
        - N/A
    - Outputs:
        - N/A

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
