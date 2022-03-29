# Readme

## **Description**

**COT** is a Python module for machine learning built on top of NumPy and Pandas, and is distributed under the MIT license. Cosine based One-sample Test (COT) is an accurate and efficient method to detect Marker Genes (MG) among many subtypes using subtype-enriched expression profiles. Basically, COT uses the cosine similarity between a molecule’s cross-subtype expression pattern and the exact mathematical definition of MG as the test statistic, and formulates the detection problem as a one-sample test. Under the assumption that a significant majority of genes are associated with the null hypothesis, COT approximates the empirical null distribution for calculating p-values. The project was developed and maintained by Virginia Tech CBIL Group.


## **Installation**

**Dependencies**

COT requires:

- python (>= 3.7.4)
- numpy (>= 1.19.5)
- scipy (>= 1.3.1)
- pandas (>= 1.2.0)
- statsmodels (>= 0.10.1)
- scikit-learn (>= 0.21.3)
- seaborn (>= 0.11.1)
- matplotlib (>= 3.1.1)

**Installation**

To install from Github, run:

    pip install git+https://github.com/MintaYLu/COT.git

To install from a local copy, please go to the main package folder and run:

    python setup.py install

## Example:
### Perform COT calculation step by step:
    1. Import the COT package
    from COT.COT import COT
    
    2. Create COT class instance and load the raw data
    cot = COT(df_raw=df_raw, normalization=False)

input:
df_raw

Note that the input of COT should be batch-corrected.

| gene | S1  | S2  | S3  | S4  |
| ---- | --- | --- | --- | --- |
| 0    | 0.5 | 0.7 | 0.7 | 0.9 |
| 1    | 1.0 | 1.0 | 0.0 | 0.0 |

output:
cot.df_raw

| gene | S1  | S2  | S3  | S4  |
| ---- | --- | --- | --- | --- |
| 0    | 0.5 | 0.7 | 0.7 | 0.9 |
| 1    | 1.0 | 1.0 | 0.0 | 0.0 |

    3. Generate the subtype mean values
    cot.generate_subtype_means(subtype_label=subtype_label)

input:
subtype_label = ["A", "A", "B", "B"]

output:
cot.subtypes {"A": ["S1", "S2"], "B": ["S3", "S4"]}

cot.df_mean

| gene | A   | B   |
| ---- | --- | --- |
| 0    | 0.6 | 0.8 |
| 1    | 1.0 | 0.0 |

    4. Generate the cosine values
    cot.generate_cos_values()

output:
cot.df_cos

| gene | cos | subtype |
| ---- | --- | ------- |
| 0    | 0.8 | B       |
| 1    | 1.0 | A       |

    5. Estimate the p-values
    cot.estimate_p_values()
Attention: too few genes may not work for predicting p-values. Please remove NaN before this step.
output:
cot.df_cos

| gene | cos | subtype | p.value | q.value |
| ---- | --- | ------- |---------|---------|
| 1    | 1.0 | A       | ?       | ?       |
| 0    | 0.8 | B       | ?       | ?       |

Cannot calculate p-values due to the limited genes numbers, please see the Example_GSE28490 for p-values computation.

    6. Obtain the subtype markers
    cot.obtain_subtype_markers()

output:
cot.markers = {"A": [1], "B": [0]}

    7. Plot the simplex
    cot. plot_simplex()

    8. Plot the heatmap
    cot.plot_heatmap()

### Or use the pipeline:
    from COT.COT import COT
    
    cot = COT(df_raw=df_raw, normalization=False)
    cot.cos_pipeline(subtype_label=subtype_label, top=2)

Then we will obtain the same output with a single step.


## **License**
This project is licensed under the MIT License - see the LICENSE.txt file for details

## **Citing**
If you have used this tool please cite:

Lu, Y., C.-T. Wu, S. J. Parker, L. Chen, G. Saylor, J. E. Van Eyk, D. M. Herrington and Y. Wang (2021). "COT: an efficient Python tool for detecting marker genes among many subtypes." bioRxiv, 2021.01.10.426146

https://doi.org/10.1101/2021.01.10.426146
