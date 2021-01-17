## Functions

- **COT(df_raw, df_mean, logarithmic_data, normalization, silent)**:  COT class constructor. During this step either the raw gene data or the subtype mean data will be imported and preprocessed, i.e., performing exponential transformation or sample-wise normalization.
    - Parameters:
        - df_raw (pandas.DataFrame): The raw gene data stored as a pandas Dataframe. df_raw has the following format: each row represents a gene, and each column represents a sample. Default value is None.
        - df_mean (pandas.DataFrame): The subtype mean data stored as a pandas Dataframe. df_mean has the following format: each row represents a gene, and each column represents a subtype. Default value is None.
        - logarithmic_data (True/False): If True, then the input data, df_raw and df_mean, is logarithmic, and we perform the exponential transformation to restore the raw gene features during the data loading process. Default value is False.
        - normalization (True/False): If True, then the feature-wise TPM normalization will be applied to the input data, df_raw and df_mean. Default value is True.
        - silent (True/False): If True, turn off the standard output. Default value is False.
    - Outputs:
        - Store the raw gene data at the class attribute, df_raw.
        - Store the subtype mean data at the class attribute, df_mean.

- **generate_subtype_means(subtype_label)**: Generate a dataframe to restore the subtype means, which are the mean values of the gene data within the same subtypes.
    - Parameters:
        - subtype_label (list(str)): The subtype label for each sample stored as a python list. Default value is None.
    - Outputs:
        - Store the subtype mean data at the class attribute, df_mean, where the index is the gene name, and the column represents the subtype name.

- **generate_cos_values()**: Compute the cos value for each subtype, then find the maximum cos value and the corresponding subtype.
    - Parameters:
        - N/A
    - Outputs:
        - Store the maximum cos value and the corresponding subtype for each gene at the class attribute, df_cos, where the index is the gene name, and the column "cos", "subtype" represent the maximum cos value and the corresponding subtype, respectively.

- **estimate_p_values()**: Estimate the p-value for the maximum cos value of each gene. During this step, the distribution function of the maximum cos value, <img src="https://render.githubusercontent.com/render/math?math=f_\text{dist}(x)">, is fitted to <img src="https://render.githubusercontent.com/render/math?math=k"> Gaussian functions, where <img src="https://render.githubusercontent.com/render/math?math=k"> is the number of subtypes. Then the p-value can be estimated as <img src="https://render.githubusercontent.com/render/math?math=p(x) = \frac{\int_{x}^{x_\text{max}}{f_\text{dist}(t) dt}}{\int_{x_\text{min}}^{x_\text{max}}{f_\text{dist}(t) dt}}">, where <img src="https://render.githubusercontent.com/render/math?math=x"> represents the maximum cos value for each gene, and <img src="https://render.githubusercontent.com/render/math?math=x_\text{min} = \frac{1}{\sqrt{k}}">, <img src="https://render.githubusercontent.com/render/math?math=x_\text{max} = 1"> are its minimum and maximum values, respectively. Furthermore, the multiple test q-values are calculated using the FDR control method. Note that the true positive instances, which appear close to <img src="https://render.githubusercontent.com/render/math?math=x = 1"> and above the fitted distribution function <img src="https://render.githubusercontent.com/render/math?math=f_\text{dist}(x)">, are removed using an iterative process.
    - Parameters:
        - N/A
    - Outputs:
        - Store the p-value for the maximum cos value of each gene at the column "p.value" of the class attribute, df_cos.
        - Store the multiple test q-value for the maximum cos value of each gene at the column "q.value" of the class attribute, df_cos.
        - Create two empty columns, "score" and "score_subtype", in the class attribute, df_cos, which are used to store the user-defined scores.

- **obtain_subtype_markers(pThre, qThre, cosThre, top, per, scoreThre)**: Obtain a dictionary showing the marker genes for each subtype, and filtered by p-value, q-value, and maximum cos value thresholds, or user-defined score threshold.
    - Parameters:
        - pThre (float, [0, 1]): If set to be a finite value, then only the instances with "p.value" <= pThre will be selected. Default value is None.
        - qThre (float, [0, 1]): If set to be a finite value, then only the instances with "q.value" <= qThre will be selected. Default value is 0.05.
        - cosThre (float, [0, 1]): If set to be a finite value, then only the instances with "cos" <= cosThre will be selected. Default value is None.
        - top (int): If set to be a finite value, then only the instances with the top n (n = parameter "top") maximum cos values will be selected. Default value is None.
        - per (int): If set to be a finite value, then only the instances with the top n (n = parameter "per") maximum cos values per subtype will be selected. Default value is None.
        - scoreThre (float): If set to be a finite value, then the class attribute, df_cos, is sorted by "score" (instead of "cos", when scoreThre is set to be None) in descending order. Moreover, only the instances with "score" <= scoreThre will be selected. Default value is None.
    - Outputs:
        - Store the marker genes for each subtype as a python dictionary, where the subtypes and corresponding marker genes are the keys and values, respectively, at the class attribute, markers.

- **plot_simplex()**: Simplex plot of the cos values for each gene, where the radius and angle represent the cos values and subtypes, respectively. The marker genes for separate subtypes, which are obtained from the class method, obtain_subtype_markers(), are labeled with different colors.
    - Parameters:
        - N/A
    - Outputs:
        - N/A

- **plot_heatmap()**: Heatmap plot of the raw gene data, where the x and y axes represent the gene and sample, respectively, for the marker genes of separate subtypes. The raw gene data is standardized and then plotted on a logarithmic scale.
    - Parameters:
        - N/A
    - Outputs:
        - N/A

- **cot_pipeline(subtype_label, pThre=None, qThre=0.05, top=None, per=None)**: The COT pipeline: generate_subtype_means(subtype_label) -> generate_cos_values() -> estimate_p_values() -> obtain_subtype_markers(pThre=pThre, qThre=qThre, top=top, per=per) -> plot_simplex() -> plot_heatmap().
    - Parameters:
        - subtype_label (list(str)): The subtype label for each sample stored as a python list. Default value is None.
        - pThre (float, [0, 1]): If set to be a finite value, then only the instances with "p.value" <= pThre will be selected. Default value is None.
        - qThre (float, [0, 1]): If set to be a finite value, then only the instances with "q.value" <= qThre will be selected. Default value is 0.05.
        - top (int): If set to be a finite value, then only the instances with the top n (n = parameter "top") maximum cos values will be selected. Default value is None.
        - per (int): If set to be a finite value, then only the instances with the top n (n = parameter "per") maximum cos values per subtype will be selected. Default value is None.
    - Outputs:
        - N/A
