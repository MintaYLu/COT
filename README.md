# COT
COT Readme


License


Copyright <2020> 


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:


The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.


THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
Readme
COT is a Python module for machine learning built on top of NumPy and Pandas, and is distributed under the MIT license. 

Installation
Dependencies
COT requires:
•	Python (>= 3.6)
•	NumPy (>= 1.13.3)
•	Pandas (>= 0.25.0)


Installation
If you already have a working installation of numpy and pandas, the easiest way to install COT is downloading the COT.py file at: <Add link here>.


Usage
Once the COT.py is downloaded, the user need to either put the COT.py file into the working directory or add it to the Python module search path by:
import sys
sys.path.append("the-path-of-the-COT.py-file")


Then, the user should be able to use the COT module. See the following code for an example:
from COT import COT


cot = COT()
cot.cos_pipeline(input_file, output_file, logarithmic_input=False, sorted_output=True, output_threshold=None)
The module will read the original data from the csv file (input_file), then generate the output csv file (output_file), which contains three columns: gene_id, cos_value, subtype_label.
If the output_threshold is set to some finite value between 0 and 1, then another csv file, which only stores the genes with a cos_value >= output_threshold, will be generated.


Description
Cosine based One-sample Test (COT) Python package is designed to detect marker genes among many subtypes using subtype-enriched expression profiles, which is applicable to multi-omics data. 


Functions
•	COT(silent): COT class constructor.
o	Parameters:
	silent (True/False): If “True”, turn off the standard output. Default value is “False”.
•	load_data(filename, logarithmic_data): Import the raw gene data from the csv file. The csv file has the following format: each row represents a gene, and each column represents the expression level of a sample; and the first column represents the gene names, and the first row represents the subtypes of the corresponding samples.
o	Parameters:
	filename (str): The filename of the gene data.
	logarithmic_data (True/False): If “True”, then the input data in the csv file is logarithmic, and we perform the exponential transformation to restore the raw gene features during the data loading process. Default value is “False”.
o	Outputs:
	Store the raw gene data at the class attribute, df_raw. In the dataframe, df_raw, the index is the gene name, and the column is the sample name, which is “subtype_name”, “subtype_name.1”, “subtype_name.2”, etc.
•	generate_subtype_means(): Generate a dataframe to restore the subtype means, which are the mean values of the gene data within the same subtypes.
o	Parameters:
	N/A
o	Ouputs:
	Store the subtype mean data at the class attribute, df_mean, where the index is the gene name, and the column represents the subtype name.
•	generate_cos_values(): Compute the cos value for each subtype, then find the maximum cos value and the corresponding subtype. 
o	Parameters:
	N/A
o	Outputs:
	Store the maximum cos value and the corresponding subtype for each gene at the class attribute, df_cos, where the index is the gene name, and the column “cos”, “subtype” represent the maximum cos value and the corresponding subtype, repectively.
•	save_cos_values(filename, threshold=None, sorted=True): Save the gene names, maximum cos values and the corresponding subtypes into a csv file.
o	Parameters:
	filename (str): The filename of the output data.
	threshold (float, [0, 1]): If threshold is set to some finite value, then another csv file will be generated, which only stores the genes with a maximum cos value >= threshold. Default value is “None”.
	sorted (True/False): If “true“, the gene will be sorted by the cos value in the descending order. Default value is “True”.
o	Outputs:
	A csv file (filename) stores all the genes with their cos values and the corresponding subtypes.
	Another csv file (filename + threshold information) when threshold is set to some finite value, which stores the data when cos value >= threshold.
•	cos_pipeline(input_file, output_file, logarithmic_input, sorted_output, output_threshold): A pipeline to process the COT algorithm.
o	Parameters:
	input_file (str): Same as filename in the class method load_data.
	output_file (str): Same as filename in the class method save_cos_values.
	logarithmic_input (True/False): Same as logarithmic_data in the class method load_data.
	sorted_output (True/False): Same as sorted in the class method save_cos_values.
	output_threshold (float, [0, 1]): Same as threshold in the class method save_cos_values.
o	Outputs:
	csv file(s) which stores the input genes with their cos values and the corresponding subtypes. Same as the outputs of the class method save_cos_values.
 








Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.
.
Contributing
State if you are open to contributions and what your requirements are for accepting them.
For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.
You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.


Authors and acknowledgment
The project was developed and maintained by Virginia Tech CBIL Group.


License
This package is distributed under the MIT license.


