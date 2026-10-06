# Targeted LC-MS/MS Quality Control and Data Integration Pipeline for PFAS Analysis
The workflow processes high-resolution LC-MS/MS data generated on a SCIEX QTOF instrument and evaluates quality control metrics following EPA 1633A guidelines.

### Workflow Overview

##### Prepare your project folder
Collect all raw input files (exported text files from SCIEX Analyst) into a dedicated project directory.
Ensure that this project directory is located in the same parent folder as the scripts (`data_analysis.ipynb`, `create_project_folder.ipynb`, and `utils.py`).

##### Make sure your instrumentation detection limits are available.
IDL values for various methods are deposited in the lab_parameters directory.
If using a new method, calculate instrument detection limits (IDLs) with: `calculate_idl.py`.

##### Generate template parameter files
Open `create_project_folder.ipynb` and set your **project_folder** variable, so that it points to the project directory containing your raw data.
This notebook detects sample and compound names in your raw data and creates input parameter templates:
   - `eis_paramters.csv`
   - `compound_parameters.csv`
   - `sample_parameters.csv`
   - `simulation_parameters.csv`
   - `mdl_parameters.csv`

Once the script run successfully you will find the input files in your project folder under **simulation_parameters**.
Update these files according to your study. Detailed instructions are included in the notebook.

##### Run the main analysis 
Execute `data_analysis.ipynb`, ensuring the `project_folder` variable is set correctly.  
The notebook will generate:
   - an Excel QA/QC report
   - a long-format `.csv` table  
Both will be accessible in your project folder under **processed_data**.

##### Dependencies & Folder Structure

Both notebooks `create_project_folder.ipynb` and `data_analysis.ipynb` rely on:
- `utils.py`
- **lab_parameters/**
Ensure these are stored together if downloaded manually.

### Input Data & Naming Conventions
Input files must come from **SCIEX Analyst** exported tables or **Waters Connect**. These are either .csv or .txt files.

Files must end with:
- `_core` → data from the core method  
- `_extended` → data from the extended method  

Files to be combined from the same batch must share the same prefix:
Example:
- batch_1_core.txt
- batch_1_extended.txt

A complete test set is available in the **test/** directory.

### Software Requirements

You may run the pipeline in:
- **JupyterLab** (via Anaconda)
- **Google Colab**
- Any environment supporting Jupyter notebooks

### Definitions

| Term | Description |
|------|-------------|
| **Extracted Internal Standards (EIS)** | Mass-labeled standards added **before extraction** (formerly *IDA*). |
| **Non-Extracted Internal Standards (NIS)** | Standards added **after extraction**, before injection (also *IPS*). |
| **Target analytes** | PFAS native compounds quantified by LC-MS/MS. |
| **HRMS channel** | High-resolution TOF channel used to confirm native ions. For triple quad (mass spec) it is the confirmation channel.|
| **MS/MS channel** | Fragmentation channel used for quantification. Also known as quantitation channel.|

### Video Tutorials
|Title|Link|
|------|-------------|
| **Overview** | [Youtube Link](https://youtu.be/pYqSxKKX_-k?si=4iWsj2EZ3F4vb5WO) |
| **What is Code?** | [Youtube link](https://youtu.be/A4K5VRTIwlY?si=dAjfHoLblOjD-Czn)|
| **QA/QC Briefing** | [Youtube Link](https://youtu.be/OOGTGDluwRA?si=0UdOD9ScEpyWxlNk)|
| **Step by step instruction (for Sciex data)** | [Youtube Link](https://youtu.be/A4K5VRTIwlY?si=dAjfHoLblOjD-Czn)|
| **Processing Waters data** | [Youtube Link](https://youtu.be/5bKdysppm-g?si=qsjoTaqhd9vVQP6-)|

### Contact
johanna.ganglbauer@uri.edu

### Ackknowledgement
This work is supported by the National Institute of Environmental Health under grant P42ES027706 in superfund research project Sources, Transport, Exposure & Effects of PFAS (STEEP).
<img src="https://web.uri.edu/wp-content/uploads/sites/1022/NIEHS_SRP_Log_horz_600.png" width="40%">
<img src="https://web.uri.edu/wp-content/uploads/sites/1022/NIEHS_SRP_Log_horz_600.png](https://web.uri.edu/wp-content/uploads/sites/1022/steep-logo.png" width="40%">
