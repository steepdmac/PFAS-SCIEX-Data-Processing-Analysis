"""
Calculates instrumentation detection limits from calibration data.
Adopt the function input in the bottom lines of the script and run the script to calculate instrumentation detection limits from raw data.
"""
from typing import Optional
import pandas as pd
import numpy as np
import os

from utils import clean_up_data, reassign_tof_nis_to_eis, get_hrms_and_msms_compounds, \
    get_sample_id_and_name, convert_waters_to_sciex

# global variable definitions
standard_identifiers = 'Avg|EIS|NIS|IDA|IPS|13C|d-|d3-|d5-|18O'

def calculate_idls(
        method_name: str, hrms_identifier: str, data_format: str,
        filepath_core: Optional[str], filepath_extended: Optional[str],
        ):
    """Reads in raw calibration data and automatically evaluates instrumentation detection limits based on signal to noise ratios of 10 for each compound and each channel.
    The results are saved as .csv in the lab_parameters subfolder.

    :param method_name: Name of the method, ideally in the format: {year}_{matrix}_{first name of researcher}. Example: 2025_serum_jingmei
    :type method_name: str
    :param hrms_identifier: Ending of compound names used to identify high-resolution mass spectrometry channels. Examples: '_TOF MS' or '_HRMS'
    :type hrms_identifier: str
    :param data_format: Data format of input data. Either 'waters' or 'sciex'.
    :type data_format: str
    :param filepath_core: Path to the raw data file of the core method containing calibration data from which IDL is calculated.
    :type filepath_core: Optional[str]
    :param filepath_extended: Path to the raw data file of the extended method containing calibration data from which IDL is calculated.
    :type filepath_extended: Optional[str]
    :raises ImportError: If the filepath_core does not point to a .csv or .txt. file.
    :raises ImportError: If the filepath_extended does not point to a .csv or .txt. file.
    :raises ValueError: If no file paths are provided for core or extended data.
    """

    # Define columns of input which are needed for further processes:
    columns_considered = [
        'Sample Index', 'Sample Name', 'Sample ID', 'Sample Type', 'Calculated Concentration',
        'Actual Concentration', 'Component Name', 'Used', 'Signal / Noise'
    ]

    if filepath_core is not None:
        # Read core data
        if filepath_core.endswith('.csv'):
            data_core = pd.read_csv(
                filepath_core, delimiter=',', encoding='utf-8', header=0,
                ).dropna(how="all")
        elif filepath_core.endswith('.txt'):
            data_core = pd.read_csv(
                filepath_core, delimiter='\t', encoding='utf-8', header=0,
                ).dropna(how="all")
        else:
            raise ImportError('Raw input file paths must either be .csv or .txt files.')
        if data_format == 'waters':
            data_core = convert_waters_to_sciex(data_core, hrms_identifier)
        data_core = data_core[columns_considered]
        mask_names = data_core['Sample Name'].str.endswith('Core')
        data_core.loc[~mask_names, 'Sample Name'] = [sample_name + ' Core' for sample_name in data_core['Sample Name'][~mask_names].to_list()]

    if filepath_extended is not None:
        # Read core data
        if filepath_extended.endswith('.csv'):
            data_extended = pd.read_csv(
                filepath_extended, delimiter=',', encoding='utf-8', low_memory=False, header=0,
                ).dropna(how="all")
        elif filepath_extended.endswith('.txt'):
            data_extended = pd.read_csv(
                filepath_extended, delimiter='\t', encoding='utf-8', header=0,
                ).dropna(how="all")
        else:
            raise ImportError('Raw input file paths must either be .csv or .txt files.')
        if data_format == 'waters':
            data_extended = convert_waters_to_sciex(data_extended, hrms_identifier)
        data_extended = data_extended[columns_considered]
        mask_names = data_extended['Sample Name'].str.endswith('Ext')
        data_extended.loc[~mask_names, 'Sample Name'] = [sample_name + ' Ext' for sample_name in data_extended['Sample Name'][~mask_names].to_list()]

    # Combine data
    if filepath_core is not None and filepath_extended is not None:
        data_extended['Sample Index'] = data_extended['Sample Index'] + data_core['Sample Index'].max() + 1
        data = pd.concat([data_core, data_extended], ignore_index=True)
    elif filepath_core is not None:
        data = data_core
    elif filepath_extended is not None:
        data = data_extended
    else:
        raise ValueError('At least one file path must be provided for core or extended data.')
    
    data['Batch Name'] = method_name  # add batch name to data
    
    # extract sample names and compound names from raw data
    sample_list = get_sample_id_and_name(data)
    
    # calls function to get complete list of samples
    data = clean_up_data(data=data, sample_list=sample_list)

    # get the assignment of NIS for HRMS channel EIS right
    data = reassign_tof_nis_to_eis(data=data)

    # get compounds dataframe and delete 'useless compounds'
    compounds, delete_compounds, compounds_available = get_hrms_and_msms_compounds(
        data=data, sample_list=sample_list, hrms_identifier=hrms_identifier, standard_identifiers=standard_identifiers
        )

    # delete detected compounds accordingly. 
    # Usualy HRMS channels from the core method have to be deleted, because they also occur in the extended method, where they are integrated with more care.
    if not delete_compounds.empty:  # check if data frame is empty
        for method in delete_compounds['from method'].unique():  # loop over 'core' method and 'extended method'
            # get all sample indices from relevant method
            if method == 'core':
                indices = [int(index) for index in sample_list['Sample Index Core'].dropna().tolist()]
            elif method == 'extended':
                indices = [int(index) for index in sample_list['Sample Index Extended'].dropna().tolist()]
            # loop over compounds to be deleted within the method and delete them accordingly
            for compound in delete_compounds.loc[delete_compounds['from method'] == method, 'Compound Name'].tolist():
                data = data.loc[~(
                    (data['Component Name'] == compound) & (data['Sample Index'].isin(indices))
                ), :]

    hrms_label = str.upper(hrms_identifier.replace('_', ''))  # get column label for HRMS channel

    idl_data = pd.DataFrame(columns=["Sample Code", "Unit"] + compounds_available)
    idl_data.loc[0, 'Sample Code'] = 'MSMS IDL'
    idl_data.loc[1, 'Sample Code'] = 'MSMS LOQ'
    idl_data.loc[2, 'Sample Code'] = f'{hrms_label} IDL'
    idl_data.loc[3, 'Sample Code'] = f'{hrms_label} LOQ'
    idl_data['Unit'] = 'ng/sample'
    
    calibration_data = data.loc[data['Sample Type'] == 'Standard', :]
    for (msms_compound, hrms_compound) in zip(
        compounds['MSMS Compound Name'].tolist(), compounds[f'{hrms_label} Compound Name'].tolist()
        ):
        msms_data = calibration_data.loc[calibration_data['Component Name'] == msms_compound, :]
        hrms_data = calibration_data.loc[calibration_data['Component Name'] == hrms_compound, :]
        min_idl = 1e-3
        for calibration_point in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            previous = msms_data.loc[(
                (msms_data['Sample ID'].str.contains(f'CS{calibration_point}')) &
                (msms_data['Used'] == True)
            ), :]
            if previous['Signal / Noise'].isna().sum() == 0:
                this_point = msms_data.loc[(
                    (msms_data['Sample ID'].str.contains(f'CS{calibration_point + 1}')) &
                    (msms_data['Used'] == True)
                ), :]
                idl = 3 * this_point['Actual Concentration'].mean() / this_point['Signal / Noise'].mean()
                idl = max(idl, min_idl)
                idl_data.loc[0, msms_compound] = round(idl, ndigits=3)
                idl_data.loc[1, msms_compound] = previous['Actual Concentration'].mean()
                break
            else:
                min_idl = previous['Actual Concentration'].mean()
        min_idl = 1e-3
        for calibration_point in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            previous = hrms_data.loc[(
                (hrms_data['Sample ID'].str.contains(f'CS{calibration_point}')) &
                (hrms_data['Used'] == True)
            ), :]
            if previous['Signal / Noise'].isna().sum() == 0:
                this_point = hrms_data.loc[(
                    (hrms_data['Sample ID'].str.contains(f'CS{calibration_point + 1}')) &
                    (hrms_data['Used'] == True)
                ), :]
                idl = 3 * this_point['Actual Concentration'].mean() / this_point['Signal / Noise'].mean()
                idl = max(idl, min_idl)
                if msms_compound is np.nan:
                    idl_data.loc[2, hrms_compound[:-1 * len(hrms_identifier)]] = round(idl, ndigits=3)
                    idl_data.loc[3, hrms_compound[:-1 * len(hrms_identifier)]] = previous['Actual Concentration'].mean()
                else:
                    idl_data.loc[2, msms_compound] = round(idl, ndigits=3)
                    idl_data.loc[3, msms_compound] = previous['Actual Concentration'].mean()
                break
            else:
                min_idl = previous['Actual Concentration'].mean() # ensure that IDL is not smaller than point where peak was not detected.

    idl_data.to_csv(
        os.path.join('lab_parameters', f'{method_name}_idl.csv'), index=False, encoding='utf-8'
    )


if __name__ == "__main__":
    calculate_idls(
        method_name='2026_watersepa_simon',
        hrms_identifier='_Qual',
        data_format='waters',
        filepath_core=r'simon/EPA_fish/260908_Lake_Trout/Lake_Trout_concentration_no_nis/wet/20260902_EPA_PFAS_Lake_Trout_no NIS_concentration_wet_core.csv',
        filepath_extended= None # r'julie/water/20251123_Water_India_extended.txt',
    )