#%%
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr


#%%
def read_netCDF(input_file:str,verbose:bool=False)->np.array:
    data_input=xr.open_dataset(input_file)
    if verbose:
        print(f"Experiment date:\t{datetime.strptime(data_input.experiment_date_time_stamp,'%Y%m%d%H%M%S%z')}")
        print(f"netCDF file created:\t{datetime.strptime(data_input.netcdf_file_date_time_stamp,'%Y%m%d%H%M%S%z')}")
        print(f"Experiment type:\t{data_input.experiment_type}")
        print(f"Number of scans:\t{data_input.number_of_scans}")
        print(f"Global m/z interval:\t [{data_input.global_mass_min}:{data_input.global_mass_max}]")
    
    global_mass_min=data_input.global_mass_min
    global_mass_max=data_input.global_mass_max

    mass_input=np.array(data_input["mass_values"][:].data)
    abundance_input=np.array(data_input["intensity_values"][:].data)
    scan_index_input=data_input["scan_index"][:].data
    mass_output=np.array(np.split(mass_input,scan_index_input[1:],axis=0))
    abundance_output=np.array(np.split(abundance_input,scan_index_input[1:],axis=0))
    if verbose:
        print(abundance_output[0][:10])
        print(mass_output[0][:10])
    return mass_output,abundance_output,int(global_mass_min),int(global_mass_max),int(data_input.number_of_scans)

#%%
def convert_netCDF_offset(input_file:str,offset:float=0.0,output:str="dataframe",verbose:bool=False):
    mass_output,abundance_output,global_mass_min,global_mass_max,num_scans=read_netCDF(input_file)

    def zip_mod(x):
        data_dict={}
        mass_values_int=np.ceil(mass_output[x]+offset).astype(np.int32)
        for i in range(len(mass_values_int)):
            if mass_values_int[i] in data_dict.keys():
                data_dict[mass_values_int[i]]+=abundance_output[x][i]
            else:
                data_dict.update({mass_values_int[i]:abundance_output[x][i]})
        return data_dict
    
    if output=="array":
        data_array=np.zeros(shape=(int(num_scans),int(global_mass_max-global_mass_min+1)),dtype=np.float64)
        scan_num=0
        for scan in map(zip_mod,range(num_scans)):
            for key,value in scan.items():
                data_array[scan_num][key-global_mass_min]=value
            scan_num+=1
        return data_array
    elif output=="dataframe":
        data_array=data_array=np.array(list(map(zip_mod,range(num_scans))))
        return pd.DataFrame.from_records(data_array,columns=[x for x in range(global_mass_min,global_mass_max+1)]).fillna(value=0)