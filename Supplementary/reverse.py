#%%
from netCDF2csv import convert_netCDF_offset
from scipy.optimize import minimize_scalar
import numpy as np
import pandas as pd

#%%
def reverse_algorithm(input_netCDF:str,input_CSV:str,program:str="",output_table:bool=True,read_csv_args:tuple=()):
    if program.lower()=="chromatof":
        original_data=pd.read_csv(input_CSV)
        original_data=original_data.drop(labels=['Sample','Time','Scan'],axis='columns')
        original_data.fillna(value=0,inplace=True)
        original_data=original_data.drop(labels=0, axis="index")
        original_data=original_data.values
        def find_offset(offset:float)->float:
            algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=offset,output="dataframe")
            algorithm_guess=algorithm_guess.drop(labels=[0,len(algorithm_guess)-1],axis="index")
            return np.linalg.norm(original_data-algorithm_guess.values)
    elif program.lower()=="chemstation":
        original_data=pd.read_csv(input_CSV)
        original_data=original_data.drop(labels=['Scan'],axis='columns')
        original_data.fillna(value=0,inplace=True)
        original_data=original_data.values
        def find_offset(offset:float)->float:
            algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=offset,output="dataframe")
            return np.linalg.norm(original_data-algorithm_guess.values)
    elif program.lower()=="amdis":
        original_data=pd.read_csv(input_CSV,delimiter='\t')
        original_data=original_data.drop(labels=['Scan','Time','TIC',original_data.columns[-1]],axis='columns')
        while sum(original_data.iloc[:,0])==0:
            original_data=original_data.drop([original_data.columns[0]],axis='columns')
        else:
            start_mz=int(original_data.columns[0])
        def find_offset(offset:float)->float:
            algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=offset,output="dataframe")
            if original_data.columns[0]!=algorithm_guess.columns[0]:
                original_data.columns=pd.Index(original_data.columns.to_numpy().astype(int)+algorithm_guess.columns[0]-start_mz).astype(int)
            algorithm_guess=algorithm_guess.drop(algorithm_guess.columns[-start_mz+algorithm_guess.columns[0]-1:],axis="columns")
            res=original_data-algorithm_guess.iloc[:-1,:]
            return np.linalg.norm(original_data.to_numpy()-algorithm_guess.to_numpy()[:-1])
    elif program.lower()=="openchrom":
        original_data=pd.read_csv(input_CSV,delimiter=';')
        original_data=original_data.drop(labels=['RT(milliseconds)','RT(minutes) - NOT USED BY IMPORT','RI'],axis='columns')
        original_data.fillna(value=0,inplace=True)
        original_data=original_data.values
        def find_offset(offset:float)->float:
            algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=offset,output="dataframe")
            return np.linalg.norm(original_data-algorithm_guess.values)
    else:
        original_data=pd.read_csv(input_CSV)
        original_data=original_data.values
        def find_offset(offset:float):
            algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=offset,output="dataframe")
            return np.linalg.norm(original_data-algorithm_guess.values)
    
    result=minimize_scalar(find_offset,bounds=[-1,0],options={"xtol":1e-4})
    if output_table:
        algorithm_guess=convert_netCDF_offset(input_file=input_netCDF,offset=result.x,output="dataframe")
        algorithm_guess.to_csv("".join([*input_netCDF.split('.')[:-1],"_",program.lower(),"_output.csv"]))
    return result