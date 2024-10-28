import os

def rename_files(directory):
    # Loop through each file in the directory
    for filename in os.listdir(directory):
        # Construct the full file path
        full_path = os.path.join(directory, filename)
        
        n = str(int(filename[4])+15)
        if 'sky' in filename:
            new_filename = "Light_" + filename[-5] + "_" + n + '.fit'                
            new_full_path = os.path.join(directory, new_filename)
        if 'dark' in filename:
            new_filename = "Dark_" + n + '.fit'                
            new_full_path = os.path.join(directory, new_filename)
        if 'bias' in filename:
            new_filename = "Bias_" + n + '.fit'                
            new_full_path = os.path.join(directory, new_filename)
        
        # Rename the file
        os.rename(full_path, new_full_path)
        print(f'Renamed: {filename} to {new_filename}')

# Specify the directory containing the files
directory_path = "C:/Users/lhung/Research/Fisheye/Data_raw/NIOB_FWS/set4"
rename_files(directory_path)