#-----------------------------------------------------------#
# iers_data.py
#
# NPS Night Skies Program
#
#
# Script that downloads and updates the IERS-A, IERS-B,
# and Leap Second tables. This script is run instead of
# allowing astropy to perform automatic updates in order
# to avoid SSL Certificate errors. Astropy uses urllib
# to download these files while this script uses requests.
# Astropy code would have to be modified to make urllib
# work, so instead we provide a custom routine here to
# perform the table downloads. We use the requests package
# to perform the download since it can be configured to 
# use SSL certficates from the default system store by 
# simply installing the pip-system-certs package.
#
# https://pypi.org/project/pip-system-certs/
#
#Input:
#    None
#Output:
#    None
#
#History:
#    Zach Vanderbosch -- Created on 11/25/2024
#-----------------------------------------------------------#

import os
import time
import stat
import pathlib
import requests
import astropy_iers_data as aid

#-----------------------------------------------------------#
#-----------------------------------------------------------#
#           Define IERS Table URLs and Filenames            #
#       https://github.com/astropy/astropy-iers-data        #
#-----------------------------------------------------------#

# Set maximum IERS Table age (in seconds).
# If an IERS Table already exists locally and is less than 
# this age, it will not be re-downloaded, saving some time.
MAX_TABLE_AGE = 86400.0  # 86400 sec = 1 day

# Path where tables are saved (astropy-iers-data package location).
# This is where all astropy routines will look for the tables when
# the astropy.utils.iers.conf parameter auto_download = False.
DATA = pathlib.Path(aid.__file__).resolve().parent / "data"

# IERS-A default file name, URL, and ReadMe with content description
IERS_A_FILE = str(DATA / "finals2000A.all")
IERS_A_URL = "https://datacenter.iers.org/data/9/finals2000A.all"
IERS_A_URL_MIRROR = "https://maia.usno.navy.mil/ser7/finals2000A.all"

# IERS-B default file name, URL, and ReadMe with content description
IERS_B_FILE = str(DATA / "eopc04.1962-now")
IERS_B_URL = "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now"

# LEAP SECONDS default file name, URL, and alternative format/URL
IERS_LEAP_SECOND_FILE = str(DATA / "Leap_Second.dat")
IERS_LEAP_SECOND_URL = "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat"
IERS_LEAP_SECOND_URL_MIRROR = "https://data.iana.org/time-zones/data/leap-seconds.list"


#-----------------------------------------------------------#

# File download function
def download_file(url, filename):

    # Use requests to retrieve file from url
    response = requests.get(url)
    status = response.status_code

    if status != 200:
        reason = response.reason
        message = f"  Could not download table at {url}\n" +\
                  f"  Status Code [{status}]: {reason}"
        return status
    else:
        print(f"  Successfully downloaded table at {url}")

    # Save response to file
    with open(filename, "w") as f:
        f.write(response.text)
    
    return status


# Function to Check File Ages if they already exist
def check_file_age(filename):

    if os.path.isfile(filename):
        time_since_mod = time.time() - os.stat(filename)[stat.ST_MTIME]
        return time_since_mod
    else:
        return np.nan



# Execute file downloads
def main():

    print('Updating IERS tables...')

    # Try downloading IERS-A file
    file_age = check_file_age(IERS_A_FILE)
    if file_age > MAX_TABLE_AGE:
        status = download_file(IERS_A_URL, IERS_A_FILE)
        if status != 200:
            status = download_file(IERS_A_URL_MIRROR, IERS_A_FILE)
    elif file_age < MAX_TABLE_AGE:
        print('  Existing %s table found with age less than %.1d day' %(
            IERS_A_FILE.split("\\")[-1], 
            MAX_TABLE_AGE/86400
        ))

    # Try downloading IERA-B file
    file_age = check_file_age(IERS_B_FILE)
    if file_age > MAX_TABLE_AGE:
        status = download_file(IERS_B_URL, IERS_B_FILE)
    elif file_age < MAX_TABLE_AGE:
        print('  Existing %s table found with age less than %.1d day' %(
            IERS_B_FILE.split("\\")[-1], 
            MAX_TABLE_AGE/86400
        ))

    # Try downloading Leap Seconds file
    file_age = check_file_age(IERS_LEAP_SECOND_FILE)
    if file_age > MAX_TABLE_AGE:
        status = download_file(IERS_LEAP_SECOND_URL, IERS_LEAP_SECOND_FILE)
        if status != 200:
            status = download_file(IERS_LEAP_SECOND_URL_MIRROR, IERS_LEAP_SECOND_FILE)
    elif file_age < MAX_TABLE_AGE:
        print('  Existing %s table found with age less than %.1d day' %(
            IERS_LEAP_SECOND_FILE.split("\\")[-1], 
            MAX_TABLE_AGE/86400
        ))


# Allow script to be run directly
if __name__ == "__main__":
    main()
