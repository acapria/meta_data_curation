# hostflags.py
import csv
def read_hostnames_from_tsv(filename):
    hostnames = []
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            hostnames.extend(row)
    return hostnames

def determine_flag(hostname):
    labhost_file = "lab_host.tsv"
    isolation_file = "isolation_source.tsv"
    
    labhost_list = read_hostnames_from_tsv(labhost_file)
    isolation_list = read_hostnames_from_tsv(isolation_file)


    # Check if hostname is in labhost_list
    if hostname in labhost_list:
        return "FLAG (Possible_labhost)"

    # Check if hostname is in isolation_list
    elif hostname in isolation_list:
        return "FLAG (Possible_IsoSource)"

    # If hostname is not in either list, return No_matching_Host flag
    else:
        return "FLAG (No_matching_Host)"

# Example usage:
# flag = determine_flag("hostname_to_check")
