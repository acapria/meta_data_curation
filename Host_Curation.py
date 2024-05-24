

#Output should be each row - accession, host field, lab host field and accessions in CSV format
# The accession numbers are from working set/genome group we make from BV-BRC
# the input is the bv-brc genome group in csv format

#developed by @AnnaCapria @RoshniBhattacharya @AbrilZuniga

#pip install biopython



#import statements
from Bio import Entrez, SeqIO
import csv
import pandas as pd

# extract the accessions from the genome group csv file
def extract_accessions(file_path):
    df = pd.read_csv(file_path)
    accession_number = df["GenBank Accessions"].tolist()
    return accession_number


#fetch the genbank record
def fetch_genbank_record(accession_number):
    Entrez.email = "bvbrcacapria@gmail.com"     #put your email
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching GenBank record: {e}")
        return None


#extract the host and labhost metadata from a genbank record
def extract_host_info(record):
    host_info = ""
    lab_host_info = ""
    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            if "host" in qualifiers:
                host_info = ', '.join(qualifiers["host"])
            if "lab_host" in qualifiers:
                lab_host_info = ', '.join(qualifiers["lab_host"])
    return host_info, lab_host_info

def main():

    file_path = "/Users/rbhattac/Desktop/Curation/misc/test_genome_group.csv" #path to your genome group csv file
    #extract the accession numbers
    accession_numbers = extract_accessions(file_path)

    #open the output file to write
    with open("/Users/rbhattac/Desktop/test_results_metadata.csv", "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(["Accession", "Host", "Lab Host"])
        for accession_number in accession_numbers:
            print(accession_number)
            record = fetch_genbank_record(accession_number)
            if record:
                accession = record.id
                # Extract host and lab host information
                host_info, lab_host_info = extract_host_info(record)
                # Write data to CSV
                csv_writer.writerow([accession, host_info, lab_host_info])


if __name__ == "__main__":
    main()
