#Code to obtain metadata from genbank and biosample records given an genbank accession number
# the input is the bv-brc genome group in csv format

#developed by @AnnaCapria @RoshniBhattacharya @AbrilZuniga

#import statements
from Bio import Entrez, SeqIO
import csv
import pandas as pd
import xml.etree.ElementTree as ET

# Provide your email to NCBI Entrez
Entrez.email = "abrilz2296@gmail.com"

# extract the accessions from the genome group csv file
def extract_accessions(file_path):
    df = pd.read_csv(file_path)
    accession_number = df["GenBank Accessions"].tolist()
    return accession_number

# extract the genbank record given an accession number
def fetch_genbank_record(accession_number):
    #Entrez.email = "bvbrcacapria@gmail.com"     #put your email
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

#extract the biosample accession number from a genbank record
def get_biosample_accession(genbank_accessions):
    handle = Entrez.elink(dbfrom="nucleotide", db="biosample", id=genbank_accessions)
    record = Entrez.read(handle)
    #print(record)
    handle.close()
    biosample_accessions = []
    if record[0]['LinkSetDb']:
        links = record[0]['LinkSetDb'][0]['Link']
        biosample_accessions = [link['Id'] for link in links]
    return biosample_accessions

#extract the biosample record given a biosample accession number
def get_biosample_record(biosample_accession):
    handle = Entrez.efetch(db="biosample", id=biosample_accession, retmode="xml")
    xml_data = handle.read()
    handle.close()
    return xml_data

#extract the metadata from a biosample record
def extract_information(xml_data, genbank_accession):
    # Parse the XML data
    root = ET.fromstring(xml_data)

    # Find the relevant sample attributes
    sample_attributes = root.findall(".//Attribute[@attribute_name]")

    # Extract information from the sample attributes
    extracted_info = {"GenBank Accession": genbank_accession}  # Include GenBank accession number
    for attr in sample_attributes:
        attribute_name = attr.attrib.get('attribute_name')
        attribute_value = attr.text
        extracted_info[attribute_name] = attribute_value

    return extracted_info

def main():

        file_path = "/Users/rbhattac/Desktop/Host_Biosample.csv" #path to your genome group csv file
        #extract the accession numbers
        accession_numbers = extract_accessions(file_path)

        #initialize lists to store the metadata
        genbank_metadata = []
        all_biosample_info = []

        for accession in accession_numbers:
            #print(accession)
            print(f"Processing genbank accession {accession}..")
            record = fetch_genbank_record(accession)
            if record:
                accession_id = accession
                # Extract host and lab host information
                host_info, lab_host_info = extract_host_info(record)
                genbank_metadata.append([accession_id, host_info, lab_host_info])

            # Extract BioSample information
            biosample_accession = get_biosample_accession(accession)

            for bioid in biosample_accession:
                biosample_xml = get_biosample_record(bioid)
                extracted_info = extract_information(biosample_xml, accession)
                all_biosample_info.append(extracted_info)

        print(f"Compiling all the extracted metadata...")
        #print(genbank_metadata)
        df1= pd.DataFrame(genbank_metadata, columns = ['GenBank Accession', 'GenBank_Host', 'GenBank_Lab Host'])
        #print(df1)
        df2 = pd.DataFrame(all_biosample_info)
        #print(df2)

        df_merged = df1.merge(df2, on='GenBank Accession', how='left')
        #print(df_merged)

        # Save the dataframe to a CSV file
        print(f"Printing the results to a csv file... ")
        df_merged.to_csv('/Users/rbhattac/Desktop/Host_biosample_results.csv', index=False)


if __name__ == "__main__":
    main()