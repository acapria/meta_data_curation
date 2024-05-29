import pandas as pd



def check_labhost(genbank_labhost, biosample_labhost):
    try:
        # If both are available
        #check if  labhost fields are blank or has NaN values
        if pd.isnull(genbank_labhost) and pd.isnull(biosample_labhost):
            return "null"
        #if one of the fields is blank or has NaN values and other is not, return the non-blank value
        elif pd.isnull(genbank_labhost):
            return biosample_labhost
        elif pd.isnull(biosample_labhost):
            return genbank_labhost
        #if both fields have values, check if they are the same. if yes return the value, if not return "LabHost_Mismatch"
        elif genbank_labhost and biosample_labhost:
            if genbank_labhost == biosample_labhost:
                return genbank_labhost
            else:
                return "LabHost_Mismatch"

    except Exception as e:
        print("Error in check_labhost:", e)
        return None

def main():
    df = pd.read_csv("/Users/rbhattac/Desktop/Curation/Pipeline/abril_Code/results_code_Abril.csv")
    #print(df['Biosample Lab Host'])
    for row in df.iterrows():
        biosample_labhost = row[1]['Biosample Lab Host']
        #print(f'Biosample Lab Host: {biosample_labhost}')
        genbank_labhost = row[1]['Genbank Lab Host']
        #print(f'Genbank Lab Host: {genbank_labhost}')
        labhost = check_labhost(genbank_labhost, biosample_labhost)
        #print(labhost)

        if labhost == "null":
            print("LabHost is null")
        elif labhost == "LabHost_Mismatch":
            print("FLAG: LabHost Mismatch")
        else:
            print(f"LabHost: {labhost}")

if __name__ == "__main__":
    main()
