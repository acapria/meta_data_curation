import pandas as pd

def check_isolation(genbank_isolation, biosample_isolation):
    try:
        # If both are available
        #check if  isolation fields are blank or has NaN values
        if pd.isnull(genbank_isolation) and pd.isnull(biosample_isolation):
            return "null"

        #if one of the fields is blank or has NaN values and other is not, return the non-blank value
        elif pd.isnull(genbank_isolation):
            return biosample_isolation
        elif pd.isnull(biosample_isolation):
            return genbank_isolation
        #if both fields have values, check if they are the same. if yes return the value, if not return "Isolation_Mismatch"
        elif genbank_isolation and biosample_isolation:
            if genbank_isolation == biosample_isolation:
                return genbank_isolation
            #check if the genbank isolation source is a substring of biosample isolation source
            elif genbank_isolation in biosample_isolation:
                return biosample_isolation
            #check if the biosample isolation source is a substring of genbank isolation source
            elif biosample_isolation in genbank_isolation:
                return genbank_isolation
            else:
                return "Isolation_Mismatch"


    except Exception as e:
        print("Error in check_isolation:", e)
        return None

def main():
    df = pd.read_csv("/Users/rbhattac/Desktop/Curation/Pipeline/abril_Code/results_code_Abril.csv")
    #print(df['Genbank Isolation Source'])

    for row in df.iterrows():
        biosample_isolation = row[1]['Biosample Isolation Source']
        print(f'Biosample Isolation Source: {biosample_isolation}')
        genbank_isolation = row[1]['Genbank Isolation Source']
        print(f'Genbank Isolation Source: {genbank_isolation}')
        isolation = check_isolation(genbank_isolation, biosample_isolation)
        print(isolation)

        # if isolation == "null":
        #     print("Isolation is null")
        # elif isolation == "Isolation_Mismatch":
        #     print("FLAG: Isolation Mismatch")
        # else:
        #     print(f"Isolation: {isolation}")

if __name__ == "__main__":
    main()

