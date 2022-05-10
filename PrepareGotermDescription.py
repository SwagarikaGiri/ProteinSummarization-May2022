import pickle
import numpy as np
import pandas as pd
import csv
from utils import gene_ontology
import string
import re

gene_onto = gene_ontology()
def load_uniprot_data():
    df = pd.read_csv('UniprotData2022.tab',sep='\t')
    df = df[['Entry','Gene ontology IDs','Function [CC]',]]
    df = df.rename(columns={"Entry": "Protein", "Gene ontology IDs":"GOTerm_IDs","Function [CC]":"Uniprot_Desc"})
    df = df.set_index('Protein')
    return df

def clean_description_string(string_):
    string_ = string_.replace('"', '')
    string_ = string_.replace('GOC:ai', '')
    return string_




def get_goterm_defination(go_term_list):
    doc_defination=""
    errored_go_term=[]
    for go in go_term_list:
        go = go.strip()
        try:
            go_details = gene_onto[go]
            print(go_details)
            defination = clean_description_string(go_details['def'])
    
            doc_defination=doc_defination+defination
        except:
            errored_go_term.append(go)
    return doc_defination

def main():
    with open('ProteinAndItsDescriptionCombinedTab.csv','w' ,newline='', encoding='utf-8') as output_csvfile:
        spamwriter = csv.writer(output_csvfile, delimiter='\t')
        list_=[]
        list_.extend(['Protein','Goterms', 'GO Description','Uniprot Description'])
        spamwriter.writerow(list_)
        df = load_uniprot_data()
        for index,row in df.iterrows():
            print(index)
            gos=row['GOTerm_IDs']
            goaDefination= get_goterm_defination(gos.split(";"))
            uniProtDefination= row['Uniprot_Desc']
            spamwriter.writerow([index,gos,goaDefination,uniProtDefination])



if __name__ == "__main__":
    main()