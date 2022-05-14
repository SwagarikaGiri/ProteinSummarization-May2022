from ast import Pass
from gc import get_objects
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
    return string_

def seperate_goterm_based_on_ontology(goterm_list,onto):
    final_goterm_list=[]
    for go in goterm_list:
        go = go.strip()
        try:
            go_details= gene_onto[go]
            if(go_details['namespace']==onto):
                final_goterm_list.append(go)
            else:
                print(go)
                print(go_details['namespace'])
        except:
            print("Did exception happened?")
    return final_goterm_list




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


def prepare_gos(goterm_list):
    return ';'.join(goterm_list)



def prepare_dataset(filename,onto):
    with open(filename+str(onto)+".tab",'w' ,newline='', encoding='utf-8') as output_csvfile:
        spamwriter = csv.writer(output_csvfile, delimiter='\t')
        list_=[]
        list_.extend(['Protein','Goterms', 'GO Description'])
        spamwriter.writerow(list_)
        df = load_uniprot_data()
        for index,row in df.iterrows():
            print(index)
            gos=row['GOTerm_IDs']
            goterm_list=seperate_goterm_based_on_ontology(gos.split(";"),onto)
            print(goterm_list)
            if(len(goterm_list)>0):
                goaDefination= get_goterm_defination(goterm_list)
                # uniProtDefination= row['Uniprot_Desc']
                spamwriter.writerow([index,prepare_gos(goterm_list),goaDefination])

def main():
    prepare_dataset('ProteinDesctiption',"biological_process")
    prepare_dataset('ProteinDesctiption',"cellular_component")
    prepare_dataset('ProteinDesctiption',"molecular_function")

    
            



if __name__ == "__main__":
    main()