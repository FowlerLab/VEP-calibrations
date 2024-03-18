#!/usr/bin/env python
# coding: utf-8

# In[12]:


#load annovar annotated csv into a dataframe in pandas
import pandas as pd

clinvar_23 = pd.read_csv("~/Downloads/clinvar_2023_annovar_37.hg19_multianno.txt", sep = '\t', header = 0)


# In[13]:


## use the parser function to store what things you want added to the original dataframe

Clnsig = []
clnvc = []
Gene = []
gene_id = []
variant_type = []
star = []
temp_list = []
temp_list2 = []

def parse (file1):
    with open(file1) as my_file:
        for line in my_file:
            reader = line.split("\t")
            significance = str(reader[145]).find("CLNSIG=")
            var_type = str(reader[145]).find("CLNVC=")
            clnvcso = str(reader[145]).find("CLNVCSO")
            origin = str(reader[145]).find("ORIGIN")             
            gene = str(reader[145]).find("GENEINFO")
            mc = str(reader[145]).find("MC=")
            one_star = str(reader[145].find("CLNREVSTAT="))
            star.append(str(reader[145][int(one_star)+11:int(significance)-1]))
            
    
            temp_list2.append(reader[145][int(mc): int(origin)])
            temp_list.append(reader[145][int(gene):int(mc)])
            Clnsig.append((reader[145][int(significance)+7:int(var_type)-1]))
            clnvc.append((reader[145][int(var_type)+6:int(clnvcso)-1]))

            
        for element in temp_list:
            colon = element.find(":")
            semicol = element.find(';')
            vertical = element.find('|')
            if int(vertical) != -1:
                gene_id.append(element[int(colon)+1:int(vertical)])
            else:
                gene_id.append(element[int(colon)+1:int(semicol)])
            Gene.append(element[9:int(colon)])
                
        for element in temp_list2:
            dash = element.find("|")
            semi = element.find(";")
            variant_type.append(element[int(dash)+1:int(semi)])
            
        print('parser finished')
    
parse("clinvar_2023_annovar_37.hg19_multianno.txt")


# In[14]:


#delete the first instance of each one of these lists, its an empty space
star.pop(0)
gene_id.pop(0)
Gene.pop(0)
Clnsig.pop(0)
clnvc.pop(0)
variant_type.pop(0)


# In[15]:


#add new columns to the dataframe from the parsed information

clinvar_23['Stars'] = star
clinvar_23['Gene'] = Gene
clinvar_23['GeneID'] = gene_id
clinvar_23['Significance'] = Clnsig
clinvar_23['clnvc'] = clnvc
clinvar_23['Variant_type'] = variant_type


# In[16]:


#rename columns to be more precise

clinvar_23 = clinvar_23.rename(columns={'AF.1': 'exome_AF', 'AF': 'genome_AF'})


# In[1]:


#save parsed file for future use

clinvar_23.to_csv(path_or_buf="~/Downloads/clinvar_2023_annovar_37_parsed_with_geneid.csv")


# In[2]:


import pandas as pd
import os.path
import csv


clinvar23_parsed = pd.read_csv("~/Downloads/clinvar_2023_annovar_37_parsed_with_geneid.csv", sep = ",", header = 0,dtype = {"Chr": str , "Start" : int, "End": int,"Ref": str, "Alt": str,"Otherinfo1": str})


# In[2]:


#create a file with only single nucleotide variants, that are missense variants and VUS (variants of uncertain significance)

dfx = clinvar23_parsed.query("`clnvc` == 'single_nucleotide_variant'")


#filter on all missense variants - missense variants have some extra words in the clinvar file so any variant type with missense in it

dfy = dfx[dfx['Variant_type'].str.contains('missense', regex = False, na = False)]


# make sure the variants are 1+ stars

dfz = dfy.query("`Stars` != 'no_assertion_criteria_provided' and `Stars` != 'no_assertion_provided'")

#write csv for a file containing any missense, one star variants 

dfz.to_csv(path_or_buf="~/Downloads/clinvar23_parsed_missense_one_star.csv")

#find any vairants that have uncertain or conflicting significance
dfa = dfz[dfz['Significance'].str.contains('Uncertain|Conflicting', regex = True, na=False)]

#write SNV, missense, VUS variants to a csv file 
dfa.to_csv(path_or_buf="~/Downloads/clinvar23_parsed_only_VUS.csv")


# In[3]:


#create file with only single nucleotide variants with greater than one star rating and missense variants

# filter dataframe based on SNV's and variants above one star
df2 = clinvar23_parsed.query("`clnvc` == 'single_nucleotide_variant'")


#filter on all missense variants - missense variants have some extra words in the clinvar file so any varianty type with missense in it
df3 = df2[df2['Variant_type'].str.contains('missense', regex = False, na = False)]

#filter on 1+ star
df4 = df3.query("`Stars` != 'no_assertion_criteria_provided' and `Stars` != 'no_assertion_provided'")

#filter on non-VUS
df5 = df4.query("`Significance` == 'Benign' or `Significance` == 'Benign/Likely_benign' or `Significance` == 'Likely_benign' or `Significance` == 'Pathogenic' or `Significance` == 'Pathogenic/Likely_pathogenic' or `Significance` == 'Likely_pathogenic'")

#write file to csv
df5.to_csv(path_or_buf="~/Downloads/clinvar_2023_annovar_37_parsed_one_star_missesne_no_vus_with_geneid.csv")


# In[4]:


#create a function to flag genes with no pathogenic variants

no_path = []


def no_pathogenic(df):
    
    if len(df) == 0:
        return

    if len(df) > 0:
            
        try: 

            pathogenic = df.query("`Significance` == 'Pathogenic'")

            #Collect all pathogenic/likely pathogenic control variants

            pathogenic_lp = df.query("`Significance` == 'Pathogenic/Likely_pathogenic'")

            #Collect all likely pathogenic control variants from the ClinVar file

            likely_p =  df.query("`Significance` == 'Likely_pathogenic'")

            #Collect all benign control variants from the ClinVar file.

            benign =  df.query("`Significance` == 'Benign'")

            benign_lb = df.query("`Significance` == 'Benign/Likely_benign'")

            #Collect all likely benign control variants from the ClinVar file. 
            likely_b =  df.query("`Significance` == 'Likely_benign'")

            #Path/Likely Path

            #Count pathogenic controls
            pathogenic_controls = len(pathogenic) + len(likely_p) + len(pathogenic_lp)

            #Count benign controls 
            benign_controls = len(benign) + len(likely_b) + len(benign_lb)
            
            if pathogenic_controls == 0:
                gene = df.iloc[0,148]
                no_path.append(gene)
                
        except ZeroDivisionError:
            pass
                
                
grouping = clinvar23_parsed.groupby(["Gene"]).apply(no_pathogenic)


# In[5]:


#take out any genes with no pathogenic variants 
df6 = df5.query("Gene != @no_path")

#write file to csv
df6.to_csv(path_or_buf="~/Downloads/clinvar_2023_annovar_37_parsed_one_star_missesne_no_vus_no_genes_w_no_path_fixed_w_gene_id.csv")


# In[ ]:




