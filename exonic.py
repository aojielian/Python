import pandas as pd
from os.path import join, isfile
from os import listdir


def get_files(func_name, exonic_func_ref_gene, input_path, output_path, pli_input_path):
    my_list = listdir(input_path)  #define my_list
    for file_name in my_list:  # for loop
        if isfile(join(input_path, file_name)):
            if file_name.lower().endswith('.csv'):  #check if file is csv
                output_file_name = file_name + '_' + 'output' + '.csv'
                get_exonic_and_stopgain(func_name, exonic_func_ref_gene, join(input_path, file_name), join(output_path, output_file_name), pli_input_path)


def get_exonic_and_stopgain(func_name, exonic_func_ref_gene, input_path, output_path, pli_input_path):
    gene_dict = pli_score_dict(pli_input_path)
    data_list = []
    data = pd.read_csv(input_path, delimiter=',', header=None, low_memory=False)
    for row_tuple in data.iterrows():
        row_index = row_tuple[0]
        row_value = row_tuple[1]
        func_ref_gene_column = row_value[5] # column 5
        exonic_func_ref_gene_column = row_value[8] # column 8
        if func_ref_gene_column == func_name and exonic_func_ref_gene_column == exonic_func_ref_gene:
            gene_names = row_value[6] # column 6
            pli_string = match_gene_pli(gene_dict, gene_names)
            gene_total_value = (row_value, pli_string)
            data_list.append(gene_total_value)
    write_organized_file_to_csv(data_list, output_path)

def match_gene_pli(gene_dict, gene_names):
    pli_string = ""
    name_list = gene_names.split(';')
    for name in name_list:
        gene_pli = gene_dict.get(name, 'not exist')
        if pli_string == "":
            pli_string = str(gene_pli)
        else:
            pli_string = pli_string + ';' + str(gene_pli)
    return pli_string

def pli_score_dict(pli_input_path):
    gene_dict = {}
    data = pd.read_csv(pli_input_path, delimiter='\t')
    for row_tuple in data.iterrows():
        row_index = row_tuple[0]
        row_value = row_tuple[1]
        gene_dict[row_value[0]] = row_value[1]
    return gene_dict

def write_organized_file_to_csv(data_list, output_path):
    heads = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene','1000g2015aug_all', 'pLI']
    with open(output_path, 'w') as csv_file_writer:
        for head in heads:
            csv_file_writer.write(head)
            csv_file_writer.write(',')
        csv_file_writer.write('\n')
        for row_value_tuple in data_list:
            row_value = row_value_tuple[0]
            pli_string = row_value_tuple[1]
            csv_file_writer.write(str(row_value[0]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[1]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[2]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[3]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[4]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[5]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[6]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[8]))
            csv_file_writer.write(',')
            csv_file_writer.write(str(row_value[55]))
            csv_file_writer.write(',')
            csv_file_writer.write(pli_string)
            csv_file_writer.write('\n')




input_path = '/oasis/projects/nsf/ddp195/a1lian/annovar/ssc_snv-VCF-annotated'
output_path = '/oasis/projects/nsf/ddp195/a1lian/annovar/ssc-exonicstopgain'
pli_input_path = '/oasis/projects/nsf/ddp195/a1lian/annovar/pli_genes.txt'
func_name = "exonic"
exonic_func_ref_gene = "stopgain"
get_files(func_name, exonic_func_ref_gene, input_path, output_path, pli_input_path)
