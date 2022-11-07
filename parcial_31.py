import os
import sys
import pandas as pd

print("\nEste programa está diseñado para traducir ADN a proteínas y " \
      "realizar un Blast contra una base de datos de E. coli")

'''Abrir el archivo fasta que se ingresa desde la terminal,
 ignorar la línea del título y eliminar los saltos de línea'''

#Cambio

ids = []
secuencia = []
string_ = ()
n = -1
num_seq = 0
with open(sys.argv[1], 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            id = line.split()[0][1:]  # para separar el ID fasta de los comentarios
            ids.append(id)  # agrega todos los ID's en una lista llamada ids
            secuencia.append('')
            n += 1
            num_seq += 1
        else:
            secuencia[n] += line



def checkprot(sequence, id):
    prot_dict = {'D': 'X', 'E': 'X', 'R': 'X', 'K': 'X',
                 'N': 'X', 'H': 'X', 'Q': 'X', 'S': 'X',
                 'T': 'A', 'A': 'T', 'G': 'C', 'V': 'X',
                 'P': 'X', 'L': 'X', 'F': 'X', 'Y': 'X',
                 'I': 'X', 'M': 'X', 'W': 'X', 'C': 'G'}
    sequence = sequence.upper()  # Para pasar a mayuscula
    protein_check = ""
    perc = 30  # Procentaje de posibles aminoacidos
    for x in sequence:
        if x in prot_dict:
            protein_check += prot_dict[x]
        else:
            continue

    if (float(protein_check.count("X")) / float(len(protein_check))) * 100 > perc:
        print(f'La secuencia {id} posiblemete corresponde a una secuencia de '
              f'aminoácidos, por favor verifique y vuelva a ingresarla')
        exit()

    return


def depuration(sequence):
    sequence = sequence.upper()  # Para pasar a mayuscula
    sequence = sequence.replace("U", "T")
    return sequence


def revcom(sequence):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                       'W': 'X', 'S': 'X', 'R': 'R', 'Y': 'X',
                       'K': 'X', 'M': 'X', 'B': 'X', 'D': 'X',
                       'H': 'X', 'V': 'X', 'N': 'X', 'E': 'X',
                       'U': 'X', 'Q': 'X', '-': '-'}  # Con indeterminaciones

    sequence = sequence.upper()  # Para pasar a mayuscula
    sequence = sequence.replace("U", "T")
    sequence = sequence[::-1]  # Sacar el reverso
    reverse_complement = ""
    for nucleotide in sequence:
        if nucleotide in complement_dict:  # revisar los caracteres en mi secuencias
            reverse_complement += complement_dict[nucleotide]
        else:
            print(f'Error: la secuencia tiene un chatarter {nucleotide} desconocido, '
                  f'por favor, revise que {nucleotide} sea una base de nucleotidos')
            exit()

    return reverse_complement


archivo = open("revcom.fasta", "w")
for id, seq in zip(ids, secuencia):
    # print("\nLa secuencia problema es:")
    # print("\n"'>' + id)
    # print(secuencia)
    checkprot(seq, id)
    reverso = revcom(seq)
    # print(f'\nEl reverso complementario de la secuencia problema es:'
    #      f'\n>{id}'
    #      f'\n{reverso}')
    archivo.write(">" + id + "\n" + reverso + "\n")
archivo.close()


def translate(sequence):
    genetic_code = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

    print("\nIniciando traducción:")
    AA_seqs = {}
    trans_seq = ""
    indeter = "X"

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if indeter in codon:
            trans_seq += "X"
        else:
            trans_seq += genetic_code[codon]
    #AA_seqs[id] = trans_seq

    return trans_seq


orf_dict = {'1': 0, '2': 1, '3': 2, '4': 0, '5': 1, '6': 2}
df_orf = pd.DataFrame(columns=['ID', 'ORF', 'length', 'seq_AA'])
min_length = 100
for id, secuencia in zip(ids, secuencia):  # Para cada fasta
    # longitud = len(sequence)
    sequence1 = depuration(secuencia)
    sequence = ''
    #print(sequence)
    for orf, pos in orf_dict.items():  # Hacer un marco
        a = 1
        seq_orf = ''
        #print(orf,pos)
        if int(orf) > 3:
            sequence = revcom(sequence1)
        else:
            sequence = sequence1

       #print(sequence)

        for x in range(pos, len(sequence) - 2, 3):  # recorrer la secuecia
            codon = sequence[x:x + 3]
            #print(codon, orf, x, a)
            if codon == 'ATG' or a == 2:
                a = 2
                seq_orf += codon
                if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':  # or x == (len(sequence) - 3 - pos):
                    a = 1
                    if len(seq_orf) > min_length:
                        seq_AA = translate(seq_orf)
                        df_orf = df_orf.append({'ID': id, 'ORF': orf, 'length': len(seq_orf),
                                                'seq_AA': seq_AA},
                                               ignore_index=True)
                    seq_orf = ''

# Se organiza el dataframe de mayor a menor longitud
df_orf = df_orf.sort_values(by="length",
                            ascending=False)

print(df_orf)

#Guardar archivos
archivo = open("sequence_query.fasta", "w")
for i in range(0, len(df_orf)):
    archivo.write(">" + df_orf.iloc[i,0] + f"_ORF{df_orf.iloc[i,1]}" + "\n" + df_orf.iloc[i,3] + "\n")
archivo.close()

# Blast con ORF's encontrados
print(f'\nSe está ejecutando blast')
comando_blast = f'blastp -db proteinas_ecoli.faa.txt -query sequence_query.fasta -outfmt "6 qseqid pident qlen slen qcovs"'
os.system(comando_blast)

# # Blast con secuencia normal
#
# print(f'\nSe está ejecutando blast con la siguiente secuencia (original): {seq_traducida}')
# comando_blast = f'blastp -db proteinas_ecoli.faa.txt -query seq_traducida.fasta -outfmt "6 qseqid pident qlen slen qcovs"'
# os.system(comando_blast)
