import os
import sys
import pandas as pd

print("\nEste programa está diseñado para traducir ADN a proteínas y realizar  "
      "un Blast contra una base de datos de Helicobacter Pylori ")

'''Abrir el archivo fasta que se ingresa desde la terminal,
 ignorar la línea del título y eliminar los saltos de línea'''

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

'''########----Verificación de secuencias ingresadas----########'''
# Esta función verifica que la secuencia ingresada no sea
# correspondiente a aminoácidos

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

# Esta función se encarga de depurar la secuencia original
# en caso de que venga con caracteres desconocidos para que
# puedan ser procesados en pasos posteriores
def depuration(sequence):
    deputation_dict = {'A': 'A', 'C': 'C', 'T': 'T', 'G': 'G',
                       'W': 'X', 'S': 'X', 'R': 'R', 'Y': 'X',
                       'K': 'X', 'M': 'X', 'B': 'X', 'D': 'X',
                       'H': 'X', 'V': 'X', 'N': 'X', 'E': 'X',
                       'U': 'X', 'Q': 'X', '-': 'X', 'X': 'X'}
    sequence = sequence.upper()
    sequence = sequence.replace("U", "T")
    seq = ""
    for nucleotide in sequence:
        # revisar que los caracteres en la secuencias se encuentren
        # en el diccionario, de lo contrario el proceso se detiene
        if nucleotide in deputation_dict:
            seq += deputation_dict[nucleotide]
        else:
            print(f'Error: la secuencia tiene un chatarter {nucleotide} desconocido, '
                  f'por favor, revise que {nucleotide} sea una base de nucleotidos')
            exit()
    return seq

# Esta función regresa el reverson complementario de la secuencia
# original y también verifica y depura caracteres especiales
def revcom(sequence):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                       'W': 'X', 'S': 'X', 'R': 'R', 'Y': 'X',
                       'K': 'X', 'M': 'X', 'B': 'X', 'D': 'X',
                       'H': 'X', 'V': 'X', 'N': 'X', 'E': 'X',
                       'U': 'X', 'Q': 'X', '-': 'X', 'X': 'X'}

dict = {'W': 'X', 'S': 'X', 'R': 'R', 'Y': 'X',
        'K': 'X', 'M': 'X', 'B': 'X', 'D': 'X',
        'H': 'X', 'V': 'X', 'N': 'X', 'E': 'X',
        'U': 'X', 'Q': 'X', '-': 'X', 'X': 'X'}

    sequence = sequence.upper()
    sequence = sequence.replace("U", "T")
    sequence = sequence[::-1]
    reverse_complement = ""
    for nucleotide in sequence:
        if nucleotide in complement_dict:
            reverse_complement += complement_dict[nucleotide]
        else:
            print(f'Error: la secuencia tiene un chatarter {nucleotide} desconocido, '
                  f'por favor, revise que {nucleotide} sea una base de nucleotidos')
            exit()
    return reverse_complement

# acá se crea un archivo revcom.fasta en donde se escriben
# los reversos complementarios de la secuencia original
# en este punto también se verifica la función checkprot
archivo = open("revcom.fasta", "w")
for id, seq in zip(ids, secuencia):
    checkprot(seq, id)
    reverso = revcom(seq)
    archivo.write(">" + id + "\n" + reverso + "\n")
archivo.close()

'''########----TRADUCCIÓN----########'''
# Esta función toma una secuencia y devuelve su traducción
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
    trans_seq = ""
    indeter = "X"

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        # En caso de que haya algún nucleótido en X
        # también se traducirá en X
        if indeter in codon:
            trans_seq += "X"
        else:
            trans_seq += genetic_code[codon]
    return trans_seq

'''########----BÚSQUEDA DE ORF'S----########'''
''' el siguiente loop toma encuentra los ORFs en los
6 marcos de lectura para el número de secuencias 
ingresadas '''

# Diccionario con los marcos de lectura
orf_dict = {'1': 0, '2': 1, '3': 2, '4': 0, '5': 1, '6': 2}
df_orf = pd.DataFrame(columns=['ID', 'ORF', 'length', 'seq_AA'])
min_length = 150
# Para cada fasta
for id, secuencia in zip(ids, secuencia):
    sequence1 = depuration(secuencia)
    sequence = ''

    # Para cada marco de lectura
    for orf, pos in orf_dict.items():
        a = 1
        seq_orf = ''

        # verifica si trabaja con secuencia original
        # o reverso complementario
        if int(orf) > 3:
            sequence = revcom(sequence1)
        else:
            sequence = sequence1

        # Recorrer la secuecia
        for x in range(pos, len(sequence) - 2, 3):
            codon = sequence[x:x + 3]

            # Buscar codon de inicio o verifica
            # si se encuentra dentro de un ORF
            if codon == 'ATG' or a == 2:
                a = 2
                seq_orf += codon

                # Buscar el codon de parada
                if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':
                    a = 1

                    # Solo admite una long mínima de 150
                    if len(seq_orf) > min_length:
                        seq_AA = translate(seq_orf)
                        df_orf = df_orf.append({'ID': id, 'ORF': orf, 'length': len(seq_orf),
                                                'seq_AA': seq_AA},
                                               ignore_index=True)
                    seq_orf = ''


# Se organiza el dataframe de mayor a menor longitud
df_orf = df_orf.sort_values(by="length",
                            ascending=False)
print(f'\n Esta tabla muestra los ORFs identificados, su'
      f'longitud, ID y el marco en el que se hallaron')
print(df_orf)

#Guardar archivos
archivo = open("SequenceAA_query.fasta", "w")
print(f'\n Se guardó los ORFs encontrados en un archivo llamado SequenceAA_query.fasta')
for i in range(0, len(df_orf)):
    archivo.write(">" + df_orf.iloc[i,0] + f"_ORF{df_orf.iloc[i,1]}"
                  + "\n" + df_orf.iloc[i,3] + "\n")
archivo.close()

# Blast con ORF's encontrados
print(f'\nSe está ejecutando blast')
comando_blast = f'blastp -db Helicobacter_pylori_prot -query SequenceAA_query.fasta -outfmt ' \
                f'"6 qseqid pident evalue bitscore stitle qlen" -max_target_seqs 5 -out salida_blast.tsv'
os.system(comando_blast)

column_names=['Seq_Query', 'Perc_ident', 'E-value', 'bitscore', 'protein', 'Query_len']
df_blast = pd.read_csv("salida_blast.tsv", sep='\t', names=column_names)
df_blast['protein']=df_blast['protein'].str.extract(r'.protein=([^\]]+).')
df_blast = df_blast[df_blast['E-value'] < 1E-100]
df_blast = df_blast[df_blast['bitscore'] > 100]
df_blast = df_blast[df_blast['Perc_ident'] > 80]
print(df_blast)

df_blast.to_csv('Resultados_blast_a_Hpylori.tsv', index=False, header=True, sep='\t')

