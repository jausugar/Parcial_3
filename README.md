# Parcial_3

Aquí se encuentra el parcial #3 de programación para ciencias biológicas.
Esta rutina puede identificar si la proteína de interés se encuentra en el organismo:
 - Helicobacter pylori strain 7.13_R1b chromosome con código de GenBank NZ_CP024072

# Intrucciones
 Posiciones Todos los archivos que se encuentran aquí junto con la(s) proteína(s) que quiere identificar
 Luego, corra el siguiente comando:

    python parcial_31.py "proteina".fasta  

Esta rutina generará 4 archivos:
 - revcom.fasta, el cual contiene un archivo fasta con las secuencias reveso complementario de la secuencia ingresada.
 - SequenceAA_query.fasta, que contiene los ORF's encontrados en aminoácidos.
 - salida_blas.tsv, que contiene todos los resultados obtenidos del blast.
 - Resultados_blast_a_Hpylori.tsv, el cual es el resultado final con las secuencias identificadas.  
