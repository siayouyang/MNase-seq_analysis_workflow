
list = []
list_file = open('overlap_1108genes.txt', 'r')
bed = open('clark_genome_2014_TSS_TTS.bed', 'r')
overlap_bed = open('overlap_1108genes.bed','w')
for a in list_file:
    list.append(a.split()[0])
for b in bed:
    if b.split()[3] in list:
        overlap_bed.write(b)
list_file.close()
overlap_bed.close()
bed.close()
print('done')
