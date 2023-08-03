# author: siayouyang
# analysis of nucleosomes paired between two samples
# python Subnucleosome\SelectWindow\SelectWindow.py --input1 Subnucleosome\SelectWindow\test.txt --input2 Subnucleosome\SelectWindow\test2.txt -o Subnucleosome\SelectWindow\test_output -g Subnucleosome\SelectWindow\test_gene.txt -b 5 -a 10
# python Subnucleosome\SelectWindow\SelectWindow.py -i Subnucleosome\SelectWindow\WT1.merge.sam.randsel.sam.smooth.wig -o Subnucleosome\SelectWindow\WT1 -g Subnucleosome\SelectWindow\science_2008_ORF
from optparse import OptionParser
import multiprocessing
import os
import re
from collections import Counter
import numpy
import pandas as pd
import math
from scipy import stats




def run():
    ##OptionParser
    parser = OptionParser()
    parser.add_option("-i", dest="input2", help="(data in WIG format")
    parser.add_option("-p", "--processors", dest="processors", type="int", default=5, help="Maximum processors to be run in parallel, will not exceed the maximum number of chromosomes . Default %default.")
    parser.add_option("-o", dest="output1", help="output name for input1")
    parser.add_option("-g", "--gene", dest="gene", help="gene data in BED format")
    parser.add_option("-b", "--before", dest="before", type="int", default=200, help="length(bp) before reference point . Default %default.")
    parser.add_option("-a", "--after", dest="after", type="int", default=800, help="length(bp) after reference point . Default %default.")
    parser.add_option('--geneLength', action='store', type='int', dest='gene_len', default=0, help='gene length filter; outputs only genes with greater length. Default %default. ')

    (options, args) = parser.parse_args()


    global  input2, max_processors, output1, gene, before, after, gene_len
    input2 = options.input2
    max_processors = options.processors
    output1 = options.output1
    gene = options.gene
    before = options.before
    after = options.after
    gene_len = options.gene_len

    #get chromosomes name
    chr_set = set()
    input2_file = open(input2, 'r')
    for a in input2_file:
        chr_set.add(a.split()[0])
    #create chromosomes name tuple
    global chr_tuple
    chr_tuple = tuple(chr_set)
    del chr_set
    input2_file.close()


def split_chr_2(chr, input2, gene, before, after, output1, gene_len):
    input2_file = open(input2, 'r')
    input2_chr_temp = open(input2 + "." + chr + ".temp", 'w')
    for a in input2_file:
        if str(a.split()[0]) == str(chr):
            input2_chr_temp.write(a)
        else:
            pass
    input2_file.close()
    input2_chr_temp.close()


def split_plus1_chr(chr, input2, gene, before, after, output1,  gene_len):
    plus1_file = open(gene, 'r')
    plus1_chr_temp = open(output1 + ".plus1.txt" + "." + chr + ".temp", 'w')
    for a in plus1_file:
        if str(a.split()[0]) == str(chr) and abs(int(a.split()[2])-int(a.split()[1])) >= gene_len :
            plus1_chr_temp.write(a)
        else:
            pass
    plus1_file.close()
    plus1_chr_temp.close()


# plus1 as reference point
def select_plus1_window(chr, input2, gene, after, before, output1, gene_len):
    plus1_chr_temp = open(output1 + ".plus1.txt" + "." + chr + ".temp", 'r')
    input2_chr_temp = open(input2 + "." + chr + ".temp", 'r')
    output_chr_temp = open(output1 + ".plus1.txt" + "." + chr + ".temp2", 'w')
    occ_list = []
    for a in input2_chr_temp:
        occ_list.append(float(a.split()[3]))
    for b in plus1_chr_temp:
        if b.split()[5] == "+":
            plus1 = int(b.split()[1])
            if int(plus1 - before) < 0:
                window = list(([0.0] * abs(plus1 - before)) + occ_list[0:plus1 + after])
            elif int(plus1 + after) > len(occ_list):
                window = list(occ_list[plus1 - before:len(occ_list)] + ([0.0] * abs(plus1 + after - len(occ_list))))
            else:
                window = list(occ_list[plus1 - before:plus1 + after])
            output_chr_temp.write(f'{b.split()[0]}\t{b.split()[1]}\t{int(b.split()[2])+1}\t{b.split()[3]}\t{b.split()[4]}\t{b.split()[5]}\t{window}\n')
        elif b.split()[5] == "-":
            plus1 = int(b.split()[2])
            if int(plus1 - after) < 0:
                window = list(([0.0] * abs(plus1 - after)) + occ_list[0:plus1 + before])
            elif int(plus1 + before) > len(occ_list):
                window = list(occ_list[plus1 - after:len(occ_list)] + ([0.0] * abs(plus1 + before - len(occ_list))))
            else:
                window = list(occ_list[plus1-after:plus1+before])
            window=list(reversed(window))
            output_chr_temp.write(f'{b.split()[0]}\t{int(b.split()[1])-1}\t{b.split()[2]}\t{b.split()[3]}\t{b.split()[4]}\t{b.split()[5]}\t{window}\n')
    del occ_list
    plus1_chr_temp.close()
    input2_chr_temp.close()
    output_chr_temp.close()


def merge_plus1_output():
    output_file = open(output1 + ".plus1.genes.txt", 'w')
    for a in chr_tuple:
        output1_chr_temp = open(output1 + ".plus1.txt" + "." + a + ".temp2", 'r')
        output_file = open(output1 + ".plus1.genes.txt", 'a')
        for b in output1_chr_temp:
            output_file.write(b)
        output1_chr_temp.close()
    output_file.close()


def average_plus1_list(before, after):
    output_file = open(output1 + ".plus1.genes.txt", 'r')
    output_average_file = open(output1 + ".plus1.average.txt", 'w')
    gene_count = 0
    sum_list = [0] * int(after+before)
    for a in output_file:
        window_list = re.sub('[\[\]\n]','',a.split('\t')[6]).split(',')
        gene_count += 1
        for b in range(0, int(after+before)):
            index = int(b)
            window_index = window_list[index]
            sum_list[index] += float(window_index)
    average_list = [0] * int(after+before)
    for c in range(0, int(after + before)):
        index = int(c)
        average_list[index] = float(sum_list[index])/float(gene_count)
    xlabel = list(range(-before,after,1))
    output_average_file.write(f'{xlabel}\n{average_list}')
    del sum_list, average_list
    output_average_file.close()
    output_file.close()


def smoothing(before, after):
    plus1_average_smooth_file = open(output1 + ".plus1.average.smoothed.txt", 'w')
    plus1_average_file = open(output1 + ".plus1.average.txt", 'r')
    plus1_average = list(re.sub('[\[\]]','', plus1_average_file.readlines()[1]).split(','))
    plus1_average_raw = []
    for b in plus1_average:
        plus1_average_raw.append(float(b))
    # smooth 21bp span
    half_span = int(21/2)
    plus1_average_smooth = []
    for e in range(0, len(plus1_average_raw)):
        raw_list = []
        n = 0
        for f in range((int(e) - half_span), (int(e) + half_span) + 1):
            if f < 0:
                raw_list.append(float(0))
            elif f >= (int(len(plus1_average_raw))):
                raw_list.append(float(0))
            else:
                raw_list.append(float(plus1_average_raw[f]))
                n += 1
        smooth_value = float(sum(raw_list) / n)
        plus1_average_smooth.append(smooth_value)
    xlabel = list(range(-before,after,1))
    plus1_average_smooth_file.write(f'{xlabel}\n{plus1_average_smooth}')
    del plus1_average_smooth
    plus1_average_file.close()


def final_os_remove(chr, input2, gene, after, before, output1, gene_len):
    os.remove(input2 + "." + chr + ".temp")
    os.remove(output1 + ".plus1.txt" + "." + chr + ".temp")
    os.remove(output1 + ".plus1.txt" + "." + chr + ".temp2")


if __name__ == "__main__":
    run()

    def multi_processes(func):
        for i in range(0, len(chr_tuple), max_processors):
            if int(i + max_processors) > int(len(chr_tuple)):
                chr_tuple_chunk = chr_tuple[i:int(len(chr_tuple))]
            else:
                chr_tuple_chunk = chr_tuple[i:i + max_processors]
            processes = []
            for chr in chr_tuple_chunk:
                processes.append(multiprocessing.Process(target=func, args=(chr,input2,gene, after, before,output1, gene_len)))
            for process in processes:
                process.start()
            for process in processes:
                process.join()

    multi_processes(split_chr_2)

    #reference point=+1
    multi_processes(split_plus1_chr)
    multi_processes(select_plus1_window)
    merge_plus1_output()
    average_plus1_list(before, after)

    #calculate average spacing using regression analysis
    smoothing(before, after)

    multi_processes(final_os_remove)

    print("Done!")
