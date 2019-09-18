import time, operator, os, sys, argparse
from subprocess import Popen, PIPE
from library_v2 import *

#======================================= parser

parser = argparse.ArgumentParser()
parser.add_argument("LTR_seq", type=str, help="full LTR sequence (5' to 3')")
parser.add_argument("target_seq", type=str, help="target sequence with PAM (5' to 3')")
parser.add_argument("input1", type=str, help="NGS result file name 1")
parser.add_argument("-input2", type=str, help="NGS result file name 2")
parser.add_argument("output", type=str, help="output file name")
parser.add_argument("-r", type=int, help="length of expected genomic sequence", default = 70)
parser.add_argument("-i", type=int, help="length of indicator sequence", default = 15)

args=parser.parse_args()

LTR_seq = args.LTR_seq.upper()
target_seq = args.target_seq.upper()
f = open(args.input1).readlines()
if (args.input2):
    f2=open(args.input2).readlines()
    f=f+f2
flen=len(f)
out_file_name=args.output
cut_len=args.r
indi_len=args.i
exin_db={}

#+======================================

(sliced_LTR, indi_1, indi_2, indi_3, indi_con, sliced_len) = find_target(LTR_seq, target_seq, indi_len)
print('find target')

indi_2mm_0=mk_indi_1mm_list(indi_1, indi_len)
indi_2mm_1=mk_indi_1mm_list(indi_2, indi_len)
indi_2mm_2=mk_indi_1mm_list(indi_3, indi_len)
indi_2mm_con=mk_indi_1mm_list(indi_con, 5)
print('table')

(sort_dic, sorted_tup)=collect_seq(indi_2mm_0, indi_2mm_1, indi_2mm_2, f, flen, cut_len, indi_len, indi_2mm_con, sliced_len)
print('collect seq')

fw=open('bwa_fastmap_input.txt','w')
for x, y in sorted_tup:
    fw.write('>'+x+'\n'+x+'\n')
fw.close()

print('Run bwa_fastmap ...')
q=Popen(['bwa','fastmap','-l',str(cut_len-20),'Genome/GRCh38','bwa_fastmap_input.txt'], stdout=PIPE).communicate()

print('bulid human Database')
exin_db=DB_modeling()

fastmap_dic=read_fastmap_results(q)
chromomap_dic=analyze_fastmap_results(out_file_name, sorted_tup, fastmap_dic, exin_db, cut_len)
make_chromomap(out_file_name, chromomap_dic)
