import time, operator, os, sys, argparse
from subprocess import Popen, PIPE

def rev_comp(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

class CrevisSeq:

    def __init__(self, args):

        self.LTR_seq = args.LTR_seq.upper()
        self.target_seq = args.target_seq.upper()
        self.input1 = args.input1
        self.input2 = args.input2
        self.output_name = args.output
        self.ref_dir = args.ref_dir
        self.cut_len = args.r
        self.indi_len = args.i
        self.exin_db = {}
        self.input_dict = {}
        self.indi_dict = {
            '0': [],
            '1': [],
            '2': [],
            'con': []
        }

    def find_target(self):

        print("Find appropriate position of target sequence in LTR sequence...")
        if self.LTR_seq.find(self.target_seq) != -1:
            self.sliced_seq = self.LTR_seq[self.LTR_seq.find(self.target_seq) + 16:]
            self.sliced_len = len(self.LTR_seq) - self.LTR_seq.find(self.target_seq) - 17
        elif self.LTR_seq.find(rev_comp(self.target_seq)) != -1:
            self.sliced_seq = self.LTR_seq[self.LTR_seq.find(rev_comp(self.target_seq)) + 5:]
            self.sliced_len = len(self.LTR_seq) - self.LTR_seq.find(rev_comp(self.target_seq)) - 6
        else:
            print("[Error] find target sequence !!")
            sys.exit()

        self.sliced_LTR = self.sliced_seq[1:]
        self.indi_0 = self.sliced_seq[:self.indi_len]
        self.indi_1 = self.sliced_seq[1:1 + self.indi_len]
        self.indi_2 = self.sliced_seq[2:2 + self.indi_len]
        self.indi_con = self.sliced_seq[2:7]

    def mk_indi_1mm_list(self):

        self.indi_dict['0'] = [self.indi_0]
        self.indi_dict['1'] = [self.indi_1]
        self.indi_dict['2'] = [self.indi_2]
        self.indi_dict['con'] = [self.indi_con]

        for i in range(self.indi_len):
            for nt1 in 'ATGC':
                self.indi_dict['0'].append(self.indi_0[:i] + nt1 + self.indi_0[i + 1:])
                self.indi_dict['1'].append(self.indi_1[:i] + nt1 + self.indi_1[i + 1:])
                self.indi_dict['2'].append(self.indi_2[:i] + nt1 + self.indi_2[i + 1:])

        for i in range(len(self.indi_con)):
            for nt1 in 'ATGC':
                self.indi_dict['con'].append(self.indi_con[:i] + nt1 + self.indi_con[i + 1:])

    def confirm_seq(self, input_name):

        amend_dic = {'0':1, '1':0, '2':2}
        t = time.time()
        pc = 0.1

        with open(input_name) as f:
            line_n = 0
            indi_cnt = 0

            for line in f:
                line_n += 1
                if line_n % 4 != 2: continue
                st = -1
                seq = line.strip()
                confirm_v = False

                for x in self.indi_dict['con']:
                    if seq.find(x) != -1:
                        confirm_v = True
                        break

                if confirm_v == False:
                    continue

                for n in ['1','0','2']:
                    for indi in self.indi_dict[n]:
                        st = seq.find(indi)
                        if st != -1:
                            indi_cnt += 1
                            amend = 1 - amend_dic[n] + self.sliced_len
                            subseq = seq[st + amend:st + amend + self.cut_len].replace('\n','')
                            if len(subseq) != self.cut_len: break
                            if subseq in self.input_dict.keys():
                                self.input_dict[subseq] += 1
                            else:
                                self.input_dict[subseq] = 1
                            break
                    if st != -1: break

        print("Find appropriate sequence : {0} of {1} in {2}".format(indi_cnt, int(line_n/4), input_name))

    def collect_seq(self):

        print("Read file 1...")
        self.confirm_seq(self.input1)
        if (self.input2):
            print("Read file 2...")
            self.confirm_seq(self.input2)

        print("Sorting sequence...")
        self.sorted_tup = sorted(self.input_dict.items(), key=operator.itemgetter(1), reverse = True)

    def run_bwa(self):

        with open('./bwa/input.txt','w') as fw:
            for x, y in self.sorted_tup:
                fw.write('>{0}\n{1}\n'.format(x,x))

        print('Run bwa fastmap...')
        self.fastmap_res = Popen(['bwa', 'fastmap', '-l', str(self.cut_len-20), self.ref_dir, './bwa/input.txt'], stdout=PIPE).communicate()

    def DB_modeling(self):

        print('Build human gtf3 database...')
        exin_db = {}
        lencnt = 0

        for i in list(range(1, 23))+['X', 'Y']:
            exin_db[str(i)] = {}
        with open('./DB/exon_intron_v2.txt') as f:
            for line in f:
                line_sp = line.strip().split()
                chro = line_sp[0]
                st = int(line_sp[7])
                ed = int(line_sp[8])
                a = st
                i_info = line_sp[0] + '\t' + line_sp[2] + '\t' + line_sp[3] + '\t' + line_sp[4] + '\t' + line_sp[5]
                if str(st)[:-3] == str(ed)[:-3]:
                    M = int(a / 1000000)
                    K = int(a / 1000)
                    if K > 1000: K = int(str(K)[-3:])
                    Bst = ('00' + str(st))[-3:]
                    Bed = ('00' + str(ed))[-3:]
                    if M not in exin_db[chro].keys():
                        exin_db[chro][M] = {}
                    if K not in exin_db[chro][M].keys():
                        exin_db[chro][M][K] = {}
                    exin_db[chro][M][K][Bst+';~;'+Bed+';'+str(lencnt)] = i_info
                else:
                    while a < ed + 1000:
                        M = int(a / 1000000)
                        K = int(a / 1000)
                        if K > 1000: K = int(str(K)[-3:])
                        B = ('00' + str(a))[-3:]
                        if M not in exin_db[chro].keys():
                            exin_db[chro][M] = {}
                        if K not in exin_db[chro][M].keys():
                            exin_db[chro][M][K] = {}
                        if a == st:
                            exin_db[chro][M][K][B + ';+;;' + str(lencnt)] = i_info
                        elif a >= ed:
                            M = int(ed / 1000000)
                            K = int(ed / 1000)
                            if K > 1000: K = int(str(K)[-3:])
                            B = ('00' + str(ed))[-3:]
                            exin_db[chro][M][K][B + ';-;;' + str(lencnt)] = i_info
                        else:
                            exin_db[chro][M][K][B + ';-;;' + str(lencnt)] = i_info
                        a += 1000
                        lencnt += 1

        self.exin_db = exin_db

    def read_fastmap_results(self):

        fastmap_list = str(self.fastmap_res[0]).split('//')

        fastmap_dict = {}
        chromomap_dict = {}
        id_chromomap = 0

        for info in fastmap_list[:-1]:
            info = info.split('\\n')
            print(info)
            if info[0] == '': info = info[1:]
            seq = info[0].split('\\t')[1]
            if len(info) == 2:
                fastmap_dict[seq] = 'Can not find position'
                continue
            sub_info_dict = {}
            for sub_info in info[1:-1]:
                sub_info_sp = sub_info.split('\\t')
                if sub_info_sp[1] == '0':
                    sub_info_dict[int(sub_info_sp[2]) - int(sub_info_sp[1])] = sub_info
            print(sub_info_dict)
            if len(sub_info_dict) == 0:
                fastmap_dict[seq] = 'Can not find position'
                continue
            info = sorted(sub_info_dict.items(), key=operator.itemgetter(0), reverse=True)[0][1]
            cut_len = int(info.split('\\t')[2])
            if info.split('\\t')[-1] == '*':
                fastmap_dict[seq] = 'fastmap find {0} sites; {1}'.format(info.split('\\t')[-2], seq[:cut_len])
            else:
                fastmap_dict[seq] = []
                chro_c = -1
                if len(info.split('\\t')[4:]) == 1:
                    i = info.split('\\t')[4].split(':')
                    chromomap_dict[str(id_chromomap)] = [i[0], str(abs(int(i[1]))), str(abs(int(i[1])) + int(info.split('\\t')[2]))]
                    id_chromomap += 1
                    chro_c = 1
                for i in info.split('\\t')[4:]:
                    i = i.split(':')
                    fastmap_dict[seq].append([i[0], i[1][:1], i[1][1:], seq[:cut_len], chro_c])

        self.fastmap_dict = fastmap_dict

    def analyze_fastmap_results(self):

        length_sorted_list = []
        for seq, info in self.fastmap_dict.items():
            yl = self.fastmap_dict[seq]
            print(yl)
            if type(yl) != str and len(yl[0][3]) == self.cut_len:
                length_sorted_list.append(yl[0][3])
        print(len(length_sorted_list))
        consolidate_dict = {}
        for x, cnt in self.sorted_tup:
            yl = self.fastmap_dict[x]
            if type(yl) != str and len(yl[0][3]) != self.cut_len:
                c = 0
                seq = yl[0][3]
                for i in length_sorted_list:
                    if i.find(seq) != -1:
                        print(seq)
                        print(yl)
                        if i in consolidate_dict.keys():
                            consolidate_dict[i] += cnt
                        else:
                            consolidate_dict[i] = cnt
                        c = 1
                        continue
                if c == 0:
                    consolidate_dict[x] = cnt
            else:
                if x in consolidate_dict.keys():
                    consolidate_dict[x] += cnt
                else:
                    consolidate_dict[x] = cnt

        chromomap_dict = {}
        chromomap_id = 0
        with open(self.output_name, 'w') as fw, open('{0}_chromomap.txt'.format(self.output_name[:self.output_name.rfind('.')]), 'w') as fwc:
            id_c = 0
            fw.write('ID\tsequence\tchromosome\tstart\tread\tgene\tclass\texon ID\n')
            for x, cnt in sorted(consolidate_dict.items(), key=operator.itemgetter(1), reverse=True):
                id_c += 1
                x = x.strip()
                yl = self.fastmap_dict[x]
                if type(yl) == str:
                    if yl.find(';') != -1:
                        yl = yl.split(';')
                        fw.write('{0}\t{1}\t\t{2}\t{3}\n'.format(id_c, y[1], cnt, yl[0]))
                    else:
                        fw.write('{0}\t{1}=t=t{2}\t{3}\n'.format(id_c, x, cnt, yl))
                    continue
                for y in yl:
                    chromomap_dict[str(chromomap_id)] = [y[0], int(y[2]), int(y[2])+len(y[3])]
                    chromomap_id += 1
                    fw.write('{0}\t{1}\t{2}\t{3}\t{4}'.format(id_c, y[3], y[0], y[2], cnt))
                    if y[4] == 1:
                        fwc.write('{0}\t{1}\t{2}\t{3}\t{4}'.format(id_c, y[3], y[0], y[2], cnt))
                    chro = y[0]
                    pos = int(y[2])
                    B = int(('00' + str(pos))[-3:])
                    M = int(pos/ 1000000)
                    K = int(('00'+str(int(pos / 1000)))[-3:])
                    try:
                        db = self.exin_db[chro][M][K]
                    except:
                        fw.write('\tNothing in database\n')
                        if y[4] == 1: fwc.write('\tNothing in database')
                        continue

                    vv = 0
                    sub_db = {}

                    for i in db.keys():
                        isp = i.split(';')
                        v = 0
                        if isp[1] == '=': v = 1
                        elif int(isp[0]) <= B and isp[1] == '+': v = 1
                        elif int(isp[0]) >= B and isp[1] == '-': v = 1
                        elif isp[1] == '~':
                            if int(isp[0]) <= B <= int(isp[3]): v = 1
                        if v == 1:
                            vv = 1
                            dbl = db[i].split('\t')
                            try:
                                ssub_db = sub_db[dbl[1]]
                                if dbl[2] not in ssub_db['chtype'].split('&'): ssub_db['chtype'] += '&' + dbl[2]
                                if ssub_db['exin'] == 'intron' and dbl[3] == 'exon': ssub_db['exin'] = 'exon'
                            except:
                                sub_db[dbl[1]] = {'chtype': dbl[2], 'exin': dbl[3]}

                    for gene, gene_list in sub_db.items():
                        fw.write('\t{0}\t{1}\t{2}'.format(gene, gene_list['chtype'], gene_list['exin']))
                        if y[4] == 1: fwc.write('\t{0}\t{1}\t{2}'.format(gene, gene_list['chtype'], gene_list['exin']))
                    if vv == 0:
                        fw.write('\tNothing in database\n')
                        if y[4] == 1: fwc.write('\tNothing in datbase\n')
                    else:
                        fw.write('\n')
                        if y[4] == 1: fwc.write('\n')

        self.chromomap_dict = chromomap_dict

    def make_chromomap(self):

        with open('{0}_chromomap.html'.format(self.output_name[:self.output_name.rfind('.')]),'w') as fw:
            s = '<head>\n\t<script src="https://unpkg.com/ideogram@1.4.1/dist/js/ideogram.min.js"></script>\n</head>'
            s += "<body>\n<script>\n\tvar ideogram = new Ideogram({ organism: 'human', chrMargin: 5, annotationsLayout: 'histogram', annotations: ["
            for x, y in self.chromomap_dict.items():
                s += "{ name: '%s', chr: '%s', start: %s, stop: %s }," %(x, y[0], y[1], y[2])
            s = s[:-1] + "] });\n</script>\n</body>"
            fw.write(s)



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("LTR_seq", type=str, help="full LTR sequence (5' to 3')")
    parser.add_argument("target_seq", type=str, help="target sequence with PAM (5' to 3')")
    parser.add_argument("input1", type=str, help="NGS result file name 1")
    parser.add_argument("-input2", type=str, help="NGS result file name 2")
    parser.add_argument("output", type=str, help="output file name")
    parser.add_argument("-r", type=int, help="length of expected genomic sequence", default=70)
    parser.add_argument("-i", type=int, help="length of indicator sequence", default=15)
    parser.add_argument("ref_dir", type=str, help="Directory of bwa index")

    return parser.parse_args()

def main():

    tt = time.time()

    c = CrevisSeq(parse_args())

    c.find_target()
    c.mk_indi_1mm_list()
    c.collect_seq()
    c.run_bwa()
    c.DB_modeling()
    c.read_fastmap_results()
    c.analyze_fastmap_results()
    c.make_chromomap()

if __name__ == '__main__':
    main()
