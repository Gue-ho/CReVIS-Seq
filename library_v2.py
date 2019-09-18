import operator, sys, time

def rev_comp(seq):
    return seq.translate(seq.maketrans('ATGC','TACG'))[::-1]

def progress_bar (i, total, str_f = '', str_ed = '', finish_str=''):
    per = "{0:.2f}".format(i*100/total)
    fill_len = int(i*50/total)
    bar = '█' * fill_len + '-' * (50-fill_len)
    print('\r%s |%s| %s %s\r' % (str_f, bar, per, str_ed), end='\r')
    if i == total: print(finish_str)

def find_target(LTR_seq, target_seq, indi_len):

    if LTR_seq.find(target_seq)!=-1:
        sliced_seq = LTR_seq[LTR_seq.find(target_seq)+16:]
        sliced_len = len(LTR_seq)-LTR_seq.find(target_seq)-17
    elif LTR_seq.find(rev_comp(target_seq))!=-1:
        sliced_seq = LTR_seq[LTR_seq.find(rev_comp(target_seq))+5:]
        sliced_len = len(LTR_seq)-LTR_seq.find(rev_comp(target_seq))-6
    else:
        print("target sequence is no located in LTR sequence!!!")
        sys.exit()

    return sliced_seq[1:], sliced_seq[:indi_len], sliced_seq[1:1+indi_len], sliced_seq[2:2+indi_len], sliced_seq[2:7], sliced_len

def mk_indi_2mm_list(indi, indi_len):


    l=[indi]

    for i in range(indi_len):
        for nt1 in 'ATGC':
            for x in range(i+1, indi_len):
                for nt2 in 'ATGC':
                    l.append(indi[:i]+nt1+indi[i+1:x]+nt2+indi[x+1:])
                   
    return l

def mk_indi_1mm_list(indi, indi_len):

    l=[indi]

    for i in range(indi_len):
        for nt1 in 'ATGC':
            l.append(indi[:i]+nt1+indi[i+1:])

    return l

def collect_seq(indi_2mm_0, indi_2mm_1, indi_2mm_2, f, flen, cut_len, indi_len, indi_2mm_con, sliced_len):
    
    dic={}
    flen=int(flen/4)
    indil_len=len(indi_2mm_0)
    indi_l=indi_2mm_1+indi_2mm_0+indi_2mm_2
    c=0
    amend_dic={0: 1, 1: 0, 2: 2}
    
    t=time.time()
    progress_bar(0, 100, 'collecting ...', '', 'data collected')
    pc=0.1
    for i in range(flen):

        seq=f[4*i+1].replace('\t','')
        st=-1
        
        v_confirm=-1

        for x in indi_2mm_con:
            if seq.find(x)!=-1:
                v_confirm=1
                break

        if v_confirm==-1: continue

        for n, x in enumerate(indi_l):
            pos=seq.find(x)
            if pos!=-1:
                amend=1-amend_dic[int(n/indil_len)]+sliced_len
                subseq=seq[pos+amend:pos+amend+cut_len].replace('\n','')
                if len(subseq)!=cut_len: break
                try: dic[subseq]+=1
                except: dic[subseq]=1
                break

        if i*100/flen > pc:
            progress_bar(pc+0.1, 100, 'collecting ...'+str(round(time.time()-t,2)), '', 'data collected')
            pc+=0.1

    sort=sorted(dic.items(), key=operator.itemgetter(1), reverse=True)

    return dic, sort

def DB_modeling():

    fdb=open('./DB/exon_intron_v2.txt').readlines()
    exin_db={}
    for i in range(1, 23):
        exin_db[str(i)]={}
    exin_db['X']={}
    exin_db['Y']={}
    lencnt=0
    for i in fdb[1:]:
        isp=i.replace('\n','').split()
        chro=isp[0]
        st=int(isp[7])
        ed=int(isp[8])
        a=st
        i_info=isp[0]+'\t'+isp[2]+'\t'+isp[3]+'\t'+isp[4]+'\t'+isp[5]
        if str(st)[:-3]==str(ed)[:-3]:
            M=int(a/1000000)
            K=int(a/1000)
            if K>1000: K=int(str(K)[-3:])
            Bst=('00'+str(st))[-3:]
            Bed=('00'+str(ed))[-3:]
            if M not in exin_db[chro].keys():
                exin_db[chro][M]={}
            if K not in exin_db[chro][M].keys():
                exin_db[chro][M][K]={}
            exin_db[chro][M][K][Bst+';~;'+Bed+';'+str(lencnt)]=i_info
            lencnt+=1
        else:
            while a<ed+1000:
                M=int(a/1000000)
                K=int(a/1000)
                if K>1000: K=int(str(K)[-3:])
                B=('00'+str(a))[-3:]
                if M not in exin_db[chro].keys():
                    exin_db[chro][M]={}
                if K not in exin_db[chro][M].keys():
                    exin_db[chro][M][K]={}
                if a==st:
                    exin_db[chro][M][K][B+';+;;'+str(lencnt)]=i_info
                elif a>=ed:
                    M=int(ed/1000000)
                    K=int(ed/1000)
                    if K>1000: K=int(str(K)[-3:])
                    B=('00'+str(ed))[-3:]
                    exin_db[chro][M][K][B+';-;;'+str(lencnt)]=i_info
                else:
                    exin_db[chro][M][K][B+';-;;'+str(lencnt)]=i_info
                a+=1000
                lencnt+=1
    return exin_db

def read_fastmap_results(q):

    fastmap_list=str(q[0]).split('//')

    fastmap_dic={}
    chromomap_dic={}
    id_chromomap=0

    for info in fastmap_list[:-1]:
        info=info.split('\\n')
        seq=info[0].split('\\t')[1]
        if len(info)==2:
            fastmap_dic[seq]='Can not find position'
            continue
        sub_info_dic={}
        for sub_info in info[1:-1]:
            sub_info_sp=sub_info.split('\\t')
            if sub_info_sp[1]=='0':
                sub_info_dic[int(sub_info_sp[2])-int(sub_info_sp[1])]=sub_info
        if len(sub_info_dic)==0:
            fastmap_dic[seq]='Can not find position'
            continue
        info=sorted(sub_info_dic.items(), key=operator.itemgetter(0), reverse=True)[0][1]
        cut_len=int(info.split('\\t')[2])
        if info.split('\\t')[-1]=='*':
            fastmap_dic[seq]='fastmap find '+info.split('\\t')[-2]+' sites;'+seq[:cut_len]
        else:
            fastmap_dic[seq]=[]
            chro_c=-1
            if len(info.split('\\t')[4:])==1:
                i=info.split('\\t')[4].split(':')
                chromomap_dic[str(id_chromomap)]=[i[0], str(abs(int(i[1]))), str(abs(int(i[1]))+int(info.split('\\t')[2]))]
                id_chromomap+=1
                chro_c=1
            for i in info.split('\\t')[4:]:
                i=i.split(':')
                fastmap_dic[seq].append([i[0], i[1][:1], i[1][1:], seq[:cut_len], chro_c])
    return fastmap_dic

def analyze_fastmap_results(filename, sorted_tup, fastmap_dic, exin_db, cut_len):
    
    length_sorted_list=[]
    for seq, info in fastmap_dic.items():
        yl=fastmap_dic[seq]
        print(yl)
        if type(yl)!=str and len(yl[0][3])==cut_len:
            length_sorted_list.append(yl[0][3])
    print(len(length_sorted_list))
    consolidate_dic={}
    for x, cnt in sorted_tup:
        yl=fastmap_dic[x]
        if type(yl)!=str and len(yl[0][3])!=cut_len:
            c=0
            seq=yl[0][3]
            for i in length_sorted_list:
                if i.find(seq)!=-1:
                    print(seq)
                    print(yl)
                    try: consolidate_dic[i]+=cnt
                    except: consolidate_dic[i]=cnt
                    c=1
                    continue
            if c==0:
                try: consolidate_dic[x]=cnt
                except: consolidate_dic[x]=cnt
        else:
            try: consolidate_dic[x]+=cnt
            except: consolidate_dic[x]=cnt
    
    chromomap_dic={}
    chromomap_id=0
    fw=open(filename,'w')
    fwc=open(filename[:filename.find('.')]+'_chromomap.txt','w')
    id_c=0
    fw.write('ID\tsequence\tchromosome\tstart\tread\tgene\tclass\texon ID\n')
    for x, cnt in sorted(consolidate_dic.items(), key=operator.itemgetter(1), reverse=True):
        id_c+=1
        x=x.replace('\n','')
        yl=fastmap_dic[x]
        if type(yl)==str:
            if yl.find(';')!=-1:
                yl=yl.split(';')
                fw.write(str(id_c)+'\t'+yl[1]+'\t\t'+str(cnt)+'\t'+yl[0]+'\n')
            else:
                fw.write(str(id_c)+'\t'+x+'\t\t'+str(cnt)+'\t'+yl+'\n')
            continue
        for y in yl:
            chromomap_dic[str(chromomap_id)]=[y[0], int(y[2]), int(y[2])+len(y[3])]
            chromomap_id+=1
            fw.write(str(id_c)+'\t'+y[3]+'\t'+y[0]+'\t'+y[2]+'\t'+str(cnt))
            if y[4]==1: fwc.write(str(id_c)+'\t'+y[3]+'\t'+y[0]+'\t'+y[2]+'\t'+str(cnt))
            chro=y[0]
            pos=int(y[2])
            B=int(('00'+str(pos))[-3:])
            M=int(pos/1000000)
            K=int(('00'+str(int(pos/1000)))[-3:])
            try:
                db=exin_db[chro][M][K]
            except:
                fw.write('\tNothing in database\n')
                if y[4]==1: fwc.write('\tNothing in database\n')
                continue

            vv=0
            sub_db={}

            for i in db.keys():
                isp=i.split(';')
                v=0
                if isp[1]=='=': v=1
                elif int(isp[0])<=B and isp[1]=='+': v=1
                elif int(isp[0])>=B and isp[1]=='-': v=1
                elif isp[1]=='~':
                    if int(isp[0])<=B<=int(isp[3]): v=1
                if v==1:
                    vv=1
                    dbl=db[i].split('\t')
                    try:
                        ssub_db=sub_db[dbl[1]]
                        if dbl[2] not in ssub_db['chtype'].split('&'): ssubdb['chtype']+='&'+dbl[2]
                        if ssub_db['exin']=='intron' and dbl[3]=='exon': ssub_db['exin']='exon'
                    except:
                        sub_db[dbl[1]]={'chtype':dbl[2], 'exin':dbl[3]}
        
            for gene, gene_list in sub_db.items():
                fw.write('\t'+gene+'\t'+gene_list['chtype']+'\t'+gene_list['exin'])
                if y[4]==1: fwc.write('\t'+gene+'\t'+gene_list['chtype']+'\t'+gene_list['exin'])
            if vv==0:
                fw.write('\tNothing in database\n')
                if y[4]==1: fwc.write('\tNothing in database\n')
            else: 
                fw.write('\n')
                if y[4]==1: fwc.write('\n')
    fw.close()
    fwc.close()
    
    return chromomap_dic

def make_chromomap(filename, chromomap_dic):
    
    fw=open(filename[:filename.find('.')]+'_chromomap.html','w')
    s='<head>\n\t<script src="https://unpkg.com/ideogram@1.4.1/dist/js/ideogram.min.js"></script>\n</head>'
    s+="<body>\n<script>\n\tvar ideogram = new Ideogram({ organism: 'human', chrMargin: 5, annotationsLayout: 'histogram', annotations: ["
    for x, y in chromomap_dic.items():
        s+="{ name: '%s', chr: '%s', start: %s, stop: %s }," %(x, y[0], y[1], y[2])
    s=s[:-1]+"] });\n</script>\n</body>"
    fw.write(s)
    fw.close()
