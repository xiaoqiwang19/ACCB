#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer
import urllib2

def convert(Strand,seq):
    strand=Strand
    if strand=="-":
        complement = {'A':'T','G':'C','C':'G','T':'A'}
        revSeqList = list(reversed(seq.upper()))
        revComSeqList = [complement[k] for k in revSeqList]     
        seq_c = ''.join(revComSeqList)   
        return seq_c
    else:
       seq_c =seq.upper() 
       return seq_c

def dup_correct(Strand,chr_num,start_pos,end_pos):
    #this query tool refer instructions http://www.oebiotech.com/Article/rhgjqswzhd.html
    url="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+chr_num+":"+start_pos+","+end_pos
    html = urllib2.urlopen(url).read()
    pat=re.compile(r'>\n([atgc]*)\n<')
    dup_seq=re.findall(pat,html)[0]
    #strand=Strand
    #dup_seq_c=convert(strand,dup_seq)    
    dup_seq_c=dup_seq.upper()
    start_pos_c=int(start_pos)-1
    end_pos_c=int(start_pos)
    return dup_seq_c,start_pos_c,end_pos_c
def del_correct(Strand,chr_num,start_pos,end_pos):
    url="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+chr_num+":"+start_pos+","+end_pos
    html = urllib2.urlopen(url).read()
    pat=re.compile(r'>\n([atgc]*)\n<')
    del_seq=re.findall(pat, html)[0]
    #strand=Strand
    #del_seq_c=convert(strand,del_seq)
    del_seq_c=del_seq.upper()
    return del_seq_c

def main(argv):
    if len(sys.argv) !=4:
        print "please check the usage:\
               python script.py  transcript_list result out"
        sys.exit()
    trans = argv[1]
    anno = argv[2]
    output = argv[3]
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    am = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37', alt_aln_method='splign',replace_reference=True)
    hn = hgvs.normalizer.Normalizer(hdp)
    my_d = {}
    my_strand={}
    with open(trans,'r') as NM :
        for i,hang in enumerate(NM,1):
            hang = hang.rstrip()
            t_version = hang.split('\t')
            my_d[t_version[0]] = t_version[1]
            my_strand[t_version[0]] = t_version[2]

    with open(anno,'r') as IN,\
        open(output,'w') as OUT :
        for n,line in enumerate(IN,1):
            line = line.rstrip()
            if line.startswith("Chr"):
                OUT.write("Chr_c\tStart_c\tEnd_c\tStrand\tAAChange.refGene\tRef_c\tAlt_c\t"+str(line) +"\n")
            else:
                content = line.split('\t')
                if content[9]!=".":
                    all_c=content[9].split(",") 
                    AAChange=[]
                    #ref=[]
                    for each_c in all_c:
                        transcrip=each_c.split(":")[1]
                        c_change=each_c.split(":")[3]
                        chr_num=content[0]
                        if transcrip in my_d :
                            Strand=my_strand[transcrip]
                            transcrip = my_d[transcrip]
                            #print transcrip
                            try:
                                c1 = hp.parse_hgvs_variant(transcrip + ':' + c_change) 
                                c1n = hn.normalize(c1)
                                g = am.c_to_g(c1n)
                                p = am.c_to_p(c1n)
                                pep = re.sub('[()]','',str(p))
                                transcrip,c_change = str(c1n).split(':')
                                np,p_change = str(pep).split(':')
                                aaChange=each_c.split(":")[0]+":"+transcrip+":"+each_c.split(":")[2]+":"+c_change+":"+p_change #updata the normalized c_change and p_change
                                AAChange.append(aaChange)                    #updata one AAChange
                                start_pos=str(g.posedit.pos.start)
                                end_pos=str(g.posedit.pos.end)               #correct for dup,end_pos need to be writen by insert format
                                pat=re.compile(r'c.\d*([ATGC]>[ATGC])')
                                ref=re.findall(pat,str(c1n))
                                if ref:  #for SNV variants          #olny treat the last c_change,because 
                                    Ref_c=convert(Strand,ref[0].split(">")[0])
                                    Alt_c=convert(Strand,ref[0].split(">")[1])
                                else:    #respectively correct dup,del,ins,etc.
                                    if "dup" in str(c1n):
                                        Ref_c="."
                                        Alt_c,start_pos,end_pos=dup_correct(Strand,chr_num,start_pos,end_pos) 
                                    elif "del" in str(c1n) and "ins" not in str(c1n):
                                        Alt_c="."
                                        Ref_c=del_correct(Strand,chr_num,start_pos,end_pos)
                                    elif "ins" in str(c1n) and "del" not in str(c1n): 
                                        Ref_c="."
                                        Alt_c=convert(Strand,str(c1n).split("ins")[1])
                                    elif "ins" in str(c1n) and "del" in str(c1n):
                                        Alt_c=convert(Strand,str(c1n).split("ins")[1])
                                        Ref_c=del_correct(Strand,chr_num,start_pos,end_pos)
                                    else:
                                        print transcrip,c_change,"this format is not in above situation,please check"   
                            except Exception,e:                                    
                                aaChange="NA"
                                log_file=open("log","a")
                                log_file.write(repr(e.message)+"\t"+transcrip+"\n")
                                c_change="NA"
                                p_change="NA"
                                start_pos="NA"
                                end_pos="NA"
                                Ref_c="NA"
                                Alt_c="NA"
                                Ref_c="NA"
                                Alt_c="NA"
                    old='\t'.join(content)
                    AAChange_c=",".join(AAChange)
                    OUT.write(content[0]+"\t"+str(start_pos)+"\t"+str(end_pos)+"\t"+Strand+"\t"+AAChange_c+"\t"+Ref_c+"\t"+Alt_c+"\t"+old+"\n")
                else:
                    old='\t'.join(content)
                    OUT.write(content[0]+"\t"+content[1]+"\t"+content[2]+"\t"+Strand+"\t"+content[9]+"\t"+content[3]+"\t"+content[4]+"\t"+old+"\n")

#    OUT.close()
#    log_file.close()                            
if __name__ == "__main__":
     main(sys.argv)
