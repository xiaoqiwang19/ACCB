#!/usr/bin/python
# -*- coding: utf-8 -*-
#modified the error about "u'No transcript definition for transcript" 
#compared with annovar_test0408.py:add a parameter :chr_list of the genome you used  
#usage: python script.py transcrip_version_strand.list list_chr input_annotion_file output
#edited on 0411,2018 
import sys
import re
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer
import urllib2
import os
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
am = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37', alt_aln_method='splign',replace_reference=True)
hn = hgvs.normalizer.Normalizer(hdp)

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
    dup_seq_c=dup_seq.upper()
    start_pos_c=int(start_pos)-1
    end_pos_c=int(start_pos)
    return dup_seq_c,start_pos_c,end_pos_c
def del_correct(Strand,chr_num,start_pos,end_pos):
    url="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+chr_num+":"+start_pos+","+end_pos
    html = urllib2.urlopen(url).read()
    pat=re.compile(r'>\n([atgc]*)\n<')
    del_seq=re.findall(pat, html)[0]
    del_seq_c=del_seq.upper()
    return del_seq_c
def transcript_correct(chr_id,transcript,start,c_change):
    #hp = hgvs.parser.Parser()
    #hdp = hgvs.dataproviders.uta.connect()
    #am = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37', alt_aln_method='splign',replace_reference=True)
    if "dup" in c_change: #I didn't consider the shorthand notation of dup,if there be, add len(aa)==0
        pat=re.compile(r'c.*(dup[ATCG]*)')
        aa=re.findall(pat,c_change)[0]
        if len(aa)==1:
            hgvs_g_creat=chr_id+":g."+str(start)+aa
            var_g = hp.parse_hgvs_variant(hgvs_g_creat)
            transcripts_list = am.relevant_transcripts(var_g)
        elif len(aa)>1:
            end=int(start)+len(aa)-1
            hgvs_g_creat=chr_id+":g."+str(start)+"_"+str(end)+aa
            var_g = hp.parse_hgvs_variant(hgvs_g_creat)
            transcripts_list = am.relevant_transcripts(var_g)
    if "del" in c_change and "ins" not in c_change: #I didn't consider the shorthand notation of del,if there be, add len(aa)==0
        pat=re.compile(r'c.*(del[ATCG]*)')
        aa=re.findall(pat,c_change)[0]
        end=int(start)+len(aa)
        hgvs_g_creat=chr_id+":g."+str(start)+"_"+str(end)+aa
        var_g = hp.parse_hgvs_variant(hgvs_g_creat)
        transcripts_list = am.relevant_transcripts(var_g)             
    elif "ins" in c_change and "del" not in c_change:
        pat=re.compile(r'c.*(ins[ATCG]*)')
        aa=re.findall(pat,c_change)[0]
        end=int(start)+1
        hgvs_g_creat=chr_id+":g."+str(start)+"_"+str(end)+aa
        var_g = hp.parse_hgvs_variant(hgvs_g_creat)
        transcripts_list = am.relevant_transcripts(var_g)
    #there is no examples for coexisting of "ins" and "del"
    elif ">" in c_change:
        pat=re.compile(r'c.\d.*\d([ATGC]{1,}>[ATGC]{1,})')
        aa=re.findall(pat,c_change)[0]
        if len(aa.split(">")[0])==1:
            hgvs_g_creat=chr_id+":g."+str(start)+aa
            var_g = hp.parse_hgvs_variant(hgvs_g_creat)
            transcripts_list = am.relevant_transcripts(var_g)
        else:
            print c_change," format is wrong,please check " # It's usually signal locus mutation,start equal end position
            transcripts_list=[]
    else:
        transcripts_list=[]
    if transcripts_list:
        transcript_orign=transcript.split(".")[0]
        for i in transcripts_list:
            if transcript_orign in i:
                transcript=i
                return transcript
                break
        return transcript
    else:
        return transcript
def main(argv):
    if len(sys.argv) !=5:
        print "please check the usage:\
               python script.py  transcrip_version_strand.list list_chr input_annotion_file output\n \
               parameter1:script.py\nparameter2:transcrip_version_strand.list\nparameter3:list_chr\nparameter4:input_annotion_file\nparameter5:output_file"
        sys.exit()
    log_file=open("log","w")
    trans = argv[1]
    chr_list=argv[2]
    anno = argv[3]
    output = argv[4]
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    am = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37', alt_aln_method='splign',replace_reference=True)
    hn = hgvs.normalizer.Normalizer(hdp)
    my_d = {}
    my_strand={}
    my_chr={}
    with open(trans,'r') as NM ,\
         open(chr_list,'r') as CHR:
        for i,hang in enumerate(NM,1):
            hang = hang.rstrip()
            t_version = hang.split('\t')
            my_d[t_version[0]] = t_version[1]
            my_strand[t_version[0]] = t_version[2]
        for eachline in CHR:
            each=eachline.strip()
            my_chr[each.split("\t")[0]]=each.split("\t")[1]
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
                    for each_c in all_c:
                        transcrip_1=each_c.split(":")[1]
                        c_change=each_c.split(":")[3]
                        chr_num=content[0]
                        if transcrip_1 in my_d :
                            Strand=my_strand[transcrip_1]
                            transcript_corrected_before = my_d[transcrip_1]
                            start=content[1]
                            transcript=transcript_corrected_before
                            chr_id=my_chr[chr_num]
                            transcrip=transcript_correct(chr_id,transcript,start,c_change)
                            #print transcrip,"sign"
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
                                log_file.write(repr(e.message)+"\t"+transcrip+"\n")
                                aaChange="NA"
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
                else:  # It mainly deals with intron variation.
                    if content[3]=="-":
                        end_pos=int(content[1])+len(content[4])
                    else:
                        end_pos=content[2] 
                    old='\t'.join(content)
                    OUT.write(content[0]+"\t"+content[1]+"\t"+str(end_pos)+"\t"+Strand+"\t"+content[9]+"\t"+content[3]+"\t"+content[4]+"\t"+old+"\n")

#    OUT.close()
#    log_file.close()                            
if __name__ == "__main__":
     main(sys.argv)
