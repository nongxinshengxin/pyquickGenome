#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq 
import argparse
import re

class QuickGenome():
    def __init__(self,fasta,output):
        self.fasta=fasta
        self.output=output
    

    def readfasta(self):
        sequences=list(SeqIO.parse(self.fasta,"fasta"))
        sequences.sort(key=lambda x:len(x),reverse=True)
        for seq in sequences:
            yield seq

    def genomeBasicInfo(self):
        genomeSize=0
        time=0
        seqG=0
        seqC=0
        n50size=0
        n90size=0
        time2=0
        time3=0
        with open(f'{self.output}/genomeBasicInfo.txt',"w") as w:
            for seq in self.readfasta():
                genomeSize+=len(seq)
                
                if time==0:
                    largestChr=len(seq)
                seqG+=seq.count("G")
                seqC+=seq.count("C")
                time+=1
            gc_content=((seqG+seqC)/genomeSize)
            for seq2 in self.readfasta():
                n50size+=len(seq2)
                time2+=1
                if n50size/genomeSize>=0.5:
                    n50=len(seq2)
                    l50=time2
                    break
            for seq3 in self.readfasta():
                n90size+=len(seq3)
                time3+=1                
                if n90size/genomeSize>=0.9:
                    n90=len(seq3)
                    l90=time3
                    break
            w.write(f"The size of genome is {genomeSize} \n")
            w.write(f"The number of contig/scaffold/chromosome is {time}\n")
            w.writelines(f"The size of largest contig/scaffold/chromosome is {largestChr}\n")
            w.writelines(f"GC content is {gc_content}\n")
            w.writelines(f"N50 is {n50}\n")
            w.writelines(f"L50 is {l50}\n")
            w.writelines(f"N90 is {n90}\n")
            w.writelines(f"L90 is {l90}\n")
        #print(genomeSize,largestChr,time,gc_content,n50)


    def findGap(self):
        with open(f'{self.output}/genomeGap.txt',"w") as w:
            for seq in self.readfasta():
                gap=re.finditer('N+',str(seq.seq),re.I)
                for pos in gap:
                    w.write(seq.id+"\t"+str(pos.span()[0])+"\t"+str(pos.span()[1])+"\n")


    def binGC(self,bin):
        with open(f'{self.output}/binGC.txt',"w") as w:

        
            for seq in self.readfasta():
                start=0
                end=start+bin
                while start<len(seq)-1:
                        
                    if end<=len(seq):
                        end=end
                    else:
                        end=len(seq)
                    perseq=seq[start:end]
                    gc_content=(perseq.count("G")+perseq.count("C"))/len(perseq)
                    w.write(seq.id+"\t"+str(start)+"\t"+str(end)+"\t"+str(gc_content)+"\n")
                    start+=bin
                    end=start+bin
                #w.write(seq.id+"\t"+str(start)+"\t"+str(end)+"\t"+str(gc_content)+"\n")





    def findTelo(self,telo,rangeArea,minRepeatNum):
        with open(f'{self.output}/telo.txt',"w") as w:
            telo=telo*minRepeatNum
            telo_re=Seq(telo).reverse_complement()*minRepeatNum
            for seq in self.readfasta():
                geneid=seq.id
                seq=seq.seq
                leftseq=seq[0:rangeArea]
                rightseq=seq[len(seq)-rangeArea:len(seq)]
                rightseq=rightseq.reverse_complement()
                if leftseq.find(telo)!=-1:
                    leftpos=leftseq.find(telo)+1
                elif leftseq.find(telo)==-1 and leftseq.find(telo_re)!=-1:
                    leftpos=leftseq.find(telo_re)+1
                else:
                    leftpos="NA"
                if rightseq.find(telo)!=-1:
                    rightpos=len(seq)-rightseq.find(telo)
                elif rightseq.find(telo_re)!=-1:
                    rightpos=len(seq)-rightseq.find(telo_re)
                elif rightseq.find(telo)==-1 and rightseq.find(telo_re)==-1:
                    rightpos="NA"
                w.write(geneid+"\t"+str(len(seq))+"\t"+str(leftpos)+"\t"+str(rightpos)+"\n")






def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('fastaFile',type=str,help='输入fasta文件路径')
    parser.add_argument('output',type=str,help='结果输出文件存放路径')
    parser.add_argument('-m','--module',type=str,default="basicInfo",help='选择模式')
    parser.add_argument('--bin',type=int,default=10000,help='滑窗长度')
    parser.add_argument('--telo',type=str,default="CCCTAA",help='端粒重复序列')
    parser.add_argument('--rangeArea',type=int,default=10000,help='染色体两端范围')
    parser.add_argument('--minRepeatNum',type=int,default=2,help='端粒重复最小次数')
    args=parser.parse_args()


    qgenome=QuickGenome(args.fastaFile,args.output)
    if args.module=="basicInfo":
        qgenome.genomeBasicInfo()
    if args.module=="binGC":
        qgenome.binGC(args.bin)
    if args.module=="findGap":
        qgenome.findGap()
    if args.module=="findTelo":
        qgenome.findTelo(args.telo,args.rangeArea,args.minRepeatNum)


if __name__=="__main__":
    main()
