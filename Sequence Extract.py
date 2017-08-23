def rc(seq):
    '''
This function recieve DNA sequence and return sequence reverse complementary 
    '''
    f=[]
    for base in seq:
        if base=="A":
            f.append('T')
        elif base=="T":
            f.append('A')
        elif base=="C":
            f.append('G')
        elif base=="G":
            f.append('C')
    f.reverse()
    r=''.join(f)
    return(r)

def gene_extract(fasta_file,Genbank_file):
        '''
This function determine the start and end of each gene according to GenBank annotation and extract each gene in an indvedual fasta file
        '''
        t=open(fasta_file,'r')
        r=open(Genbank_file,'r')
        z=list(r)
        k=list(t)
        del k[0]
        e=[]
        for base in k:
                mm=base.replace('\n','')
                e.append(mm)

        j=''.join(e)
        c=0
        x=[]
        while c<=len(z)-1:
                if 'gene   ' in z[c]:
                        rr=z[c]+z[c+1]
                        x.append(rr)
                        c=c+1
                else:
                        c=c+1
        v=[]
        b=[]
        for base in x:
                rr=base.split('               ')
                b.append(rr[0])
                v.append(rr[1])
        m=[]
        for base in v:
            bhm=base.replace('"','')
            bh=bhm.replace('\n','')
            m.append(bh)

        d=[]
        for base in m:
            bhm=base.replace('/gene=','')
            bh=bhm.replace(' ','')
            bhm=bh.replace('/','')
            d.append(bhm)

            
        i=[]
        for base in b:
            bhm=base.replace('complement','')
            bh=bhm.replace('\n','')
            bhm=bh.replace(' ','')
            bh=bhm.replace('gene','')
            bhm=bh.replace('join','')
            bh=bhm.replace('(','')
            bhm=bh.replace(')','')
            i.append(bhm)

        aw,u,g,ch=[],[],[],0
        for base in i:
                if ','in base:
                        ii=base.split(',')
                        u.append(ii[0])
                        u.append(ii[1])
                        g.append(d[ch])
                        g.append(d[ch]+'(1)')
                        ch=ch+1
                else:
                        u.append(base)
                        g.append(d[ch])
                        ch=ch+1

        vb,cc,hh,zz,hn=[],[],[],[],1
        for base in u:
                bn=base.split('..')
                bb=int(bn[0])-1
                qq=int(bn[1])-1
                xc=(qq-bb)+1
                hh.append(xc)
                while bb<=qq:
                        vb.append(j[bb])
                        bb=bb+1
                        hn=hn+1
                        if hn==78:
                                vb.append('\n ')
                                hn=1
                        if bb==qq:
                                sc=''.join(vb)
                                we=len(sc)
                                zz.append(we)
                                vb=[]
                                cc.append(sc)
                                hn=1
        uu=0
        while uu<=len(cc)-1:
                op=open(g[uu]+'.fas','w')
                op.writelines('>'+str(hh[uu])+g[uu]+'  macha\n')
                op.writelines(cc[uu])
                uu=uu+1
        op.close()
        return "Number of return Genes= "+ str(uu-1)

def intron_extract(fasta_file,Genbank_file):
        '''
This function determine the start and end of each gene according to GenBank annotation and extract the whole entergenic regions in a single fasta file
        '''
        t=open(fasta_file,'r')
        r=open(Genbank_file,'r')
        z=list(r)
        k=list(t)
        del k[0]
        e=[]
        for base in k:
                mm=base.replace('\n','')
                e.append(mm)

        j=''.join(e)
        c=0
        x=[]
        while c<=len(z)-1:
                if 'gene   ' in z[c]:
                        rr=z[c]+z[c+1]
                        x.append(rr)
                        c=c+1
                else:
                        c=c+1
        v=[]
        b=[]
        for base in x:
                rr=base.split('               ')
                b.append(rr[0])
                v.append(rr[1])
        m=[]
        for base in v:
            bhm=base.replace('"','')
            bh=bhm.replace('\n','')
            m.append(bh)

        d=[]
        for base in m:
            bhm=base.replace('/gene=','')
            bh=bhm.replace(' ','')
            bhm=bh.replace('/','')
            d.append(bhm)

            
        i=[]
        for base in b:
            bhm=base.replace('complement','')
            bh=bhm.replace('\n','')
            bhm=bh.replace(' ','')
            bh=bhm.replace('gene','')
            bhm=bh.replace('join','')
            bh=bhm.replace('(','')
            bhm=bh.replace(')','')
            i.append(bhm)

        aw,u,g,ch=[],[],[],0
        for base in i:
                if ','in base:
                        ii=base.split(',')
                        u.append(ii[0])
                        u.append(ii[1])
                        g.append(d[ch])
                        g.append(d[ch]+'(1)')
                        ch=ch+1
                else:
                        u.append(base)
                        g.append(d[ch])
                        ch=ch+1

        vb,cc,hh,zz,hn=[],[],[],[],1
        for base in u:
                bn=base.split('..')
                bb=int(bn[0])-1
                qq=int(bn[1])-1
                xc=(qq-bb)+1
                hh.append(xc)
                while bb<=qq:
                        vb.append(j[bb])
                        bb=bb+1
                        hn=hn+1
                        if hn==78:
                                vb.append('\n ')
                                hn=1
                        if bb==qq:
                                sc=''.join(vb)
                                we=len(sc)
                                zz.append(we)
                                vb=[]
                                cc.append(sc)
                                hn=1

        qqq=[]
        for base in cc:
                mer=base.replace('\n','')
                sem=mer.replace(' ','')
                qqq.append(sem)
        aa=1
        aaa=[]
        ddd=[]
        for base in qqq:
                while aa<=len(base)-2:
                        aaa.append(base[aa])
                        aa=aa+1
                cccc=''.join(aaa)
                aa=1
                aaa=[]
                ddd.append(cccc)
        rrr=[j]
        for base in ddd:
                if base in j:
                        vvv=rrr[0].replace(base,'NNNNN')
                        rrr.append(vvv)
                        del rrr[0]
                else:
                        print('not found')
        doc=rrr[0]
        op=open('intron.fasta','w')
        op.write(doc)
        op.close()
        return op


def primer_pick_up(Forward_primer,Reverse_primer, numer_of_fasta_files):
    '''
This function return the sequence between forward and reverse primers
    '''
    a=Forward_primer
    b=Reverse_primer
    c=rc(Reverse_primer)
    d=rc(Forward_primer)
    pm=1
    while pm<=numer_of_fasta_files: 
        op=open(str(pm)+'.fas','r')
        z=list(op)
        xx=z[0].replace('\n','')
        doc=z[0]
        del z[0]
        r=''.join(z)
        x=r.replace('\n','')
        e=x.replace(a,'-'+a)
        f=e.replace(b,b+'-')
        g=f.replace(c,'-'+c)
        h=g.replace(d,d+'-')
        i=h.split('-')
        j=[]
        for base in i:
            if (c and d) in base:
                    j.append(base)
            elif (a and b) in base:
                    j.append(base)
        try:
            yy=str(pm)
            op=open('Gene'+yy+'.fas','w')
            op.writelines(yy+doc+'\n')
            op.writelines(j[0])
            op.close()
        pm=pm+1
    return 
    

