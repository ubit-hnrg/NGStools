
# coding: utf-8

# In[40]:


import argparse
import vcf
import os
import glob
import time


# In[20]:


## parse args
description =['reduce a bam file acording to vcf file']
if True:
    parser = argparse.ArgumentParser(description= description[0])
    parser.add_argument('-b','--bamfile',required=True,help='input bam file')
    parser.add_argument('-v','--vcffile', help='input vcf file to be usde for filtering coordinates', required=True)
    parser.add_argument('-o','--outpath',required=True)
    parser.add_argument('-p','--prefix',dest = 'prefix',default = '',required = False)
    parser.add_argument('-d','--delta',required = False, default = 500)

    args = parser.parse_args()
    bamfile = args.bamfile
    vcffile = args.vcffile
    outpath = args.outpath
    prefix = args.prefix
    delta = args.delta
else:
    bamfile = './high-homozygosis.vcf'
    vcfile = '1A-mindp20'
    outpath = './'  
    prefix = 'shortBam-CEDIE1A'
    delta = 500    

#bamfile = '/home/ariel/Projects/Gutierrez/CEDIE/201803/bamfiles/chr2-190-192k.bam'
#vcfile = '/home/ariel/Projects/Gutierrez/CEDIE/201803/annot/CEDIE1A_minDP20_final_annot.vcf'
#vcfile = '/home/ariel/test.vcf'




# In[27]:


# create output dir

if not os.path.isdir(outpath):
    os.system('mkdir %s'%outpath)

t1 = time.time()
k = 0 
for rec in vcf.Reader(open(vcffile)):
    ch,pos =  rec.CHROM, rec.POS
    #if ch!=2:
    #    continue
    start=pos - delta
    end = pos + delta
    tempfile = '%s%s_%s.temp'%(outpath,prefix,k)
    runsamtools = 'samtools view %s -b %s:%s-%s > %s'%(bamfile,ch,start,end,tempfile)
    #print runsamtools
    os.system(runsamtools)
    k = k+1



# In[39]:


pattern = outpath+prefix+'*.temp'    
tempfiles = glob.glob(pattern)
tojoin = ' '.join(tempfiles)

mergebamfile = outpath+prefix+'_merged.bam'
runmerge = 'samtools merge %s %s'%(mergebamfile,tojoin)
print runmerge
os.system(runmerge)
os.system('rm %s'%tojoin)

deltatime= time.time()-t1
print '%s variats proceced'%k
print 'ranged extracted: +/- %s bases around each variant'%delta
print 'elapsed time : %s'%deltatime


# In[43]:


time.clock()


# In[48]:


time.time()

