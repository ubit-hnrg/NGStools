
# coding: utf-8

# In[3]:


import pandas as pd
import re
import subprocess
import os
import argparse


# In[4]:


description = ['This script convert vcf files into tsv by using vcf2tsv script in vcflib \n Vcf2tsv fails to conver flags features that dont have explicetly equal the symbol \n This script first add a number 1 to these flags .Also a report in tsv format is \n incuded for describing all features in INFO field']

if (True):
    parser = argparse.ArgumentParser(description= description[0])
    parser.add_argument('-i','--inputfile',required=True,help='file.vcf as input')
    parser.add_argument('-o','--outprefix', help='prefit to be used in output and feature description tsv files', required=True)
    parser.add_argument('-x','--xls', dest='xls', action='store_true')
    parser.set_defaults(xls=False)
    
    #parser.add_argument('--no-feature', dest='feature', action='store_false')

    #parser.add_argument('-x','--xls',required=False,,default = False, help='Boolean if you want to save outputs in xls format too. ssconvert of gnumeric package is required. Default False')
    
                                 
    args = parser.parse_args()
    ifile = args.inputfile
    prefix = args.outprefix
    xls=args.xls
    
else:
    ifile = './high-homozygosis.vcf'
    prefix = '1A-mindp20'                                 
    xls = False              
        
outfile = prefix +'.tsv'
infofile = 'INFO'+prefix+'.tsv' 


# In[5]:


print 'PARSING HEADER'
nl = 0
with open(ifile) as f:
    Ids =[];Number=[];Type=[];Description=[]
    for line in f:
        nl = nl+1    
        if (nl> 1000):
            break
        if re.match(r'(##INFO)', line):
            field = line.split(',')
            fields = field[0:3]+[''.join(field[3:]).strip('\n>')]
            if(len(fields)!=4):
                print 'warning!'

            Ids.append(fields[0].split('=')[-1])
            Number.append(fields[1].split('=')[-1])
            Type.append(fields[2].split('=')[-1])
            Description.append(fields[3].split('=')[-1])

# you may also want to remove whitespace characters like `\n` at the end of each line
#content = [x.strip() for x in content] 


# In[ ]:


info = pd.DataFrame({'Id':Ids,'Number':Number,'Type':Type,'Description':Description})
info = info[['Id','Number','Type','Description']]
info.to_csv(infofile,index = False,sep ='\t')
print 'info fields parsed'
print 'there are %s INFO features'%info.shape[0]



# In[97]:


print 'backup of input vcf'

flagedfile = 'flaged_%s.vcf'%prefix

#make a sequtity copy
os.system('cp %s %s'%(ifile,flagedfile))


#define flag features (boolean columns in the final tsv file)
flags = info[info.Type=='Flag'].Id
print 'fixing vcf file (FLAG fields) before to change format'

for fl in flags:
    ##reemplazo flags dentro de INFO    

    os.system("sed -i -e 's/;%s;/;%s=1;/g' %s"%(fl,fl,flagedfile)) 
    #reemplazo flags que son el último campo de INFO
    os.system("sed -i -e 's/;%s\t/;%s=1\t/g' %s"%(fl,fl,flagedfile)) 
   

print 'excecuting vcf2tsv - vcflib'
#ejecuto vcf2tsv de vcflib
os.system('vcf2tsv %s > %s'%(flagedfile,outfile))
    


# In[ ]:


# CONVERT OUTPUTS TO xls in two sheets


# In[ ]:


print 'Saving into xls file'
if xls: 
    os.system('ssconvert --merge-to %s.xls %s %s'%(prefix,outfile,infofile))


# In[ ]:


#"vcffilter -f 'SNPEFF_IMPACT = HIGH & HOM' CEDIE1A_minDP20_final_annot.vcf > high-homozygosis.vcf "
#bcftools query --print-header -f %s /home/ariel/Projects/Gutierrez/CEDIE/201803/annot/high-homozygosis.vcf > all_info.table"%tabIds

#subprocess.call(['sed','-i','-e','%s'%pattern,'%s'%flagedfile]) # no me esta funcando este por ahora
    
#tabIds="'[\\t%"+'\\t%'.join(flags)+"]\\n'"
#print tabIds

