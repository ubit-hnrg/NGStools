samples_file = 'samples.list'
ubams_file = 'ubams.list'

with open(samples_file) as sf:
    samples = sf.readlines()
samples = [x.strip() for x in samples] 

with open(ubams_file) as ubf:
    ubams = ubf.readlines()
ubams = [x.strip() for x in ubams] 


open_files = []
for i in range(len(samples)):
    sample = samples[i]
    ubam = ubams[i]
    filename ='%s.txt'%sample
    if sample not in open_files:
        with open(filename,'w') as f:
            f.write("%s\n"%ubam)
            open_files.append(sample)
    else:
        with open(filename,'a') as f:
            f.write("%s\n"%ubam)

        

