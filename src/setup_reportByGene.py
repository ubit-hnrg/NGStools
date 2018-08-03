import subprocess
requirements = ['pandas','pyensembl','mygene']

for package  in requirements:
	process = subprocess.Popen(['pip','install',package],stdout=subprocess.PIPE)
	print process.communicate()
