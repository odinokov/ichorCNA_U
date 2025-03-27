import os,sys
scriptHome = os.path.realpath(__file__).split('/')[:-1]
scriptHome = '/'.join(scriptHome)
sys.path.append('{0}/src/'.format(scriptHome))

from multiprocessing import Pool
import subprocess
import time

usage = """
<yaml input config file>
"""

argv = sys.argv[1:]

try:
	yamlInput = argv.pop(0)
except:
	sys.exit(usage)


#
#	Helper functions for job managment
#
def runCmd(cmd):
	import os
	os.system(cmd)

def runScript(cmd):
	import os
	os.system('bash {0}'.format(cmd))

def multithread(cmds,func,numCores):
	from multiprocessing import Pool
	p=Pool(int(numCores))
	for c in cmds:
		p.apply_async(func, args=(c,))
	p.close()
	p.join()
	return

def launchJobs(inFolder):
	files = os.listdir(inFolder)
	for f in files:
		time.sleep(0.25)
		os.system('sbatch {0}/{1}'.format(inFolder,f))
	return

def run_jobs(inFolder,jobName="none",slurm=True,username="none"):
	if slurm == False:
		cmds = os.listdir(inFolder)
		scriptPathCmds = ['{0}/{1}'.format(inFolder,x) for x in cmds]
		if sData.cores > 1:
			multithread(scriptPathCmds,runScript,sData.cores)
		else:
			for c in scriptPathCmds:
				os.system('bash {0}'.format(c))
	else:
		launchJobs(inFolder)
		while areJobsRunning(jobName,username):
			time.sleep(60)
	return

def areJobsRunning(jobName,username):
	cmd = 'squeue -u {0} -O name:50'.format(username)
	jobs = subprocess.check_output(cmd,shell=True).decode('utf-8')
	jobs = jobs.split( )[1:]
	for j in jobs:
		if jobName in j:
			return True
	return False

#
#	Input data structure
#

class inData(object):
	def __init__(self,f,scriptHome):
		import yaml
		with open(f) as fData:
			yam = yaml.safe_load(fData)

		self.runDepthCalc = yam['runDepthCalc']
		self.runIchor = yam['runIchorCNA']
		self.runStats = yam['calcStats']
		self.queue = yam['slurmQueue']
		self.username = yam['slurmUser']

		self.chrs = yam['chrs']
		self.binSize = yam['binSize']
		self.normalStart = yam['initial_normal_proportion']
		self.ploidy = yam['ploidy'] 
		self.scStates = yam['scStates']
		self.maxCN = yam['maxCN']
		self.includeHOMD = yam['includeHOMD']
		self.txnE = yam['txnE']
		self.txnStrength = yam['txnStrength']
		self.plotExt = yam['plotFileType']
		self.plotYlim = yam['plotYlim']
		self.cores = yam['cores']

		self.scriptHome = scriptHome
		self.outFolder = yam['outputDir']
		self.rscriptPath = yam['RscriptPath']
		self.bedtoolsPath = yam['bedtoolsPath']
		self.juliaPath = yam["juliaPath"]
		self.bedDepthPath = '{0}/src/getBedDepth.jl'.format(scriptHome)
		if yam['ichorCNAPath'] == "default":
			self.ichorCNAPath = '{0}/src_ichorCNA/'.format(scriptHome)
		else:
			self.ichorCNAPath = yam['ichorCNAPath']

		self.dataType = yam["dataType"].lower()
		self.chromSizeFile = yam["chromSizes"]
		self.refGenome = yam["referenceGenome"]
		self.normal = yam["normalCtrl"]
		self.gcWig = yam["gcWig"]
		self.mapWig = yam["mapWig"]

		self.estimateNormal = yam['estimateNormal']
		self.estimatePloidy = yam['estimatePloidy']
		self.estimateClonality = yam['estimateClonality']

		self.ploidy = [str(x) for x in self.ploidy]
		self.scStates = [str(x) for x in self.scStates]
		self.normalStart = [str(x) for x in self.normalStart]
		self.maxSubCNA = yam['maxSubcloneCNA']
		self.maxSub = yam['maxSubcloneFrac']
		self.blacklist = yam["blacklist"]

		self.samples = {}
		files = os.listdir(yam["inputDir"])
		files = [x for x in files if ".bai" not in x]
		files = [x for x in files if ".crai" not in x]
		files = [x for x in files if ".DS_Store" not in x]
		files = [x for x in files if ".tbi" not in x]
		for f in files:
			if self.dataType == "xam":
				n = f.replace(".bam","")
				n = n.replace(".cram","")
				self.samples[n] = os.path.abspath("{0}/{1}".format(yam["inputDir"],f))
			elif self.dataType == "bed":
				n = f.replace(".bed","")
				n = n.replace(".tsv","")
				n = n.replace(".bgz","")
				n = n.replace(".gz","")
				n = n.replace(".zip","")
				self.samples[n] = os.path.abspath("{0}/{1}".format(yam["inputDir"],f))

class Result(object):
	def __init__(self):
		self.name = ''
		self.purity = 0.0
		self.ploidy = 0.0
		self.scFrac = 0.0
		self.coverage = 0.0
		self.scGenome = 0.0
		self.scCNA = 0.0
		return

def write_beddepth(sData):
	chroms = [str(x) for x in sData.chrs]
	for s in sData.samples:
		inFile = sData.samples[s]
		mem = 50000
		samCmd1 = ""
		samCmd2 = ""
		args = [sData.bedDepthPath,sData.binSize,','.join(chroms),inFile,sData.outFolder,s,
			sData.bedtoolsPath,sData.chromSizeFile,sData.juliaPath]
		hcCmd = """{8} {0} \\
	{3} \\
	{1} \\
	{7} \\
	{4}/wigFiles/{5}.wig \\
	--chroms {2} \\
	--toolpath {6} \n\n""".format(*args)

		out = open("{0}/jobScripts/binDepth/{1}.depth.bed.slurm".format(sData.outFolder,s),"w")

		out.write("""#!/bin/bash
#
#SBATCH -p {2}		# partition (queue)
#SBATCH -N 1		# number of nodes
#SBATCH -n 1		# number of cores
#SBATCH --mem {3}		# memory (in mb)
#SBATCH -t 0-05:00		# time (D-HH:MM)
#SBATCH -o {1}/logs/bedDepth_{0}.log
#SBATCH --job-name=bedDepth_{0}
#SBATCH --chdir={1}\n\n""".format(s,sData.outFolder,sData.queue,mem))

		out.write(hcCmd)
		out.close()
	return


def write_ichorCNA(sData):
	scriptPaths = []
	if sData.normal == "":
		if sData.binSize == 1000:
			normalPanel = '{0}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds'.format(sData.ichorCNAPath)
		else:
			normalPanel = '{0}/inst/extdata/HD_ULP_PoN_{1}kb_median_normAutosome_mapScoreFiltered_median.rds'.format(sData.ichorCNAPath,sData.binSize)
	else:
		normalPanel = sData.normal

	allChr = [str(x) for x in sData.chrs]
	if sData.gcWig == "" and sData.ichorCNAPath != "None":
		allChr = [x.replace("chr","") for x in allChr]

	for s in sData.samples:
		if sData.gcWig == "" and sData.ichorCNAPath != "None":
			gcWig = "{0}/inst/extdata/gc_hg19_{1}kb.wig".format(sData.ichorCNAPath,sData.binSize)
		else:
			gcWig = sData.gcWig
		if sData.mapWig == "" and sData.ichorCNAPath != "None":
			mapWig = "{0}/inst/extdata/map_hg19_{1}kb.wig".format(sData.ichorCNAPath,sData.binSize)
		else:
			mapWig = sData.mapWig
		if sData.blacklist == "" or sData.blacklist == False:
			blacklist = "NULL"
		else:
			blacklist = sData.blacklist
		allChrSet = "\\\"" + '\\\",\\\"'.join(allChr) + "\\\""

		args = [s,sData.outFolder,sData.binSize,allChrSet,max(sData.chrs),
						sData.txnE,sData.txnStrength,','.join(sData.normalStart),','.join(sData.ploidy),
						sData.maxCN,sData.estimateNormal,sData.estimatePloidy,sData.estimateClonality,
						','.join(sData.scStates),sData.includeHOMD,sData.ichorCNAPath,sData.scriptHome,
						normalPanel,gcWig,mapWig,sData.maxSubCNA,sData.maxSub,blacklist,sData.rscriptPath]

		scriptPaths.append("{0}/jobScripts/ichorCNA/{1}.slurm".format(sData.outFolder,s))
		out = open("{0}/jobScripts/ichorCNA/{1}.slurm".format(sData.outFolder,s),"w")


		out.write("""#!/bin/bash
#
#SBATCH -p {2}		# partition (queue)
#SBATCH -N 1		# number of nodes
#SBATCH -n 1		# number of cores
#SBATCH --mem 5000		# memory (in mb)
#SBATCH -t 0-05:00		# time (D-HH:MM)
#SBATCH -o {1}/logs/ichorCNA_{0}.log	#Logfile
#SBATCH --job-name=ichorCNA_{0}
#SBATCH --mail-type=FAIL
#SBATCH --chdir={1}\n\n""".format(s,sData.outFolder,sData.queue))

		out.write("""{23} {16}/src/runIchorCNA.R \\
	--id {0} \\
	--WIG {1}/wigFiles/{0}.wig \\
	--ploidy "c({8})" \\
	--normal "c({7})" \\
	--maxCN {9} \\
	--gcWig  {18} \\
	--mapWig {19}  \\
	--normalPanel {17} \\
	--libdir {15} \\
	--chrs "c({3})" \\
	--chrTrain "c({3})" \\
	--chrNormalize "c({3})" \\
	--estimateNormal {10} \\
	--estimatePloidy {11} \\
	--estimateScPrevalence {12} \\
	--scStates "c({13})" \\
	--txnE {5} \\
	--txnStrength {6} \\
	--outDir {1}/ichorOut \\
	--includeHOMD {14} \\
	--maxFracCNASubclone {20} \\
	--maxFracGenomeSubclone {21} \\
	--centromere {22}\n""".format(*args))
		out.close()
	return

def read_outputs(sData):
	outFolder = "{0}/ichorOut".format(sData.outFolder)
	outputs = os.listdir(outFolder)
	RES = {}
	for f in outputs:
		if f.endswith(".params.txt"):
			sample = f.replace('.params.txt','')
			RES[sample] = Result()
			RES[sample].name = sample
			with open('{0}/{1}'.format(outFolder,f),'r') as results:
				for ln in results:
					ln = ln.rstrip()
					if "Tumor Fraction" in ln:
						ln = ln.split(':')[-1]
						ln = ln.replace(" ","")
						ln = ln.replace("\t","")
						RES[sample].purity = ln
					if "Ploidy" in ln:
						ln = ln.split(':')[-1]
						ln = ln.replace(" ","")
						ln = ln.replace("\t","")
						RES[sample].ploidy = ln
					if "Subclone Fraction" in ln:
						ln = ln.split(':')[-1]
						ln = ln.replace(" ","")
						ln = ln.replace("\t","")
						RES[sample].scFrac = ln
					if "Fraction Genome Subclonal" in ln:
						ln = ln.split(':')[-1]
						ln = ln.replace(" ","")
						ln = ln.replace("\t","")
						RES[sample].scGenome = ln
					if "Fraction CNA Subclonal" in ln:
						ln = ln.split(':')[-1]
						ln = ln.replace(" ","")
						ln = ln.replace("\t","")
						RES[sample].scCNA = ln

	out = open('{0}/summary_table.tsv'.format(sData.outFolder),'w')
	header = ["SampleName","TumorFraction","Ploidy","SubcloneFraction","SubcloneGenomeFraction",
		"FractionCNASubclonal"]
	out.write('\t'.join(header) + '\n')
	for r in sorted(RES):
		args = [RES[r].name,RES[r].purity,RES[r].ploidy,RES[r].scFrac,
			RES[r].scGenome]
		args = [str(x) for x in args]
		out.write('\t'.join(args) + '\n')
	out.close()
	return

def build_gisic_inputs(inFolder,summFile,markerCount,outFile,minTumor=0.001):
	import os
	SAMPLE = {}
	with open(summFile,'r') as f:
		for ln in f:
			if 'SampleName' in ln:
				continue
			ln = ln.rstrip()
			ln = ln.split()
			SAMPLE[ln[0]] = float(ln[1])
	files = os.listdir(inFolder)
	files = [x for x in files if '.cna.seg' in x]
	DATA = []
	markerCount = str(markerCount)
	for f in files:
		sName = f.replace(".cna.seg","")
		if sName in SAMPLE:
			if SAMPLE[sName] >= minTumor:
				with open('{0}/{1}'.format(inFolder,f),'r') as data:
					for ln in data:
						if 'chr	start' in ln:
							continue
						ln = ln.rstrip()
						ln = ln.split()
						lnData = [sName,ln[0],ln[1],ln[2],markerCount,ln[5]]
						DATA.append(lnData)
	out = open(outFile,'w')
	for i in DATA:
		ln = '\t'.join(i)
		out.write(ln + '\n')
	out.close()
	return

####################################################################################################

sData = inData(yamlInput,scriptHome)
os.system('mkdir -p {0}'.format(sData.outFolder))
os.system('mkdir -p {0}/wigFiles'.format(sData.outFolder))
os.system('mkdir -p {0}/stats'.format(sData.outFolder))
os.system('mkdir -p {0}/temp'.format(sData.outFolder))
os.system('mkdir -p {0}/ichorOut'.format(sData.outFolder))
os.system('mkdir -p {0}/jobScripts'.format(sData.outFolder))
os.system('mkdir -p {0}/jobScripts/ichorCNA'.format(sData.outFolder))
os.system('mkdir -p {0}/jobScripts/binDepth'.format(sData.outFolder))
os.system('mkdir -p {0}/logs'.format(sData.outFolder))
os.system('cp {0} {1}/'.format(yamlInput,sData.outFolder))

print("\n")

write_beddepth(sData)
write_ichorCNA(sData)

if sData.runDepthCalc == True:
	print("Running bin depth counters...")
	if sData.queue.upper() == "LOCAL":
		run_jobs("{0}/jobScripts/binDepth".format(sData.outFolder),slurm=False)
	elif sData.queue != False:
		run_jobs("{0}/jobScripts/binDepth".format(sData.outFolder),
			jobName="bedDepth",username=sData.username)
	if sData.gcWig == "" and sData.ichorCNAPath != "None":
		os.system("grep -rl 'chrom=chr' {0}/wigFiles | xargs sed -i 's/chrom=chr/chrom=/g'".format(sData.outFolder))
	os.system("rm {0}/*.ichorCNAtempFixedBins.bed".format(sData.outFolder))

if sData.runIchor == True:
	print("Running ichorCNA...")
	if sData.queue.upper() == "LOCAL":
		run_jobs("{0}/jobScripts/ichorCNA".format(sData.outFolder),slurm=False)
	else:
		run_jobs("{0}/jobScripts/ichorCNA".format(sData.outFolder),
			jobName="ichorCNA",username=sData.username)
if sData.runStats == True:
	read_outputs(sData)
