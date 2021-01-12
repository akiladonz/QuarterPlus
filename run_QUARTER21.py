

''' 
created by Sina Ghadermarzi 

2020 July 6

gets a fasta file that (for now) has just one entry and runs the whole QUARTER2 process up to the final result.csv or html file
the results will be saved in the same folder as the input file but with names result.csv and result.html

'''





import sys
import os
from pywrappers.py_asaquick import *
from pywrappers.py_disembl import *
from pywrappers.py_disopred3 import *
from pywrappers.py_espritz import *
from pywrappers.py_globpipe import *
from pywrappers.py_iupred import *
from pywrappers.py_spot import *
from pywrappers.py_vsl2b import *
from pywrappers.py_asaquick import *

predictor = {
	"0":	disembl_rem465,
	"1":	disembl_hotloops,
	"2":	espritz_disorder,
	"3":	espritz_nmr,
	"4":	espritz_xray,
	"5":	globpipe,
	"6":	iupred_long,
	"7":	iupred_short,
	"8":	spot,
	"9":	vsl2b,
	"10":   disopred3
}

thresholds = {
	"0":	0.873,
	"1":	0.713,
	"2":	0.593,
	"3":	0.819,
	"4":	0.783,
	"5":	0.708,
	"6":	0.831,
	"7":	0.875,
	"8":	None,
	"9":	0.865,
	"10":	None	
}


def gen_taskid(prog_name):
	random_code= (''.join(choice(digits) for i in range(8)))
	time_string=  datetime.now().strftime("%Y%m%d%H%M%S")
	taskid = time_string+"_"+prog_name+"_"+ random_code
	return taskid


###first create a temporary folder for everything
this_script_directory = os.path.dirname(os.path.realpath(__file__))
taskid = gen_taskid("QUARTER2")
work_folder = this_script_directory+"/tmp"+taskid
os.mkdir(work_folder)
fastafile_address = work_folder+"/seqs.fasta"

with open(work_folder+"/command.txt","w") as commandf:
	commandf.writelines(" ".join(sys.argv))

if len(sys.argv)==3:
	method_idx = sys.argv[1]
	infile_address  = sys.argv[2]
	with open(infile_address) as inf, open(work_folder+"/seqs.fasta","w") as fastaf_fixed:
		lines = inf.read().replace("\r","")
		fastaf_fixed.writelines("\n".join(lines[1:]))
elif len(sys.argv)==2:
	input_from_frontend_txtfile = sys.argv[1]
	with open(input_from_frontend_txtfile) as inf, open(work_folder+"/seqs.fasta","w") as fastaf_fixed:
		lines = inf.read().replace("\r","").strip("\n").split("\n")
		method_idx = lines[0]
		fastaf_fixed.writelines("\n".join(lines[1:]))
	fastafile_address = work_folder+"/seqs.fasta"
else:
	print("error")
	sys.exit("Invalid input arguments.")



##############################################






# prepare inputs .rsa .seg and .seq files for the test.sh script
####first read the fasta file

with open(fastafile_address) as fastaf:
	fasta_spl = fastaf.read().strip("\n").split("\n")
	pid= fasta_spl[0][1:]
	seq = "".join(fasta_spl[1:])
###for now, the above code will just read one sequence from the fasta and everything will be put in the work_folder
##########################################################

### write the .seq file
with open(work_folder + "/seq.seq","w") as seqf:
	seqf.writelines(pid+"\t"+seq+"\n")
##########################################################





###run the seg tool, read the output, and put it in the right .seg format
cmd = "/home/biomine/programs/Wootton/seg " + seq + " -x > " + work_folder + "/seq.cmplx"
# print(cmd)
os.system(cmd)
with open(work_folder + "/seq.cmplx") as tmp_cmplx_file, open(work_folder + "/seq.seg","w") as segf:
	res = tmp_cmplx_file.read().strip("\n")
	res  = "".join(res.split("\n"))
	segf.writelines(pid+"\t"+res+"\n")
os.system("rm "+work_folder + "/seq.cmplx")
###########################################################



### run asaquick and put the result in the .rsa file

with  open(work_folder + "/seq.rsa","w") as rsaf:
	res = asaquick(seq)
	res  = " ".join([str(s) for s in res])
	rsaf.writelines(pid+" "+res+"\n")

#########################################################


###now all the input files are ready and we can run the script
if method_idx != "11":
	##calculate disorder prediction by the one given predictor 
	pred_values  = predictor[method_idx](seq)
	pred_strs = [str(x) for x in pred_values]
	with open (work_folder + "/"+method_idx+".dis","w") as disf:
		disf.writelines(pid+"\t"+ " ".join(pred_strs))

	cmd = this_script_directory +"/calculate_QAs.sh"
	cmd+= " "+ work_folder + "/seq.seq"
	cmd+= " "+ work_folder + "/"+method_idx+".dis"
	cmd+= " "+ work_folder + "/seq.rsa"
	cmd+= " "+ work_folder + "/seq.seg"
	cmd+= " "+ method_idx
	cmd+= " "+ work_folder
	os.system(cmd)

	#read the result and put them into one file
	out_str=">"+pid+"\n"+seq+"\n"
	# read the disorder prediction for method number i into an array
	with open(work_folder + "/" +method_idx+".dis") as disf:
		dis_spl = disf.read().rstrip("\n").split("\t")[1].split()
		out_str  += ",".join(dis_spl)+"\n"
	# read the QA results for method i in to an array
	with open(work_folder + "/" +method_idx+".result") as qaf:
		qa_spl = qaf.read().rstrip("\n").split("\n")
		out_str  += ",".join(qa_spl)+"\n"
	#the result will be in method_idx.result: we can read and present it together with predictions
else:
	##calculate the disorder prediction for all predictors and run QA script for all of them
	for i in range(11):
		curr_idx  = str(i)
		pred_values  = predictor[curr_idx](seq)
		pred_strs = [str(x) for x in pred_values]
		with open (work_folder + "/"+curr_idx+".dis","w") as disf:
			disf.writelines(pid+"\t"+ " ".join(pred_strs))
		cmd = this_script_directory +"/calculate_QAs.sh"
		cmd+= " "+ work_folder + "/seq.seq"
		cmd+= " "+ work_folder + "/"+curr_idx+".dis"
		cmd+= " "+ work_folder + "/seq.rsa"
		cmd+= " "+ work_folder + "/seq.seg"
		cmd+= " "+ curr_idx
		cmd+= " "+ work_folder
		os.system(cmd)
	## after calculating disorder and QA predictions for all predictors, we combine the result into one file
	out_str=">"+pid+"\n"+seq+"\n"
	for i in range(11):
		# read the disorder prediction for method number i into an array
		with open(work_folder + "/" +str(i)+".dis") as disf:
			dis_spl = disf.read().rstrip("\n").split("\t")[1].split()
			out_str  += ",".join(dis_spl)+"\n"
		# read the QA results for method i in to an array
		with open(work_folder + "/" +str(i)+".result") as qaf:
			qa_spl = qaf.read().rstrip("\n").split("\n")
			out_str  += ",".join(qa_spl)+"\n"
		# out_str = out_str.rstrip("\n")


with open(work_folder + "/disqa.result","w") as outf:
	outf.writelines(out_str)