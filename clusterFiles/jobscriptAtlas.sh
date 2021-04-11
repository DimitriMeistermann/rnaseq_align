#$ -cwd
#$ -V
#$ -q hugemem.q
#$ -l mem_free=40G
#$ -e  ./log/untrackedErr.txt
#$ -o  ./log/untrackedOut.txt

# Usually, all log are stored in the output folder, logs specific to SGE are stored in the Snakefile directory in ./log

{exec_job}
