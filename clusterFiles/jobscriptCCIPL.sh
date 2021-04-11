#SBATCH -o ./log/untrackedOut.txt
#SBATCH -e ./log/untrackedErr.txt
#SBATCH -p Bird 
#SBATCH -w bigmem001

{exec_job}