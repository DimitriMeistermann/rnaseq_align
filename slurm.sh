#SBATCH -o ./log/SGEuntrackedOut.txt
#SBATCH -e ./log/SGEuntrackedErr.txt
#SBATCH -p Bird 
#SBATCH -w bigmem001

{exec_job}