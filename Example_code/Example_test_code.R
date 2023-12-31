# Test code

##############################

# in R install the libraries

install.packages("devtools")

library(devtools)

install_github("celphin/RepeatOBserverV1") #to install the package

library(RepeatOBserverV1) # to load the package


#############################
# on the terminal 

# change to directory you want to work in
cd ~/scratch/repeats/auto_script/

# download the Setup_Run_Repeats.sh
wget https://github.com/celphin/RepeatOBserverV1/blob/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

#-----------------------------------
# download genome, Citrus limon
wget https://download.cncb.ac.cn/gwh/Plants/Citrus_limon_lemon_GWHCBFQ00000000.1/GWHCBFQ00000000.1.genome.fasta.gz

gunzip GWHCBFQ00000000.1.genome.fasta.gz

mv GWHCBFQ00000000.1.genome.fasta Lemon.fasta

#---------------------------
# make slurm script

cat << EOF > Auto_Lemon.sh
#!/bin/bash
#SBATCH --account=<your-account>
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Lemon -f Lemon.fasta -h H0 -c 15 -m 128000M

EOF

sbatch Auto_Lemon.sh

###########################
# look at the output in:
cd ~/scratch/repeats/auto_script/output_chromosomes/Lemon_H0/Summary_output