# 01_load_data
setwd('/projectnb/modislc/users/mkmoon/KoreaPhenology/data/')
for(tt in 1:3){
  system(paste('qsub -V -pe omp 2 -l mem_per_core=8G -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/KoreaPhenology/run_01.sh ',tt,sep=''))      
}  



