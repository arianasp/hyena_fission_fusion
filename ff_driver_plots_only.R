## Driver for updating plots only

#----------------------------- SET UP DIRECTORIES---------------------------------------

user <- Sys.info()['user']
if(user == 'strau'){
  raw.data.directory <- '~/../Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/raw_data/'
  processed.data.directory <- '~/../Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/processed_data/'
  results.directory <- '~/../Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/results/'
  code.directory <- '~/../Documents/code/hyena_fission_fusion/'
}else if(user == 'straussed'){
  raw.data.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/raw_data/'
  processed.data.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/processed_data/'
  results.directory <- '~/Dropbox/Documents/Research/Full_projects/2021 Fission fusion social hubs/results/'
  code.directory <- '~/Documents/code/hyena_fission_fusion/'
}else{
  remote.stem <- '/Volumes/EAS_shared/'
  code.stem <- '~/Dropbox/code_ari/'
}

#----------------------------- SOURCE FUNCTIONS ---------------------------------------
source(paste0(code.directory, 'ff_functions_library.R'))


print('--------------------------- MAIN RESULTS ---------------------------------')
#Adjust subdirectory to be where to store extracted data based on parameters. 
R.fusion <- 100
R.fission <- 200
if(R.fusion == 100 & R.fission == 200){
  results.directory.param <- paste0(results.directory, '1_main_output_R100_200/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
} else {
  results.directory.param <- paste0(results.directory, 'stability_check_R', R.fusion, '_', R.fission, '/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
}
output.dirs <- c(data.outdir, plots.outdir)
generate_figures(data.outdir = output.dirs[1], plots.outdir = output.dirs[2], code.directory)

print('--------------------------- CHECK SMALLER THRESHOLD ---------------------------------')
#Adjust subdirectory to be where to store extracted data based on parameters. 
R.fusion <- 50
R.fission <- 100
if(R.fusion == 100 & R.fission == 200){
  results.directory.param <- paste0(results.directory, '1_main_output_R100_200/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
} else {
  results.directory.param <- paste0(results.directory, 'stability_check_R', R.fusion, '_', R.fission, '/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
}
output.dirs <- c(data.outdir, plots.outdir)
generate_figures(output.dirs[1], output.dirs[2], code.directory)

print('--------------------------- CHECK LARGER THRESHOLD ---------------------------------')
#Adjust subdirectory to be where to store extracted data based on parameters. 
R.fusion <- 200
R.fission <- 300
if(R.fusion == 100 & R.fission == 200){
  results.directory.param <- paste0(results.directory, '1_main_output_R100_200/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
} else {
  results.directory.param <- paste0(results.directory, 'stability_check_R', R.fusion, '_', R.fission, '/')
  data.outdir <- paste0(results.directory.param, 'data/')
  plots.outdir <- paste0(results.directory.param, 'plots/')
}
output.dirs <- c(data.outdir, plots.outdir)
generate_figures(output.dirs[1], output.dirs[2], code.directory)
