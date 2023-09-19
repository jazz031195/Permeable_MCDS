#!/bin/bash -l

process=""
loc=""
path=""
straight=""
concentration=""
time_steps=""
name=""
input=""
permeability=""

print_usage() {
  printf "Usage: 
  t : number of time steps
  c : concentration per mm3
  p : number of processes
  l : location (intra or extra)
  P : path to folder 
  S : if axons are straight : true, else false 
  n : name of file 
  I : name of input swc file 
  k : permeability "
}

while getopts t:c:n:p:l:P:S:I:k: opts; do
   case ${opts} in
      t) time_steps=${OPTARG} ;;
      c) concentration=${OPTARG} ;;
      p) process=${OPTARG} ;;
      l) loc=${OPTARG} ;;
      P) path=${OPTARG} ;;
      S) straight=${OPTARG} ;;
      n) name=${OPTARG} ;;
      I) input=${OPTARG} ;;
      k) permeability=${OPTARG} ;;
      *) print_usage
        exit 1 ;;
   esac
done

path_to_conf="MCDC_Simulator_public/instructions/axons/configuration_run.conf"
path_to_model="MCDC_Simulator_public/instructions/axons/configuration_model.conf"

cp "$path_to_model" "$path_to_conf"
path_to_scheme="$path/docs/scheme_files/PGSE_sample_scheme_new.scheme"

if [ "$straight" == "true" ]; then
  #path_to_axons="$path/instructions/demos/output/axons/straight/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}.swc"
  #path_to_exp="$path/instructions/demos/output/axons/straight/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}_${loc}"
  sed -i s:replace_tortuous:"false":g "${path_to_conf}"
else
  #path_to_axons="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}.swc"
  #path_to_exp="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}_${loc}"
  sed -i s:replace_tortuous:"true":g "${path_to_conf}"
fi
path_to_exp="$path/instructions/${name}"
path_to_axons="$path/instructions/${input}"
sed -i s:replace_exp_prefix:"${path_to_exp}":g "${path_to_conf}"
sed -i s:replace_scheme_file:"${path_to_scheme}":g "${path_to_conf}"
sed -i s:replace_axon_list:"${path_to_axons}":g "${path_to_conf}"
sed -i s:replace_process:"${process}":g "${path_to_conf}"
sed -i s:replace_loc:"${loc}":g "${path_to_conf}"
sed -i s:replace_concentration:"${concentration}":g "${path_to_conf}"
sed -i s:replace_timesteps:"${time_steps}":g "${path_to_conf}"
sed -i s:replace_permeability:"${permeability}":g "${path_to_conf}"