#!/usr/bin/csh
set start_time=`date +%s`

module restore postdoc_work

set rsno_params  = ( -1.1662 -1.5390 -1.5601 -1.3936 -1.2739 -1.7979 -1.6640 -1.9532 -1.7008 -1.4012 \
-1.2438 -1.6094 -1.3161 -1.5490 -1.4130 -1.7756 -1.6207 -1.3761 -1.5467 -1.8103 \
-1.1176 -1.1707 -1.8719 -1.1091 -1.5732 -1.9218 -1.7375 -1.1296 -1.1492 -1.5762 \
-1.3516 -1.3200 -1.0940 -1.2564 -1.8513 -1.8130 -1.6986 -1.1607 -1.5495 -1.6511 \
-1.8117 -1.4098 -1.8445 -1.5602 -1.2460 -1.4585 -1.4374 -1.7378 -1.7524 -1.6781 \
-1.3135 -1.6909 -1.8998 -1.8430 -1.3014 -1.9527 -1.5665 -1.5061 -1.9653 -1.3262 \
-1.5849 -1.8442 -1.5210 -1.7941 -1.9308 -1.1359 -1.5627 -1.5851 -1.0821 -1.2352 \
-1.1719 -1.7196 -1.5995 -1.9995 -1.4933 -1.3038 -1.1022 -1.4972 -1.5062 -1.2231 )
#########################
set ksno_params = ( 0.2418 0.2797 0.2871 0.2892 0.2274 0.3015 0.2640 0.2647 0.2345 0.3129 \
0.2914 0.2982 0.2993 0.3193 0.2781 0.3000 0.2967 0.2550 0.2988 0.2175 \
0.3119 0.3119 0.2560 0.3004 0.2719 0.2886 0.2267 0.2838 0.2139 0.3038 \
0.2696 0.2630 0.2690 0.2325 0.3110 0.2864 0.3021 0.2178 0.3172 0.2401 \
0.2910 0.2851 0.2354 0.2511 0.2847 0.2890 0.2520 0.2577 0.3165 0.2512 \
0.2104 0.3177 0.2373 0.2767 0.2604 0.2353 0.2367 0.2286 0.2736 0.2422 \
0.3087 0.2890 0.3072 0.2542 0.2809 0.2893 0.2250 0.2195 0.2424 0.3190 \
0.2977 0.2629 0.2587 0.3003 0.2807 0.2950 0.2111 0.3045 0.2447 0.2760 )


set datea     = ${1}
###
echo "Starting Post Assimilation Processes for $datea"
###
set year = `echo $datea | cut -b1-4`
set month = `echo $datea | cut -b5-6`
set day = `echo $datea | cut -b7-8`
set hour = `echo $datea | cut -b9-10`
###
set fore_date = `echo $datea +24 | ./advance_time`
set fyear = `echo $fore_date | cut -b1-4`
set fmonth = `echo $fore_date | cut -b5-6`
set fday = `echo $fore_date | cut -b7-8`
set fhour = `echo $fore_date | cut -b9-10`
###
set past_date = `echo $datea -24 | ./advance_time`
###
set RUN_ASSIM = True
set RUN_POST_ASSIM = True

if ( ${RUN_ASSIM} == True) then
  echo "Running Assimilation..."
endif
if ( ${RUN_POST_ASSIM} == True) then
  echo "Running Post Assimilation processes..."
endif

set run_dir = /glade/scratch/criedel/ICEPACK_RUNS/rundir
cd ${run_dir}
if ( ${RUN_ASSIM} == True) then
  set n = 1
  while ($n <= 80)
    set icnum = `echo $n + 10000 | bc | cut -b2-5`
    if ($n == 64) then
      cd mem$icnum
      echo "mv iced.${year}-${month}-${day}-00000.nc truth_data.nc"
      mv iced.${year}-${month}-${day}-00000.nc truth_data.nc
      cd ../
      @ n++
      continue
    endif
    cd mem$icnum
    echo "mv iced.${year}-${month}-${day}-00000.nc dart_restart.nc"
    mv iced.${year}-${month}-${day}-00000.nc dart_restart.nc
    cd ../
    @ n++
  end
  rm obs_seq.out
  ln -s /glade/scratch/criedel/ICEPACK_RUNS/OBS_SEQ_FILES/obs_seq.out_${year}${month}${day} obs_seq.out 
  ln -s /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${past_date}/output_priorinf_mean.nc input_priorinf_mean.nc
  ln -s /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${past_date}/output_priorinf_sd.nc input_priorinf_sd.nc

  echo "Running the filter..."
  ./filter >& filter.out 
  
  set check_finished = `grep -q 'Finished' filter.out  && echo $?`
  if ( ${check_finished} != 0) then
    echo "Filter did not finish running"
    touch ${run_dir}/STOP
    exit
  endif
  rm input_priorinf_*.nc
  echo "Filter is done running."
endif
##
if ( ${RUN_POST_ASSIM} == True) then
  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}
  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses
  mkdir /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts
  mv filter.out /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/
  ############
  cp input_*.nc output_*.nc obs_seq.final /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/
  sleep 2
  set check_old = `du -csh --apparent-size input_*.nc output_*.nc obs_seq.final | grep 'total' | cut -b1-2`
  set check_new = `du -csh --apparent-size /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/$datea/*.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/$datea/obs_seq.final | grep 'total' | cut -b1-2`
  echo ${check_old}
  echo ${check_new}
  if ( ${check_old} != ${check_new} ) then
    echo "Files did not copy correctly"
    touch ${run_dir}/STOP
    exit
  else
    rm input_*.nc output_*.nc obs_seq.final
  endif
  ############
  set n = 1
  while ($n <= 80)
    set icnum = `echo $n + 10000 | bc | cut -b2-5`
    if ($n == 64) then
      echo "Formatt truth member file...."
      cd mem$icnum
      echo "mv truth_data.nc iced.${year}-${month}-${day}-00000.nc"
      mv truth_data.nc iced.${year}-${month}-${day}-00000.nc
      cp iced.${year}-${month}-${day}-00000.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses/iced.${year}-${month}-${day}-00000.nc_$icnum
      #####
      rm icepack_in
      cat >! script.sed << EOF
      /RESTART_FILENAME/c\
      ice_ic         = iced.${year}-${month}-${day}-00000.nc
      /KSNOW/c\
      ksno              = ${ksno_params[$n]}
      /RSNOW/c\
      R_snw             = ${rsno_params[$n]}
      /MEMA/c\
      atm_data_file   = 'ATM_FORCING_${icnum}.txt'
      /MEMO/c\
      ocn_data_file   = 'OCN_FORCING_${icnum}.txt'
EOF
      sed -f script.sed ../icepack_in.template >! icepack_in
      echo "./icepack >& icepack.out"
      ./icepack >& icepack.out
      sleep 1
      ##
      set check_finished = `grep -q 'ICEPACK COMPLETED SUCCESSFULLY' icepack.out && echo $?`
      if ($check_finished != 0) then
        echo "ICEPACK DID NOT FINISH....STOP!"
        touch ${run_dir}/STOP
        exit
      endif
      if (! -e restart/iced.${fyear}-${fmonth}-${fday}-00000.nc) then
        echo "ICEPACK DID NOT CREATE RESTART FILE....STOP!"
        touch ${run_dir}/STOP
        exit
      endif
      mv icepack.out /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/icepack.out_$icnum
      echo "mv history/icepack.h.${year}${month}${day}.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/icepack.h.${year}${month}${day}.nc_$icnum"
      mv history/icepack.h.${year}${month}${day}.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/icepack.h.${year}${month}${day}.nc_$icnum
      mv ice_diag.full_ITD /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/ice_diag.full_ITD_$icnum
      rm ice_diag.*
      rm iced.${year}-${month}-${day}-00000.nc
      echo "mv restart/iced.${fyear}-${fmonth}-${fday}-00000.nc ."
      mv restart/iced.${fyear}-${fmonth}-${fday}-00000.nc .
      ########################
      cd ../
      @ n++
      continue
    endif
    cd mem$icnum 
    cp dart_restart.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses/dart_restart.nc_$icnum
    cp restart_state.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses/restart_state.nc_$icnum
    cp ../input.nml .
    rm update_state.out
    echo "/glade/work/criedel/DAtests/CICE-SCM/DART_JEFF/DART_JEFF/DART/models/cice-scm/work/dart_to_cice >& update_state.out"
    /glade/work/criedel/DAtests/CICE-SCM/DART_JEFF/DART_JEFF/DART/models/cice-scm/work/dart_to_cice >& update_state.out
    set check_finished2 = `grep -q 'Finished' update_state.out  && echo $?`
    if ( ${check_finished2} != 0) then
      echo "dart_to_cice did not finish...STOP!"
      touch ${run_dir}/STOP
      exit
    endif
    echo "mv dart_restart.nc iced.${year}-${month}-${day}-00000.nc"
    mv dart_restart.nc iced.${year}-${month}-${day}-00000.nc
    echo "cp iced.${year}-${month}-${day}-00000.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses/iced.${year}-${month}-${day}-00000.nc_$icnum"
    cp iced.${year}-${month}-${day}-00000.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/analyses/iced.${year}-${month}-${day}-00000.nc_$icnum 
    rm restart_state.nc
    rm icepack_in
    cat >! script.sed << EOF
    /RESTART_FILENAME/c\
    ice_ic         = iced.${year}-${month}-${day}-00000.nc   
    /KSNOW/c\
    ksno              = ${ksno_params[$n]}
    /RSNOW/c\
    R_snw             = ${rsno_params[$n]}
    /MEMA/c\
    atm_data_file   = 'ATM_FORCING_${icnum}.txt'
    /MEMO/c\
    ocn_data_file   = 'OCN_FORCING_${icnum}.txt'
EOF
  sed -f script.sed ../icepack_in.template >! icepack_in 
    echo "./icepack >& icepack.out"
    ./icepack >& icepack.out
    sleep 1
    ##
    set check_finished = `grep -q 'ICEPACK COMPLETED SUCCESSFULLY' icepack.out && echo $?`
    if ($check_finished != 0) then
      echo "ICEPACK DID NOT FINISH....STOP!"
      touch ${run_dir}/STOP
      exit
    endif
    if (! -e restart/iced.${fyear}-${fmonth}-${fday}-00000.nc) then
      echo "ICEPACK DID NOT CREATE RESTART FILE....STOP!"
      touch ${run_dir}/STOP
      exit
    endif 
    mv icepack.out /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/icepack.out_$icnum
    echo "mv history/icepack.h.${year}${month}${day}.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/icepack.h.${year}${month}${day}.nc_$icnum"
    mv history/icepack.h.${year}${month}${day}.nc /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/icepack.h.${year}${month}${day}.nc_$icnum
    mv ice_diag.full_ITD /glade/scratch/criedel/ICEPACK_RUNS/ARCHIVE/${datea}/forecasts/ice_diag.full_ITD_$icnum
    rm ice_diag.*
    rm iced.${year}-${month}-${day}-00000.nc
    echo "mv restart/iced.${fyear}-${fmonth}-${fday}-00000.nc ."
    mv restart/iced.${fyear}-${fmonth}-${fday}-00000.nc .
    cd ../
    @ n++
  end
  echo "Done with post assimilation step!"
endif
echo "Done cycling for $datea"
set end_time=`date +%s`
echo "execution time was `expr $end_time - $start_time` s."




