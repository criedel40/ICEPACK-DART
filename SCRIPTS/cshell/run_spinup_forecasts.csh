#!/usr/bin/csh

set storage_loc = /glade/scratch/criedel/ICEPACK_RUNS/SPINUP_FORE/VARYATM_CONSOCN 
set sim_length = 78840

set rsno_params = ( -1.3196 -1.5510 -1.1303 -1.5012 -1.5114 -1.7008 -1.0834 -1.4158 -1.5610 -1.5157 \
                    -1.7532 -1.5099 -1.4833 -1.4996 -1.7797 -1.5947 -1.6612 -1.5463 -1.4264 -1.6242 \
                    -1.4833 -1.2607 -1.0970 -1.3577 -1.6344 -1.2499 -1.5865 -1.3744 -1.4356 -1.4785 \
                    -1.6796 -1.8481 -1.6594 -1.4105 -1.6255 -1.9362 -1.4469 -1.6827 -1.6328 -1.5340 \
                    -1.5606 -1.5870 -1.3229 -1.5983 -1.2586 -1.7795 -1.3200 -1.5444 -1.7076 -1.5537 \
                    -1.4274 -1.4003 -1.2811 -1.3838 -1.4599 -1.5082 -1.7577 -1.7889 -1.5483 -1.2722 \
                    -1.4947 -1.5546 -1.2938 -1.3353 -1.4075 -1.0184 -1.3902 -1.6125 -1.5433 -1.8346 \
                    -1.7459 -1.4164 -1.3006 -1.4901 -1.2707 -1.6546 -1.4071 -1.5818 -1.6313 -1.7968 )

set ksno_params = ( 0.2208 0.2346 0.2751 0.3183 0.3099 0.2786 0.3225 0.2519 0.2706 0.2204 \
                    0.2214 0.2581 0.2514 0.1997 0.2517 0.1879 0.2166 0.2487 0.2273 0.3338 \
                    0.2809 0.2271 0.2031 0.2402 0.2116 0.2188 0.3326 0.2507 0.3112 0.2942 \
                    0.3110 0.2583 0.1806 0.2980 0.3298 0.2917 0.2569 0.2532 0.2248 0.2561 \
                    0.2604 0.2541 0.2645 0.2531 0.2441 0.2065 0.1903 0.2089 0.2108 0.3043 \
                    0.2915 0.2255 0.2034 0.2981 0.2928 0.1964 0.3155 0.2012 0.1982 0.2702 \
                    0.2268 0.2154 0.2126 0.2510 0.3117 0.2206 0.1871 0.2359 0.2887 0.3055 \
                    0.2528 0.2372 0.2397 0.2354 0.3238 0.3297 0.2611 0.1664 0.2644 0.1926 )

set sim_num = 26

while ($sim_num <= 26)
  echo "-----------------------------------"
  echo "Working on simulation => ${sim_num}"
  ###########
  #PREP WORK
  echo "mkdir ${storage_loc}/spinupyear${sim_num}"
  mkdir ${storage_loc}/spinupyear${sim_num}
  ##########################################
  set mem = 1
  while ($mem <= 80)
    echo "running member $mem"
    set inst_string = `printf %04d $mem`
    if ($sim_num == 1) then
      echo "mkdir mem${inst_string}"
      mkdir mem${inst_string}
      echo "mkdir mem${inst_string}/history"
      mkdir mem${inst_string}/history
      echo "mkdir mem${inst_string}/restart"
      mkdir mem${inst_string}/restart
      echo "ln -s /glade/scratch/criedel/ICEPACK_RUNS/TEST_SETUP/icepack mem${inst_string}/"
      ln -s /glade/scratch/criedel/ICEPACK_RUNS/TEST_SETUP/icepack mem${inst_string}/
    endif
    cd mem${inst_string}
    if ($sim_num == 1) then
      set restart_flag = .false.
    else
      set restart_flag = .true.
    endif
    cat >! script.sed << EOF
    /SIM_LEN/c\
    npt              = $sim_length
    /RESTART_LOGIC/c\
    restart           = $restart_flag
    /KSNOW/c\
    ksno              = ${ksno_params[$mem]}
    /RSNOW/c\
    R_snw             = ${rsno_params[$mem]}
EOF
  sed -f script.sed ../icepack_in.template >! icepack_in

     if ($sim_num > 1) then
       rm iced.2012-01-01-00000.nc 
       set prev_run = `expr ${sim_num} - 1`
       echo "cp ${storage_loc}/spinupyear${prev_run}/mem${inst_string}/iced.2012-01-01-00000.nc ."
       cp ${storage_loc}/spinupyear${prev_run}/mem${inst_string}/iced.2012-01-01-00000.nc .
     endif
     rm icepack.out script.sed
     ./icepack >& icepack.out
     set check_finished = `grep -q 'ICEPACK COMPLETED SUCCESSFULLY' icepack.out && echo $?`
     if ($check_finished != 0) then
       echo "ICEPACK DID NOT FINISH....STOP!"
       exit
     endif   
     if (! -e restart/iced.2012-01-01-00000.nc) then
       echo "ICEPACK DID NOT FINISH....STOP!"
       exit
     endif
    
     echo "mkdir ${storage_loc}/spinupyear${sim_num}/mem${inst_string}"
     mkdir ${storage_loc}/spinupyear${sim_num}/mem${inst_string}
     echo "mv ice_diag.* ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/"
     mv ice_diag.* ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/
     echo "mv history/icepack.h.*.nc ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/"
     mv history/icepack.h.*.nc ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/
     echo "mv restart/iced.*.nc ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/"
     mv restart/iced.*.nc ${storage_loc}/spinupyear${sim_num}/mem${inst_string}/
     cd ../
    @ mem = $mem + 1
  end

  @ sim_num = $sim_num + 1
end





