  ÛG  §   k820309    ?          19.1        Iíb                                                                                                          
       ../../../models/utilities/default_model_mod.f90 DEFAULT_MODEL_MOD              GET_MODEL_SIZE ADV_1STEP GET_STATE_META_DATA MODEL_INTERPOLATE SHORTEST_TIME_BETWEEN_ASSIMILATIONS END_MODEL STATIC_INIT_MODEL INIT_TIME FAIL_INIT_TIME INIT_CONDITIONS FAIL_INIT_CONDITIONS NC_WRITE_MODEL_ATTS NC_WRITE_MODEL_VARS PERT_MODEL_COPIES GET_CLOSE_OBS GET_CLOSE_STATE CONVERT_VERTICAL_OBS CONVERT_VERTICAL_STATE READ_MODEL_TIME WRITE_MODEL_TIME                      @                              
       R8 I8 I4 MISSING_R8                      @                              
       TIME_TYPE SET_TIME                      @                              
       LOCATION_TYPE SET_LOCATION SET_LOCATION_MISSING GET_CLOSE_TYPE GET_CLOSE_OBS GET_CLOSE_STATE CONVERT_VERTICAL_OBS CONVERT_VERTICAL_STATE                      @                              
  	     ERROR_HANDLER E_ERR E_MSG NMLFILEUNIT DO_OUTPUT FIND_NAMELIST_IN_FILE CHECK_NAMELIST_READ DO_NML_FILE DO_NML_TERM                      @                              
       NC_CHECK          @          @                              
       ENSEMBLE_TYPE                      @                              
       READ_MODEL_TIME WRITE_MODEL_TIME               @                                       u #SET_LOCATION_SINGLE    #SET_LOCATION_ARRAY    &         @   @                                                        #LON 	   #LAT 
   #VERT_LOC    #WHICH_VERT    #LOCATION_TYPE              
                                 	     
                
                                 
     
                
                                      
                
                                             &         @   @                                                        #LIST    #LOCATION_TYPE              
                                                   
              &                                                         À  @                                '                    #SECONDS    #DAYS                  D                                                               D                                                                  @                                '                     #LON    #LAT    #VLOC    #WHICH_VERT                  D                                             
                 D                                            
                 D                                            
                 D                                                               @  @               À                '                    #NT    #TYPE_TO_CUTOFF_MAP    #GTT                  $                                                            $                                                                       &                                                       $                                          P       °            #GET_CLOSE_TYPE_BY_TYPE              &                                                         À  @              D                '°                   #NUM    #MAXDIST    #LON_OFFSET    #LOC_BOX    #COUNT     #START !   #BOT_LAT "   #TOP_LAT #   #BOT_LON $   #TOP_LON %   #LON_WIDTH &   #LAT_WIDTH '   #LON_CYCLIC (                 D                                                               D                                            
              D                                                                       &                   &                                                      D                                          p                             &                                                      D                                           ¸                             &                   &                                                      D                              !                                        &                   &                                                         D                             "     x         
                 D                             #              
                 D                             $           	   
                 D                             %           
   
                 D                             &              
                 D                             '               
                 D                              (     ¨                          @  @               @           )     'h                   #NUM_VARS *   #NUM_COPIES +   #MY_NUM_COPIES ,   #MY_NUM_VARS -   #MY_COPIES .   #MY_VARS /   #COPIES 0   #VARS 1   #TIME 2   #DISTRIBUTION_TYPE 3   #VALID 4   #ID_NUM 5   #TASK_TO_PE_LIST 6   #PE_TO_TASK_LIST 7   #MY_PE 8   #LAYOUT_TYPE 9   #TRANSPOSE_TYPE :   #NUM_EXTRAS ;   #CURRENT_TIME <                 $                             *                                 $                              +                                $                              ,                                $                              -                              $                              .                                         &                                                       $                             /            `                             &                                                       $                             0            ¨                 
            &                   &                                                       $                             1                            
            &                   &                                                       $                              2            h             	      #TIME_TYPE              &                                                         $                              3     °      
                    $                              4     ´                          $                              5     ¸                        $                              6            À                            &                                                      $                              7                                        &                                                         $                              8     P                          $                              9     T                          $                              :     X                          $                              ;     \                          $                              <            `             #TIME_TYPE                  @  @               D           =     'h                   #NUM_VARS >   #NUM_COPIES ?   #MY_NUM_COPIES @   #MY_NUM_VARS A   #MY_COPIES B   #MY_VARS C   #COPIES D   #VARS E   #TIME F   #DISTRIBUTION_TYPE G   #VALID H   #ID_NUM I   #TASK_TO_PE_LIST J   #PE_TO_TASK_LIST K   #MY_PE L   #LAYOUT_TYPE M   #TRANSPOSE_TYPE N   #NUM_EXTRAS O   #CURRENT_TIME P                 $                             >                                 $                              ?                                $                              @                                $                              A                              $                              B                                         &                                                       $                             C            `                             &                                                       $                             D            ¨                 
            &                   &                                                       $                             E                            
            &                   &                                                       $                              F            h             	      #TIME_TYPE              &                                                         $                              G     °      
                    $                              H     ´                          $                              I     ¸                        $                              J            À                            &                                                      $                              K                                        &                                                         $                              L     P                          $                              M     T                          $                              N     X                          $                              O     \                          $                              P            `             #TIME_TYPE    #        @                                   Q                 
   #GC R   #BASE_LOC S   #BASE_TYPE T   #LOCS U   #LOC_QTYS V   #LOC_TYPES W   #NUM_CLOSE X   #CLOSE_IND Y   #DIST Z   #ENS_HANDLE [             
                                  R                   #GET_CLOSE_TYPE              
                                 S                     #LOCATION_TYPE              
                                  T                     
                                 U                    "                &                                           #LOCATION_TYPE              
                                  V                    #             &                                                     
                                  W                    $             &                                                                                      X                                                       Y                    %              &                                                                                     Z                   
 &              &                                                     
                                 [     h             #ENSEMBLE_TYPE =   #        @                                   \                 
   #GC ]   #BASE_LOC ^   #BASE_TYPE _   #LOCS `   #LOC_QTYS a   #LOC_INDX b   #NUM_CLOSE c   #CLOSE_IND d   #DIST e   #ENS_HANDLE f             
                                  ]                   #GET_CLOSE_TYPE              
                                 ^                     #LOCATION_TYPE              
                                  _                     
                                 `                    '                &                                           #LOCATION_TYPE              
                                  a                    (             &                                                     
                                 b                    )             &                                                                                      c                                                       d                    *              &                                                                                     e                   
 +              &                                                     
                                 f     h             #ENSEMBLE_TYPE =   #        @                                   g                    #ENS_HANDLE h   #NUM i   #LOCS j   #LOC_QTYS k   #LOC_TYPES l   #WHICH_VERT m   #STATUS n             
                                  h     h             #ENSEMBLE_TYPE =             
                                  i                     
                                 j                                    &                                           #LOCATION_TYPE              
                                  k                                 &                                                     
                                  l                                 &                                                     
                                  m                                                      n                                  &                                           #        @                                   o                    #ENS_HANDLE p   #NUM q   #LOCS r   #LOC_QTYS s   #LOC_INDX t   #WHICH_VERT u   #ISTATUS v             
                                  p     h             #ENSEMBLE_TYPE =             
                                  q                     
                                 r                                    &                                           #LOCATION_TYPE              
                                  s                                 &                                                     
                                 t                                 &                                                     
                                  u                                                      v            &        @                                 w                           #FILENAME x   #TIME_TYPE              
                                x                    1 #        @                                   y                    #NCID z   #DART_TIME {             
                                  z                     
                                  {                   #TIME_TYPE    %         @                                 |                            #         @                                   }                    #X ~   #TIME              
                                ~                   
               &                                                     
                                                     #TIME_TYPE    #         @                                                       #STATE_HANDLE    #INDEX_IN    #LOCATION    #VAR_TYPE              
                                       h             #ENSEMBLE_TYPE )             
                                                      D                                                      #LOCATION_TYPE              F @                                           #         @                                                       #STATE_HANDLE    #ENS_SIZE    #LOCATION    #OBS_QUANTITY    #EXPECTED_OBS    #ISTATUS              
                                       h             #ENSEMBLE_TYPE )             
                                                       
                                                      #LOCATION_TYPE              
                                                      D                                                    
     p          5  p        r        5  p        r                               D                                                          p          5  p        r        5  p        r                      &         @                                                             #TIME_TYPE    #         @                                                        #         @                                                        #         @                                                       #TIME              D                                                     #TIME_TYPE    #         @                                                       #TIME              D                                                     #TIME_TYPE    #         @                                                       #X              D                                                   
               &                                           #         @                                                       #X              D                                                   
               &                                           #         @                                                       #NCID    #DOMAIN_ID              
                                                       
                                             #         @                                                       #NCID    #DOMAIN_ID    #STATE_ENS_HANDLE    #MEMBERINDEX    #TIMEINDEX              
                                                       
                                                       
                                       h             #ENSEMBLE_TYPE )             
                                                      
                                            #         @                                                        #STATE_ENS_HANDLE ¡   #ENS_SIZE ¢   #PERT_AMP £   #INTERF_PROVIDED ¤             
                                 ¡     h              #ENSEMBLE_TYPE )             
                                  ¢                     
                                 £     
                D                                 ¤                   J      fn#fn '   ê   r  b   uapp(DEFAULT_MODEL_MOD    \  T   J  TYPES_MOD !   °  S   J  TIME_MANAGER_MOD      É   J  LOCATION_MOD    Ì  ²   J  UTILITIES_MOD %   ~  I   J  NETCDF_UTILITIES_MOD %   Ç  N   J  ENSEMBLE_MANAGER_MOD !     a   J  DART_TIME_IO_MOD .   v  q       gen@SET_LOCATION+LOCATION_MOD 1   ç        SET_LOCATION_SINGLE+LOCATION_MOD 5   z  @   a   SET_LOCATION_SINGLE%LON+LOCATION_MOD 5   º  @   a   SET_LOCATION_SINGLE%LAT+LOCATION_MOD :   ú  @   a   SET_LOCATION_SINGLE%VERT_LOC+LOCATION_MOD <   :  @   a   SET_LOCATION_SINGLE%WHICH_VERT+LOCATION_MOD 0   z  m      SET_LOCATION_ARRAY+LOCATION_MOD 5   ç     a   SET_LOCATION_ARRAY%LIST+LOCATION_MOD +   s  g      TIME_TYPE+TIME_MANAGER_MOD ;   Ú  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5   "	  H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS +   j	  |       LOCATION_TYPE+LOCATION_MOD 3   æ	  H   %   LOCATION_TYPE%LON+LOCATION_MOD=LON 3   .
  H   %   LOCATION_TYPE%LAT+LOCATION_MOD=LAT 5   v
  H   %   LOCATION_TYPE%VLOC+LOCATION_MOD=VLOC A   ¾
  H   %   LOCATION_TYPE%WHICH_VERT+LOCATION_MOD=WHICH_VERT ,     y      GET_CLOSE_TYPE+LOCATION_MOD /     H   a   GET_CLOSE_TYPE%NT+LOCATION_MOD ?   Ç     a   GET_CLOSE_TYPE%TYPE_TO_CUTOFF_MAP+LOCATION_MOD 0   [  °   a   GET_CLOSE_TYPE%GTT+LOCATION_MOD 4     û      GET_CLOSE_TYPE_BY_TYPE+LOCATION_MOD <     H   %   GET_CLOSE_TYPE_BY_TYPE%NUM+LOCATION_MOD=NUM D   N  H   %   GET_CLOSE_TYPE_BY_TYPE%MAXDIST+LOCATION_MOD=MAXDIST J     ¬   %   GET_CLOSE_TYPE_BY_TYPE%LON_OFFSET+LOCATION_MOD=LON_OFFSET D   B     %   GET_CLOSE_TYPE_BY_TYPE%LOC_BOX+LOCATION_MOD=LOC_BOX @   Ö  ¬   %   GET_CLOSE_TYPE_BY_TYPE%COUNT+LOCATION_MOD=COUNT @     ¬   %   GET_CLOSE_TYPE_BY_TYPE%START+LOCATION_MOD=START D   .  H   %   GET_CLOSE_TYPE_BY_TYPE%BOT_LAT+LOCATION_MOD=BOT_LAT D   v  H   %   GET_CLOSE_TYPE_BY_TYPE%TOP_LAT+LOCATION_MOD=TOP_LAT D   ¾  H   %   GET_CLOSE_TYPE_BY_TYPE%BOT_LON+LOCATION_MOD=BOT_LON D     H   %   GET_CLOSE_TYPE_BY_TYPE%TOP_LON+LOCATION_MOD=TOP_LON H   N  H   %   GET_CLOSE_TYPE_BY_TYPE%LON_WIDTH+LOCATION_MOD=LON_WIDTH H     H   %   GET_CLOSE_TYPE_BY_TYPE%LAT_WIDTH+LOCATION_MOD=LAT_WIDTH J   Þ  H   %   GET_CLOSE_TYPE_BY_TYPE%LON_CYCLIC+LOCATION_MOD=LON_CYCLIC 3   &  x     ENSEMBLE_TYPE+ENSEMBLE_MANAGER_MOD <     H   a   ENSEMBLE_TYPE%NUM_VARS+ENSEMBLE_MANAGER_MOD >   æ  H   a   ENSEMBLE_TYPE%NUM_COPIES+ENSEMBLE_MANAGER_MOD A   .  H   a   ENSEMBLE_TYPE%MY_NUM_COPIES+ENSEMBLE_MANAGER_MOD ?   v  H   a   ENSEMBLE_TYPE%MY_NUM_VARS+ENSEMBLE_MANAGER_MOD =   ¾     a   ENSEMBLE_TYPE%MY_COPIES+ENSEMBLE_MANAGER_MOD ;   R     a   ENSEMBLE_TYPE%MY_VARS+ENSEMBLE_MANAGER_MOD :   æ  ¬   a   ENSEMBLE_TYPE%COPIES+ENSEMBLE_MANAGER_MOD 8     ¬   a   ENSEMBLE_TYPE%VARS+ENSEMBLE_MANAGER_MOD 8   >  £   a   ENSEMBLE_TYPE%TIME+ENSEMBLE_MANAGER_MOD E   á  H   a   ENSEMBLE_TYPE%DISTRIBUTION_TYPE+ENSEMBLE_MANAGER_MOD 9   )  H   a   ENSEMBLE_TYPE%VALID+ENSEMBLE_MANAGER_MOD :   q  H   a   ENSEMBLE_TYPE%ID_NUM+ENSEMBLE_MANAGER_MOD C   ¹     a   ENSEMBLE_TYPE%TASK_TO_PE_LIST+ENSEMBLE_MANAGER_MOD C   M     a   ENSEMBLE_TYPE%PE_TO_TASK_LIST+ENSEMBLE_MANAGER_MOD 9   á  H   a   ENSEMBLE_TYPE%MY_PE+ENSEMBLE_MANAGER_MOD ?   )  H   a   ENSEMBLE_TYPE%LAYOUT_TYPE+ENSEMBLE_MANAGER_MOD B   q  H   a   ENSEMBLE_TYPE%TRANSPOSE_TYPE+ENSEMBLE_MANAGER_MOD >   ¹  H   a   ENSEMBLE_TYPE%NUM_EXTRAS+ENSEMBLE_MANAGER_MOD @     _   a   ENSEMBLE_TYPE%CURRENT_TIME+ENSEMBLE_MANAGER_MOD 3   `  x     ENSEMBLE_TYPE+ENSEMBLE_MANAGER_MOD <   Ø  H   a   ENSEMBLE_TYPE%NUM_VARS+ENSEMBLE_MANAGER_MOD >      H   a   ENSEMBLE_TYPE%NUM_COPIES+ENSEMBLE_MANAGER_MOD A   h  H   a   ENSEMBLE_TYPE%MY_NUM_COPIES+ENSEMBLE_MANAGER_MOD ?   °  H   a   ENSEMBLE_TYPE%MY_NUM_VARS+ENSEMBLE_MANAGER_MOD =   ø     a   ENSEMBLE_TYPE%MY_COPIES+ENSEMBLE_MANAGER_MOD ;        a   ENSEMBLE_TYPE%MY_VARS+ENSEMBLE_MANAGER_MOD :       ¬   a   ENSEMBLE_TYPE%COPIES+ENSEMBLE_MANAGER_MOD 8   Ì   ¬   a   ENSEMBLE_TYPE%VARS+ENSEMBLE_MANAGER_MOD 8   x!  £   a   ENSEMBLE_TYPE%TIME+ENSEMBLE_MANAGER_MOD E   "  H   a   ENSEMBLE_TYPE%DISTRIBUTION_TYPE+ENSEMBLE_MANAGER_MOD 9   c"  H   a   ENSEMBLE_TYPE%VALID+ENSEMBLE_MANAGER_MOD :   «"  H   a   ENSEMBLE_TYPE%ID_NUM+ENSEMBLE_MANAGER_MOD C   ó"     a   ENSEMBLE_TYPE%TASK_TO_PE_LIST+ENSEMBLE_MANAGER_MOD C   #     a   ENSEMBLE_TYPE%PE_TO_TASK_LIST+ENSEMBLE_MANAGER_MOD 9   $  H   a   ENSEMBLE_TYPE%MY_PE+ENSEMBLE_MANAGER_MOD ?   c$  H   a   ENSEMBLE_TYPE%LAYOUT_TYPE+ENSEMBLE_MANAGER_MOD B   «$  H   a   ENSEMBLE_TYPE%TRANSPOSE_TYPE+ENSEMBLE_MANAGER_MOD >   ó$  H   a   ENSEMBLE_TYPE%NUM_EXTRAS+ENSEMBLE_MANAGER_MOD @   ;%  _   a   ENSEMBLE_TYPE%CURRENT_TIME+ENSEMBLE_MANAGER_MOD +   %  Ì       GET_CLOSE_OBS+LOCATION_MOD .   f&  \   a   GET_CLOSE_OBS%GC+LOCATION_MOD 4   Â&  [   a   GET_CLOSE_OBS%BASE_LOC+LOCATION_MOD 5   '  @   a   GET_CLOSE_OBS%BASE_TYPE+LOCATION_MOD 0   ]'     a   GET_CLOSE_OBS%LOCS+LOCATION_MOD 4   ü'     a   GET_CLOSE_OBS%LOC_QTYS+LOCATION_MOD 5   (     a   GET_CLOSE_OBS%LOC_TYPES+LOCATION_MOD 5   )  @   a   GET_CLOSE_OBS%NUM_CLOSE+LOCATION_MOD 5   T)     a   GET_CLOSE_OBS%CLOSE_IND+LOCATION_MOD 0   à)     a   GET_CLOSE_OBS%DIST+LOCATION_MOD 6   l*  [   a   GET_CLOSE_OBS%ENS_HANDLE+LOCATION_MOD -   Ç*  Ë       GET_CLOSE_STATE+LOCATION_MOD 0   +  \   a   GET_CLOSE_STATE%GC+LOCATION_MOD 6   î+  [   a   GET_CLOSE_STATE%BASE_LOC+LOCATION_MOD 7   I,  @   a   GET_CLOSE_STATE%BASE_TYPE+LOCATION_MOD 2   ,     a   GET_CLOSE_STATE%LOCS+LOCATION_MOD 6   (-     a   GET_CLOSE_STATE%LOC_QTYS+LOCATION_MOD 6   ´-     a   GET_CLOSE_STATE%LOC_INDX+LOCATION_MOD 7   @.  @   a   GET_CLOSE_STATE%NUM_CLOSE+LOCATION_MOD 7   .     a   GET_CLOSE_STATE%CLOSE_IND+LOCATION_MOD 2   /     a   GET_CLOSE_STATE%DIST+LOCATION_MOD 8   /  [   a   GET_CLOSE_STATE%ENS_HANDLE+LOCATION_MOD 2   ó/  ¤       CONVERT_VERTICAL_OBS+LOCATION_MOD =   0  [   a   CONVERT_VERTICAL_OBS%ENS_HANDLE+LOCATION_MOD 6   ò0  @   a   CONVERT_VERTICAL_OBS%NUM+LOCATION_MOD 7   21     a   CONVERT_VERTICAL_OBS%LOCS+LOCATION_MOD ;   Ñ1     a   CONVERT_VERTICAL_OBS%LOC_QTYS+LOCATION_MOD <   ]2     a   CONVERT_VERTICAL_OBS%LOC_TYPES+LOCATION_MOD =   é2  @   a   CONVERT_VERTICAL_OBS%WHICH_VERT+LOCATION_MOD 9   )3     a   CONVERT_VERTICAL_OBS%STATUS+LOCATION_MOD 4   µ3  ¤       CONVERT_VERTICAL_STATE+LOCATION_MOD ?   Y4  [   a   CONVERT_VERTICAL_STATE%ENS_HANDLE+LOCATION_MOD 8   ´4  @   a   CONVERT_VERTICAL_STATE%NUM+LOCATION_MOD 9   ô4     a   CONVERT_VERTICAL_STATE%LOCS+LOCATION_MOD =   5     a   CONVERT_VERTICAL_STATE%LOC_QTYS+LOCATION_MOD =   6     a   CONVERT_VERTICAL_STATE%LOC_INDX+LOCATION_MOD ?   «6  @   a   CONVERT_VERTICAL_STATE%WHICH_VERT+LOCATION_MOD <   ë6  @   a   CONVERT_VERTICAL_STATE%ISTATUS+LOCATION_MOD 1   +7  m       READ_MODEL_TIME+DART_TIME_IO_MOD :   7  L   a   READ_MODEL_TIME%FILENAME+DART_TIME_IO_MOD 2   ä7  a       WRITE_MODEL_TIME+DART_TIME_IO_MOD 7   E8  @   a   WRITE_MODEL_TIME%NCID+DART_TIME_IO_MOD <   8  W   a   WRITE_MODEL_TIME%DART_TIME+DART_TIME_IO_MOD    Ü8  P       GET_MODEL_SIZE    ,9  Y       ADV_1STEP    9     a   ADV_1STEP%X    :  W   a   ADV_1STEP%TIME $   h:         GET_STATE_META_DATA 1   ì:  [   a   GET_STATE_META_DATA%STATE_HANDLE -   G;  @   a   GET_STATE_META_DATA%INDEX_IN -   ;  [   a   GET_STATE_META_DATA%LOCATION -   â;  @   a   GET_STATE_META_DATA%VAR_TYPE "   "<  §       MODEL_INTERPOLATE /   É<  [   a   MODEL_INTERPOLATE%STATE_HANDLE +   $=  @   a   MODEL_INTERPOLATE%ENS_SIZE +   d=  [   a   MODEL_INTERPOLATE%LOCATION /   ¿=  @   a   MODEL_INTERPOLATE%OBS_QUANTITY /   ÿ=  ´   a   MODEL_INTERPOLATE%EXPECTED_OBS *   ³>  ´   a   MODEL_INTERPOLATE%ISTATUS 4   g?  _       SHORTEST_TIME_BETWEEN_ASSIMILATIONS    Æ?  H       END_MODEL "   @  H       STATIC_INIT_MODEL    V@  R       INIT_TIME    ¨@  W   a   INIT_TIME%TIME    ÿ@  R       FAIL_INIT_TIME $   QA  W   a   FAIL_INIT_TIME%TIME     ¨A  O       INIT_CONDITIONS "   ÷A     a   INIT_CONDITIONS%X %   B  O       FAIL_INIT_CONDITIONS '   ÒB     a   FAIL_INIT_CONDITIONS%X $   ^C  a       NC_WRITE_MODEL_ATTS )   ¿C  @   a   NC_WRITE_MODEL_ATTS%NCID .   ÿC  @   a   NC_WRITE_MODEL_ATTS%DOMAIN_ID $   ?D         NC_WRITE_MODEL_VARS )   ÖD  @   a   NC_WRITE_MODEL_VARS%NCID .   E  @   a   NC_WRITE_MODEL_VARS%DOMAIN_ID 5   VE  [   a   NC_WRITE_MODEL_VARS%STATE_ENS_HANDLE 0   ±E  @   a   NC_WRITE_MODEL_VARS%MEMBERINDEX .   ñE  @   a   NC_WRITE_MODEL_VARS%TIMEINDEX "   1F         PERT_MODEL_COPIES 3   ÀF  [   a   PERT_MODEL_COPIES%STATE_ENS_HANDLE +   G  @   a   PERT_MODEL_COPIES%ENS_SIZE +   [G  @   a   PERT_MODEL_COPIES%PERT_AMP 2   G  @   a   PERT_MODEL_COPIES%INTERF_PROVIDED 