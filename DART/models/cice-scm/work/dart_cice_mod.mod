  �  !   k820309    ?          19.1        J�b                                                                                                          
       ../../../models/cice-scm/dart_cice_mod.f90 DART_CICE_MOD              SET_MODEL_TIME_STEP GET_HORIZ_GRID_DIMS GET_NCAT_DIM READ_HORIZ_GRID                                                     
                @       �   @                              
                            @                              
       R8 RAD2DEG PI SECPERDAY DIGITS12                      @                              
  	    TIME_TYPE GET_DATE SET_DATE GET_TIME SET_TIME SET_CALENDAR_TYPE GET_CALENDAR_STRING PRINT_DATE PRINT_TIME i@ i@                      @                              
       GET_UNIT OPEN_FILE CLOSE_FILE FILE_EXIST REGISTER_MODULE ERROR_HANDLER FIND_NAMELIST_IN_FILE CHECK_NAMELIST_READ E_ERR E_MSG FIND_TEXTFILE_DIMS                      @                              
       NC_CHECK                                                            #TIME_MINUS    &         @   @                                                       #TIME1    #TIME2 
   #TIME_TYPE 	             
                                                     #TIME_TYPE 	             
                                  
                   #TIME_TYPE 	                                                              #TIME_EQ    %         @   @                                                       #TIME1    #TIME2              
                                                     #TIME_TYPE 	             
                                                     #TIME_TYPE 	                 @�                                      u #SET_CALENDAR_TYPE_INTEGER    #SET_CALENDAR_TYPE_STRING    #         @     @                                                #MYTYPE              
                                             #         @     @                                               #CALSTRING              
                                                    1                �  @                           	     '                    #SECONDS    #DAYS                 � D                                                              � D                                                &         @                                                             #TIME_TYPE 	   #         @                                                       #NX              D @                                           #         @                                                       #NCAT              D @                                           #         @                                                       #NX    #TLAT    #TLON              
                                                      D @                                                  
     p          5 � p        r        5 � p        r                               D @                                                  
     p          5 � p        r        5 � p        r                         �   A      fn#fn #   �   U   b   uapp(DART_CICE_MOD    6  @   J  TYPESIZES    v  @   J  NETCDF    �  a   J  TYPES_MOD !     �   J  TIME_MANAGER_MOD    �  �   J  UTILITIES_MOD %   �  I   J  NETCDF_UTILITIES_MOD &   �  P      i@+TIME_MANAGER_MOD ,   4  u      TIME_MINUS+TIME_MANAGER_MOD 2   �  W   a   TIME_MINUS%TIME1+TIME_MANAGER_MOD 2      W   a   TIME_MINUS%TIME2+TIME_MANAGER_MOD &   W  M      i@+TIME_MANAGER_MOD )   �  f      TIME_EQ+TIME_MANAGER_MOD /   
  W   a   TIME_EQ%TIME1+TIME_MANAGER_MOD /   a  W   a   TIME_EQ%TIME2+TIME_MANAGER_MOD 7   �  }       gen@SET_CALENDAR_TYPE+TIME_MANAGER_MOD ;   5  T      SET_CALENDAR_TYPE_INTEGER+TIME_MANAGER_MOD B   �  @   a   SET_CALENDAR_TYPE_INTEGER%MYTYPE+TIME_MANAGER_MOD :   �  W      SET_CALENDAR_TYPE_STRING+TIME_MANAGER_MOD D      L   a   SET_CALENDAR_TYPE_STRING%CALSTRING+TIME_MANAGER_MOD +   l  g       TIME_TYPE+TIME_MANAGER_MOD ;   �  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5   	  H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS $   c	  _       SET_MODEL_TIME_STEP $   �	  P       GET_HORIZ_GRID_DIMS '   
  @   a   GET_HORIZ_GRID_DIMS%NX    R
  R       GET_NCAT_DIM "   �
  @   a   GET_NCAT_DIM%NCAT     �
  d       READ_HORIZ_GRID #   H  @   a   READ_HORIZ_GRID%NX %   �  �   a   READ_HORIZ_GRID%TLAT %   <  �   a   READ_HORIZ_GRID%TLON 