  ÉI  º   k820309    ?          19.1        èOb                                                                                                          
       ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 ENSEMBLE_MANAGER_MOD       )       COPIES_IN_WINDOW MEAN_ROW SET_NUM_EXTRA_COPIES GET_ALLOW_TRANSPOSE INIT_ENSEMBLE_MANAGER END_ENSEMBLE_MANAGER GET_ENSEMBLE_TIME ENSEMBLE_TYPE DUPLICATE_ENS GET_VAR_OWNER_INDEX GET_MY_NUM_COPIES GET_MY_COPIES GET_MY_NUM_VARS GET_MY_VARS COMPUTE_COPY_MEAN COMPUTE_COPY_MEAN_SD GET_COPY PUT_COPY ALL_VARS_TO_ALL_COPIES ALL_COPIES_TO_ALL_VARS ALLOCATE_VARS DEALLOCATE_VARS COMPUTE_COPY_MEAN_VAR GET_COPY_OWNER_INDEX SET_ENSEMBLE_TIME BROADCAST_COPY PREPARE_TO_WRITE_TO_VARS PREPARE_TO_WRITE_TO_COPIES PREPARE_TO_READ_FROM_VARS PREPARE_TO_READ_FROM_COPIES PREPARE_TO_UPDATE_VARS PREPARE_TO_UPDATE_COPIES PRINT_ENS_HANDLE SET_CURRENT_TIME MAP_TASK_TO_PE MAP_PE_TO_TASK GET_CURRENT_TIME ALLOCATE_SINGLE_COPY PUT_SINGLE_COPY GET_SINGLE_COPY DEALLOCATE_SINGLE_COPY                      @                              
       R8 I4 I8 MISSING_R8                      @                              
       DO_NML_FILE DO_NML_TERM ERROR_HANDLER E_ERR E_MSG DO_OUTPUT NMLFILEUNIT FIND_NAMELIST_IN_FILE CHECK_NAMELIST_READ TIMESTAMP SET_OUTPUT                      @                              
       TIME_TYPE SET_TIME                      @                              
       TASK_COUNT MY_TASK_ID SEND_TO RECEIVE_FROM TASK_SYNC BROADCAST_SEND BROADCAST_RECV                      @                              
       INDEX_SORT               @                                      u #INDEX_SORT_REAL_INT    #INDEX_SORT_REAL_I8 
   #INDEX_SORT_INT_INT    #INDEX_SORT_USER    #         @     @                                                #X    #INDX    #NUM 	            
                                                     
    p          5 O p            5 O p                                                                                             p          5 O p            5 O p                                    
                                  	           #         @     @                            
                    #X    #INDX    #NUM             
                                                     
 	   p          4 5 O p            4 5 O p                                                                                        
    p          4 5 O p            4 5 O p                                    
                                            #         @     @                                               #X    #INDX    #NUM             
                                                          p          5 O p            5 O p                                                                                             p          5 O p            5 O p                                    
                                             #         @     @                                                #INDX    #NUM    #COMPAREFUNC                                                                       p          5 O p            5 O p                                    
                                             %         @                                                           #A    #B              
                                                      
                                                          À  @                               '                    #SECONDS    #DAYS                  D                                                               D                                                %         @                                                            #ENS_HANDLE              
                                       h             #ENSEMBLE_TYPE    %         @                                                            #ENS_HANDLE              
                                       h             #ENSEMBLE_TYPE    #         @                                                        #ENS_HANDLE !   #N "             
D                                 !     h              #ENSEMBLE_TYPE              
                                  "           %         @                                 #                           #ENS_HANDLE $             
                                  $     h             #ENSEMBLE_TYPE    #         @                                   %                    #ENS_HANDLE &   #NUM_COPIES '   #NUM_VARS (   #DISTRIBUTION_TYPE_IN )   #LAYOUT_TYPE *   #TRANSPOSE_TYPE_IN +             D @                               &     h              #ENSEMBLE_TYPE              
  @                               '                     
                                 (                     
 @                               )                     
 @                               *                     
 @                               +           #         @                                   ,                    #ENS_HANDLE -             
D                                 -     h              #ENSEMBLE_TYPE    #         @                                   .                    #ENS_HANDLE /   #INDX 0   #MTIME 1             
                                  /     h             #ENSEMBLE_TYPE              
                                  0                     D                                 1                    #TIME_TYPE                      @               @                'h                   #NUM_VARS 2   #NUM_COPIES 3   #MY_NUM_COPIES 4   #MY_NUM_VARS 5   #MY_COPIES 6   #MY_VARS 7   #COPIES 8   #VARS 9   #TIME :   #DISTRIBUTION_TYPE ;   #VALID <   #ID_NUM =   #TASK_TO_PE_LIST >   #PE_TO_TASK_LIST ?   #MY_PE @   #LAYOUT_TYPE A   #TRANSPOSE_TYPE B   #NUM_EXTRAS C   #CURRENT_TIME D                 $                             2                                 $                              3                                $                              4                                $                              5                              $                              6                                         &                                                       $                             7            `                             &                                                       $                             8            ¨                 
            &                   &                                                       $                             9                            
            &                   &                                                       $                              :            h             	      #TIME_TYPE              &                                                         $                              ;     °      
                    $                              <     ´                          $                              =     ¸                        $                              >            À                            &                                                      $                              ?                                        &                                                         $                              @     P                          $                              A     T                          $                              B     X                          $                              C     \                          $                              D            `             #TIME_TYPE    #         @                                   E                    #ENS1 F   #ENS2 G   #DUPLICATE_TIME H             
                                  F     h             #ENSEMBLE_TYPE              
D                                 G     h              #ENSEMBLE_TYPE              
                                  H           #         @                                   I                    #ENS_HANDLE J   #VAR_NUMBER K   #OWNER L   #OWNERS_INDEX M             
                                  J     h             #ENSEMBLE_TYPE              
                                 K                     D                                 L                      D                                 M            %         @                                 N                           #ENS_HANDLE O             
                                  O     h             #ENSEMBLE_TYPE    #         @                                   P                    #ENS_HANDLE Q   #COPIES R             
                                  Q     h             #ENSEMBLE_TYPE              D@                               R                                  &                                           %         @                                 S                           #ENS_HANDLE T             
                                  T     h             #ENSEMBLE_TYPE    #         @                                   U                    #ENS_HANDLE V   #VARS W             
                                  V     h             #ENSEMBLE_TYPE              D@                              W                                  &                                           #         @                                   X                    #ENS_HANDLE Y   #START_COPY Z   #END_COPY [   #MEAN_COPY \             
D                                 Y     h              #ENSEMBLE_TYPE              
                                  Z                     
                                  [                     
                                  \           #         @                                   ]                    #ENS_HANDLE ^   #START_COPY _   #END_COPY `   #MEAN_COPY a   #SD_COPY b             
D                                 ^     h              #ENSEMBLE_TYPE              
                                  _                     
                                  `                     
                                  a                     
                                  b           #         @                                   c                    #RECEIVING_PE d   #ENS_HANDLE e   #COPY f   #VARS g   #MTIME h             
  @                               d                     
  @                               e     h             #ENSEMBLE_TYPE              
  @                               f                     D@                              g                   
               &                                                     F @                               h                    #TIME_TYPE    #         @                                   i                    #SENDING_PE j   #ENS_HANDLE k   #COPY l   #VARS m   #MTIME n             
  @                               j                     
D @                               k     h              #ENSEMBLE_TYPE              
  @                               l                     
 @                              m                   
 	             &                                                     
 @                               n                   #TIME_TYPE    #         @                                   o                    #ENS_HANDLE p   #LABEL q             
D @                               p     h              #ENSEMBLE_TYPE              
 @                             q                    1 #         @                                   r                    #ENS_HANDLE s   #LABEL t             
D @                               s     h              #ENSEMBLE_TYPE              
 @                             t                    1 #         @                                   u                    #ENS_HANDLE v             
D                                 v     h              #ENSEMBLE_TYPE    #         @                                   w                    #ENS_HANDLE x             
D                                 x     h              #ENSEMBLE_TYPE    #         @                                   y                    #ENS_HANDLE z   #START_COPY {   #END_COPY |   #MEAN_COPY }   #VAR_COPY ~             
D                                 z     h              #ENSEMBLE_TYPE              
                                  {                     
                                  |                     
                                  }                     
                                  ~           #         @                                                      #ENS_HANDLE    #COPY_NUMBER    #OWNER    #OWNERS_INDEX              
                                       h             #ENSEMBLE_TYPE              
                                                       D                                                       D                                             #         @                                                       #ENS_HANDLE    #INDX    #MTIME              
D                                      h              #ENSEMBLE_TYPE              
                                                       
                                                     #TIME_TYPE    #         @                                                      #ENS_HANDLE    #COPY    #ARRAYDATA              
  @                                    h             #ENSEMBLE_TYPE              
  @                                                    D@                                                 
 
              &                                           #         @                                                       #ENS_HANDLE              
                                      h              #ENSEMBLE_TYPE    #         @                                                       #ENS_HANDLE              
                                      h              #ENSEMBLE_TYPE    #         @                                                       #ENS_HANDLE              
                                       h             #ENSEMBLE_TYPE    #         @                                                       #ENS_HANDLE              
                                       h             #ENSEMBLE_TYPE    #         @                                                       #ENS_HANDLE              
                                      h              #ENSEMBLE_TYPE    #         @                                                       #ENS_HANDLE              
                                      h              #ENSEMBLE_TYPE    #         @                                                      #ENS_HANDLE    #FORCE    #LABEL    #CONTENTS    #LIMIT              
                                       h             #ENSEMBLE_TYPE              
 @                                                    
 @                                                 1           
 @                                                    
 @                                          #         @                                                       #ENS_HANDLE    #T               
D                                      h              #ENSEMBLE_TYPE              
                                                      #TIME_TYPE    %         @                                ¡                           #ENS_HANDLE ¢   #T £             
                                  ¢     h             #ENSEMBLE_TYPE              
                                  £           %         @                                ¤                           #ENS_HANDLE ¥   #P ¦             
                                  ¥     h             #ENSEMBLE_TYPE              
                                  ¦           #         @                                   §                    #ENS_HANDLE ¨   #T ©             
                                  ¨     h             #ENSEMBLE_TYPE              D                                 ©                    #TIME_TYPE    #         @                                   ª                    #ENS_HANDLE «   #X ¬             
                                  «     h             #ENSEMBLE_TYPE            
D                                ¬                   
               &                                           #         @                                   ­                    #ENS_HANDLE ®   #COPY ¯   #X °             
D                                 ®     h              #ENSEMBLE_TYPE              
                                  ¯                     
                                 °                   
              &                                           #         @                                   ±                    #ENS_HANDLE ²   #COPY ³   #X ´             
                                  ²     h             #ENSEMBLE_TYPE              
                                  ³                     
D                                ´                   
               &                                           #         @                                   µ                    #ENS_HANDLE ¶   #X ·             
                                  ¶     h             #ENSEMBLE_TYPE            
D @                              ·                   
               &                                                  c      fn#fn *       b   uapp(ENSEMBLE_MANAGER_MOD      T   J  TYPES_MOD    [  Ç   J  UTILITIES_MOD !   "  S   J  TIME_MANAGER_MOD "   u     J  MPI_UTILITIES_MOD      K   J  SORT_MOD (   S         gen@INDEX_SORT+SORT_MOD -   ñ  b      INDEX_SORT_REAL_INT+SORT_MOD /   S  ¤   a   INDEX_SORT_REAL_INT%X+SORT_MOD 2   ÷  ¤   a   INDEX_SORT_REAL_INT%INDX+SORT_MOD 1     @   a   INDEX_SORT_REAL_INT%NUM+SORT_MOD ,   Û  b      INDEX_SORT_REAL_I8+SORT_MOD .   =	  ¬   a   INDEX_SORT_REAL_I8%X+SORT_MOD 1   é	  ¬   a   INDEX_SORT_REAL_I8%INDX+SORT_MOD 0   
  @   a   INDEX_SORT_REAL_I8%NUM+SORT_MOD ,   Õ
  b      INDEX_SORT_INT_INT+SORT_MOD .   7  ¤   a   INDEX_SORT_INT_INT%X+SORT_MOD 1   Û  ¤   a   INDEX_SORT_INT_INT%INDX+SORT_MOD 0     @   a   INDEX_SORT_INT_INT%NUM+SORT_MOD )   ¿  l      INDEX_SORT_USER+SORT_MOD .   +  ¤   a   INDEX_SORT_USER%INDX+SORT_MOD -   Ï  @   a   INDEX_SORT_USER%NUM+SORT_MOD A     ^      INDEX_SORT_USER%COMPAREFUNC+SORT_MOD=COMPAREFUNC 7   m  @   a   INDEX_SORT_USER%COMPAREFUNC%A+SORT_MOD 7   ­  @   a   INDEX_SORT_USER%COMPAREFUNC%B+SORT_MOD +   í  g      TIME_TYPE+TIME_MANAGER_MOD ;   T  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5     H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS !   ä  `       COPIES_IN_WINDOW ,   D  [   a   COPIES_IN_WINDOW%ENS_HANDLE      `       MEAN_ROW $   ÿ  [   a   MEAN_ROW%ENS_HANDLE %   Z  _       SET_NUM_EXTRA_COPIES 0   ¹  [   a   SET_NUM_EXTRA_COPIES%ENS_HANDLE '     @   a   SET_NUM_EXTRA_COPIES%N $   T  `       GET_ALLOW_TRANSPOSE /   ´  [   a   GET_ALLOW_TRANSPOSE%ENS_HANDLE &     ¸       INIT_ENSEMBLE_MANAGER 1   Ç  [   a   INIT_ENSEMBLE_MANAGER%ENS_HANDLE 1   "  @   a   INIT_ENSEMBLE_MANAGER%NUM_COPIES /   b  @   a   INIT_ENSEMBLE_MANAGER%NUM_VARS ;   ¢  @   a   INIT_ENSEMBLE_MANAGER%DISTRIBUTION_TYPE_IN 2   â  @   a   INIT_ENSEMBLE_MANAGER%LAYOUT_TYPE 8   "  @   a   INIT_ENSEMBLE_MANAGER%TRANSPOSE_TYPE_IN %   b  X       END_ENSEMBLE_MANAGER 0   º  [   a   END_ENSEMBLE_MANAGER%ENS_HANDLE "     m       GET_ENSEMBLE_TIME -     [   a   GET_ENSEMBLE_TIME%ENS_HANDLE '   Ý  @   a   GET_ENSEMBLE_TIME%INDX (     W   a   GET_ENSEMBLE_TIME%MTIME    t  x      ENSEMBLE_TYPE '   ì  H   a   ENSEMBLE_TYPE%NUM_VARS )   4  H   a   ENSEMBLE_TYPE%NUM_COPIES ,   |  H   a   ENSEMBLE_TYPE%MY_NUM_COPIES *   Ä  H   a   ENSEMBLE_TYPE%MY_NUM_VARS (        a   ENSEMBLE_TYPE%MY_COPIES &         a   ENSEMBLE_TYPE%MY_VARS %   4  ¬   a   ENSEMBLE_TYPE%COPIES #   à  ¬   a   ENSEMBLE_TYPE%VARS #     £   a   ENSEMBLE_TYPE%TIME 0   /  H   a   ENSEMBLE_TYPE%DISTRIBUTION_TYPE $   w  H   a   ENSEMBLE_TYPE%VALID %   ¿  H   a   ENSEMBLE_TYPE%ID_NUM .        a   ENSEMBLE_TYPE%TASK_TO_PE_LIST .        a   ENSEMBLE_TYPE%PE_TO_TASK_LIST $   /  H   a   ENSEMBLE_TYPE%MY_PE *   w  H   a   ENSEMBLE_TYPE%LAYOUT_TYPE -   ¿  H   a   ENSEMBLE_TYPE%TRANSPOSE_TYPE )      H   a   ENSEMBLE_TYPE%NUM_EXTRAS +   O   _   a   ENSEMBLE_TYPE%CURRENT_TIME    ®   p       DUPLICATE_ENS #   !  [   a   DUPLICATE_ENS%ENS1 #   y!  [   a   DUPLICATE_ENS%ENS2 -   Ô!  @   a   DUPLICATE_ENS%DUPLICATE_TIME $   "         GET_VAR_OWNER_INDEX /   "  [   a   GET_VAR_OWNER_INDEX%ENS_HANDLE /   ô"  @   a   GET_VAR_OWNER_INDEX%VAR_NUMBER *   4#  @   a   GET_VAR_OWNER_INDEX%OWNER 1   t#  @   a   GET_VAR_OWNER_INDEX%OWNERS_INDEX "   ´#  `       GET_MY_NUM_COPIES -   $  [   a   GET_MY_NUM_COPIES%ENS_HANDLE    o$  d       GET_MY_COPIES )   Ó$  [   a   GET_MY_COPIES%ENS_HANDLE %   .%     a   GET_MY_COPIES%COPIES     º%  `       GET_MY_NUM_VARS +   &  [   a   GET_MY_NUM_VARS%ENS_HANDLE    u&  b       GET_MY_VARS '   ×&  [   a   GET_MY_VARS%ENS_HANDLE !   2'     a   GET_MY_VARS%VARS "   ¾'         COMPUTE_COPY_MEAN -   C(  [   a   COMPUTE_COPY_MEAN%ENS_HANDLE -   (  @   a   COMPUTE_COPY_MEAN%START_COPY +   Þ(  @   a   COMPUTE_COPY_MEAN%END_COPY ,   )  @   a   COMPUTE_COPY_MEAN%MEAN_COPY %   ^)         COMPUTE_COPY_MEAN_SD 0   ð)  [   a   COMPUTE_COPY_MEAN_SD%ENS_HANDLE 0   K*  @   a   COMPUTE_COPY_MEAN_SD%START_COPY .   *  @   a   COMPUTE_COPY_MEAN_SD%END_COPY /   Ë*  @   a   COMPUTE_COPY_MEAN_SD%MEAN_COPY -   +  @   a   COMPUTE_COPY_MEAN_SD%SD_COPY    K+         GET_COPY &   Ô+  @   a   GET_COPY%RECEIVING_PE $   ,  [   a   GET_COPY%ENS_HANDLE    o,  @   a   GET_COPY%COPY    ¯,     a   GET_COPY%VARS    ;-  W   a   GET_COPY%MTIME    -         PUT_COPY $   .  @   a   PUT_COPY%SENDING_PE $   Y.  [   a   PUT_COPY%ENS_HANDLE    ´.  @   a   PUT_COPY%COPY    ô.     a   PUT_COPY%VARS    /  W   a   PUT_COPY%MTIME '   ×/  c       ALL_VARS_TO_ALL_COPIES 2   :0  [   a   ALL_VARS_TO_ALL_COPIES%ENS_HANDLE -   0  L   a   ALL_VARS_TO_ALL_COPIES%LABEL '   á0  c       ALL_COPIES_TO_ALL_VARS 2   D1  [   a   ALL_COPIES_TO_ALL_VARS%ENS_HANDLE -   1  L   a   ALL_COPIES_TO_ALL_VARS%LABEL    ë1  X       ALLOCATE_VARS )   C2  [   a   ALLOCATE_VARS%ENS_HANDLE     2  X       DEALLOCATE_VARS +   ö2  [   a   DEALLOCATE_VARS%ENS_HANDLE &   Q3         COMPUTE_COPY_MEAN_VAR 1   ä3  [   a   COMPUTE_COPY_MEAN_VAR%ENS_HANDLE 1   ?4  @   a   COMPUTE_COPY_MEAN_VAR%START_COPY /   4  @   a   COMPUTE_COPY_MEAN_VAR%END_COPY 0   ¿4  @   a   COMPUTE_COPY_MEAN_VAR%MEAN_COPY /   ÿ4  @   a   COMPUTE_COPY_MEAN_VAR%VAR_COPY %   ?5         GET_COPY_OWNER_INDEX 0   Å5  [   a   GET_COPY_OWNER_INDEX%ENS_HANDLE 1    6  @   a   GET_COPY_OWNER_INDEX%COPY_NUMBER +   `6  @   a   GET_COPY_OWNER_INDEX%OWNER 2    6  @   a   GET_COPY_OWNER_INDEX%OWNERS_INDEX "   à6  m       SET_ENSEMBLE_TIME -   M7  [   a   SET_ENSEMBLE_TIME%ENS_HANDLE '   ¨7  @   a   SET_ENSEMBLE_TIME%INDX (   è7  W   a   SET_ENSEMBLE_TIME%MTIME    ?8  q       BROADCAST_COPY *   °8  [   a   BROADCAST_COPY%ENS_HANDLE $   9  @   a   BROADCAST_COPY%COPY )   K9     a   BROADCAST_COPY%ARRAYDATA )   ×9  X       PREPARE_TO_WRITE_TO_VARS 4   /:  [   a   PREPARE_TO_WRITE_TO_VARS%ENS_HANDLE +   :  X       PREPARE_TO_WRITE_TO_COPIES 6   â:  [   a   PREPARE_TO_WRITE_TO_COPIES%ENS_HANDLE *   =;  X       PREPARE_TO_READ_FROM_VARS 5   ;  [   a   PREPARE_TO_READ_FROM_VARS%ENS_HANDLE ,   ð;  X       PREPARE_TO_READ_FROM_COPIES 7   H<  [   a   PREPARE_TO_READ_FROM_COPIES%ENS_HANDLE '   £<  X       PREPARE_TO_UPDATE_VARS 2   û<  [   a   PREPARE_TO_UPDATE_VARS%ENS_HANDLE )   V=  X       PREPARE_TO_UPDATE_COPIES 4   ®=  [   a   PREPARE_TO_UPDATE_COPIES%ENS_HANDLE !   	>         PRINT_ENS_HANDLE ,   >  [   a   PRINT_ENS_HANDLE%ENS_HANDLE '   ë>  @   a   PRINT_ENS_HANDLE%FORCE '   +?  L   a   PRINT_ENS_HANDLE%LABEL *   w?  @   a   PRINT_ENS_HANDLE%CONTENTS '   ·?  @   a   PRINT_ENS_HANDLE%LIMIT !   ÷?  _       SET_CURRENT_TIME ,   V@  [   a   SET_CURRENT_TIME%ENS_HANDLE #   ±@  W   a   SET_CURRENT_TIME%T    A  g       MAP_TASK_TO_PE *   oA  [   a   MAP_TASK_TO_PE%ENS_HANDLE !   ÊA  @   a   MAP_TASK_TO_PE%T    
B  g       MAP_PE_TO_TASK *   qB  [   a   MAP_PE_TO_TASK%ENS_HANDLE !   ÌB  @   a   MAP_PE_TO_TASK%P !   C  _       GET_CURRENT_TIME ,   kC  [   a   GET_CURRENT_TIME%ENS_HANDLE #   ÆC  W   a   GET_CURRENT_TIME%T %   D  _       ALLOCATE_SINGLE_COPY 0   |D  [   a   ALLOCATE_SINGLE_COPY%ENS_HANDLE '   ×D     a   ALLOCATE_SINGLE_COPY%X     cE  i       PUT_SINGLE_COPY +   ÌE  [   a   PUT_SINGLE_COPY%ENS_HANDLE %   'F  @   a   PUT_SINGLE_COPY%COPY "   gF     a   PUT_SINGLE_COPY%X     óF  i       GET_SINGLE_COPY +   \G  [   a   GET_SINGLE_COPY%ENS_HANDLE %   ·G  @   a   GET_SINGLE_COPY%COPY "   ÷G     a   GET_SINGLE_COPY%X '   H  _       DEALLOCATE_SINGLE_COPY 2   âH  [   a   DEALLOCATE_SINGLE_COPY%ENS_HANDLE )   =I     a   DEALLOCATE_SINGLE_COPY%X 