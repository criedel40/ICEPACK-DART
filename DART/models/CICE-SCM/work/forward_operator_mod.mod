  �:  ~   k820309    ?          19.1        �3�b                                                                                                          
       ../../../assimilation_code/modules/observations/forward_operator_mod.f90 FORWARD_OPERATOR_MOD              GET_OBS_ENS_DISTRIB_STATE GET_EXPECTED_OBS_DISTRIB_STATE                      @                              
       R8 I8 MISSING_R8                      @                              
       TIME_TYPE                      @                              
       ERROR_HANDLER E_ERR                      @                              
       MY_TASK_ID                      @                              
       INIT_OBS DESTROY_OBS OBS_SEQUENCE_TYPE OBS_TYPE GET_OBS_VALUES GET_QC GET_OBS_DEF GET_OBS_FROM_KEY          @       �   @                              
       OBS_DEF_TYPE GET_OBS_DEF_ERROR_VARIANCE GET_EXPECTED_OBS_FROM_DEF_DISTRIB_STATE GET_OBS_DEF_TYPE_OF_OBS                      @                              
       ASSIMILATE_THIS_TYPE_OF_OBS EVALUATE_THIS_TYPE_OF_OBS          @       �   @                              
       ENSEMBLE_TYPE COMPUTE_COPY_MEAN_VAR PREPARE_TO_READ_FROM_VARS PREPARE_TO_WRITE_TO_VARS GET_MY_NUM_COPIES COPIES_IN_WINDOW GET_ALLOW_TRANSPOSE ALL_VARS_TO_ALL_COPIES ALL_COPIES_TO_ALL_VARS ALLOCATE_SINGLE_COPY GET_SINGLE_COPY PUT_SINGLE_COPY DEALLOCATE_SINGLE_COPY                      @                         	     
       CREATE_STATE_WINDOW FREE_STATE_WINDOW GET_STATE                      @                         
     
       CHECK_OUTLIER_THRESHOLD GET_DART_QC INPUT_QC_OK GOOD_DART_QC DARTQC_BAD_INCOMING_QC DARTQC_ASSIM_GOOD_FOP DARTQC_EVAL_GOOD_FOP DARTQC_NOT_IN_NAMELIST               �  @                                '                    #SECONDS    #DAYS                 � D                                                              � D                                                              �  @                               '�              	      #NUM_COPIES    #NUM_QC    #NUM_OBS    #MAX_NUM_OBS    #COPY_META_DATA    #QC_META_DATA    #FIRST_TIME    #LAST_TIME    #OBS                 � D                                                              � D                                                             � D                                                             � D                                                .           �D                                                @                     &                                                                           �              y                                               .           �D                                         X       @                     &                                                                           �              y                                                            � D                                   �                          � D                                   �                        �D                                          �       (      	     #OBS_TYPE              &                                                                   �              y#OBS_TYPE                                                                  �  @              �                '(                   #KEY    #DEF    #VALUES (   #QC )   #PREV_TIME *   #NEXT_TIME +   #COV_GROUP ,                � D                                                              � D                                   �                     #OBS_DEF_TYPE                   �  @              E                '�              
      #LOCATION    #KIND    #TIME     #ERROR_VARIANCE !   #KEY "   #WRITE_EXTERNAL_FO #   #HAS_EXTERNAL_FO $   #EXTERNAL_FO %   #EXTERNAL_FO_KEY &   #ENS_SIZE '                � D                                                        #LOCATION_TYPE                  �  @                              '                    #X                 � D                                            
                � D                                                            � D                                                        #TIME_TYPE                 � D                            !               
                � D                             "                               � D                             #     $                                    �                                                                         � D                             $     (                                    �                                                                        � D                            %            0                 
            &                                                        � D                             &     x       	                   � D                             '     |       
                �D                             (            �                
            &                                                                   �              y
                                                         �D                             )            �                
            &                                                                   �              y
                                                            � D                              *                              � D                              +                              � D                              ,                                 �  @               A           -     '�              
      #LOCATION .   #KIND /   #TIME 0   #ERROR_VARIANCE 1   #KEY 2   #WRITE_EXTERNAL_FO 3   #HAS_EXTERNAL_FO 4   #EXTERNAL_FO 5   #EXTERNAL_FO_KEY 6   #ENS_SIZE 7                � D                              .                           #LOCATION_TYPE                 � D                              /                               � D                              0                          #TIME_TYPE                 � D                             1               
                � D                              2                               � D                              3     $                                    �                                                                         � D                              4     (                                    �                                                                        � D                             5            0                 
            &                                                        � D                              6     x       	                   � D                              7     |       
                        @               @           8     'h                   #NUM_VARS 9   #NUM_COPIES :   #MY_NUM_COPIES ;   #MY_NUM_VARS <   #MY_COPIES =   #MY_VARS >   #COPIES ?   #VARS @   #TIME A   #DISTRIBUTION_TYPE B   #VALID C   #ID_NUM D   #TASK_TO_PE_LIST E   #PE_TO_TASK_LIST F   #MY_PE G   #LAYOUT_TYPE H   #TRANSPOSE_TYPE I   #NUM_EXTRAS J   #CURRENT_TIME K                � $                             9                                � $                              :                               � $                              ;                               � $                              <                             � $                              =                                         &                                                      � $                             >            `                             &                                                      � $                             ?            �                 
            &                   &                                                      � $                             @                            
            &                   &                                                      � $                              A            h             	      #TIME_TYPE              &                                                        � $                              B     �      
                   � $                              C     �                         � $                              D     �                       � $                              E            �                            &                                                     � $                              F                                        &                                                        � $                              G     P                         � $                              H     T                         � $                              I     X                         � $                              J     \                         � $                              K            `             #TIME_TYPE                  @  @               D           L     'h                   #NUM_VARS M   #NUM_COPIES N   #MY_NUM_COPIES O   #MY_NUM_VARS P   #MY_COPIES Q   #MY_VARS R   #COPIES S   #VARS T   #TIME U   #DISTRIBUTION_TYPE V   #VALID W   #ID_NUM X   #TASK_TO_PE_LIST Y   #PE_TO_TASK_LIST Z   #MY_PE [   #LAYOUT_TYPE \   #TRANSPOSE_TYPE ]   #NUM_EXTRAS ^   #CURRENT_TIME _                � $                             M                                � $                              N                               � $                              O                               � $                              P                             � $                              Q                                         &                                                      � $                             R            `                             &                                                      � $                             S            �                 
            &                   &                                                      � $                             T                            
            &                   &                                                      � $                              U            h             	      #TIME_TYPE              &                                                        � $                              V     �      
                   � $                              W     �                         � $                              X     �                       � $                              Y            �                            &                                                     � $                              Z                                        &                                                        � $                              [     P                         � $                              \     T                         � $                              ]     X                         � $                              ^     \                         � $                              _            `             #TIME_TYPE    #         @                                   `                    #ENS_HANDLE a   #OBS_FWD_OP_ENS_HANDLE b   #QC_ENS_HANDLE c   #SEQ d   #KEYS e   #OBS_VAL_INDEX f   #INPUT_QC_INDEX g   #OBS_ERR_VAR_COPY h   #OBS_VAL_COPY i   #OBS_KEY_COPY j   #OBS_GLOBAL_QC_COPY k   #OBS_EXTRA_QC_COPY l   #OBS_MEAN_COPY m   #OBS_VAR_COPY n   #ISPRIOR o   #PRIOR_QC_COPY p             
D @                               a     h              #ENSEMBLE_TYPE 8             
D @                               b     h              #ENSEMBLE_TYPE 8             
D @                               c     h              #ENSEMBLE_TYPE 8             
  @                               d     �              #OBS_SEQUENCE_TYPE              
  @                               e                                 &                                                     
  @                               f                     
  @                               g                     
                                  h                     
                                  i                     
                                  j                     
  @                               k                     
                                  l                     
  @                               m                     
  @                               n                     
  @                               o                     
D @                              p                   
               &                                           #         @                                  q                    #SEQ r   #KEYS s   #STATE_TIME t   #ISPRIOR u   #ISTATUS v   #ASSIMILATE_THIS_OB x   #EVALUATE_THIS_OB y   #STATE_ENS_HANDLE z   #NUM_ENS w   #COPY_INDICES {   #EXPECTED_OBS |             
  @                               r     �              #OBS_SEQUENCE_TYPE              
 @                               s                    
             &                                                     
  @                               t                   #TIME_TYPE              
  @                               u                    D @                               v                         p          5 � p 	       r w       5 � p 	       r w                               D @                               x                      D @                               y                      
  @                               z     h             #ENSEMBLE_TYPE 8             
  @                               w                    
  @                               {                        p          5 � p 	       r w       5 � p 	       r w                              
D @                              |                    
     p          5 � p 	       r w       5 � p 	       r w                        �   f      fn#fn *     I   b   uapp(FORWARD_OPERATOR_MOD    O  Q   J  TYPES_MOD !   �  J   J  TIME_MANAGER_MOD    �  T   J  UTILITIES_MOD "   >  K   J  MPI_UTILITIES_MOD !   �  �   J  OBS_SEQUENCE_MOD    ,  �   J  OBS_DEF_MOD    �  v   J  OBS_KIND_MOD %   J  H  J  ENSEMBLE_MANAGER_MOD &   �  p   J  DISTRIBUTED_STATE_MOD $     �   J  QUALITY_CONTROL_MOD +   �  g      TIME_TYPE+TIME_MANAGER_MOD ;   ?  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5   �  H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS 3   �  �      OBS_SEQUENCE_TYPE+OBS_SEQUENCE_MOD I   �  H   %   OBS_SEQUENCE_TYPE%NUM_COPIES+OBS_SEQUENCE_MOD=NUM_COPIES A   �  H   %   OBS_SEQUENCE_TYPE%NUM_QC+OBS_SEQUENCE_MOD=NUM_QC C   7	  H   %   OBS_SEQUENCE_TYPE%NUM_OBS+OBS_SEQUENCE_MOD=NUM_OBS K   	  H   %   OBS_SEQUENCE_TYPE%MAX_NUM_OBS+OBS_SEQUENCE_MOD=MAX_NUM_OBS Q   �	  �   %   OBS_SEQUENCE_TYPE%COPY_META_DATA+OBS_SEQUENCE_MOD=COPY_META_DATA M   �
  �   %   OBS_SEQUENCE_TYPE%QC_META_DATA+OBS_SEQUENCE_MOD=QC_META_DATA I   �  H   %   OBS_SEQUENCE_TYPE%FIRST_TIME+OBS_SEQUENCE_MOD=FIRST_TIME G     H   %   OBS_SEQUENCE_TYPE%LAST_TIME+OBS_SEQUENCE_MOD=LAST_TIME ;   O    %   OBS_SEQUENCE_TYPE%OBS+OBS_SEQUENCE_MOD=OBS *   _  �      OBS_TYPE+OBS_SEQUENCE_MOD 2     H   %   OBS_TYPE%KEY+OBS_SEQUENCE_MOD=KEY 2   J  b   %   OBS_TYPE%DEF+OBS_SEQUENCE_MOD=DEF )   �  �      OBS_DEF_TYPE+OBS_DEF_MOD ;   �  c   %   OBS_DEF_TYPE%LOCATION+OBS_DEF_MOD=LOCATION +   �  W      LOCATION_TYPE+LOCATION_MOD /   U  H   %   LOCATION_TYPE%X+LOCATION_MOD=X 3   �  H   %   OBS_DEF_TYPE%KIND+OBS_DEF_MOD=KIND 3   �  _   %   OBS_DEF_TYPE%TIME+OBS_DEF_MOD=TIME G   D  H   %   OBS_DEF_TYPE%ERROR_VARIANCE+OBS_DEF_MOD=ERROR_VARIANCE 1   �  H   %   OBS_DEF_TYPE%KEY+OBS_DEF_MOD=KEY M   �  �   %   OBS_DEF_TYPE%WRITE_EXTERNAL_FO+OBS_DEF_MOD=WRITE_EXTERNAL_FO I   x  �   %   OBS_DEF_TYPE%HAS_EXTERNAL_FO+OBS_DEF_MOD=HAS_EXTERNAL_FO A     �   %   OBS_DEF_TYPE%EXTERNAL_FO+OBS_DEF_MOD=EXTERNAL_FO I   �  H   %   OBS_DEF_TYPE%EXTERNAL_FO_KEY+OBS_DEF_MOD=EXTERNAL_FO_KEY ;   �  H   %   OBS_DEF_TYPE%ENS_SIZE+OBS_DEF_MOD=ENS_SIZE 8   @  �   %   OBS_TYPE%VALUES+OBS_SEQUENCE_MOD=VALUES 0   4  �   %   OBS_TYPE%QC+OBS_SEQUENCE_MOD=QC >   (  H   %   OBS_TYPE%PREV_TIME+OBS_SEQUENCE_MOD=PREV_TIME >   p  H   %   OBS_TYPE%NEXT_TIME+OBS_SEQUENCE_MOD=NEXT_TIME >   �  H   %   OBS_TYPE%COV_GROUP+OBS_SEQUENCE_MOD=COV_GROUP )      �       OBS_DEF_TYPE+OBS_DEF_MOD ;   �  c   %   OBS_DEF_TYPE%LOCATION+OBS_DEF_MOD=LOCATION 3   R  H   %   OBS_DEF_TYPE%KIND+OBS_DEF_MOD=KIND 3   �  _   %   OBS_DEF_TYPE%TIME+OBS_DEF_MOD=TIME G   �  H   %   OBS_DEF_TYPE%ERROR_VARIANCE+OBS_DEF_MOD=ERROR_VARIANCE 1   A  H   %   OBS_DEF_TYPE%KEY+OBS_DEF_MOD=KEY M   �  �   %   OBS_DEF_TYPE%WRITE_EXTERNAL_FO+OBS_DEF_MOD=WRITE_EXTERNAL_FO I   -  �   %   OBS_DEF_TYPE%HAS_EXTERNAL_FO+OBS_DEF_MOD=HAS_EXTERNAL_FO A   �  �   %   OBS_DEF_TYPE%EXTERNAL_FO+OBS_DEF_MOD=EXTERNAL_FO I   e  H   %   OBS_DEF_TYPE%EXTERNAL_FO_KEY+OBS_DEF_MOD=EXTERNAL_FO_KEY ;   �  H   %   OBS_DEF_TYPE%ENS_SIZE+OBS_DEF_MOD=ENS_SIZE 3   �  x      ENSEMBLE_TYPE+ENSEMBLE_MANAGER_MOD <   m  H   a   ENSEMBLE_TYPE%NUM_VARS+ENSEMBLE_MANAGER_MOD >   �  H   a   ENSEMBLE_TYPE%NUM_COPIES+ENSEMBLE_MANAGER_MOD A   �  H   a   ENSEMBLE_TYPE%MY_NUM_COPIES+ENSEMBLE_MANAGER_MOD ?   E  H   a   ENSEMBLE_TYPE%MY_NUM_VARS+ENSEMBLE_MANAGER_MOD =   �  �   a   ENSEMBLE_TYPE%MY_COPIES+ENSEMBLE_MANAGER_MOD ;   !  �   a   ENSEMBLE_TYPE%MY_VARS+ENSEMBLE_MANAGER_MOD :   �  �   a   ENSEMBLE_TYPE%COPIES+ENSEMBLE_MANAGER_MOD 8   a   �   a   ENSEMBLE_TYPE%VARS+ENSEMBLE_MANAGER_MOD 8   !  �   a   ENSEMBLE_TYPE%TIME+ENSEMBLE_MANAGER_MOD E   �!  H   a   ENSEMBLE_TYPE%DISTRIBUTION_TYPE+ENSEMBLE_MANAGER_MOD 9   �!  H   a   ENSEMBLE_TYPE%VALID+ENSEMBLE_MANAGER_MOD :   @"  H   a   ENSEMBLE_TYPE%ID_NUM+ENSEMBLE_MANAGER_MOD C   �"  �   a   ENSEMBLE_TYPE%TASK_TO_PE_LIST+ENSEMBLE_MANAGER_MOD C   #  �   a   ENSEMBLE_TYPE%PE_TO_TASK_LIST+ENSEMBLE_MANAGER_MOD 9   �#  H   a   ENSEMBLE_TYPE%MY_PE+ENSEMBLE_MANAGER_MOD ?   �#  H   a   ENSEMBLE_TYPE%LAYOUT_TYPE+ENSEMBLE_MANAGER_MOD B   @$  H   a   ENSEMBLE_TYPE%TRANSPOSE_TYPE+ENSEMBLE_MANAGER_MOD >   �$  H   a   ENSEMBLE_TYPE%NUM_EXTRAS+ENSEMBLE_MANAGER_MOD @   �$  _   a   ENSEMBLE_TYPE%CURRENT_TIME+ENSEMBLE_MANAGER_MOD 3   /%  x     ENSEMBLE_TYPE+ENSEMBLE_MANAGER_MOD <   �&  H   a   ENSEMBLE_TYPE%NUM_VARS+ENSEMBLE_MANAGER_MOD >   �&  H   a   ENSEMBLE_TYPE%NUM_COPIES+ENSEMBLE_MANAGER_MOD A   7'  H   a   ENSEMBLE_TYPE%MY_NUM_COPIES+ENSEMBLE_MANAGER_MOD ?   '  H   a   ENSEMBLE_TYPE%MY_NUM_VARS+ENSEMBLE_MANAGER_MOD =   �'  �   a   ENSEMBLE_TYPE%MY_COPIES+ENSEMBLE_MANAGER_MOD ;   [(  �   a   ENSEMBLE_TYPE%MY_VARS+ENSEMBLE_MANAGER_MOD :   �(  �   a   ENSEMBLE_TYPE%COPIES+ENSEMBLE_MANAGER_MOD 8   �)  �   a   ENSEMBLE_TYPE%VARS+ENSEMBLE_MANAGER_MOD 8   G*  �   a   ENSEMBLE_TYPE%TIME+ENSEMBLE_MANAGER_MOD E   �*  H   a   ENSEMBLE_TYPE%DISTRIBUTION_TYPE+ENSEMBLE_MANAGER_MOD 9   2+  H   a   ENSEMBLE_TYPE%VALID+ENSEMBLE_MANAGER_MOD :   z+  H   a   ENSEMBLE_TYPE%ID_NUM+ENSEMBLE_MANAGER_MOD C   �+  �   a   ENSEMBLE_TYPE%TASK_TO_PE_LIST+ENSEMBLE_MANAGER_MOD C   V,  �   a   ENSEMBLE_TYPE%PE_TO_TASK_LIST+ENSEMBLE_MANAGER_MOD 9   �,  H   a   ENSEMBLE_TYPE%MY_PE+ENSEMBLE_MANAGER_MOD ?   2-  H   a   ENSEMBLE_TYPE%LAYOUT_TYPE+ENSEMBLE_MANAGER_MOD B   z-  H   a   ENSEMBLE_TYPE%TRANSPOSE_TYPE+ENSEMBLE_MANAGER_MOD >   �-  H   a   ENSEMBLE_TYPE%NUM_EXTRAS+ENSEMBLE_MANAGER_MOD @   
.  _   a   ENSEMBLE_TYPE%CURRENT_TIME+ENSEMBLE_MANAGER_MOD *   i.  n      GET_OBS_ENS_DISTRIB_STATE 5   �/  [   a   GET_OBS_ENS_DISTRIB_STATE%ENS_HANDLE @   20  [   a   GET_OBS_ENS_DISTRIB_STATE%OBS_FWD_OP_ENS_HANDLE 8   �0  [   a   GET_OBS_ENS_DISTRIB_STATE%QC_ENS_HANDLE .   �0  _   a   GET_OBS_ENS_DISTRIB_STATE%SEQ /   G1  �   a   GET_OBS_ENS_DISTRIB_STATE%KEYS 8   �1  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_VAL_INDEX 9   2  @   a   GET_OBS_ENS_DISTRIB_STATE%INPUT_QC_INDEX ;   S2  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_ERR_VAR_COPY 7   �2  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_VAL_COPY 7   �2  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_KEY_COPY =   3  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_GLOBAL_QC_COPY <   S3  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_EXTRA_QC_COPY 8   �3  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_MEAN_COPY 7   �3  @   a   GET_OBS_ENS_DISTRIB_STATE%OBS_VAR_COPY 2   4  @   a   GET_OBS_ENS_DISTRIB_STATE%ISPRIOR 8   S4  �   a   GET_OBS_ENS_DISTRIB_STATE%PRIOR_QC_COPY /   �4  �       GET_EXPECTED_OBS_DISTRIB_STATE 3   �5  _   a   GET_EXPECTED_OBS_DISTRIB_STATE%SEQ 4   86  �   a   GET_EXPECTED_OBS_DISTRIB_STATE%KEYS :   �6  W   a   GET_EXPECTED_OBS_DISTRIB_STATE%STATE_TIME 7   7  @   a   GET_EXPECTED_OBS_DISTRIB_STATE%ISPRIOR 7   [7  �   a   GET_EXPECTED_OBS_DISTRIB_STATE%ISTATUS B   8  @   a   GET_EXPECTED_OBS_DISTRIB_STATE%ASSIMILATE_THIS_OB @   O8  @   a   GET_EXPECTED_OBS_DISTRIB_STATE%EVALUATE_THIS_OB @   �8  [   a   GET_EXPECTED_OBS_DISTRIB_STATE%STATE_ENS_HANDLE 7   �8  @   a   GET_EXPECTED_OBS_DISTRIB_STATE%NUM_ENS <   *9  �   a   GET_EXPECTED_OBS_DISTRIB_STATE%COPY_INDICES <   �9  �   a   GET_EXPECTED_OBS_DISTRIB_STATE%EXPECTED_OBS 