  �=  �   k820309    ?          19.1        F�b                                                                                                          
       ../../../assimilation_code/modules/observations/obs_kind_mod.f90 OBS_KIND_MOD       d       GET_NAME_FOR_TYPE_OF_OBS ASSIMILATE_THIS_TYPE_OF_OBS EVALUATE_THIS_TYPE_OF_OBS GET_QUANTITY_FOR_TYPE_OF_OBS GET_INDEX_FOR_TYPE_OF_OBS WRITE_TYPE_OF_OBS_TABLE READ_TYPE_OF_OBS_TABLE GET_TYPE_OF_OBS_FROM_MENU MAP_TYPE_OF_OBS_TABLE USE_EXT_PRIOR_THIS_TYPE_OF_OBS GET_NAME_FOR_QUANTITY GET_INDEX_FOR_QUANTITY SET_NAMEVALUE_FOR_QUANTITY GET_NUM_ITEMS_FOR_QUANTITY GET_ITEMNAME_FOR_QUANTITY GET_ITEMVALUE_FOR_QUANTITY HAS_BOUNDS_FOR_QUANTITY GET_NUM_TYPES_OF_OBS GET_NUM_QUANTITIES QTY_STATE_VARIABLE QTY_SEAICE_AGREG_CONCENTR QTY_SEAICE_AGREG_FREEBOARD QTY_SEAICE_AGREG_VOLUME QTY_SEAICE_AGREG_SNOWVOLUME QTY_SEAICE_AGREG_THICKNESS QTY_SEAICE_AGREG_SNOWDEPTH QTY_U_SEAICE_COMPONENT QTY_V_SEAICE_COMPONENT QTY_SEAICE_ALBEDODIRVIZ QTY_SEAICE_ALBEDODIRNIR QTY_SEAICE_ALBEDOINDVIZ QTY_SEAICE_ALBEDOINDNIR QTY_SEAICE_CATEGORY QTY_SEAICE_CONCENTR QTY_SEAICE_VOLUME QTY_SEAICE_SNOWVOLUME QTY_SEAICE_SURFACETEMP QTY_SEAICE_FIRSTYEARAREA QTY_SEAICE_ICEAGE QTY_SEAICE_LEVELAREA QTY_SEAICE_LEVELVOLUME QTY_SEAICE_MELTPONDAREA QTY_SEAICE_MELTPONDDEPTH QTY_SEAICE_MELTPONDLID QTY_SEAICE_MELTPONDSNOW QTY_SEAICE_SALINITY001 QTY_SEAICE_SALINITY002 QTY_SEAICE_SALINITY003 QTY_SEAICE_SALINITY004 QTY_SEAICE_SALINITY005 QTY_SEAICE_SALINITY006 QTY_SEAICE_SALINITY007 QTY_SEAICE_SALINITY008 QTY_SEAICE_ICEENTHALPY001 QTY_SEAICE_ICEENTHALPY002 QTY_SEAICE_ICEENTHALPY003 QTY_SEAICE_ICEENTHALPY004 QTY_SEAICE_ICEENTHALPY005 QTY_SEAICE_ICEENTHALPY006 QTY_SEAICE_ICEENTHALPY007 QTY_SEAICE_ICEENTHALPY008 QTY_SEAICE_SNOWENTHALPY001 QTY_SEAICE_SNOWENTHALPY002 QTY_SEAICE_SNOWENTHALPY003 QTY_SOM_TEMPERATURE QTY_SEAICE_FY QTY_SEAICE_AGREG_FY QTY_SEAICE_AGREG_SURFACETEMP QTY_SALINITY QTY_TEMPERATURE QTY_POTENTIAL_TEMPERATURE QTY_PRESSURE QTY_VELOCITY QTY_U_CURRENT_COMPONENT QTY_V_CURRENT_COMPONENT QTY_W_CURRENT_COMPONENT QTY_DRY_LAND QTY_SEA_SURFACE_PRESSURE QTY_SEA_SURFACE_HEIGHT QTY_SEA_SURFACE_ANOMALY QTY_U_WIND_COMPONENT QTY_V_WIND_COMPONENT MAX_DEFINED_QUANTITIES SYN_SEAICE_CONCENTR SAT_U_SEAICE_COMPONENT SAT_V_SEAICE_COMPONENT SAT_SEAICE_CONCENTR SAT_SEAICE_VOLUME SAT_SEAICE_SNOWVOLUME SAT_SEAICE_SURFACETEMP SAT_SEAICE_FY SAT_SEAICE_AGREG_FY SAT_SEAICE_AGREG_SURFACETEMP SAT_SEAICE_AGREG_FREEBOARD SAT_SEAICE_AGREG_CONCENTR SAT_SEAICE_AGREG_VOLUME SAT_SEAICE_AGREG_SNOWVOLUME SAT_SEAICE_AGREG_THICKNESS SAT_SEAICE_AGREG_SNOWDEPTH MAX_DEFINED_TYPES_OF_OBS                      @                              
       OBSTYPELENGTH R8 MISSING_R8                      @                              
  
     ERROR_HANDLER E_ERR E_WARN LOGFILEUNIT FIND_NAMELIST_IN_FILE LOG_IT CHECK_NAMELIST_READ DO_OUTPUT ASCII_FILE_FORMAT STRING_TO_REAL $         @                                                           #OBS_TYPE_IND                      
  @                                          %         @                                                           #OBS_TYPE_IND              
  @                                          %         @                                                           #OBS_TYPE_IND              
  @                                          %         @                                 	                           #OBS_TYPE_IND 
             
  @                               
           %         @                                                           #OBS_TYPE_NAME              
                                                    1 #         @                                                       #IFILE    #FFORM    #USE_LIST              
                                                       
 @                                                 1           
@                                                   
             &                                           #         @                                                       #IFILE    #PRE_I_FORMAT    #FFORM              
                                                       
                                                       
 @                                                 1 %         @                                                             %         @                                                            #OBS_DEF_INDEX              
                                             %         @                                                            #OBS_TYPE_IND              
  @                                          $         @                                                          #OBS_QTY_IND                      
  @                                          %         @                                                            #OBS_QTY_NAME              
                                                    1 #         @                                                       #OBS_QTY_IND    #ITEMNAME     #ITEMVALUE !             
  @                                                    
  @                                                  1           
  @                             !                    1 %         @                                 "                           #OBS_QTY_IND #             
  @                               #           $         @                                $                           #OBS_QTY_IND %   #ITEM_INDEX &                     
  @                               %                     
  @                               &           $         @                                '                           #OBS_QTY_IND (   #ITEMNAME )                     
  @                               (                     
                                )                    1 %         @                                 *                           #OBS_QTY_IND +   #MINBOUND ,   #MAXBOUND -             
  @                               +                     D                                ,     
                 D                                -     
       %         @                                 .                            %         @                                 /                                                                         0                                                       0                                             1                                                      1                                             2                                                      2                                             3                                                      3                                             4                                                      4                                             5                                                      5                                             6                                                      6                                             7                                                      7                                             8                                                      8                                             9                                       	               9                                             :                                       
               10                                             ;                                                      11                                             <                                                      12                                             =                                                      13                                             >                                                      14                                             ?                                                      15                                             @                                                      16                                             A                                                      17                                             B                                                      18                                             C                                                      19                                             D                                                      20                                             E                                                      21                                             F                                                      22                                             G                                                      23                                             H                                                      24                                             I                                                      25                                             J                                                      26                                             K                                                      27                                             L                                                      28                                             M                                                      29                                             N                                                      30                                             O                                                      31                                             P                                                       32                                             Q                                       !               33                                             R                                       "               34                                             S                                       #               35                                             T                                       $               36                                             U                                       %               37                                             V                                       &               38                                             W                                       '               39                                             X                                       (               40                                             Y                                       )               41                                             Z                                       *               42                                             [                                       +               43                                             \                                       ,               44                                             ]                                       -               45                                             ^                                       .               46                                             _                                       /               47                                             `                                       0               48                                             a                                       1               49                                             b                                       2               50                                             c                                       3               51                                             d                                       4               52                                             e                                       5               53                                             f                                       6               54                                             g                                       7               55                                             h                                       8               56                                             i                                       9               57                                             j                                       :               58                                             k                                       ;               59                                             l                                       <               60                                             m                                       =               61                                             n                                       >               62                                             o                                       >               62                                             p                                                      1                                             q                                                      2                                             r                                                      3                                             s                                                      4                                             t                                                      5                                             u                                                      6                                             v                                                      7                                             w                                                      8                                             x                                       	               9                                             y                                       
               10                                             z                                                      11                                             {                                                      12                                             |                                                      13                                             }                                                      14                                             ~                                                      15                                                                                                   16                                             �                                                      16   �   V      fn#fn "   �   1	  b   uapp(OBS_KIND_MOD    '
  \   J  TYPES_MOD    �
  �   J  UTILITIES_MOD )   F  j       GET_NAME_FOR_TYPE_OF_OBS 6   �  @   a   GET_NAME_FOR_TYPE_OF_OBS%OBS_TYPE_IND ,   �  b       ASSIMILATE_THIS_TYPE_OF_OBS 9   R  @   a   ASSIMILATE_THIS_TYPE_OF_OBS%OBS_TYPE_IND *   �  b       EVALUATE_THIS_TYPE_OF_OBS 7   �  @   a   EVALUATE_THIS_TYPE_OF_OBS%OBS_TYPE_IND -   4  b       GET_QUANTITY_FOR_TYPE_OF_OBS :   �  @   a   GET_QUANTITY_FOR_TYPE_OF_OBS%OBS_TYPE_IND *   �  c       GET_INDEX_FOR_TYPE_OF_OBS 8   9  L   a   GET_INDEX_FOR_TYPE_OF_OBS%OBS_TYPE_NAME (   �  l       WRITE_TYPE_OF_OBS_TABLE .   �  @   a   WRITE_TYPE_OF_OBS_TABLE%IFILE .   1  L   a   WRITE_TYPE_OF_OBS_TABLE%FFORM 1   }  �   a   WRITE_TYPE_OF_OBS_TABLE%USE_LIST '   	  p       READ_TYPE_OF_OBS_TABLE -   y  @   a   READ_TYPE_OF_OBS_TABLE%IFILE 4   �  @   a   READ_TYPE_OF_OBS_TABLE%PRE_I_FORMAT -   �  L   a   READ_TYPE_OF_OBS_TABLE%FFORM *   E  P       GET_TYPE_OF_OBS_FROM_MENU &   �  c       MAP_TYPE_OF_OBS_TABLE 4   �  @   a   MAP_TYPE_OF_OBS_TABLE%OBS_DEF_INDEX /   8  b       USE_EXT_PRIOR_THIS_TYPE_OF_OBS <   �  @   a   USE_EXT_PRIOR_THIS_TYPE_OF_OBS%OBS_TYPE_IND &   �  i       GET_NAME_FOR_QUANTITY 2   C  @   a   GET_NAME_FOR_QUANTITY%OBS_QTY_IND '   �  b       GET_INDEX_FOR_QUANTITY 4   �  L   a   GET_INDEX_FOR_QUANTITY%OBS_QTY_NAME +   1  v       SET_NAMEVALUE_FOR_QUANTITY 7   �  @   a   SET_NAMEVALUE_FOR_QUANTITY%OBS_QTY_IND 4   �  L   a   SET_NAMEVALUE_FOR_QUANTITY%ITEMNAME 5   3  L   a   SET_NAMEVALUE_FOR_QUANTITY%ITEMVALUE +     a       GET_NUM_ITEMS_FOR_QUANTITY 7   �  @   a   GET_NUM_ITEMS_FOR_QUANTITY%OBS_QTY_IND *      y       GET_ITEMNAME_FOR_QUANTITY 6   �  @   a   GET_ITEMNAME_FOR_QUANTITY%OBS_QTY_IND 5   �  @   a   GET_ITEMNAME_FOR_QUANTITY%ITEM_INDEX +     w       GET_ITEMVALUE_FOR_QUANTITY 7   �  @   a   GET_ITEMVALUE_FOR_QUANTITY%OBS_QTY_IND 4   �  L   a   GET_ITEMVALUE_FOR_QUANTITY%ITEMNAME (     }       HAS_BOUNDS_FOR_QUANTITY 4   �  @   a   HAS_BOUNDS_FOR_QUANTITY%OBS_QTY_IND 1   �  @   a   HAS_BOUNDS_FOR_QUANTITY%MINBOUND 1     @   a   HAS_BOUNDS_FOR_QUANTITY%MAXBOUND %   Y  P       GET_NUM_TYPES_OF_OBS #   �  P       GET_NUM_QUANTITIES #   �  q       QTY_STATE_VARIABLE *   j  q       QTY_SEAICE_AGREG_CONCENTR +   �  q       QTY_SEAICE_AGREG_FREEBOARD (   L  q       QTY_SEAICE_AGREG_VOLUME ,   �  q       QTY_SEAICE_AGREG_SNOWVOLUME +   .  q       QTY_SEAICE_AGREG_THICKNESS +   �  q       QTY_SEAICE_AGREG_SNOWDEPTH '     q       QTY_U_SEAICE_COMPONENT '   �  q       QTY_V_SEAICE_COMPONENT (   �  q       QTY_SEAICE_ALBEDODIRVIZ (   c  r       QTY_SEAICE_ALBEDODIRNIR (   �  r       QTY_SEAICE_ALBEDOINDVIZ (   G  r       QTY_SEAICE_ALBEDOINDNIR $   �  r       QTY_SEAICE_CATEGORY $   +   r       QTY_SEAICE_CONCENTR "   �   r       QTY_SEAICE_VOLUME &   !  r       QTY_SEAICE_SNOWVOLUME '   �!  r       QTY_SEAICE_SURFACETEMP )   �!  r       QTY_SEAICE_FIRSTYEARAREA "   e"  r       QTY_SEAICE_ICEAGE %   �"  r       QTY_SEAICE_LEVELAREA '   I#  r       QTY_SEAICE_LEVELVOLUME (   �#  r       QTY_SEAICE_MELTPONDAREA )   -$  r       QTY_SEAICE_MELTPONDDEPTH '   �$  r       QTY_SEAICE_MELTPONDLID (   %  r       QTY_SEAICE_MELTPONDSNOW '   �%  r       QTY_SEAICE_SALINITY001 '   �%  r       QTY_SEAICE_SALINITY002 '   g&  r       QTY_SEAICE_SALINITY003 '   �&  r       QTY_SEAICE_SALINITY004 '   K'  r       QTY_SEAICE_SALINITY005 '   �'  r       QTY_SEAICE_SALINITY006 '   /(  r       QTY_SEAICE_SALINITY007 '   �(  r       QTY_SEAICE_SALINITY008 *   )  r       QTY_SEAICE_ICEENTHALPY001 *   �)  r       QTY_SEAICE_ICEENTHALPY002 *   �)  r       QTY_SEAICE_ICEENTHALPY003 *   i*  r       QTY_SEAICE_ICEENTHALPY004 *   �*  r       QTY_SEAICE_ICEENTHALPY005 *   M+  r       QTY_SEAICE_ICEENTHALPY006 *   �+  r       QTY_SEAICE_ICEENTHALPY007 *   1,  r       QTY_SEAICE_ICEENTHALPY008 +   �,  r       QTY_SEAICE_SNOWENTHALPY001 +   -  r       QTY_SEAICE_SNOWENTHALPY002 +   �-  r       QTY_SEAICE_SNOWENTHALPY003 $   �-  r       QTY_SOM_TEMPERATURE    k.  r       QTY_SEAICE_FY $   �.  r       QTY_SEAICE_AGREG_FY -   O/  r       QTY_SEAICE_AGREG_SURFACETEMP    �/  r       QTY_SALINITY     30  r       QTY_TEMPERATURE *   �0  r       QTY_POTENTIAL_TEMPERATURE    1  r       QTY_PRESSURE    �1  r       QTY_VELOCITY (   �1  r       QTY_U_CURRENT_COMPONENT (   m2  r       QTY_V_CURRENT_COMPONENT (   �2  r       QTY_W_CURRENT_COMPONENT    Q3  r       QTY_DRY_LAND )   �3  r       QTY_SEA_SURFACE_PRESSURE '   54  r       QTY_SEA_SURFACE_HEIGHT (   �4  r       QTY_SEA_SURFACE_ANOMALY %   5  r       QTY_U_WIND_COMPONENT %   �5  r       QTY_V_WIND_COMPONENT '   �5  r       MAX_DEFINED_QUANTITIES $   o6  q       SYN_SEAICE_CONCENTR '   �6  q       SAT_U_SEAICE_COMPONENT '   Q7  q       SAT_V_SEAICE_COMPONENT $   �7  q       SAT_SEAICE_CONCENTR "   38  q       SAT_SEAICE_VOLUME &   �8  q       SAT_SEAICE_SNOWVOLUME '   9  q       SAT_SEAICE_SURFACETEMP    �9  q       SAT_SEAICE_FY $   �9  q       SAT_SEAICE_AGREG_FY -   h:  r       SAT_SEAICE_AGREG_SURFACETEMP +   �:  r       SAT_SEAICE_AGREG_FREEBOARD *   L;  r       SAT_SEAICE_AGREG_CONCENTR (   �;  r       SAT_SEAICE_AGREG_VOLUME ,   0<  r       SAT_SEAICE_AGREG_SNOWVOLUME +   �<  r       SAT_SEAICE_AGREG_THICKNESS +   =  r       SAT_SEAICE_AGREG_SNOWDEPTH )   �=  r       MAX_DEFINED_TYPES_OF_OBS 