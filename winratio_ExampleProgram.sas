/******************************************************************************/
/* PROD PROGRAM:    /home/reprieve/POOLED/final/f2024_winratio/publication/macros/winratio_ExampleProgram.sas
/* WORK PROGRAM:    /home/reprieve/POOLED/final/f2024_winratio/publication/macros/winratio_ExampleProgram.sas
/* 
/* PURPOSE:         Example program for the WINRATIO suite of macros supplied with Davies-Smith (CID 2025 :)
/* 
/* INPUT:           /home/reprieve/POOLED/final/f2024_winratio/sas_datasets/derived/adttwinratio.sas7bdat
/*
/* MACROS USED:     /home/reprieve/POOLED/f2024_winratio/publication/macros/winratio.sas
/*                  /home/reprieve/POOLED/f2024_winratio/publication/macros/StratifiedWR.sas
/*
/* AUTHOR:          Heather Ribaudo
/* CREATION DATE:   09JUN2025
/* 
/* NOTES:           See Pocock Eu J Cardiology, 2012
/******************************************************************************/

%include "winratio.sas" ;
%include "StratifiedWR.sas" ;


/* ************************************************************************************* */
/* OVERALL WIN RATIO */
%winratio( idsn    = ADTTWR 
         , odsn    = winratio_unstrat 
         , trtvar  = TRT01P 
         , trtval  = Pitavastatin
         , prntlst = /* %str('5002' '5003' '5022' '5184' '5275' '5444') */ 
         , debug   = 
         /* NO VALUE GIVEN FOR DEBUG WILL GIVE STEP BY STEP PRINTOUT OF PROCESSING */
         /* STRING OF IDENTIFIERS IN PRNTLST IDENTIFIES A SMALL SUBSET OF USUBJIDs TO PRINT */
         ) ;
title  "*********************************************************************" ;
title2 "Basic output from winratio macro - Dataset: WINRATIO_UNSTRAT";
proc print data=winratio_unstrat heading=horizontal noobs label;
run;


/* ************************************************************************************* */
/*  STRATIFIED WR - ONE STRATIFICATION FACTOR */
%StratifiedWR( idsn      = ADTTWR 
             , odsn      = winratio_SEX 
             , trtvar    = TRT01P 
             , trtval    = Pitavastatin
             , stratvars = NATALSEX 
             ) ; 
title  "*********************************************************************" ;
title2 "Basic output from StratifiedWR macro - Dataset: WINRATIO_SEX";
proc print data=winratio_SEX heading=horizontal noobs label;
run;


/* ************************************************************************************* */
/*  STRATIFIED WR - TWO STRATIFICATION FACTORS */
%StratifiedWR( idsn      = ADTTWR 
             , odsn      = winratio_SEXGDB 
             , trtvar    = TRT01P 
             , trtval    = Pitavastatin
             , stratvars = NATALSEX GBDGP1
             ) ; 
title  "*************************************************************************" ;
title2 "Basic output from StratifiedWR macro - Dataset: WINRATIO_SEXGBD";
proc print data=winratio_SEXGDB heading=horizontal noobs label;
run;


ENDSAS ; /*******************************************************/


/**************************************************************************/
                                                      
/* First obs in ADTTWR - key variables */

/* USUBJID  TRT01P        PARAM                             PARAMCD  PARAMN    AVAL   CNSR  EVNT    GBDGP1     natalsex */

/*  5001    Placebo       Tier 1: CV or Undetermined Death   EVT1       1    7.93429    1     0   High Income   Male    */
/*  5001    Placebo       Tier 2: Stroke                     EVT2       2    7.93429    1     0   High Income   Male    */
/*  5001    Placebo       Tier 3: MI                         EVT3       3    7.93429    1     0   High Income   Male    */
/*  5001    Placebo       Tier 4: Other CV Event             EVT4       4    7.93429    1     0   High Income   Male    */
/*  5001    Placebo       Tier 5: CV Procedure               EVT5       5    7.93429    1     0   High Income   Male    */
/*  5003    Placebo       Tier 1: CV or Undetermined Death   EVT1       1    6.31348    1     0   High Income   Male    */
/*  5003    Placebo       Tier 2: Stroke                     EVT2       2    6.31348    1     0   High Income   Male    */
/*  5003    Placebo       Tier 3: MI                         EVT3       3    5.22382    0     1   High Income   Male    */
/*  5003    Placebo       Tier 4: Other CV Event             EVT4       4    6.31348    1     0   High Income   Male    */
/*  5003    Placebo       Tier 5: CV Procedure               EVT5       5    6.31348    1     0   High Income   Male    */
/*  5184    Pitavastatin  Tier 1: CV or Undetermined Death   EVT1       1    5.22382    0     1   High Income   Male    */
/*  5184    Pitavastatin  Tier 2: Stroke                     EVT2       2    5.22382    2     0   High Income   Male    */
/*  5184    Pitavastatin  Tier 3: MI                         EVT3       3    5.22382    2     0   High Income   Male    */
/*  5184    Pitavastatin  Tier 4: Other CV Event             EVT4       4    5.22382    2     0   High Income   Male    */
/*  5184    Pitavastatin  Tier 5: CV Procedure               EVT5       5    5.22382    2     0   High Income   Male    */
/*  5275    Pitavastatin  Tier 1: CV or Undetermined Death   EVT1       1    8.01369    1     0   High Income   Male    */
/*  5275    Pitavastatin  Tier 2: Stroke                     EVT2       2    8.01369    1     0   High Income   Male    */
/*  5275    Pitavastatin  Tier 3: MI                         EVT3       3    8.01369    1     0   High Income   Male    */
/*  5275    Pitavastatin  Tier 4: Other CV Event             EVT4       4    7.24709    0     1   High Income   Male    */
/*  5275    Pitavastatin  Tier 5: CV Procedure               EVT5       5    8.01369    1     0   High Income   Male    */
 
