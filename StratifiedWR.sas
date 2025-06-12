/***********************************************************************************************************************/
/* PROD PROGRAM:    /home/reprieve/POOLED/f2024_winratio/publication/macros/StratifiedWR.sas
/* WORK PROGRAM:    /home/reprieve/POOLED/f2024_winratio/publication/macros/StratifiedWR.sas
/* 
/* PURPOSE:         Wrapper macro to perform stratified analyses with repeated running of the %winratio macro across
/*                    the individual strata and then combining the results per Dong 2018.
/*
/* SOURCE PRGM:     /home/reprieve/POOLED/macros/StratifiedWR.sas
/*
/* INPUT:           Has same input dataset requirements as %winratio macro.
/*                  Dataset (&IDSN) in the form of CDISC BDS TTE. Dataset needs a parameter defined for each event tier
/*                    including the following:
/*                   AVAL = Time of event or censoring for that tier
/*                   EVNT = Event indicator (1=vent/0=No event)
/*                  PARAM = Description of the event for the tier
/*                 PARAMN = Numeric index for the tier from 1-n where n>1
/*                  Plus USUBJID and treatment code variable (defined by macro parameter TRTVAR)
/*
/* MACRO INPUT PARAMETERS:
/*                    IDSN : Name of input dataset in the form of a BDS TTE
/*                    ODSN : Desired name for the output dataset
/*                  TRTVAR : Name of treatment group variable
/*                  TRTVAL : Value of TRTVAR representing treatment
/*               STRATVARS : Variable(s) defining stratification
/*                   DEBUG : Debug mode (DEFAULT=CANCEL - no debug  - all LST output will be suppressed)
/*
/* OUTPUT:          Dataset (&ODSN) giving the Treatment and control win counts along with the Win Ratio
/*                     and its variance/covariance.
/*
/* MACROS USED:     Calls %winratio macro and its submacros defined for variance estimation.
/* 
/* AUTHOR:          Heather Ribaudo
/* CREATION DATE:   06JUN2025 (This version)
/***********************************************************************************************************************/

%macro StratifiedWR( IDSN= , ODSN= , TRTVAR= , TRTVAL= , STRATVARS= , DEBUG=CANCEL ) ;

    title1 "*** STRATIFIED WR MACRO EXECUTION *** " ;

    %local ntiers nstratvars stratfreq stratcomma strnum I ;

    %let nstratvars = %sysfunc(countw(&STRATVARS)) ; /* HOW MANY STRATIFICATION VARIABLES */

    /* [HJR20250109] DEFINE HOW MANY TIERS - USED AT THE END OF THE MACRO IN RESULTS ASSEMBLY */
    title2 "How many event tiers?";
    proc sql %if %length(&DEBUG)>0 %then noprint ;;
        SELECT max(PARAMN) INTO :ntiers TRIMMED FROM &IDSN ;
    quit ;
    %if %eval(&NTIERS<2) %then %do ;
       %put ERROR: Dataset needs at least 2 tiers ;
       %return ;
    %end ;

    /* , AND * DELIMITED STRINGS OF THE STRATIFICATION VARIABLES - FOR PROC FREQ AND SQL CALLS */
    %let stratfreq  = %sysfunc(tranwrd(&STRATVARS,%str( ),%str(*))) ;
    %let stratcomma = %sysfunc(tranwrd(&STRATVARS,%str( ),%str(,))) ;
    %put RESOLVED MACRO VARIABLES DEFINING STRATIFICATION FACTORS: &=NSTRATVARS &=STRATFREQ &=STRATCOMMA ;

    /* CREATE AN IDENTIFIER FOR EACH STRATUM 
      - BASICALLY ALL OBSERVED COMBINATIONS OF STRATA VARIABLES COMBINED INTO SINGLE STRINGS IN ONE VARIABLE */
    ODS output %if &NSTRATVARS=1 %then OneWayFreqs=STRATA ;
               %else                          List=STRATA ;;
    title2 "Stratification levels: Frequency of &STRATFREQ FROM &IDSN(where=(PARAMN=1))" ;
    proc freq data=&IDSN(where=(PARAMN=1)) ;
        tables &STRATFREQ / list ;
    run ;
    data STRATA ;
        set STRATA  end=eof;
        attrib stratnum  length=3   label="Numeric Stratum counter" ;
        attrib stratdesc length=$60 label="Stratum description" ;
        stratnum = _N_ ;
        stratdesc = catx("",&STRATCOMMA) ;
        keep  &STRATVARS stratnum stratdesc frequency ;
        if eof then call symput ("Nlevels",compress(stratnum)) ;
    run ;
    title2 "Dataset: STRATA - All stratum levels";
    proc print data=STRATA  noobs ;
    run &DEBUG ;
    %put NUMBER OF STRATIFICATION LEVELS: &=NLEVELS ;

    /* PUT THE STRATA DETAILS BACK TO THE ORIGINAL DATASET */
    proc sort data=&IDSN ;
        by &STRATVARS ;
    run ;
    data IDSNPLUS ;
        merge &IDSN STRATA ;
        by &STRATVARS ;
    run ;
    title2 "First 5 obs in IDSNPLUS";
    proc print data=IDSNPLUS (obs=5) noobs ;
        var USUBJID &TRTVAR param paramcd aval &STRATVARS stratnum stratdesc ;
    run &DEBUG;



    /* ***************************************************************************************** */
    /*  CORE OF WHAT THE MACRO IS DOING - RUN WIN RATIO MACRO OVER THE STRATA */

    title1 "*** STRATIFIED WR MACRO EXECUTION *** " ;
    %winratio( idsn=IDSNPLUS, odsn= STRATRESULTS , trtvar=&TRTVAR , trtval=&TRTVAL , byvar= STRATNUM , debug=&DEBUG ) ;

    /* ***************************************************************************************** */


    /* MERGE TOGETHER THE STRATUM SPECIFIC RESULTS AND CALCULATE A STRATIFIED ESTIMATE AND VARIANCE
                 FROM STRATUM SPECIFIC RESULTS */
    data StratResults Sum_1 ;
        merge STRATA StratResults ;
        by stratnum ;

        /* ADDITIONAL PROCESSING TO GET THE STRATIFIED STATISTIC */
        /* PER DONG et al. 2018 WEIGHTED USING INVERSE OF THE STRATUM SIZE */
        if frequency>0 then Wm = (1/frequency) ;
        WRNum    = Wm*win_trt ;
        WRDenom  = Wm*win_con ;

        varA =sig2_trt * Wm * Wm  ;
        varB =sig2_con * Wm * Wm  ;
        varC =sig_trt_con * Wm * Wm  ;
    run ;

    proc means data=Sum_1 noprint ;
        output out=OverallResults 
                sum(Frequency n_trt n_con d c: t: win_trt win_con WRNum WRDenom v:)= 
                mean(alpha)= ; /* mean alpha will return the value of alpha to the output dataset */
    run ;

    data &ODSN.;
        length stratnum 3 ;
        set OverallResults (in=in1 drop=theta_KL0) 
            StratResults 
            ;
        if in1 then do ; /* If Overall results ... */
            *stratnum= 0 ;
            stratdesc = "Overall (stratified)" ;
            if WRnum>0 and WRDenom>0 then do ;
                /* SEE DONG G 2018 (p.782) */
                WinRatio = WRNum / WRDenom ;
                lambdahat  = sum(WRNum,WRDenom) / 2 ; /* UNDER Ho */
                lambdahat2 = lambdahat*lambdahat ; 
                sig2_log_wr = (VarA+VarB-(2*VarC)) / lambdahat2 ;
                *** 100*(1-alpha)% CI of the win ratio;
                WR_L = exp(log(WinRatio) - probit(1-ALPHA/2)*sqrt(sig2_log_wr));
                WR_U = exp(log(WinRatio) + probit(1-ALPHA/2)*sqrt(sig2_log_wr));
            end ;
        end ;
        length WRCI $20 ;
        WRCI = catx("",put(WinRatio,4.2),cats(" (",put(WR_L,4.2),","),cats(put(WR_U,4.2),")")) ;
        label frequency = "Number of participants" 
              WRCI      = "Win Ratio (95% CI)" ;
        drop WRNum WRDenom v: _type_ _freq_ lambdahat ;
    run ;

    title2 "DATASET: &ODSN - Stratified Results...  ";
    proc print data=&ODSN noobs ;
    run &DEBUG ;

%mend ;
