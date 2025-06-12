/***********************************************************************************************************************/
/* PROD PROGRAM:    /home/reprieve/POOLED/f2024_winratio/publication/macros/winratio.sas
/* WORK PROGRAM:    /home/reprieve/POOLED/f2024_winratio/publication/macros/winratio.sas
/* 
/* PURPOSE:         CALCULATES THE WIN RATIO AND VAR FROM INPUT DATA
/* 
/* SOURCE PRGM:     /home/reprieve/POOLED/macros/winratio.sas
/*                  Original code taken from the online appendix from Dong G Pharaceutical Statistics 2016
/*
/* INPUT:           Dataset (&IDSN) in the form of CDISC BDS TTE. 
/*
/*                  Dataset needs a parameter defined for each event tier including the following:
/*                   AVAL = Time of event or censoring for that tier
/*                   EVNT = Event indicator (1=vent/0=No event)
/*                  PARAM = Description of the event for the tier
/*                 PARAMN = Numeric index for the tier from 1-n where n is 2 or higher 
/*                  BYVAR = Any BY variables to run analysis by 
/*
/*                  Plus USUBJID and treatment code variable (defined by macro parameter &TRTVAR)
/*            E.g., 
/*            USUBJID    TRT01P          PARAM                               PARAMCD    PARAMN      AVAL     CNSR    EVNT  
/*            
/*             5001      Placebo         Tier 1: CV or Undetermined Death     EVT1         1      7.93429      1       0   
/*             5001      Placebo         Tier 2: Stroke                       EVT2         2      7.93429      1       0   
/*             5001      Placebo         Tier 3: MI                           EVT3         3      7.93429      1       0   
/*             5001      Placebo         Tier 4: Other CV Event               EVT4         4      7.93429      1       0   
/*             5001      Placebo         Tier 5: CV Procedure                 EVT5         5      7.93429      1       0   
/*             5002      Placebo         Tier 1: CV or Undetermined Death     EVT1         1      3.58111      1       0   
/*             5002      Placebo         Tier 2: Stroke                       EVT2         2      3.58111      1       0   
/*             5002      Placebo         Tier 3: MI                           EVT3         3      3.58111      1       0   
/*             5002      Placebo         Tier 4: Other CV Event               EVT4         4      3.58111      1       0   
/*             5002      Placebo         Tier 5: CV Procedure                 EVT5         5      3.58111      1       0   
/*             5184      Pitavastatin    Tier 1: CV or Undetermined Death     EVT1         1      5.22382      0       1   
/*             5184      Pitavastatin    Tier 2: Stroke                       EVT2         2      5.22382      2       0   
/*             5184      Pitavastatin    Tier 3: MI                           EVT3         3      5.22382      2       0   
/*             5184      Pitavastatin    Tier 4: Other CV Event               EVT4         4      5.22382      2       0   
/*             5184      Pitavastatin    Tier 5: CV Procedure                 EVT5         5      5.22382      2       0   
/*             5275      Pitavastatin    Tier 1: CV or Undetermined Death     EVT1         1      8.01369      1       0   
/*             5275      Pitavastatin    Tier 2: Stroke                       EVT2         2      8.01369      1       0   
/*             5275      Pitavastatin    Tier 3: MI                           EVT3         3      8.01369      1       0   
/*             5275      Pitavastatin    Tier 4: Other CV Event               EVT4         4      7.24709      0       1   
/*             5275      Pitavastatin    Tier 5: CV Procedure                 EVT5         5      8.01369      1       0   
/*                                                                                                                         
/*
/* MACRO INPUT PARAMETERS:
/*                    IDSN : Name of input dataset in the form of a BDS TTE
/*                    ODSN : Desired name for the output dataset
/*                  TRTVAR : Name of treatment group variable
/*                  TRTVAL : Value of TRTVAR representing treatment (case-sensitive, not in quotes)
/*                   BYVAR : Space delimited list of BY variables
/*                   ALPHA : Alpha level for CI (DEFAULT=0.05)  
/*                      SE : Is the variance/SE under the null or the alternative
/*                 PRNTLST : List of participant identifiers to print to listing (DEFAULT= (no USUBJID selected))  
/*                   DEBUG : Debug mode (DEFAULT=CANCEL - no debug most LST output will be suppressed)
/*                            
/* OUTPUT:          Dataset (&ODSN) giving the Treatment and control win counts along with the Win Ratio
/*                     and its variance/covariance.
/*
/* AUTHOR:          Heather Ribaudo
/* CREATION DATE:   01June2025 (this version for distribution)
/* 
/* MODIFICATIONS:   Updates to source macro from Dong
/*                  - Added pre-processing of the input dataset to allow data source to be a CDISC BDS TTE
/*                  - Extended from 3 fixed tiers to a flexible number of tiers (2 or more)
/*                  - Changed variable and dataset names throughout the variance and covariance calculations to be more
/*                      transparent as to what is happenining
/*                  - Changed processing to avoid sorting and dataset merges to reduce processing time
/*                  - Standard error for the CI calculated based on the alternative (estimated treatment win prob rather 
/*                      than theta_0)
/***********************************************************************************************************************/



/* **************************************************************************************************** */                                                                                                                

%macro winratio( IDSN= , ODSN= , TRTVAR = , TRTVAL= , BYVAR= , ALPHA=0.05 , SE=NULL , PRNTLST = '99999' , DEBUG=CANCEL ) ;

    %local i NTIERS NBY BYSQL SQLON ;
    
    title1 "**** WIN RATIO MACRO EXECUTION ****" ;

    /* CREATE MACRO VARIABLES FOR SQL PROCESSING OF BY VARIABLES */
    %if %length(&BYVAR)>0 %then %do ;
        /* CREATE MACRO VARIABLE THAT IS A COUNT OF THE NUMBER OF BY VARIABLES */
        %let NBY = %sysfunc(countw(&BYVAR)) ;
        %let BYSQL = a.%scan(&BYVAR,1) ;
        %let SQLON = a.%scan(&BYVAR,1) = b.%scan(&BYVAR,1) ;
        
        %if &NBY > 1 %then %do ;
            %do i =1 %to %sysfunc(countw(&BYVAR)) ;
                %let BYSQL = &BYSQL , a.%scan(&BYVAR,1) ;
                %let SQLON = &SQLON AND a.%scan(&BYVAR,&I) = b.%scan(&BYVAR,&I) ;
            %end ;
        %end ;
    %end ;

    /* IF PRNTLST IS EMPTY - PROVIDE A NULL VALUE */
    %if %length(&PRNTLST)=0 %then %let prntlst = ' ' ;;
    
    /* INDEX DATASET RATHER THAN SORTING IT */
    data IDSN (index=(ID=(&BYVAR USUBJID &TRTVAR))) ;
        set &IDSN ;
    run ;

    /* TRANSPOSE THE INPUT BDS TTE DATASET TO GET THE DESIRED STRUCTURE FOR THE SOURCE MACRO (from Dong 2016) */
    /* FIRST THE TIME OF THE EVENT OR CENSORING */
    proc transpose data=IDSN out=EVENTTIMES prefix=TM ;
        by      &BYVAR USUBJID &TRTVAR ;
        var     AVAL ;
        id      PARAMN ;
        idlabel PARAM ;
    run ;
    title2 "Dataset: EVENTTIMES(where=(USUBJID in (&PRNTLST)))" ;
    proc print data=EVENTTIMES (where=(USUBJID in (&PRNTLST))) noobs ;
    run &DEBUG ;

    /* THEN THE EVENT/CENSORING INDICATOR */
    proc transpose data=IDSN out=EVNT prefix=EVT ;
        by      &BYVAR USUBJID &TRTVAR ;
        var     EVNT ;
        id      PARAMN ;
        idlabel PARAM ;
    run ;
    title2 "Dataset: EVNT(where=(USUBJID in (&PRNTLST)))" ;
    proc print data=EVNT (where=(USUBJID in (&PRNTLST))) noobs ;
    run &DEBUG ;
    
    /* TOTAL NUMBER OF TIERS */
    title2 "How many event tiers?";
    proc sql %if %length(&DEBUG)>0 %then noprint ;;
        SELECT max(PARAMN) INTO :ntiers TRIMMED FROM &IDSN ;
    quit ;
    %if %eval(&NTIERS<2) %then %do ;
       %put ERROR: Dataset needs at least 2 tiers ;
       %return ;
    %end ;

    /*****************************************************************************************************/
    /********* Rename patient ID, event, event day and censoring day                                   ***/
    /*****************************************************************************************************/
    /*** [HJR20241031] Updates to previous code through this section to process tier specific parameters
                       with do loops stopping at NTIERS */
    /*** Treatment Group ***/                   
    /*** Naming convention: _trt for treatment group ***/
    /***   For example, evt1_trt = 'Event of Endpoint 1 (1=Yes, 0 = No) in Treatment group' ***/
    data trt (keep = pid_trt &BYVAR %do t=1 %to &NTIERS ;
                                          evt&T._trt day&T._trt 
                                   %end ;
                                   )  ;
        merge EVNT EVENTTIMES ;
        by &BYVAR USUBJID ;
        if &TRTVAR = "&TRTVAL" ;
  
        rename USUBJID=pid_trt 
               %do t=1 %to &NTIERS ; evt&T=evt&T._trt tm&T=day&T._trt  %end ;;
    run ;
    title2 "Dataset: TRT(where=(pid_trt in (&PRNTLST)))" ;
    proc print data=TRT (where=(pid_trt in (&PRNTLST))) noobs ;
    run &DEBUG ;

    /*** Control Group ***/
    /*** Naming convention: _con for Control group ***/
    /***   For example, evt1_con = 'Event of Endpoint 1 (1=Yes, 0 = No) in Control group' ***/
    data con (keep = pid_con &BYVAR %do t=1 %to &NTIERS ;
                                          evt&T._con day&T._con 
                                   %end ;
                                   ) ;
        merge EVNT EVENTTIMES ;
        by &BYVAR USUBJID ;
        if &TRTVAR ne "&TRTVAL" ;
  
        rename USUBJID=pid_con  
               %do t=1 %to &NTIERS ; evt&T=evt&T._con tm&T=day&T._con %end ;;
    run ;
    title2 "Dataset: CON(where=(pid_trt in (&PRNTLST)))" ;
    proc print data=CON (where=(pid_con in (&PRNTLST))) noobs ;
    run &DEBUG ;
    

    /*****************************************************************************************************/
    /********* Compare and determine winners/losers/ties                                               ***/
    /*****************************************************************************************************/
    /* CREATES ALL PAIRWISE COMPARISONS */
    /* CHANGE BYVAR STRING TO A COMMA DELIMITED STRING */
    proc sql;
        CREATE TABLE final AS
        SELECT   a.* /* SAME VARIABLES AS BELOW PLUS &BYVAR IF IT THERE */
               , b.pid_con %do T=1 %to &NTIERS ; , b.evt&T._con , b.day&T._con %end ;
        FROM     trt AS a
               , con AS b 
        %if %length(&BYVAR)>0 %then 
           WHERE &SQLON  ;; /* THIS ENSURES THAT PAIRS ARE ONLY CREATED WITHIN EACH BYVAR */
    quit;
    
    data compare %if &NBY = 1 %then (index=(&BYVAR)) ;
                 %if &NBY > 1 %then (index=(BY=(&BYVAR))) ;; 
        set final ;
     
        /* win = 0: winner has not been determined, win = 1: winner has been determined */
        win = 0 ;
    
        /* wincat: winning category */
        /* wincat = 'd' if a pair of patients are tied */
        /* Otherwise indexed by t or c (whether treatment or control wins) with event tier of the win */
        length wincat $2 ;
    
        /*** [HJR 20230130]: REVISED FROM THE ORIGINAL CODE FOR MORE INTUITIVE VARIABLE NAMES, VARIABLE NUMBER
                              OF TIERS AND TO ACCOUNT FOR CENSORING - OTHERWISE, THE LOGIC IS THE SAME                ***/
        /****************************************************************************************************************/
        /*** FIRST NEED TO HANDLE CENSORING **/
        /*** FROM Pocock 2012 p177 (item 2): If one patient had a Tier 1 (CV death), the other must be followed for longer
                                            in order to know definitely who had CV death first ***/
        /*** Code logic: 
               if evt1_trt > evt1_con   =>   TRT had an event - CONTROL is censored
               if day1_con < day1_trt   =>   CONTROL follow-up ended before TRT event 
                         evt1_trt = 0   =>   Censor the TRT event to induce the required tie ***/

        %do t=1 %to &NTIERS ;
            if evt&T._trt > evt&T._con AND day&T._con < day&T._trt then evt&T._trt=0 ;
        %end ;
        /* ASSIGN OVERALL EVENT INDICATOR FOR PARTICIPANT - i.e., did the participant have any events */
        evt_trt = max(evt1_trt %do t=2 %to &NTIERS ; , evt&T._trt  %end ;) ;
  
        /* if event in control but not in treatment but censoring time in treatment is before event time for control
           - censor control event (same logic as above) */
        %do t=1 %to &NTIERS ;
            if evt&T._con > evt&T._trt AND day&T._trt < day&T._con  then evt&T._con=0 ;
        %end ;
        /* ASSIGN OVERALL EVENT INDICATOR FOR PARTICIPANT - i.e., did they have any events */
        evt_con = max(evt1_con %do t=2 %to &NTIERS ; , evt&T._con  %end ;) ;
      
        /* THIS CODE IS FROM THE DONG MACRO BUT PUT INTO A DO LOOP FOR FLEXIBLE NUMBER OF TIERS */
        if evt1_trt = 1 then do ;
            if evt1_con = 1 and      day1_trt < day1_con then wincat='c1' ;  /*** Control   patient won, per Endpoint 1 */
            else if evt1_con = 1 and day1_trt > day1_con then wincat='t1' ;  /*** Treatment patient won, per Endpoint 1 */
                                                          /*** Note: If day1_trt = day1_con then winner is not declared */
            else if evt1_con < 1                         then wincat='c1' ;
  
            if wincat ne '' then win = 1 ;  /*** Winner has been determined */
        end ;
  
        %do t=2 %to &NTIERS ;    
            if win = 0 and evt&T._trt = 1 then do ;  /*** If a Winner has not been yet been declared ... */
                if evt1_con = 1 and evt1_trt=0                         then wincat='t1';   /*** Control has Endpoint 1 but tx. does not */
                %do st=2 %to %eval(&T-1) ;
                    %if &T>2 %then %do ;
                        else if evt&ST._con = 1 and evt&ST._trt=0      then wincat="t&ST"; /*** Same logic for tx. win for subseq. tiers */
                    %end ;
                %end ;
                    else if evt&T._con = 1 and day&T._trt < day&T._con then wincat="c&T" ;  
                    else if evt&T._con = 1 and day&T._trt > day&T._con then wincat="t&T" ;  
                    else if evt&T._con < 1                             then wincat="c&T" ;
      
                if wincat ne '' then win = 1 ;
            end ;
        %end ;

        if win = 0 and evt_trt <1 then do ;
            if evt1_con = 1 and evt1_trt=0                             then wincat='t1';
            %do t=2 %to &NTIERS ; 
                else if evt&T._con = 1 and evt&T._trt=0                then wincat="t&T";
            %end ;
            
            if wincat ne '' then win = 1 ;
        end ;
  
        if wincat='' then wincat='d' ;  /*** The pair of patients are tied (result is a DRAW - hence d) */
    run ;
    title2 "Check censoring of events and assignment of winners are happening as planned" ;
    title3 "Dataset: COMPARE(where=(pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST)))" ;
    proc print data=COMPARE(where=((pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST))) ) noobs ;
        format day: 6.2 ;
    run &DEBUG ;
    
    /*****************************************************************************************************/
    /******** calculate win ratio and theta_K0 (theta_L0)                                              ***/
    /*****************************************************************************************************/
    /* OBSERVED NUMBERS OF WINS BY TREATMENT GROUP AND TIER */
    title2 'Dataset: WIN - Gives numbers of wins by treatment group and tier' ;
    proc freq data=compare noprint ;
        %if %length(&BYVAR)>0 %then by &BYVAR ;;
        tables wincat / missing out=WIN ;
    run ;


    /* MAKE SURE ALL COMBINATIONS OF WINS ARE INCLUDED IN RESULTS */
    /* THIS IS DONE BY CREATING A BASIC DATASET (_win) THAT HAS ALL COMBINATIONS OF BY VARIABLES AND WINS AND LOSSES 
        THAT THE OBSERVED WINS DATASET (win) IS MERGED ON TO */

    /* FIRST CREATE A DATASET THAT IS JUST THE BY VARIABLE  */
    %if %length(&BYVAR)>0 %then %do ;
        proc sort data=win(keep=&BYVAR) out=byvars nodupkey ;
            by &BYVAR ;
        run ;
        title2 'Dataset: byvars' ;
        proc print data=byvars noobs ;
        run &DEBUG ;
    %end ;

    /* CREATE THE EMPTY DATASET ADDED ON TO THE BYVARIABLE - IF APPLICABLE */
    data _win ;  
        %if %length(&BYVAR)>0 %then set byvars ;;
        length wincat $2 wincatlabel $25 ;
        %do T = 1 %to &NTIERS ;
            wincat = "c&T" ; wincatlabel = "Control wins on Tier &T" ;   output ;
            wincat = "t&T" ; wincatlabel = "Treatment wins on Tier &T" ; output ;
        %end ;
        wincat = 'd' ; wincatlabel = "Total ties" ; output ;
    run ;
    proc sort data=_win ; by &BYVAR wincat ; run ;
    title2 'Dataset: _WIN - Empty dataset to make sure final output has combinations of by variables and wins' ;
    proc print data=_win noobs ;
    run &DEBUG ;

    
    /* MERGE TOGETHER THE RESULTS DATASET (win) WITH THE EMPTY STRUCTURAL DATASET (_win) */
    proc sort data=win ; by &BYVAR wincat ; run ;
    data win ;
        merge _win win(in=in2) ;
        by &BYVAR wincat ;
        if not in2 then count=0 ;
    run ;
    proc sort data=win ;
        by &BYVAR wincat ;
    run ;
    title2 'Dataset: WIN - Basic counts of all comparisons' ;
    proc print data=win noobs ;
    run &DEBUG ;

    /* FINAL RESULTS IN WIDE FORM TO FACILITATE PROCESSING */    
    proc transpose data=win out=win2 ;
        %if %length(&BYVAR)>0 %then by &BYVAR ;;
        var count ;
        id  wincat ;
        idlabel wincatlabel ;
    run ;
    
    /* SUMMARIZE THE RESULTS */
    data win2 ;
        set win2 ;

        win_trt = sum(of t1-t&NTIERS) ;
        win_con = sum(of c1-c&NTIERS) ;
        total = sum(of t1-t&NTIERS, c1-c&NTIERS, d) ;

        /* WIN RATIO NOT CALCULATED IF ZERO WINS IN EITHER GROUP */
        if win_con > 0 and win_trt >0 then do ;
            WinRatio = win_trt/win_con ;
            theta_KL0 = (win_trt + win_con)/(2*total) ;
            theta_KA  = win_trt/total ;
        end ;
        mergevar = 1 ; /* FACILITATES MERGING WITH VARIANCE ESTIMATION WHEN NO BYVAR */
        
        label win_trt   = "Total wins in Treatment group (&TRTVAL)"
              win_con   = "Total wins in Control group"
              total     = 'Total number of pairs'
              WinRatio  = 'Win ratio'
              theta_KL0 = 'theta K0/L0' 
              theta_KA  = 'theta KA' 
              ;

        drop _name_ _label_ ;
    run ;
    title2 'Number of winners and win ratio' ;
    title3 'Dataset: WIN2 - Summarized counts and WR statistic' ;
    proc print data=win2 label noobs ;
        var &BYVAR theta_KL0 theta_KA win_trt win_con WinRatio ;
    run &DEBUG ;
    
    /* ************************************************************************************************ */
    /* [HJR 20230131]: THE LOGIC OF THE CODE BELOW MOSTLY UNCHANGED FROM THE ORIGINAL SOURCE CODE       */
    /* UPDATES INCLUDE:                                                                                 */
    /*  0) PROCESSING OF theta_kl0 win_trt win_con winratio AS A DATASET RATHER THAN CONSTANT MACRO     */
    /*         VARIABLES TO ALLOW FOR BY VARIABLE PROCESSING                                            */            
    /*  1) ADDITION OF A MERGEBY VARIABLE (=1) TO AVOID NO BY STATEMENT MERGE WARNINGS                  */
    /*  2) SUMMING DONE IN SQL RATHER THAN MEANS TO AVOID DATASET SORTING FOR SPEED                     */
    /*  3) UPDATED DATASET AND VARIABLE NAMING CONVENTIONS                                              */
    /*  4) ADDITION OF TITLE STATEMENTS THROUGHOUT AND PUT STATEMENT TO THE LOG FOR MACRO VARIABLES     */
    /* ************************************************************************************************ */

    /* ADD ESTIMATES OF theta_kl0 win_trt win_con winratio TO THE BASE OF THE KERNAL FUNCTION DATASETS */
    /*   - THESE ARE NEEDED IN THE VARIANCE/COVARIANCE ESTIMATION */
    proc sql ;
        CREATE table KLBASE AS
        SELECT A.pid_trt , A.pid_con , A.wincat , B.* 
                , count(distinct A.pid_trt) AS n_trt  /* ADDED TO FACILITATE DATASET PROCESSING (SEE 0 ABOVE) */
                , count(distinct A.pid_con) AS n_con  /* ADDED TO FACILITATE DATASET PROCESSING (SEE 0 ABOVE) */
        FROM      compare  AS A
        %if %length(&BYVAR)=0 %then %do ;
                  , win2 (keep=theta_kl0 theta_kA win_trt win_con winratio)          AS B
        %end ;
        %if %length(&BYVAR)>0 %then %do ;
            LEFT JOIN win2 (keep=&BYVAR theta_kl0 theta_kA win_trt win_con winratio) AS B
            ON       &SQLON 
            GROUP BY &BYSQL 
        %end ;
        ;
    quit ;

    /* IS STANDARD ERROR UNDER THE NULL OR ALTERNATIVE */
    %if %length(&SE)>0 and %index(&SE,ALT)  %then %let theta = theta_kA ;
    %else                                         %let theta = theta_kl0 ;;


    /*****************************************************************************************************/
    /******* Construct kernel functions for variance estimate (SEE DONG 2016. SECTION 3)               ***/
    /*****************************************************************************************************/
    /* KERNEL FUNCTION K
    /* Kij = 1, if the ith patient in the treatment group won over the jth patient in the control group  */ 
    /*     = 0, if the ith patient in the treatment group lost to the jth patient in the control group   */ 
    /*     = 0, if the ith patient in the treatment group tied with the jth patient in the control group */
    data K ;
        set KLBASE ;   
        if indexc(wincat,'t')  then k=1 ;  /*** if Treatment won */
        else                        k=0 ;
    run ;
    title2 'Kernal function K' ;
    title3 'Dataset: K (where=((pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST))) )' ;
    proc print data=K (where=((pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST))) ) noobs ;
    run &DEBUG ;
    
    /* KERNEL FUNCTION L DERIVED IN THE SAME WAY */
    data L ;
        set KLBASE ;
        if indexc(wincat,'c')  then k=1 ;  /*** if Control won */
        else                        k=0 ;
    run ;
    title2 'Kernal function L' ;
    title3 'Dataset: L (where=((pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST))) )' ;
    proc print data=L (where=((pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST))) ) noobs ;
    run &DEBUG ;
    

    /*****************************************************************************************************/
    /******* Calculate variances - MACRO                                                               ***/
    /*****************************************************************************************************/
    /** THIS IS PERFORMING THE TRIPLE SUM IN EQs 8b and 8c (p432) of Dong 2016).
    It is called 4 times to get sig2_10 (8b) and sig2_20 (8c) for the treatment and control groups in turn
    It is mostly straightforward to follow, the most obscure part is ....
          DoubleSum = (k - THETA_KL0) * (sum_k1 - k - (&N2-1) * THETA_KL0)
      which is the sum over j and j prime ne j of (kij-theta)(kijprime-theta) (8b)
      where sum_k1 is the total number of wins for participant i
    
    Very roughly, 
          if Kij = 1 then there are (sum_k1 - 1) instances where (kijprime - theta)=(1-theta)
          if Kij = 0 then there are (sum_k1)     instances where (kijprime - theta)=(1-theta)
          otherwise,                                             (kijprime - theta)=  -theta
      Since there are N2-1 instances of j ne jprime
        summing the expression over j ne jprime = (sum_k1 - k - (&N2-1) * THETA_KL0)
      Hence, DoubleSum = (k - THETA_KL0) * (sum_k1 - k - (&N2-1) * THETA_KL0) !! **/

    /* [HJR 20240916]: THETA_KL0 CHANGED TO &THETA BELOW TO HANDLE SE CALCULATION UNDER THE ALTERNATIVE */
    /****************************************************************************************************/
    %macro sig(DSIN=, PID1=, PID2=, N1=, N2=, BYVAR=, BYSQL=, DSOUT=, VAROUT= , THETAHAT= , DEBUG=CANCEL);
        
        /*** DSIN   : input kernel function dataset (K/L)                  ***/
        /*** PID1   : patient ID in the 1st group (pid_trt/pid_con)        ***/
        /*** PID2   : patient ID in the 2nd group (pid_trt/pid_con)        ***/
        /*** N1     : number of patients in the 1st group (n_trt/n_con)    ***/
        /*** N2     : number of patients in the 2nd group (n_trt/n_con)    ***/
        /*** DSOUT  : output dataset                                       ***/
        /*** VAROUT : variable name of variance component being calculated ***/

        title2 "Calculating variances: Kernal=&DSIN Parameter=&VAROUT";
        title3 "Dataset: &DSIN" '(where=(pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST)) )' ;
        proc print data=&DSIN (where=(pid_trt in (&PRNTLST) AND pid_con in (&PRNTLST)) ) noobs ;
        run &DEBUG ;
 
        /* FOR PROCESS SPEED, SOURCE CODE PROC MEANS FOLLOWED BY SQL JOIN IS ALL DONE IN ONE JOIN */
        proc sql ;
            CREATE TABLE DSINplus AS
            SELECT    * , sum(K) AS sum_k1 /* sum_k1 = sum of k=1 (or L=1 depending on &DSIN) */
                                           /* THIS IS A PROC MEANS IN THE DONG MACRO BUT DOING IT IN SQL 
                                                      AVOIDS THE NEED FOR DATASET SORTING */
            FROM      &DSIN AS a
            GROUP BY  %if %length(&BYSQL)>0 %then &BYSQL ,;  a.&PID1
            ORDER BY  %if %length(&BYSQL)>0 %then &BYSQL ,;  a.&PID1, a.&PID2
            ;
        quit ;
        title3 "Dataset: DSINplus" '(where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) )' ;
        proc print data=DSINplus (where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) ) noobs ;
        run &DEBUG ;
        
        data DoubleSum (index=(IDX=(&BYVAR &N1 &N2))) ; /* INDEXING AVOIDS THE NEED TO SORT */
            set DSINplus ;
            /* THIS IS THE SUM OVER j AND jprime IN EQUATION 8B/8C */
            DoubleSum  = (k - &THETAHAT) * (sum_k1 - k - (&N2-1) * &THETAHAT) ; 
        run ;
        title3 "Dataset: DoubleSum" '(where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) )' ;
        proc print data=DoubleSum (where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) ) noobs ;
        run &DEBUG ;
       
        /* THIS IS NOW THE OUTER SUM (SUM OVER I) OF 8b/8c */
        proc means data=DoubleSum noprint ;
            by &BYVAR &N1 &N2 ; /* N1 and N2 ARE IN BY STATEMENT SO THEY ARE RETAINED IN DATASET FOR FINAL CALCULATION */
            var DoubleSum ;
            output out=&DSOUT sum=TripleSum ;
        run ;
       
        data &DSOUT(index=(IDX=(&BYVAR n_trt n_con))) ; /* n_trt AND n_con INCLUDED IN INDEX TO SIMPLIFY CODE TO STILL
                                                           WORK WHEN NO &BYVAR - THEY ARE SAME CONSTANTS WITHIN DATASET */
            set &DSOUT(drop=_TYPE_ _FREQ_);      
            if &N2>0 then /* ADDED IN CASE OF SMALL STRATUM SIZES */
                  &VAROUT = TripleSum  * &N1 * &N2 / (&N2 - 1) ;
        run ;
        title3 "Dataset: &DSOUT" ;
        proc print data=&DSOUT label noobs ;
        run &DEBUG ;

    %mend sig ;

    /* Sigma square hat - Result for 8b */
    %sig(dsin=K ,pid1=pid_trt ,pid2=pid_con ,n1=N_TRT ,n2=N_CON ,byvar=&BYVAR ,bysql=&BYSQL ,dsout=sig2_trt1 ,varout=sig2_trt1 ,thetahat=&THETA ,debug=&DEBUG);
    /* Sigma square hat - Result for 8c */
    %sig(dsin=K ,pid1=pid_con ,pid2=pid_trt ,n1=N_CON ,n2=N_TRT ,byvar=&BYVAR ,bysql=&BYSQL ,dsout=sig2_trt2 ,varout=sig2_trt2 ,thetahat=&THETA ,debug=&DEBUG);
    
    /* Sigma square hat - 8b for control */
    %sig(dsin=L ,pid1=pid_con ,pid2=pid_trt ,n1=N_CON ,n2=N_TRT ,byvar=&BYVAR ,bysql=&BYSQL ,dsout=sig2_con1 ,varout=sig2_con1 ,thetahat=&THETA ,debug=&DEBUG);
    /* Sigma square hat - 8c for control */
    %sig(dsin=L ,pid1=pid_trt ,pid2=pid_con ,n1=N_TRT ,n2=N_CON ,byvar=&BYVAR ,bysql=&BYSQL ,dsout=sig2_con2 ,varout=sig2_con2 ,thetahat=&THETA ,debug=&DEBUG);
    

    /*****************************************************************************************************/
    /******* Calculate covariance - MACRO                                                              ***/
    /******* Follows the same principle as above for equations represented in 11a and 11b              ***/
    /*****************************************************************************************************/
    %macro sig_cov(DSIN1=, DSIN2=, PID1=, PID2=, N1=, N2=, BYVAR= , BYSQL= , SQLON= , DSOUT=, VAROUT= , THETAHAT= , DEBUG=CANCEL );
        
        /*** DSIN1  : dataset of kernel function 1 (K/L)                                               ***/
        /*** DSIN2  : dataset of kernel function 2 (K/L)                                               ***/
        /*** PID1   : patient ID in the 1st group (pid_trt/pid_con)                                    ***/
        /*** PID2   : patient ID in the 2nd group (pid_trt/pid_con)                                    ***/
        /*** N1     : number of patients in the 1st group (&n_trt/&n_con)                              ***/
        /*** N2     : number of patients in the 2nd group (&n_trt/&n_con)                              ***/
        /*** DSOUT  : output dataset                                                                   ***/
        /*** VAROUT : variable name of covariance component being calculated                           ***/
        /*** THETA  : variable to use for theta_KL0 in SE calculation (theta_kl0 under NULL, winratio  ***/
        /***            under alternative)                                                             ***/

        title2 "Calculating covariances: Parameter=&VAROUT" ;
        /* BRING sum_K1 (from DSIN1) and K (aka L) from DSIN2 INTO DSIN1 */
        /* DONE IN A SQL JOIN FOR PROCESSING SPEED WHEN DATASETS ARE BIG */    
        proc sql ; 
            /* THIS IS A PROC MEANS IN THE DONG MACRO BUT DOING IT IN SQL AVOIDS THE NEED FOR DATASET SORTING */
            CREATE TABLE sum_k1 AS
            SELECT    %if %length(&BYSQL)>0 %then &BYSQL ,; &PID1 
                          , sum(K) AS sum_k1 /* sum_k1 = sum of k=1 (or L=1 depending on &DSIN) */
            FROM      &DSIN2 AS a   
            GROUP BY  %if %length(&BYSQL)>0 %then &BYSQL ,; &PID1 
            ;  /* END CREATE TABLE */
            CREATE TABLE DSIN1plus  AS
            SELECT    a.*, b.sum_k1
            FROM      &DSIN1  AS a 
            LEFT JOIN sum_k1  AS b 
               ON     %if %length(&SQLON)>0 %then &SQLON AND ;  a.&PID1 = b.&PID1
            ORDER BY %if %length(&BYSQL)>0 %then &BYSQL ,; &PID1 , &PID2
            ;  /* END CREATE TABLE */
            CREATE TABLE DSIN1plusplus AS
            SELECT    a.*, b.k AS L
            FROM      DSIN1plus AS a 
            LEFT JOIN &DSIN2     AS b 
               ON     %if %length(&SQLON)>0 %then &SQLON AND ;  a.&PID1 = b.&PID1 AND a.&PID2 = b.&PID2
            ORDER BY %if %length(&BYSQL)>0 %then &BYSQL ,; &PID1, &PID2 
            ;  /* END CREATE TABLE */
        quit ;
        title3 "Dataset: DSIN1plus(where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) )" ;
        proc print data=DSIN1plus(where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) ) noobs ;
        run ;
        title3 "Dataset: DSIN1plusplus (where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) )" ;
        proc print data=DSIN1plusplus(where=(&PID1 in (&PRNTLST) AND &PID2 in (&PRNTLST)) ) noobs ;
        run ;
        
        data DoubleSum (index=(IDX=(&BYVAR &N1 &N2))) ; /* INDEXING DATASET SO WE DO NOT NEED SORT FOR NEXT DATA STEP */
            set DSIN1plusplus ;
            DoubleSum = (k-&THETAHAT) * (sum_k1 - L - (&N2-1) * &THETAHAT) ;
        run ;
        
        proc means data=DoubleSum noprint ;
            by  &BYVAR &N1 &N2 ;
            var DoubleSum ;
            output out=TripleSum sum=TripleSum ;
        run ;
        
        data &DSOUT(index=(IDX=(&BYVAR n_trt n_con))) ; /* n_trt AND n_con ARE INCLUDED IN THE INDEX TO FACILITATE 
                                                                  MERGING WHEN NO &BYVAR */
            set TripleSum(drop=_TYPE_ _FREQ_) ;
            if &N2>0 then /* ADDED IN CASE OF SMALL STRATUM SIZES */
                  &VAROUT = TripleSum * &N1 * &N2 / (&N2 - 1) ;
        run ;
        title3 "Dataset: &DSOUT" ;
        proc print data=&DSOUT label noobs ;
        run &DEBUG ;
    %mend sig_cov ;
    /* Sigma hat - tc10 */
    %sig_cov(dsin1=K ,dsin2=L ,pid1=pid_trt ,pid2=pid_con ,n1=N_TRT ,n2=N_CON ,byvar=&BYVAR ,bysql=&BYSQL ,sqlon=&SQLON ,dsout=sig_trt_con1 ,varout=sig_trt_con1 ,thetahat=&THETA ,debug=&DEBUG);
    /* Sigma hat - tc20 */
    %sig_cov(dsin1=K ,dsin2=L ,pid1=pid_con ,pid2=pid_trt ,n1=N_CON ,n2=N_TRT ,byvar=&BYVAR ,bysql=&BYSQL ,sqlon=&SQLON ,dsout=sig_trt_con2 ,varout=sig_trt_con2 ,thetahat=&THETA ,debug=&DEBUG);
    

    /*****************************************************************************************************/
    /******* VARIANCE AND COVARIANCE ESTIMATION PUTS ALL OF THIS TOGETHER                                */
    /*****************************************************************************************************/
    data SIG ;
        merge sig2_trt1
              sig2_trt2
              sig2_con1
              sig2_con2
              sig_trt_con1
              sig_trt_con2
              ;
        by &BYVAR n_trt n_con ; /* FACILITATES THE INDEXING - INDEXING STATEMENT FORMAT IS SIMPLIFIED IF ONLY ONE BYVAR */
        mergevar = 1 ;          /* FACILITATES MERGING WITH WIN RATIO RESULTS WHEN NO BYVAR */

        if n_con > 0 and n_trt>0 then do ; /* AVOID NOTES FOR 0 DENOMS. AND ENSURES MISSING VALUES FOR ALL PARAMETERS */
            /* Sigma square hat - t0 */
            sig2_trt = (sig2_trt1/N_TRT) + (sig2_trt2/N_CON) ;
    
            /* Sigma square hat - c0 */
            sig2_con = (sig2_con1/N_CON) + (sig2_con2/N_TRT) ;
    
            /* Sigma hat - tc0 */
            sig_trt_con = (sig_trt_con1/N_TRT) + (sig_trt_con2/N_CON) ;

            _theta = "&SE" ;
        end ;
    
        label sig2_trt    = 'Sigma square hat - t0'
              sig2_con    = 'Sigma square hat - c0'
              sig_trt_con = 'Sigma hat - tc0'
              _theta      = "Theta used in variance estimation"
              ;
    
        keep &BYVAR mergevar n_trt n_con sig2_trt sig2_con sig_trt_con _theta ;
    run ;
    proc sort data=SIG ;
        by &BYVAR mergevar ;
    run ;
    title2 'Dataset: SIG - Variance and covariance' ;
    proc print data=SIG label noobs ;
        var &BYVAR n_trt n_con sig2_trt sig2_con sig_trt_con _theta ;
    run &DEBUG ;
    
    /* FINAL OUTPUT DATA SET */
    data &ODSN;
        merge WIN2 SIG;
        by &BYVAR mergevar ;
      
        if winratio>0 then do ;
            /* SIGMA SQUARE HAT OF LOG(WIN RATIO) - MOST CLEARLY DEFINED IN DONG 2018 (bottom of page 782) */
            sig2_log_wr = (sig2_trt + sig2_con - 2*sig_trt_con)/((WIN_TRT + WIN_CON)*(WIN_TRT + WIN_CON)/4) ;
        
            /* 100*(1-alpha)% CI of the win ratio */
            WR_L = exp(log(WINRATIO) - probit(1-&ALPHA/2)*sqrt(sig2_log_wr)) ;
            WR_U = exp(log(WINRATIO) + probit(1-&ALPHA/2)*sqrt(sig2_log_wr)) ;
        end ;
      
        alpha = &ALPHA ;
      
        label sig2_log_wr = "Sigma square hat of log(win ratio) under &SE"
              WR_L        = "Lower limit of 100*(1-alpha)% CI of the Win ratio"
              WR_U        = "Upper limit of 100*(1-alpha)% CI of the Win ratio"
              WinRatio    = 'Win ratio'
              n_trt       = "Treatment group sample size"
              n_con       = "Control group sample size"
              alpha       = 'Alpha'
              ;
  
        rename total=TotalPairs ;
        drop mergevar ;
    run ;
    title2 "Dataset: %upcase(&ODSN) - Win ratio and its 100*(1-alpha)% CI" ;
    proc print data=&ODSN label noobs ;
    run &DEBUG ;

%mend winratio ;
