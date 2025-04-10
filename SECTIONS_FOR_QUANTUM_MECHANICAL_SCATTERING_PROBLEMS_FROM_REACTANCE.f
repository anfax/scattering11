! ACRLDCS. PROGRAM FOR CALCULATING DIFFERENTIAL AND INTEGRAL CROSS        ACRL0000
! 1   SECTIONS FOR QUANTUM MECHANICAL SCATTERING PROBLEMS FROM REACTANCE  ACRL0000
! 2   OR TRANSITION MATRICES.  BRANDT, M.A., TRUHLAR, D.G., SMITH, R.L.   ACRL0000
! REF. IN COMP. PHYS. COMMUN. 5 (1973) 456                                ACRL0000
! SEE ERRATUM COMP. PHYS. COMMUN. 7(1974)177                              ACRL0000
! DCS,T16,CM61000.21206065                                                ACRL0001
! FUN(S,,,,,,40000)                                                       ACRL0002
! LGO.                                                                    ACRL0003
! $                                                                       ACRL0004
      PROGRAM DRIVER (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE9)      ACRL0005
C     THIS IS THE DRIVER FOR THE TEST RUN.                              ACRL0006
C                                                                       ACRL0007
C         THE COMMON VARIABLES NCALL,NFACT,IREAD,IWRITE, AND ITAPE ARE  ACRL0008
C    INITIALIZED HERE.  NCALL IS USED ONLY IN THE COMMON BLOCK /FACT/   ACRL0009
C    AND MAY BE GIVEN ANY VALUE OTHER THAN -1867.  THE ARRAY FL IS USED ACRL0010
C    TO STORE THE LOGARITHMS OF THE FACTORIALS NEEDED BY SUBPROGRAM Z.  ACRL0011
C    NFACT SHOULD BE SET TO A VALUE GREATER THAN OR EQUAL TO 4*(MAXIMUM ACRL0012
C    OF ANY ANGULAR MOMENTUM OCCURRING IN THIS RUN) + 2.  IF NFACT IS   ACRL0013
C    SET LARGER THAN 322, THE DIMENSION OF FL MUST BE CHANGED.  THE PRO-ACRL0014
C    GRAM HAS NOT BEEN CHECKED OUT FOR VALUES OF NFACT GREATER THAN 322.ACRL0015
C                                                                       ACRL0016
      COMMON/FACT/FL(322),NCALL,NFACT                                   ACRL0017
C                                                                       ACRL0018
C         THE COMMON BLOCK /IOSEND/ IS USED TO SPECIFY THE TAPE NUMBERS ACRL0019
C    USED BY THE COMPUTER STATION FOR INPUT/OUTPUT.  IREAD IS THE TAPE  ACRL0020
C    NUMBER OF THE DEVICE FROM WHICH THE DATA IS READ.  IWRITE IS THE   ACRL0021
C    TAPE NUMBER OF THE DEVICE TO WHICH THE ANSWERS ARE WRITTEN.  ITAPE ACRL0022
C    IS THE TAPE NUMBER OF THE EXTERNAL DEVICE USED TO STORE THE        ACRL0023
C    T MATRIX IF NECESSARY.  THE COMMON BLOCK /IOSEND/ IS USED BY SUB-  ACRL0024
C    PROGRAMS CROSS, INOUT, MXLNEQ, AND F6J.  (SEE CARD NO.41)          ACRL0025
C                                                                       ACRL0026
      COMMON /IOSEND/ IREAD,IWRITE,ITAPE                                ACRL0027
C                                                                       ACRL0028
C         THE VARIABLE DIMENSIONS USED BY SUBPROGRAM CROSS ARE SET HERE.ACRL0029
C                                                                       ACRL0030
      DATA LA,LMAXP,JATMAX,JAXMAX,INMX,JMAXP,JTAT,NMAX,KMAX,MMAX/181,11,ACRL0031
     13,2,183,8,6,14,5,4/                                               ACRL0032
C                                                                       ACRL0033
      DIMENSION PSTOR(181,11),LEVEL(3),LEVELP(3),JJ(3),JJP(3),          ACRL0034
     1          B(11),SIGMA(183,2),INDEX(183),                          ACRL0035
     2          PROBL(8),PROBLO(8),OPAC(8),OPACO(8),                    ACRL0036
     3          JINDEX(14,6),K1MIN(6),L3MIN(6),NSAVE(6),KSTOR(14),      ACRL0037
     4          R(14,14),G(14,14),H(14,14),TEMP(14),                    ACRL0038
     5          TREAL(5,4,6),TIMAG(5,4,6)                               ACRL0039
C                                                                       ACRL0040
      DATA NCALL,NFACT,IREAD,IWRITE,ITAPE/6,70,5,6,9/                   ACRL0041
C                                                                       ACRL0042
      CALL CROSS (LA,LMAXP,JATMAX,JAXMAX,INMX,JMAXP,JTAT,NMAX,KMAX,MMAX,ACRL0043
     1PSTOR,LEVEL,LEVELP,JJ,JJP,B,SIGMA,INDEX,PROBL,PROBLO,OPAC,OPACO,  ACRL0044
     2JINDEX,K1MIN,L3MIN,NSAVE,KSTOR,R,G,H,TEMP,TREAL,TIMAG)            ACRL0045
      STOP                                                              ACRL0046
      END                                                               ACRL0047
      SUBROUTINE CROSS (LA,LMAXP,JATMAX,JAXMAX,INMX,JMAXP,JTAT,NMAX,    ACRL0048
     1KMAX,MMAX,PSTOR,LEVEL,LEVELP,JJ,JJP,B,SIGMA,INDEX,PROBL,PROBLO,   ACRL0049
     2OPAC,OPACO,JINDEX,K1MIN,L3MIN,NSAVE,KSTOR,R,G,H,TEMP,TREAL,TIMAG) ACRL0050
C                                                                       ACRL0051
C         SUBPROGRAM CROSS CALCULATES THE INTEGRAL AND DIFFERENTIAL     ACRL0052
C    CROSS SECTIONS FROM THE SCATTERING T MATRIX USING FORMULAS GIVEN   ACRL0053
C    IN J. BLATT AND L. BIEDENHARN, REV. MOD. PHYS. 24, 258 (1952) AND  ACRL0054
C    CORRECTED BY R. HUBY, PROC. PHYS. SOC. (LONDON) 67A, 1103 (1954).  ACRL0055
C         THE INITIAL CHANNEL ANGULAR MOMENTUM IS J AND THE FINAL       ACRL0056
C    CHANNEL ANGULAR MOMENTUM IS JP.  J1, WHEN USED AS A VARI-          ACRL0057
C    ABLE, IS THE TOTAL ANGULAR MOMENTUM, AND WHEN J1 IS USED AS AN     ACRL0058
C    ARRAY SUBSCRIPT, (J1-1) IS THE TOTAL ANGULAR MOMENTUM.  (JTOT-1) ISACRL0059
C    THE MAXIMUM TOTAL ANGULAR MOMENTUM IN A DATA SET.  (JMIN - 1) IS   ACRL0060
C    THE MINIMUM TOTAL ANGULAR MOMENTUM IN A DATA SET.  A DATA SET CON- ACRL0061
C    TAINS A SET OF R OR T J-BLOCKS(SEE CARD NO.104) AND THE REMAINING  ACRL0062
C    DATA NEEDED TO FIND A CROSS SECTION (ALL THE DATA READ IN THE PARTSACRL0063
C    OF CROSS FROM CARD NO.315 TO CARD NO.541 AND FROM CARD NO.912 TO   ACRL0064
C    CARD NO.920).  THE DESCRIPTION OF THE COMMON VARIABLES STARTS WITH ACRL0065
C    CARD NO.94 AND THE DESCRIPTION OF THE VARIABLES IN THE PARAMETER   ACRL0066
C    LIST STARTS WITH CARD NO.171.  OTHER VARIABLES FOR WHICH THE USER  ACRL0067
C    MUST SUPPLY THE VALUES ARE DESCRIBED IN THE PART FROM CARD NO.267  ACRL0068
C    TO CARD NO.465.                                                    ACRL0069
C         BECAUSE DO LOOPS WERE USED TO PERFORM THE BL COEFFICIENTS SUMSACRL0070
C    OVER RELATIVE ORBITAL ANGULAR MOMENTUM, L, CROSS HAS CERTAIN       ACRL0071
C    REQUIREMENTS.  CASE 1 -- THE CHANNEL PARITY, SEE TEXT, IS NOT CON- ACRL0072
C                       SERVED.  EACH J TO JP CROSS SECTION TO BE FOUND ACRL0073
C                       FROM A DATA SET REQUIRES A 2*J+1 BY 2*JP+1 SUB- ACRL0074
C                       MATRIX OF THE T MATRIX FOR EACH J1.             ACRL0075
C                   CASE 2 -- THE CHANNEL PARITY IS CONSERVED.          ACRL0076
C                       ONLY A J+1 BY JP+1 SUBMATRIX IS NEEDED FOR EACH ACRL0077
C                       J1.  J AND JP MUST HAVE THE SAME PARITY IN      ACRL0078
C                       THE DATA SET(SEE DESCRIPTION OF THE VARIABLE    ACRL0079
C                       JPAR1 STARTING AT CARD NO.319).                 ACRL0080
C    IN EITHER CASE, THE SUBMATRIX ROWS AND COLUMNS MUST BE IN THE STAN-ACRL0081
C    DARD ORDER OF INCREASING L.  AT LOW J1, THE TRIANGLE RULE MAY NOT  ACRL0082
C    ALLOW 2*J+1 OR 2*JP+1 VALUES OF L.  IN SUCH CASES THE SUBMATRIX    ACRL0083
C    REQUIRED WILL BE SMALLER THAN STATED ABOVE BUT THE ROWS AND COLUMNSACRL0084
C    MUST STILL BE IN THE STANDARD ORDER REQUIRED BY CROSS.             ACRL0085
C         IF THE DATA IS NOT USABLE IN ITS PRESENT FORM BECAUSE THE     ACRL0086
C    DATA ROWS ARE NOT IN THE STANDARD ORDER, SUBPROGRAM TRNSFR ALLOWS  ACRL0087
C    THE USER THE OPPORTUNITY TO REORDER THE DATA ROWS.                 ACRL0088
C                                                                       ACRL0089
      LOGICAL LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LT    ACRL0090
      COMMON /IOSEND/ IREAD,IWRITE,ITAPE                                ACRL0091
      COMMON/LCOM/LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LTACRL0092
C                                                                       ACRL0093
C                        COMMON VARIABLES                               ACRL0094
C                                                                       ACRL0095
C       THE EXPLANATION OF THE /IOSEND/ VARIABLES BEGINS ON CARD NO.19. ACRL0096
C                                                                       ACRL0097
C       LPR1 DETERMINES IF THE T MATRIX HAS ALREADY BEEN STORED FOR THE ACRL0098
C         NEXT INTERMEDIATE CROSS SECTION CALCULATION. (SEE CARD NO.179)ACRL0099
C                                                                       ACRL0100
C       LPR2 MAY BE SET .TRUE.(.FALSE.) IF MORE THAN(ONLY) ONE DATA SET ACRL0101
C        WILL BE READ IN THIS CALL TO CROSS (SEE CARDS NOS.912 TO 919)  ACRL0102
C                                                                       ACRL0103
C       LPR3 IS .TRUE. ONLY IF A T J-BLOCK IS TO BE READ. A T J-BLOCK ISACRL0104
C    THE BLOCK OF TOTAL ANG. MOMENTUM J IN THE BLOCK DIAGONAL T MATRIX. ACRL0105
C                                                                       ACRL0106
C       LOREAD SHOULD BE SET .TRUE. FOR J-BLOCKS WITH NOPT.LT.2 (SEE    ACRL0107
C          CARD NO.450 FOR NOPT DESCRIPTION). LOREAD MUST BE SET .FALSE.ACRL0108
C          FOR THE FIRST J-BLOCK WITH NOPT.GE.2.  THEREAFTER, LOREAD    ACRL0109
C          SHOULD BE SET .FALSE. ONLY WHEN A NEW REORDERING SEQUENCE IS ACRL0110
C          REQUIRED BY SUBPROGRAM TRNSFR (SEE CARDS NOS.141,142,1026,   ACRL0111
C          1275) AND SHOULD BE SET .TRUE. OTHERWISE. IF THE SIZE OF THE ACRL0112
C          J-BLOCKS CHANGES, THE ROW ORDER OF THE J-BLOCKS ALSO CHANGES.ACRL0113
C                                                                       ACRL0114
C       LSTOR IS .TRUE. IF MORE THAN ONE CROSS SECTION IS TO BE FOUND   ACRL0115
C          AND ENOUGH STORAGE IS AVAILABLE FOR PSTOR TO BE USED.        ACRL0116
C           LSTOR DETERMINES IF THE ARRAY PSTOR IS TO BE USED.          ACRL0117
C                                                                       ACRL0118
C       LEGEND IS USED TO SAVE COMPUTER TIME IN CALCULATING THE         ACRL0119
C          LEGENDRE POLYNOMIALS AND IS EXPLAINED IN FUNCTION P STARTING ACRL0120
C                            WITH CARD NO.1537.                         ACRL0121
C                                                                       ACRL0122
C       LZ DETERMINES PRINTING OPTIONS IN SUBPROGRAM INOUT AND IS USED  ACRL0123
C    IN EVALUATING EQ.(4.7). ALL EQ. NOS. REFER TO BLATT AND BIEDENHARN.ACRL0124
C                                                                       ACRL0125
C       LCHECK IS .TRUE. ONLY IF THE J-BLOCK READ IS TO BE SEARCHED BY  ACRL0126
C                          SUBPROGRAM INOUT.                            ACRL0127
C                                                                       ACRL0128
C       LABORT IS AN OUTPUT VARIABLE FOR INOUT WHICH HAS THE VALUE      ACRL0129
C    .FALSE. THE FIRST TIME SUBPROGRAM INOUT IS CALLED WHILE A DATA SET ACRL0130
C    IS BEING READ.  IF THE J-BLOCK IS FOUND NOT TO BE SYMMETRIC, LABORTACRL0131
C    IS SET .TRUE. AND ALL THE INCORRECT ENTRIES ARE PRINTED.  THE RE-  ACRL0132
C    MAINING DATA IN THE PRESENT DATA SET ARE READ AND THE J-BLOCKS FOR ACRL0133
C    WHICH LCHECK IS .TRUE. ARE CHECKED FOR SYMMETRY. (SEE CARD NO.1053)ACRL0134
C                                                                       ACRL0135
C       LT IS EXPLAINED IN SUBPROGRAM PARITY BEGINNING AT CARD NO.1514. ACRL0136
C                                                                       ACRL0137
C                  PROVIDED EXTERNALS USED BY CROSS                     ACRL0138
C                                                                       ACRL0139
C   SUBROUTINE INOUT TAKES CARE OF THE INPUT/OUTPUT REQUESTS OF CROSS.  ACRL0140
C   SUBROUTINE TRNSFR REORDERS THE R OR T J-BLOCK ROWS AND COLUMNS      ACRL0141
C                        TO PROVIDE THE ORDER USED BY CROSS.            ACRL0142
C   SUBROUTINE MXMPLY DOES MATRIX MULTIPLICATION.                       ACRL0143
C   SUBROUTINE MXLNEQ DOES MATRIX INVERSION.                            ACRL0144
C   SUBROUTINE OPACIT CALCULATES THE OPACITY FOR EACH L.                ACRL0145
C   SUBROUTINE PARITY (1) DETERMINES WHETHER A GIVEN TERM IN EQ. (4.7)  ACRL0146
C                         OF BLATT AND BIEDENHARN SATISFIES THEIR EQS.  ACRL0147
C                         (4.8) AND (4.9), AND                          ACRL0148
C                     (2) DETERMINES WHETHER A GIVEN CONTRIBUTION TO    ACRL0149
C                         OPAC AND PROBL IS AN EVEN OR ODD PARITY ONE.  ACRL0150
C   FUNCTION P FINDS THE LEGENDRE POLYNOMIALS.                          ACRL0151
C   FUNCTION Z CALCULATES THE Z COEFFICIENT OF J. BLATT, L. BIEDENHARN  ACRL0152
C                 AND M. ROSE, REV.MOD.PHYS. 24, 249 (1952) WITH THE    ACRL0153
C                 PHASE FACTOR REMOVED AS RECOMMENDED BY R. HUBY.       ACRL0154
C                                                                       ACRL0155
C         ALL ARRAYS EXCEPT LABEL ARE DIMENSIONED IN THE DRIVER AND THE ACRL0156
C    VARIABLE DIMENSIONS ARE SET WITH A DATA STATEMENT AT CARD NO.31.   ACRL0157
C    IF THE PARAMETER LIST IS TOO LONG, SUBPROGRAM CROSS CAN BE CHANGED ACRL0158
C    TO A MAIN PROGRAM BY DROPPING CARDS NOS. 43 - 50, 91, AND REPLACINGACRL0159
C    CARDS NOS. 163 - 168 WITH CARDS NOS. 34 - 39 AND CARD NO. 925 WITH ACRL0160
C    CARD NO. 46.                                                       ACRL0161
C                                                                       ACRL0162
      DIMENSION PSTOR(LA,LMAXP),LEVEL(JATMAX),LEVELP(JATMAX),JJ(JATMAX),ACRL0163
     1          JJP(JATMAX),B(LMAXP),SIGMA(INMX,JAXMAX),INDEX(INMX),    ACRL0164
     2          PROBL(JMAXP),PROBLO(JMAXP),OPAC(JMAXP),OPACO(JMAXP),    ACRL0165
     3          JINDEX(NMAX,JTAT),K1MIN(JTAT),L3MIN(JTAT),NSAVE(JTAT),  ACRL0166
     4          KSTOR(NMAX),R(NMAX,NMAX),G(NMAX,NMAX),H(NMAX,NMAX),     ACRL0167
     5          TEMP(NMAX),TREAL(KMAX,MMAX,JTAT),TIMAG(KMAX,MMAX,JTAT)  ACRL0168
     6         ,LABEL(8)                                                ACRL0169
C                                                                       ACRL0170
C                         PARAMETER LIST                                ACRL0171
C                                                                       ACRL0172
C  LA  -  A VARIABLE DIMENSION  -  IF(LSTOR) LA=181, IF(.NOT.LSTOR) LA=1ACRL0173
C                                                                       ACRL0174
C  LMAXP  -  A VARIABLE DIMENSION  -  2*MAXIMUM(JTOT - 1) + 1           ACRL0175
C                                                                       ACRL0176
C  JATMAX  -  A VARIABLE DIMENSION  -  MAX(JAPMAX), WHERE JAPMAX IS THE ACRL0177
C     NUMBER OF INTERMEDIATE CROSS SECTIONS TO BE FOUND FROM A DATA SET.ACRL0178
C            AN INTERMEDIATE CROSS SECTION IS ONE FOUND BY EQ.(4.5).    ACRL0179
C                                                                       ACRL0180
C  JAXMAX  -  A VARIABLE DIMENSION  -  MAX(JAPMAX + 1) FOR ANY DATA SET ACRL0181
C       WITH JDEG.NE.0(SEE CARD NO.366 FOR JDEG).  IF JDEG.EQ.0 FOR ALL ACRL0182
C          THE DATA SETS IN THIS CALL TO CROSS, JAXMAX=1.               ACRL0183
C                                                                       ACRL0184
C  INMX  -  A VARIABLE DIMENSION  -  INMX=183                           ACRL0185
C                                                                       ACRL0186
C  JMAXP  -  A VARIABLE DIMENSION  -  MAXIMUM(JTOT + J) CONSIDERED      ACRL0187
C                                                                       ACRL0188
C  JTAT  -  A VARIABLE DIMENSION  -  MAXIMUM (JTOT-JMIN) + 1  CONSIDEREDACRL0189
C                                                                       ACRL0190
C  NMAX  -  A VARIABLE DIMENSION  -  MAXIMUM NUMBER OF ROWS OF ANY      ACRL0191
C                                      R OR T J-BLOCK TO BE READ        ACRL0192
C                                                                       ACRL0193
C  KMAX  -  A VARIABLE DIMENSION  -  KMAX=NMAX IF THE T MATRIX IS NOT TOACRL0194
C          BE STORED ON TAPE ITAPE.  IF TAPE ITAPE IS TO BE USED,       ACRL0195
C     KMAX=MAX(JPSEM).  (SEE CARDS NOS. 411 AND 412 FOR JSEM AND JPSEM) ACRL0196
C                                                                       ACRL0197
C  MMAX  -  A VARIABLE DIMENSION  -  IF KMAX=NMAX, THEN MMAX=NMAX --    ACRL0198
C                                     OTHERWISE MMAX=MAX(JSEM)          ACRL0199
C                                                                       ACRL0200
C  PSTOR  -  A VARIABLY DIMENSIONED ARRAY  -  LEGENDRE POLYNOMIALS USED ACRL0201
C                              IN FINDING THE DIFFERENTIAL CROSS SECTIONACRL0202
C                                                                       ACRL0203
C  LEVEL  -  A VARIABLY DIMENSIONED ARRAY  -  THE JAPMAX VALUES OF      ACRL0204
C      (THE NON-ANGULAR QUANTUM NUMBER)*(THE SIGN OF THE CHANNEL PARITY)ACRL0205
C            FOR THE INITIAL STATES OF THE JAPMAX CALCULATIONS          ACRL0206
C THE NON-ANGULAR QUANTUM NUMBER REPRESENTS ALL THE QUANTUM NUMBERS OF AACRL0207
C      QUANTUM STATE EXCLUDING CHANNEL ANGULAR MOMENTUM, RELATIVE       ACRL0208
C        ORBITAL ANGULAR MOMENTUM, AND TOTAL ANGULAR MOMENTUM.          ACRL0209
C                                                                       ACRL0210
C  LEVELP  -  A VARIABLY DIMENSIONED ARRAY  -  THE SAME AS LEVEL BUT    ACRL0211
C                FOR THE FINAL STATES INSTEAD OF THE INITIAL STATES     ACRL0212
C                                                                       ACRL0213
C  JJ  -  A VARIABLY DIMENSIONED ARRAY  -  CHANNEL ANGULAR MOMENTA OF   ACRL0214
C                       THE INITIAL STATES FOR THE JAPMAX CALCULATIONS  ACRL0215
C                                                                       ACRL0216
C  JJP  -  A VARIABLY DIMENSIONED ARRAY  -  CHANNEL ANGULAR MOMENTA OF  ACRL0217
C                         THE FINAL STATES FOR THE JAPMAX CALCULATIONS  ACRL0218
C                                                                       ACRL0219
C  B  -  A VARIABLY DIMENSIONED ARRAY  -  BL COEFFICIENTS OF BLATT AND  ACRL0220
C                                       BIEDENHARN AS CORRECTED BY HUBY ACRL0221
C                                                                       ACRL0222
C  SIGMA  -  A VARIABLY DIMENSIONED ARRAY  -  DIFFERENTIAL CROSS SECTIONACRL0223
C                                                                       ACRL0224
C  INDEX  -  A VARIABLY DIMENSIONED ARRAY  -  INTEGERS 0-180            ACRL0225
C                                                                       ACRL0226
C  PROBL  -  A VARIABLY DIMENSIONED ARRAY  -  TRANSITION PROBABILITIES, ACRL0227
C         STORED IN ORDER OF INCREASING L, FROM THE INITIAL TO FINAL    ACRL0228
C         STATE--EVEN PARITY CONTRIBUTION                               ACRL0229
C                                                                       ACRL0230
C  PROBLO  -  A VARIABLY DIMENSIONED ARRAY  -  ODD PARITY PART OF PROBL ACRL0231
C                                                                       ACRL0232
C  OPAC  -  A VARIABLY DIMENSIONED ARRAY  -  TRANSITION PROBABILITIES,  ACRL0233
C         STORED IN ORDER OF INCREASING L, FOR ALL TRANSITIONS OUT OF   ACRL0234
C         THE INITIAL STATE--EVEN PARITY CONTRIBUTION                   ACRL0235
C                                                                       ACRL0236
C  OPACO  -  A VARIABLY DIMENSIONED ARRAY  -  ODD PARITY PART OF OPAC   ACRL0237
C                                                                       ACRL0238
C  JINDEX  -  A VARIABLY DIMENSIONED ARRAY  -  QUANTUM STATE LABELS FOR ACRL0239
C                                        THE ROWS OF THE T J-BLOCK      ACRL0240
C                                                                       ACRL0241
C  K1MIN  -  A VARIABLY DIMENSIONED ARRAY  -  THE INITIAL STATE COLUMN  ACRL0242
C                                           INDEX FOR EACH J-BLOCK      ACRL0243
C                                                                       ACRL0244
C  L3MIN  -  A VARIABLY DIMENSIONED ARRAY  -  THE FINAL STATE ROW       ACRL0245
C                                           INDEX FOR EACH J-BLOCK      ACRL0246
C                                                                       ACRL0247
C  NSAVE  -  A VARIABLY DIMENSIONED ARRAY  -  T J-BLOCK SIZE FOR EACH J1ACRL0248
C                                                                       ACRL0249
C  KSTOR  -  A VARIABLY DIMENSIONED ARRAY  -  EXPLAINED IN TRNSFR ON    ACRL0250
C                                                CARD NO.1285.          ACRL0251
C                                                                       ACRL0252
C  R  -  A VARIABLY DIMENSIONED ARRAY  -  REACTANCE(R) MATRIX J-BLOCK ORACRL0253
C                                         A TEMPORARY STORAGE ARRAY     ACRL0254
C                                                                       ACRL0255
C  G  -  A VARIABLY DIMENSIONED ARRAY  -  A TEMPORARY STORAGE ARRAY     ACRL0256
C                                                                       ACRL0257
C  H  -  A VARIABLY DIMENSIONED ARRAY  -  A TEMPORARY STORAGE ARRAY     ACRL0258
C                                                                       ACRL0259
C  TEMP  -  A VARIABLY DIMENSIONED ARRAY  -  A TEMPORARY STORAGE ARRAY  ACRL0260
C                                                                       ACRL0261
C  TREAL  -  A VARIABLY DIMENSIONED ARRAY  -  REAL PART OF THE T MATRIX ACRL0262
C                                                                       ACRL0263
C  TIMAG  -  A VARIABLY DIMENSIONED ARRAY  -  IMAG. PART OF THE T MATRIXACRL0264
C                                                                       ACRL0265
C                                                                       ACRL0266
C ********************************************************************  ACRL0267
C * EPSD SHOULD BE SET BY THE USER TO 1.E-M WHERE M IS THE NUMBER OF *  ACRL0268
C *   SIGNIFICANT DIGITS FOR THE MACHINE IN USE. (SEE CARD NO.274)   *  ACRL0269
C * THE FUNCTION CPTIME SHOULD BE CHANGED TO A TIMING ROUTINE USED   *  ACRL0270
C *             AT THE USER COMPUTER STATION OR ELIMINATED.          *  ACRL0271
C *                (IT IS USED IN CARD NOS.277 AND 922)              *  ACRL0272
C ********************************************************************  ACRL0273
      DATA PI,EPSD,PI4,ID,IE,IG/3.1415926535897,1.E-14,12.5663706143588,ACRL0274
     1  1,1,1/                                                          ACRL0275
      LEGEND=.TRUE.                                                     ACRL0276
      TIME1=CPTIME(.1)                                                  ACRL0277
C       ANG2 IS THE FACTOR USED IN CHANGING FROM BOHR**2 TO ANGSTROM**2 ACRL0278
      ANG2=(0.529)**2                                                   ACRL0279
C              DATA RECORD NUMBER ONE                                   ACRL0280
C                                                                       ACRL0281
C       READ IN LABEL FOR THIS CALL TO CROSS                            ACRL0282
C                                                                       ACRL0283
C              DATA RECORD NUMBER TWO                                   ACRL0284
C                                                                       ACRL0285
C       READ IN THE INITIAL VALUES OF LPR2 AND LSTOR                    ACRL0286
C                                                                       ACRL0287
C              DATA RECORD NUMBER THREE                                 ACRL0288
C                                                                       ACRL0289
C      READ IN IAPMAX AND THE INITIAL VALUE OF JAPMAX.                  ACRL0290
C                                                                       ACRL0291
C         IAPMAX IS THE NUMBER OF DATA SETS IN THIS CALL TO CROSS AND   ACRL0292
C                      JAPMAX IS DEFINED ABOVE.                         ACRL0293
C                                                                       ACRL0294
      CALL INOUT (1,LABEL,TMP2,TMP3,TMP4,IAPMAX,JAPMAX,IC,8,IE,IG)      ACRL0295
C                                                                       ACRL0296
      LPR1=.FALSE.                                                      ACRL0297
      LABORT=.FALSE.                                                    ACRL0298
C           SET DELTA EQUAL TO AN ANGLE OF ONE DEGREE IN RADIANS.       ACRL0299
      DELTA=PI/180.                                                     ACRL0300
      FACTOR=0.6*DELTA*PI                                               ACRL0301
      IF (.NOT.LSTOR) GO TO 3                                           ACRL0302
C             STORE ALL THE LEGENDRE POLYNOMIALS NEEDED LATER           ACRL0303
      X=0.0                                                             ACRL0304
      DO 2 I=1,LA                                                       ACRL0305
      Y=COS(X)                                                          ACRL0306
      DO 1 K=1,LMAXP                                                    ACRL0307
      L=K-1                                                             ACRL0308
   1  PSTOR(I,K)=P(L,Y)                                                 ACRL0309
   2  X=X+DELTA                                                         ACRL0310
   3  CONTINUE                                                          ACRL0311
C                                                                       ACRL0312
C              DATA RECORD NUMBER FOUR                                  ACRL0313
C                                                                       ACRL0314
C         READ IN THE VALUE OF JPAR1 AND IEPS FOR THE FIRST DATA SET ANDACRL0315
C    THE JAPMAX VALUES OF THE INITIAL AND FINAL STATES OF THE           ACRL0316
C    INTERMEDIATE CROSS SECTIONS TO BE FOUND FROM THE FIRST DATA SET.   ACRL0317
C                                                                       ACRL0318
C        JPAR1=2(3) OR 5(6) IF L CHANGES BY TWO(ONE) IN GOING FROM ONE  ACRL0319
C                ROW OF THE T J-BLOCK TO THE NEXT ROW WHILE KEEPING ALL ACRL0320
C                   THE OTHER QUANTUM NUMBERS CONSTANT.  FOR EXAMPLE,   ACRL0321
C        JPAR1=2 IF THE CHANNEL PARITY IS CONSERVED AND ALL T MATRIX    ACRL0322
C           ELEMENTS HAVE THE SAME CHANNEL PARITY FOR THIS DATA SET.    ACRL0323
C        JPAR1=3 IF CHANNEL PARITY IS NOT CONSERVED.                    ACRL0324
C               ALSO TAPE ITAPE WILL NOT BE USED IF JPAR1.LE.3.         ACRL0325
C        JPAR1= 5, 6 HAVE THE SAME MEANING AS JPAR1-3 EXCEPT THAT TAPE  ACRL0326
C          ITAPE IS USED FOR STORING THE CURRENT SET OF T J-BLOCKS.     ACRL0327
C                                                                       ACRL0328
C       IEPS IS USED TO DETERMINE EPS, EPS = (10.)**(-IEPS), WHERE EPS  ACRL0329
C                  IS THE SYMMETRY TOLERANCE USED BY SUBPROGRAM INOUT.  ACRL0330
C                                                                       ACRL0331
C      LEVEL, LEVELP, JJ AND JJP ARE EXPLAINED IN THE PARAMETER LIST.   ACRL0332
C                                                                       ACRL0333
      CALL INOUT (2,LEVEL,JJ,LEVELP,JJP,JAPMAX,JPAR1,IEPS,JATMAX,JATMAX,ACRL0334
     1  IG)                                                             ACRL0335
      DO 68 IAP=1,IAPMAX                                                ACRL0336
C       JPAR IS THE INCREMENT USED BY THE DO LOOPS WHICH DO THE CROSS   ACRL0337
C                           SECTION SUMS OVER L.                        ACRL0338
      JPAR=1                                                            ACRL0339
      IF (JPAR1.EQ.2) JPAR=2                                            ACRL0340
      IF (JPAR1.EQ.5) JPAR=2                                            ACRL0341
      EPS=10.**(-IEPS)                                                  ACRL0342
      DO 66 JAP=1,JAPMAX                                                ACRL0343
      IF (JPAR1.GT.3) REWIND ITAPE                                      ACRL0344
      LZ = .FALSE.                                                      ACRL0345
      IF ((IAP.NE.1.OR.JAP.NE.1).AND.(JAP.NE.1.OR.(.NOT.LPR2)))LZ=.TRUE.ACRL0346
      JS=1                                                              ACRL0347
C              DATA RECORD NUMBERS FIVE OR TEN                          ACRL0348
C                                                                       ACRL0349
C      READ THE LABEL FOR THE NEXT CROSS SECTION TO BE CALCULATED       ACRL0350
C                                                                       ACRL0351
C              DATA RECORD NUMBERS SIX OR ELEVEN                        ACRL0352
C                                                                       ACRL0353
C      READ THE DATA NECESSARY FOR THE CALCULATION                      ACRL0354
C                                                                       ACRL0355
C        JTOT IS THE MAXIMUM TOTAL ANGULAR MOMENTUM+1 FOR ANY J-BLOCK   ACRL0356
C                   READ IN THIS DATA SET.  AFTER THE T MATRIX HAS BEEN ACRL0357
C                   STORED, JTOT IS THE MAXIMUM TOTAL ANGULAR MOMENTUM+1ACRL0358
C                   USED IN THIS CROSS SECTION CALCULATION.             ACRL0359
C                                                                       ACRL0360
C        JMIN IS THE MINIMUM TOTAL ANGULAR MOMENTUM+1 FOR ANY J-BLOCK   ACRL0361
C                   READ IN THIS DATA SET.  AFTER THE T MATRIX HAS BEEN ACRL0362
C                   STORED, JMIN IS THE MINIMUM TOTAL ANGULAR MOMENTUM+1ACRL0363
C                   USED IN THIS CROSS SECTION CALCULATION.             ACRL0364
C                                                                       ACRL0365
C        JDEG IS USED IN THE STORING OF SIGMA --                        ACRL0366
C                                                                       ACRL0367
C         JDEG=0 IF EQ.(3.16) OF BLATT AND BIEDENHARN IS NOT NEEDED FOR ACRL0368
C                 ANY OF THE JAPMAX CALCULATIONS.                       ACRL0369
C                                                                       ACRL0370
C         IF EQ.(3.16) IS TO BE USED TO ADD INTERMEDIATE CROSS SECTIONS,ACRL0371
C            JDEG=3 FOR THE CALCULATION OF THE FIRST TERM OF EQ.(3.16), ACRL0372
C            JDEG=1 FOR THE CALCULATION OF THE LAST TERM OF (3.16), AND ACRL0373
C            JDEG=2 OTHERWISE.                                          ACRL0374
C         IF ANY OF THE CURRENT JAPMAX CALCULATIONS WILL NOT NEED TO USEACRL0375
C            EQ.(3.16), JDEG=4 FOR THAT CALCULATION.  OF THE INTER-     ACRL0376
C            MEDIATE CROSS SECTION CALCULATIONS WHICH USE EQ.(3.16), THEACRL0377
C            ONES FOR WHICH EQ.(3.16) WILL COMPLETED FIRST SHOULD BE    ACRL0378
C            DONE LAST, I.E. THE JAPMAX CALCULATIONS SHOULD BE SEQUENCEDACRL0379
C            IN ORDER OF DECREASING JDEG.  THIS IS BECAUSE THE PARTIAL  ACRL0380
C            SUM IN EQ.(3.16) IS STORED IN SIGMA IN THE ORDER WHICH THE ACRL0381
C            INTERMEDIATE CROSS SECTIONS ARE FOUND, AND THE ABOVE ORDER ACRL0382
C            KEEPS THE PARTIAL SUMS INTACT.                             ACRL0383
C                                                                       ACRL0384
C        TK IS THE INITIAL STATE WAVE VECTOR SQUARED OF UNITS BOHR**-2. ACRL0385
C                                                                       ACRL0386
C        TMASS IS THE REDUCED MASS OF THE SCATTERING SYSTEM IN UNITS OF ACRL0387
C                             ELECTRON MASS = 1.                        ACRL0388
C        WQ IS THE STATISTICAL FACTOR IN EQUATION 3.16 OF BLATT AND     ACRL0389
C          BIEDENHARN.  FOR EXAMPLE WQ=0.25 FOR SINGLET ELECTRON        ACRL0390
C          SCATTERING AND WQ=0.75 FOR TRIPLET ELECTRON SCATTERING.      ACRL0391
C        IF EQ.(3.16) IS NOT TO BE USED, THE DEFAULT VALUE OF WQ IS 1.0.ACRL0392
C                                                                       ACRL0393
      CALL INOUT (3,LABEL,TK,TMASS,WQ,JTOT,JMIN,JDEG,8,IE,IG)           ACRL0394
C                                                                       ACRL0395
      IF ((JDEG.EQ.0).OR.(JDEG.EQ.4)) WQ=1.0                            ACRL0396
      IF (LEVEL(JAP).GE.0) JS=2                                         ACRL0397
      PIK=PI/TK                                                         ACRL0398
      PIK2=2.0*PIK                                                      ACRL0399
C      SET THE VALUE OF J AND JP TO BE USED FOR THIS CROSS SECTION.     ACRL0400
      J=JJ(JAP)                                                         ACRL0401
      JP=JJP(JAP)                                                       ACRL0402
C      NVALU AND NPVALU ARE USED TO DECODE THE DATA IN JINDEX.          ACRL0403
      NVALU=100*IABS(LEVEL(JAP))                                        ACRL0404
      NPVALU=100*IABS(LEVELP(JAP))                                      ACRL0405
C       IF TAPE ITAPE IS TO BE USED FOR STORING THE T MATRIX--          ACRL0406
C             SINCE ONLY A JPSEM BY JSEM SUBMATRIX OF EACH J-BLOCK WILL ACRL0407
C                   BE USED TO FIND THE BL COEFFICIENTS,                ACRL0408
C                   MAX(JPSEM) AND MAX(JSEM) DETERMINE WHAT VALUES      ACRL0409
C                              SHOULD BE USED FOR KMAX AND MMAX.        ACRL0410
      JSEM=2*J/JPAR+1                                                   ACRL0411
      JPSEM=2*JP/JPAR+1                                                 ACRL0412
C             ZERO THE NECESSARY ARRAYS                                 ACRL0413
      JMAX=JTOT+J                                                       ACRL0414
      DO 4 L1=1,JMAX                                                    ACRL0415
      OPAC(L1)=0.0                                                      ACRL0416
      OPACO(L1)=0.0                                                     ACRL0417
      PROBLO(L1)=0.0                                                    ACRL0418
   4  PROBL(L1)=0.0                                                     ACRL0419
      JSUM=2*J+1                                                        ACRL0420
      JPSUM=2*JP+1                                                      ACRL0421
      JEX=1                                                             ACRL0422
      IF (LEVEL(JAP)*LEVELP(JAP).GE.0) JEX=2                            ACRL0423
      FACJEX=(-1.0)**(JP-J)                                             ACRL0424
      FUDGE=1.0/(FLOAT(JSUM)*TK)                                        ACRL0425
      IALFA=IABS(LEVEL(JAP))                                            ACRL0426
      IALFAP=IABS(LEVELP(JAP))                                          ACRL0427
C                                                                       ACRL0428
      IF (LPR1) GO TO 13                                                ACRL0429
      LPR1=.TRUE.                                                       ACRL0430
C                                                                       ACRL0431
C       THIS IS ENTERED ONLY FOR READING THE R OR T MATRIX DATA SET AND ACRL0432
C                             FOR STORING THE T MATRIX.                 ACRL0433
   5  LZ=.FALSE.                                                        ACRL0434
      JTEMP=JMIN                                                        ACRL0435
C              DATA RECORD NUMBER SEVEN                                 ACRL0436
C                                                                       ACRL0437
C        READ IN THE J-BLOCK SPECIFICATIONS (LPR3,LCHECK,LOREAD,N,J1,   ACRL0438
C                         NOPT, AND JINDEX)                             ACRL0439
C                                                                       ACRL0440
C        N IS THE NUMBER OF ROWS IN THE R OR T MATRIX J-BLOCK TO BE READACRL0441
C                                                                       ACRL0442
C        J1 IS THE TOTAL ANGULAR MOMENTUM+1 OF THE J-BLOCK TO BE READ.  ACRL0443
C                                                                       ACRL0444
C     *********************************************************         ACRL0445
C     *  THE R OR T J-BLOCKS SHOULD BE READ WITH INCREASING   *         ACRL0446
C     *             J1 FROM J1 = JMIN TO J1 = JTOT.           *         ACRL0447
C     *********************************************************         ACRL0448
C                                                                       ACRL0449
C      NOPT=0--THE J-BLOCK IS TO BE READ WITH ITS ROWS IN THE ORDER MEN-ACRL0450
C              TIONED ABOVE WITHOUT EACH ROW NECESSARILY STARTING ON A  ACRL0451
C              NEW CARD. (SEE CARD NO.1043)                             ACRL0452
C      NOPT=1--THE J-BLOCK IS TO BE READ BY INDIVIDUAL ROWS IN THE ORDERACRL0453
C              MENTIONED ABOVE.                                         ACRL0454
C      NOPT=2--THE J-BLOCK IS TO BE READ BY INDIVIDUAL ROWS OUT OF ORDERACRL0455
C      NOPT=3--THE J-BLOCK IS TO BE READ WITH THE ROWS OUT OF ORDER     ACRL0456
C              WITHOUT EACH ROW NECESSARILY STARTING ON A NEW CARD.     ACRL0457
C                                                                       ACRL0458
C         JINDEX(I,JT) IS READ WITH THE VALUE --                        ACRL0459
C                     JINDEX(I,JT) = S + 100*ALPHA                      ACRL0460
C                       WHERE S IS THE CHANNEL ANGULAR MOMENTUM CORRES- ACRL0461
C                       PONDING TO THE ITH ROW OF THE J-BLOCK WITH TOTALACRL0462
C                       ANGULAR MOMENTUM (J1-1), AND ALPHA IS THE NON-  ACRL0463
C                       ANGULAR QUANTUM NUMBER OF THAT ROW.             ACRL0464
C                          JT = J1 - JMIN + 1                           ACRL0465
C                                                                       ACRL0466
      CALL INOUT (4,TEMP,TMP2,TMP3,JINDEX,N,J1,NOPT,NMAX,NMAX,JTAT)     ACRL0467
C                                                                       ACRL0468
C             SAVE THE SIZE OF THE J-BLOCK.                             ACRL0469
      JT=J1-JMIN+1                                                      ACRL0470
      NSAVE(JT)=N                                                       ACRL0471
      JX=J1-1                                                           ACRL0472
C              DATA RECORD NUMBERS EIGHT AND NINE                       ACRL0473
      IF (LPR3) GO TO 10                                                ACRL0474
C                                                                       ACRL0475
C       THIS IS REACHED ONLY WHEN AN R J-BLOCK IS TO BE READ.           ACRL0476
      CALL INOUT (5,KSTOR,EPS,H,R,N,NOPT,JX,NMAX,NMAX,NMAX)             ACRL0477
C                                                                       ACRL0478
      IF ((J1.EQ.JTOT).AND.LABORT) GO TO 67                             ACRL0479
C                                                                       ACRL0480
C       IF LABORT IS .TRUE. THE T MATRIX IS NOT NEEDED SINCE NO CROSS   ACRL0481
C          SECTIONS WILL BE FOUND WITH THIS DATA SET. (SEE CARD NO.129) ACRL0482
      IF (LABORT) GO TO 5                                               ACRL0483
C                                                                       ACRL0484
C       IF LABORT IS .FALSE., CONTINUE BY FINDING THE T J-BLOCK.        ACRL0485
C       PUT R*R INTO G                                                  ACRL0486
      CALL MXMPLY (R,R,G,N,N,N,NMAX,NMAX,NMAX,TEMP)                     ACRL0487
C       PUT 2*R INTO R, -2*R*R INTO H, AND 1 + R*R INTO G               ACRL0488
      DO 7 I=1,N                                                        ACRL0489
      DO 6 JK=1,N                                                       ACRL0490
      R(I,JK)=2.0*R(I,JK)                                               ACRL0491
   6  H(I,JK)=-2.0*G(I,JK)                                              ACRL0492
   7  G(I,I)=G(I,I)+1.0                                                 ACRL0493
C       PUT THE INVERSE OF G INTO G.  G WILL CONTAIN 1/(1 + R*R)        ACRL0494
      CALL MXLNEQ (G,N,NMAX,D,N,EPSD,TEMP)                              ACRL0495
      CALL MXMPLY (R,G,R,N,N,N,NMAX,NMAX,NMAX,TEMP)                     ACRL0496
C      R NOW CONTAINS 2*R/(1 + R*R), THE IMAGINARY PART OF THE T J-BLOCKACRL0497
      CALL MXMPLY (G,H,G,N,N,N,NMAX,NMAX,NMAX,TEMP)                     ACRL0498
C       G NOW CONTAINS -2*R*R/(1 + R*R), THE REAL PART OF THE T J-BLOCK.ACRL0499
      IF (JPAR1.GT.3) GO TO 12                                          ACRL0500
C                                                                       ACRL0501
C  THIS IS REACHED IF TAPE ITAPE WILL NOT BE USED TO STORE THE T MATRIX ACRL0502
      DO 9 IK=1,N                                                       ACRL0503
      DO 8 IL=1,N                                                       ACRL0504
      TREAL(IL,IK,JT)=G(IL,IK)                                          ACRL0505
   8  TIMAG(IL,IK,JT)=R(IL,IK)                                          ACRL0506
   9  CONTINUE                                                          ACRL0507
      IF (J1.EQ.JTOT) GO TO 13                                          ACRL0508
C                                                                       ACRL0509
C       THIS IS REACHED ONLY IF J1.NE.JTOT SO GO BACK AND READ THE NEXT ACRL0510
C        MATRIX DATA SPECIFICATION USED TO READ THE NEXT R OR T J-BLOCK.ACRL0511
      GO TO 5                                                           ACRL0512
  10  IF (JPAR1.GT.3) GO TO 11                                          ACRL0513
C                                                                       ACRL0514
C       THIS IS REACHED ONLY IF A T J-BLOCK IS TO BE READ INTO CORE.    ACRL0515
C                                                                       ACRL0516
C       FIRST READ THE REAL PART OF THE T J-BLOCK INTO THE ARRAY TREAL. ACRL0517
C                                                                       ACRL0518
      CALL INOUT (5,KSTOR,EPS,H,TREAL(1,1,JT),N,NOPT,JX,NMAX,KMAX,MMAX) ACRL0519
      LOREAD=.TRUE.                                                     ACRL0520
C                                                                       ACRL0521
C       LOREAD IS SET .TRUE. BECAUSE THE REAL AND IMAGINARY PARTS OF THEACRL0522
C          T J-BLOCK SHOULD HAVE THE SAME ROW ORDER.                    ACRL0523
C                                                                       ACRL0524
C      NOW READ THE IMAGINARY PART OF THE T J-BLOCK INTO THE ARRAY TIMAGACRL0525
C                                                                       ACRL0526
      CALL INOUT (5,KSTOR,EPS,H,TIMAG(1,1,JT),N,NOPT,JX,NMAX,KMAX,MMAX) ACRL0527
      IF ((J1.EQ.JTOT).AND.LABORT) GO TO 67                             ACRL0528
      IF (J1.EQ.JTOT) GO TO 13                                          ACRL0529
C                                                                       ACRL0530
C       THIS IS REACHED ONLY IF J1.NE.JTOT SO GO BACK AND READ THE NEXT ACRL0531
C        MATRIX DATA SPECIFICATION USED TO READ THE NEXT R OR T J-BLOCK.ACRL0532
      GO TO 5                                                           ACRL0533
C                                                                       ACRL0534
C   THIS IS REACHED ONLY IF THE T MATRIX IS TO BE STORED ON TAPE ITAPE. ACRL0535
C       FIRST READ THE REAL PART OF THE T J-BLOCK INTO THE ARRAY G.     ACRL0536
  11  CALL INOUT (5,KSTOR,EPS,H,G,N,NOPT,JX,NMAX,NMAX,NMAX)             ACRL0537
C                                                                       ACRL0538
      LOREAD=.TRUE.                                                     ACRL0539
C       NOW READ THE IMAGINARY PART OF THE T J-BLOCK INTO THE ARRAY R.  ACRL0540
      CALL INOUT (5,KSTOR,EPS,H,R,N,NOPT,JX,NMAX,NMAX,NMAX)             ACRL0541
      IF (LABORT.AND.(J1.EQ.JTOT)) GO TO 67                             ACRL0542
      IF (LABORT) GO TO 5                                               ACRL0543
C                                                                       ACRL0544
C       THIS IS REACHED WITH THE REAL PART OF THE T J-BLOCK IN G AND THEACRL0545
C         IMAGINARY PART IN R READY TO BE STORED ON TAPE ITAPE.         ACRL0546
  12  CALL INOUT (6,TMP1,TMP2,R,G,N,J1,IC,NMAX,NMAX,NMAX)               ACRL0547
      IF (J1.EQ.JTOT) REWIND ITAPE                                      ACRL0548
      IF (J1.EQ.JTOT) GO TO 13                                          ACRL0549
C                                                                       ACRL0550
C       THIS IS REACHED ONLY IF J1.NE.JTOT SO GO BACK AND READ THE NEXT ACRL0551
C        MATRIX DATA SPECIFICATION USED TO READ THE NEXT R OR T J-BLOCK.ACRL0552
      GO TO 5                                                           ACRL0553
C                                                                       ACRL0554
C       THIS IS REACHED AFTER THE T MATRIX HAS BEEN STORED SO BEGIN THE ACRL0555
C            THE NEXT INTERMEDIATE CROSS SECTION CALCULATION TO BE      ACRL0556
C                        DONE WITH THIS DATA SET.                       ACRL0557
C                                                                       ACRL0558
C       INFORM SIGMA WHERE THIS CROSS SECTION IS TO BE STORED.          ACRL0559
  13  JAX=JAP                                                           ACRL0560
      IF (JDEG.LE.2) JAX=JAPMAX+1                                       ACRL0561
      IF (JDEG.EQ.0) JAX=1                                              ACRL0562
C                                                                       ACRL0563
C    IN THIS SECTION THE OPACITY IS FOUND AND IF TAPE ITAPE WAS USED THEACRL0564
C         REQUIRED PART OF THE T MATRIX IS STORED IN TREAL AND TIMAG.   ACRL0565
      DO 25 J1=JTEMP,JTOT                                               ACRL0566
      IF (JPAR1.GT.3) GO TO 17                                          ACRL0567
      IF (J1.LT.JMIN) GO TO 25                                          ACRL0568
      JT=J1-JTEMP+1                                                     ACRL0569
      JX=J1-1                                                           ACRL0570
      N=NSAVE(JT)                                                       ACRL0571
C       STORE THE INITIAL STATE ROW INDEX IN K1MIN AND THE FINAL STATE  ACRL0572
C         COLUMN INDEX IN L3MIN FOR EACH TOTAL ANGULAR MOMENTUM (J1-1). ACRL0573
      K1MIN(JT)=0                                                       ACRL0574
      DO 14 KK=1,N                                                      ACRL0575
      K1MIN(JT)=K1MIN(JT)+1                                             ACRL0576
      IF ((JINDEX(KK,JT)-NVALU).EQ.J) GO TO 15                          ACRL0577
  14  CONTINUE                                                          ACRL0578
  15  KT=KK                                                             ACRL0579
      CALL OPACIT (J,JX,JPAR,NVALU,N,KK,JS,JINDEX(1,JT),TREAL(1,1,JT),  ACRL0580
     1  TIMAG(1,1,JT),OPAC,OPACO,NMAX,KMAX,MMAX,JMAXP)                  ACRL0581
      L3MIN(JT)=0                                                       ACRL0582
      DO 16 LK=1,N                                                      ACRL0583
      LJ=LK                                                             ACRL0584
      L3MIN(JT)=L3MIN(JT)+1                                             ACRL0585
      IF ((JINDEX(LK,JT)-NPVALU).EQ.JP) GO TO 24                        ACRL0586
  16  CONTINUE                                                          ACRL0587
C                                                                       ACRL0588
C   READ BACK FROM ITAPE THE T J-BLOCK OF TOTAL ANGULAR MOMENTUM (J2-1) ACRL0589
C                                                                       ACRL0590
  17  CALL INOUT (7,TMP1,TMP2,R,G,N,J2,IC,NMAX,NMAX,NMAX)               ACRL0591
      IF (J2.LT.JMIN) GO TO 25                                          ACRL0592
      JX=J2-1                                                           ACRL0593
      JT=J2-JTEMP+1                                                     ACRL0594
      K1MIN(JT)=1                                                       ACRL0595
      L3MIN(JT)=1                                                       ACRL0596
      K2MIN=0                                                           ACRL0597
      DO 18 KK=1,N                                                      ACRL0598
      K2MIN=K2MIN+1                                                     ACRL0599
      IF ((JINDEX(KK,JT)-NVALU).EQ.J) GO TO 19                          ACRL0600
  18  CONTINUE                                                          ACRL0601
  19  KT=KK                                                             ACRL0602
      CALL OPACIT (J,JX,JPAR,NVALU,N,KK,JS,JINDEX(1,JT),G,R,OPAC,OPACO, ACRL0603
     1  NMAX,NMAX,NMAX,JMAXP)                                           ACRL0604
      L4MIN=0                                                           ACRL0605
      DO 20 LK=1,N                                                      ACRL0606
      LJ=LK                                                             ACRL0607
      L4MIN=L4MIN+1                                                     ACRL0608
      IF ((JINDEX(LK,JT)-NPVALU).EQ.JP) GO TO 21                        ACRL0609
  20  CONTINUE                                                          ACRL0610
C        IN THIS SECTION THE NEEDED T MATRIX ELEMENTS ARE STORED.       ACRL0611
  21  K1MAX=K2MIN+JSEM-1                                                ACRL0612
      IF (JX.LT.J) K1MAX=K2MIN+2*JX/JPAR                                ACRL0613
      L4MAX=L4MIN+JPSEM-1                                               ACRL0614
      IF (JX.LT.JP) L4MAX=L4MIN+2*JX/JPAR                               ACRL0615
      K=1                                                               ACRL0616
      DO 23 IK=K2MIN,K1MAX                                              ACRL0617
      L=1                                                               ACRL0618
      DO 22 IL=L4MIN,L4MAX                                              ACRL0619
      TREAL(L,K,JT)=G(IL,IK)                                            ACRL0620
      TIMAG(L,K,JT)=R(IL,IK)                                            ACRL0621
  22  L=L+1                                                             ACRL0622
  23  K=K+1                                                             ACRL0623
C                  CHECK THE JINDEX VALUES                              ACRL0624
  24  K=KT+JSEM-1                                                       ACRL0625
      IF (JX.LT.J) K=KT+2*JX/JPAR                                       ACRL0626
      L=LJ+JPSEM-1                                                      ACRL0627
      IF (JX.LT.JP) L=LJ+2*JX/JPAR                                      ACRL0628
      I=JINDEX(LJ,JT)                                                   ACRL0629
      II=JINDEX(L,JT)                                                   ACRL0630
      III=JINDEX(KT,JT)                                                 ACRL0631
      IIII=JINDEX(K,JT)                                                 ACRL0632
      IF ((I.EQ.II).AND.(III.EQ.IIII)) GO TO 25                         ACRL0633
C                                                                       ACRL0634
C               THIS IS REACHED WHEN AN ERROR OCCURS                    ACRL0635
C                                                                       ACRL0636
      CALL INOUT (8,TMP1,TMP2,TMP3,TMP4,JX,IB,IC,ID,IE,IG)              ACRL0637
      GO TO 67                                                          ACRL0638
  25  CONTINUE                                                          ACRL0639
C                                                                       ACRL0640
C  THIS SECTION CALCULATES THE BL COEFFICIENTS .                        ACRL0641
      K1=1                                                              ACRL0642
      LMAX=2*JTOT-1                                                     ACRL0643
      DO 58 LAMBDA=1,LMAX                                               ACRL0644
      K1=K1+1                                                           ACRL0645
      B(LAMBDA)=0.0                                                     ACRL0646
      L=LAMBDA-1                                                        ACRL0647
C          BEGIN THE SUM OVER J1 IN EQN. (4.7) OF BLATT AND BIEDENHARN  ACRL0648
      DO 57 JV=JMIN,JTOT                                                ACRL0649
      JT=JV-JTEMP+1                                                     ACRL0650
      J1=JV-1                                                           ACRL0651
      J1T2=2*J1                                                         ACRL0652
C       USE THE TRIANGLE RULE TO FIND THE LIMITS ON THE SUMS            ACRL0653
      IF (J1-J) 26,27,26                                                ACRL0654
  26  L1MIN=IABS(J1-J)+1                                                ACRL0655
      GO TO 28                                                          ACRL0656
  27  L1MIN=1                                                           ACRL0657
  28  IF (J1-JP) 29,30,29                                               ACRL0658
  29  L1PMIN=IABS(J1-JP)+1                                              ACRL0659
      GO TO 31                                                          ACRL0660
  30  L1PMIN=1                                                          ACRL0661
  31  L1MAX=J1+J+1                                                      ACRL0662
      L1PMAX=J1+JP+1                                                    ACRL0663
      L1C=K1MIN(JT)                                                     ACRL0664
C             BEGIN THE INNER SUM OVER L1 AND J1                        ACRL0665
      DO 56 L1Q=L1MIN,L1MAX,JPAR                                        ACRL0666
      LX=L1Q                                                            ACRL0667
      L1=L1Q-1                                                          ACRL0668
      L1PL=L1+L                                                         ACRL0669
      L1ML=IABS(L1-L)                                                   ACRL0670
      LS=JEX+L1                                                         ACRL0671
      LR=JEX+L1PL                                                       ACRL0672
      LZ=.TRUE.                                                         ACRL0673
      IF ((2*L1.LT.L).OR.(J1T2.LT.L).OR.(K1.LT.2)) GO TO 32             ACRL0674
      LZ=.FALSE.                                                        ACRL0675
      YC=Z(L1,J1,L1,J1,J,L)                                             ACRL0676
  32  L1PC=L3MIN(JT)                                                    ACRL0677
C                                                                       ACRL0678
C         SUM1 IS THE THREEFOLD SUM OVER J,L AND LP OF (4.7)            ACRL0679
      SUM1=0.0                                                          ACRL0680
C                                                                       ACRL0681
C             BEGIN THE INNER SUM OVER L1P, L1 AND J1                   ACRL0682
      DO 53 L1R=L1PMIN,L1PMAX,JPAR                                      ACRL0683
      L1P=L1R-1                                                         ACRL0684
      CALL PARITY (LS,L1P)                                              ACRL0685
      IF (LT) GO TO 53                                                  ACRL0686
      L1PPL=L1P+L                                                       ACRL0687
      L1PML=IABS(L1P-L)                                                 ACRL0688
      IF (J1T2.LT.L) GO TO 40                                           ACRL0689
      IF (LZ) GO TO 36                                                  ACRL0690
C   SUM1 AND SUM2 ARE BOTH ZERO FOR ODD VALUES OF THE VARIABLE L.       ACRL0691
      IF (2*L1P.LT.L) GO TO 33                                          ACRL0692
      YA=Z(L1P,J1,L1P,J1,JP,L)                                          ACRL0693
      SUM1=SUM1+YA*(TREAL(L1PC,L1C,JT)**2+TIMAG(L1PC,L1C,JT)**2)        ACRL0694
  33  IF (L1PMAX-L1R) 36,36,34                                          ACRL0695
  34  L2PMIN=L1P+JPAR+1                                                 ACRL0696
      L2PMAX=J1+JP+1                                                    ACRL0697
      L2PC=L1PC+1                                                       ACRL0698
      SUM2=0.0                                                          ACRL0699
C                                                                       ACRL0700
C  BEGIN THE SUM OVER L2P CONTAINED IN THE CURLY BRACKETS(C.B.) OF (4.7)ACRL0701
C                 THIS SUM IS THE VARIABLE SUM2.                        ACRL0702
C                                                                       ACRL0703
      DO 35 L2R=L2PMIN,L2PMAX,JPAR                                      ACRL0704
      L2P=L2R-1                                                         ACRL0705
      CALL PARITY (L1P,L2P)                                             ACRL0706
      IF ((L2P.LT.L1PML).OR.(L2P.GT.L1PPL).OR.LT) GO TO 35              ACRL0707
      YB=Z(L1P,J1,L2P,J1,JP,L)                                          ACRL0708
      SUM2=SUM2+YB*(TREAL(L1PC,L1C,JT)*TREAL(L2PC,L1C,JT)+TIMAG(L1PC,L1CACRL0709
     1,JT)*TIMAG(L2PC,L1C,JT))                                          ACRL0710
  35  L2PC=L2PC+1                                                       ACRL0711
      B(LAMBDA)=B(LAMBDA)+YC*2.0*SUM2                                   ACRL0712
  36  IF (L1MAX-LX) 40,40,37                                            ACRL0713
  37  L2MIN=L1+JPAR+1                                                   ACRL0714
      L2MAX=J1+J+1                                                      ACRL0715
      L2C=L1C+1                                                         ACRL0716
C                                                                       ACRL0717
C             BEGIN THE SUM OVER L2 CONTAINED IN THE C.B.               ACRL0718
      DO 39 L2Q=L2MIN,L2MAX,JPAR                                        ACRL0719
      L2=L2Q-1                                                          ACRL0720
      CALL PARITY (L2,L1PL)                                             ACRL0721
      IF ((L2.LT.L1ML).OR.(L2.GT.L1PL).OR.LT) GO TO 39                  ACRL0722
      Y=Z(L1,J1,L2,J1,J,L)                                              ACRL0723
      L2PMIN=L1PMIN                                                     ACRL0724
      L2PMAX=L1PMAX                                                     ACRL0725
      L2PC=L3MIN(JT)                                                    ACRL0726
      SUM3=0.0                                                          ACRL0727
C                                                                       ACRL0728
C      BEGIN THE INNER SUM, SUM3, OVER L2P AND L2 CONTAINED IN THE C.B. ACRL0729
      DO 38 L2R=L2PMIN,L2PMAX,JPAR                                      ACRL0730
      L2P=L2R-1                                                         ACRL0731
      CALL PARITY (LR,L2P)                                              ACRL0732
      IF (.NOT.LT) CALL PARITY (L2P,L1PPL)                              ACRL0733
      IF ((L2P.LT.L1PML).OR.(L2P.GT.L1PPL).OR.LT) GO TO 38              ACRL0734
      YB=Z(L1P,J1,L2P,J1,JP,L)                                          ACRL0735
      SUM3=SUM3+YB*(TREAL(L1PC,L1C,JT)*TREAL(L2PC,L2C,JT)+TIMAG(L1PC,L1CACRL0736
     1,JT)*TIMAG(L2PC,L2C,JT))                                          ACRL0737
  38  L2PC=L2PC+1                                                       ACRL0738
      B(LAMBDA)=B(LAMBDA)+2.0*Y*SUM3                                    ACRL0739
  39  L2C=L2C+1                                                         ACRL0740
  40  IF (J1-L) 41,42,41                                                ACRL0741
  41  J2MIN=MAX0(J1+1,IABS(J1-L))                                       ACRL0742
      GO TO 43                                                          ACRL0743
  42  J2MIN=J1+1                                                        ACRL0744
  43  J2MAX=MIN0(J1+L,JTOT-1)                                           ACRL0745
      IF (J2MIN.GT.J2MAX) GO TO 53                                      ACRL0746
C                                                                       ACRL0747
C              BEGIN THE SUM OVER J2 CONTAINED IN THE C.B.              ACRL0748
      DO 52 J2=J2MIN,J2MAX                                              ACRL0749
      JY=J2-JTEMP+2                                                     ACRL0750
      IF (J2-J) 44,45,44                                                ACRL0751
  44  L2MIN=IABS(J2-J)+1                                                ACRL0752
      GO TO 46                                                          ACRL0753
  45  L2MIN=1                                                           ACRL0754
  46  IF (J2-JP) 47,48,47                                               ACRL0755
  47  L2PMIN=IABS(J2-JP)+1                                              ACRL0756
      GO TO 49                                                          ACRL0757
  48  L2PMIN=1                                                          ACRL0758
  49  L2MAX=J2+J+1                                                      ACRL0759
      L2PMAX=J2+JP+1                                                    ACRL0760
      L2C=K1MIN(JY)                                                     ACRL0761
C                                                                       ACRL0762
C         BEGIN THE INNER SUM OVER L2 AND J2 CONTAINED IN THE C.B.      ACRL0763
      DO 51 L2Q=L2MIN,L2MAX,JPAR                                        ACRL0764
      L2=L2Q-1                                                          ACRL0765
      CALL PARITY (L2,L1PL)                                             ACRL0766
      IF ((L2.LT.L1ML).OR.(L2.GT.L1PL).OR.LT) GO TO 51                  ACRL0767
      Y=Z(L1,J1,L2,J2,J,L)                                              ACRL0768
      L2PC=L3MIN(JY)                                                    ACRL0769
      SUM4=0.0                                                          ACRL0770
C                                                                       ACRL0771
C   BEGIN THE INNER SUM, SUM4, OVER L2P, L2 AND J2 CONTAINED IN THE C.B.ACRL0772
      DO 50 L2R=L2PMIN,L2PMAX,JPAR                                      ACRL0773
      L2P=L2R-1                                                         ACRL0774
      CALL PARITY (LR,L2P)                                              ACRL0775
      IF (.NOT.LT) CALL PARITY (L2P,L1PPL)                              ACRL0776
      IF ((L2P.LT.L1PML).OR.(L2P.GT.L1PPL).OR.LT) GO TO 50              ACRL0777
      YB=Z(L1P,J1,L2P,J2,JP,L)                                          ACRL0778
      SUM4=SUM4+YB*(TREAL(L1PC,L1C,JT)*TREAL(L2PC,L2C,JY)+TIMAG(L1PC,L1CACRL0779
     1,JT)*TIMAG(L2PC,L2C,JY))                                          ACRL0780
  50  L2PC=L2PC+1                                                       ACRL0781
      B(LAMBDA)=B(LAMBDA)+2.0*Y*SUM4                                    ACRL0782
  51  L2C=L2C+1                                                         ACRL0783
  52  CONTINUE                                                          ACRL0784
  53  L1PC=L1PC+1                                                       ACRL0785
      IF (.NOT.LZ) B(LAMBDA)=B(LAMBDA)+YC*SUM1                          ACRL0786
C                                                                       ACRL0787
C         HERE THE TRANSITION PROBABILITIES ARE STORED BY THEIR INITIAL ACRL0788
C    RELATIVE ORBITAL ANGULAR MOMENTUM (LX - 1).                        ACRL0789
      IF (L.EQ.0) GO TO 54                                              ACRL0790
      GO TO 56                                                          ACRL0791
  54  CALL PARITY (JS,L1)                                               ACRL0792
      FJT=FLOAT(2*L1+1)                                                 ACRL0793
      IF (.NOT.LT) GO TO 55                                             ACRL0794
      PROBLO(LX)=PROBLO(LX)+YC*SUM1/FJT*FACJEX                          ACRL0795
      GO TO 56                                                          ACRL0796
  55  PROBL(LX)=PROBL(LX)+YC*SUM1/FJT*FACJEX                            ACRL0797
  56  L1C=L1C+1                                                         ACRL0798
  57  CONTINUE                                                          ACRL0799
      B(LAMBDA)=0.25*FACJEX*B(LAMBDA)*WQ                                ACRL0800
      IF (K1.GE.2) K1=0                                                 ACRL0801
  58  CONTINUE                                                          ACRL0802
C                                                                       ACRL0803
C  THIS SECTION CALCULATES THE DIFFERENTIAL CROSS SECTION.              ACRL0804
      THETA=0.0                                                         ACRL0805
      DO 60 I=1,181                                                     ACRL0806
      INDEX(I)=I-1                                                      ACRL0807
      IF (.NOT.LSTOR) X=COS(THETA)                                      ACRL0808
      SIGMA(I,JAX)=0.0                                                  ACRL0809
      DO 59 LAMBDA=1,LMAX                                               ACRL0810
      L=LAMBDA-1                                                        ACRL0811
      IF (LSTOR) SIGMA(I,JAX)=SIGMA(I,JAX)+PSTOR(I,LAMBDA)*B(LAMBDA)    ACRL0812
      IF (.NOT.LSTOR) SIGMA(I,JAX)=SIGMA(I,JAX)+P(L,X)*B(LAMBDA)        ACRL0813
  59  CONTINUE                                                          ACRL0814
      SIGMA(I,JAX)=FUDGE*SIGMA(I,JAX)                                   ACRL0815
  60  THETA=THETA+DELTA                                                 ACRL0816
C                                                                       ACRL0817
C       WRITE J AND JP                                                  ACRL0818
      CALL INOUT (9,IALFA,IALFAP,JJ,JJP,JAP,IB,IC,JATMAX,JATMAX,IG)     ACRL0819
C  HERE THE INITIAL K.E. IS CALCULATED IN VARIOUS SETS OF UNITS.        ACRL0820
      TK=TK/(2.*TMASS)                                                  ACRL0821
      TL=TK*27.210                                                      ACRL0822
      TJ=1.602E-12*TL                                                   ACRL0823
      VEL=SQRT(2.0*TJ/(9.11E-28*TMASS))                                 ACRL0824
      TM=1.4397E+13*TJ                                                  ACRL0825
      TN=TJ/(6.627*2.99793E-17)                                         ACRL0826
      CALL INOUT (10,TK,TL,TJ,TM,TN,FUDGE,VEL,ID,IE,IG)                 ACRL0827
C                                                                       ACRL0828
C                                                                       ACRL0829
C       WRITE OUT THE OPACITY AND TRANSITION PROBABLITY                 ACRL0830
      X=WQ/(FLOAT(JSUM))                                                ACRL0831
      DO 61 L1=1,JMAX                                                   ACRL0832
      OPAC(L1)=OPAC(L1)*X                                               ACRL0833
      OPACO(L1)=OPACO(L1)*X                                             ACRL0834
      PROBL(L1)=PROBL(L1)*X                                             ACRL0835
  61  PROBLO(L1)=PROBLO(L1)*X                                           ACRL0836
      CALL INOUT (11,PROBL,PROBLO,OPAC,OPACO,JMAX,IB,IC,JMAXP,JMAXP,IG) ACRL0837
C                                                                       ACRL0838
C  THIS SECTION CALCULATES THE PARTIAL INTEGRAL CROSS SECTION.          ACRL0839
      QOM=0.0                                                           ACRL0840
      QM=0.0                                                            ACRL0841
      XL=-PIK                                                           ACRL0842
      DO 62 L1=1,JMAX                                                   ACRL0843
      XL=XL+PIK2                                                        ACRL0844
      QM=QM+XL*PROBL(L1)                                                ACRL0845
      QOM=QOM+XL*PROBLO(L1)                                             ACRL0846
      OPAC(L1)=OPAC(L1)+OPACO(L1)                                       ACRL0847
      PROBL(L1)=PROBL(L1)+PROBLO(L1)                                    ACRL0848
      PROBLO(L1)=QM                                                     ACRL0849
  62  OPACO(L1)=QOM                                                     ACRL0850
      CALL INOUT (12,PROBL,OPAC,PROBLO,OPACO,JMAX,JDEG,IC,JMAXP,JMAXP,  ACRL0851
     1  IG)                                                             ACRL0852
      X=DELTA                                                           ACRL0853
      TOTSIG=0.0                                                        ACRL0854
      K=2                                                               ACRL0855
C                                                                       ACRL0856
C  THIS SECTION CALCULATES THE INTEGRAL CROSS SECTION BY WEDDLE'S RULE. ACRL0857
      DO 63 I=1,30                                                      ACRL0858
      TOTSIG=TOTSIG+5.0*SIGMA(K,JAX)*SIN(X)                             ACRL0859
      X=X+DELTA                                                         ACRL0860
      K=K+1                                                             ACRL0861
      TOTSIG=TOTSIG+SIGMA(K,JAX)*SIN(X)                                 ACRL0862
      X=X+DELTA                                                         ACRL0863
      K=K+1                                                             ACRL0864
      TOTSIG=TOTSIG+6.0*SIGMA(K,JAX)*SIN(X)                             ACRL0865
      X=X+DELTA                                                         ACRL0866
      K=K+1                                                             ACRL0867
      TOTSIG=TOTSIG+SIGMA(K,JAX)*SIN(X)                                 ACRL0868
      X=X+DELTA                                                         ACRL0869
      K=K+1                                                             ACRL0870
      TOTSIG=TOTSIG+5.0*SIGMA(K,JAX)*SIN(X)                             ACRL0871
      X=X+DELTA                                                         ACRL0872
      K=K+1                                                             ACRL0873
      IF (K.EQ.181) GO TO 64                                            ACRL0874
      TOTSIG=TOTSIG+2.0*SIGMA(K,JAX)*SIN(X)                             ACRL0875
      X=X+DELTA                                                         ACRL0876
  63  K=K+1                                                             ACRL0877
C                              FACTOR = PI**2/300.                      ACRL0878
  64  TOTSG1=FACTOR*TOTSIG                                              ACRL0879
      TOTSG2=PI4*B(1)*FUDGE                                             ACRL0880
      TOTSG4=TOTSG2*ANG2                                                ACRL0881
      SIGMA(182,JAX)=TOTSG1                                             ACRL0882
      SIGMA(183,JAX)=TOTSG4                                             ACRL0883
C                                                                       ACRL0884
C       WRITE OUT THE BL COEFFICIENTS.                                  ACRL0885
      CALL INOUT (13,INDEX,B,WQ,TMP4,LMAX,IB,IC,INMX,LMAXP,IG)          ACRL0886
C                                                                       ACRL0887
C       WRITE OUT THE CROSS SECTION                                     ACRL0888
      CALL INOUT (14,TOTSG1,INDEX,TOTSG2,SIGMA,TOTSG4,JDEG,JAX,ID,INMX, ACRL0889
     1  JAXMAX)                                                         ACRL0890
      IF ((JDEG.NE.2).AND.(JDEG.NE.1)) GO TO 66                         ACRL0891
C                                                                       ACRL0892
C          THIS IS REACHED WHEN ANOTHER TERM IS TO BE ADDED TO (3.16)   ACRL0893
      DO 65 I=1,183                                                     ACRL0894
  65  SIGMA(I,JAP)=SIGMA(I,JAP)+SIGMA(I,JAX)                            ACRL0895
      IF (JDEG.NE.1) GO TO 66                                           ACRL0896
C                                                                       ACRL0897
C                  WRITE OUT THE CROSS SECTION FOUND BY (3.16)          ACRL0898
      JDEG=9                                                            ACRL0899
      TOTSG1=SIGMA(182,JAP)                                             ACRL0900
      TOTSG4=SIGMA(183,JAP)                                             ACRL0901
      CALL INOUT (14,TOTSG1,INDEX,TOTSG1,SIGMA,TOTSG4,JDEG,JAP,ID,INMX, ACRL0902
     1  JAXMAX)                                                         ACRL0903
  66  CONTINUE                                                          ACRL0904
  67  LABORT=.FALSE.                                                    ACRL0905
C                                                                       ACRL0906
C      THIS IS REACHED WHEN CROSS IS FINISHED WITH THE CURRENT DATA SET.ACRL0907
      IF (IAP.EQ.IAPMAX) GO TO 68                                       ACRL0908
C                                                                       ACRL0909
C       THIS IS REACHED IF MORE DATA SETS ARE TO BE READ.               ACRL0910
      LPR1=.FALSE.                                                      ACRL0911
C              DATA RECORD NUMBERS TWELVE AND FOUR                      ACRL0912
C                                                                       ACRL0913
C            READ THE NEXT VALUE OF JAPMAX                              ACRL0914
C                                                                       ACRL0915
C         READ THE NEXT VALUE OF JPAR1, IEPS, AND THE SET OF INITIAL ANDACRL0916
C    FINAL QUANTUM STATES FOR THE NEXT SET OF CROSS SECTIONS.           ACRL0917
C                                                                       ACRL0918
      IF (LPR2) CALL INOUT (15,LEVEL,JJ,LEVELP,JJP,JAPMAX,JPAR1,IEPS,   ACRL0919
     1  JATMAX,JATMAX,IG)                                               ACRL0920
  68  CONTINUE                                                          ACRL0921
      TIME2=CPTIME(.1)                                                  ACRL0922
      DIFF=TIME2-TIME1                                                  ACRL0923
      CALL INOUT (16,DIFF,TMP2,TMP3,TMP4,IA,IB,IC,ID,IE,IG)             ACRL0924
      RETURN                                                            ACRL0925
      END                                                               ACRL0926
      SUBROUTINE INOUT (N,TMP1,TMP2,TMP3,TMP4,IA,IB,IC,ID,IE,IG)        ACRL0927
C                                                                       ACRL0928
C         THIS IS THE INPUT/OUTPUT PART OF THE PROGRAM.  INOUT          ACRL0929
C    CONTAINS ALL THE READ AND WRITE INSTRUCTIONS EXCEPT THE TWO IN     ACRL0930
C    SUBPROGRAM MXLNEQ AND THE ONE IN SUBPROGRAM F6J.                   ACRL0931
C                                                                       ACRL0932
      LOGICAL LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LT    ACRL0933
      COMMON/LCOM/LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LTACRL0934
C                                                                       ACRL0935
C         THE INITIAL VALUES OF THE /LCOM/ VARIABLES LPR2,LSTOR,LCHECK, ACRL0936
C    AND LOREAD ARE READ AT STATEMENTS 1 AND 4.  THE VALUES OF LZ AND   ACRL0937
C    LABORT MAY BE CHANGED IN STATEMENTS 9 TO 19.                       ACRL0938
C       THE EXPLANATION OF THE /IOSEND/ VARIABLES BEGINS ON CARD NO.19. ACRL0939
      COMMON /IOSEND/ K,L,M                                             ACRL0940
C                                                                       ACRL0941
C     K = IREAD                                                         ACRL0942
C     L = IWRITE                                                        ACRL0943
C     M = ITAPE                                                         ACRL0944
C    ****************************************************************   ACRL0945
C    *    STATEMENTS 21 AND 22 REQUIRE THE EXTERNAL DEVICE ITAPE.   *   ACRL0946
C    ****************************************************************   ACRL0947
C                                                                       ACRL0948
      DIMENSION TMP1(ID),TMP2(IE),TMP3(ID,ID),TMP4(IE,IG)               ACRL0949
C                                                                       ACRL0950
C                         PARAMETER LIST                                ACRL0951
C                                                                       ACRL0952
C     N  -  AN INPUT VARIABLE  -  USED IN THE COMPUTED GO TO STATEMENT  ACRL0953
C         THE OTHER VARIABLES IN THE PARAMETER LIST ARE INPUT OR OUTPUT ACRL0954
C    VARIABLES WHICH ARE WRITTEN OR READ IN THE SUBPROGRAM.             ACRL0955
C         AT STATEMENT NOS. 2, 4-8, 9, 24-26, 28, AND 30-33 INTEGER     ACRL0956
C    (REAL) VARIABLES AND/OR ARRAYS ARE READ AND/OR WRITTEN USING REAL  ACRL0957
C    (INTEGER) FORMAT SPECIFICATIONS.  IN EVERY CASE THESE VARIABLE     ACRL0958
C    AND/OR ARRAY NAMES ARE USED ONLY TO TRANSFER NUMBERS AND ARE NOT   ACRL0959
C    USED IN ARITHMETIC STATEMENTS WHICH WOULD CHANGE THE VARIABLE AND/ ACRL0960
C    OR ARRAY TYPE.  THE VARIABLE AND/OR ARRAY NAMES IN THE PARAMETER   ACRL0961
C    LISTS OF THE CALLS TO SUBPROGRAM INOUT ARE MATCHED TO THE FORMAT   ACRL0962
C    SPECIFICATIONS SO THAT THE INPUT DATA IS EVENTUALLY STORED IN THE  ACRL0963
C    PROPER TYPE LOCATION BEFORE IT IS USED AND THE OUTPUT IS PRINTED   ACRL0964
C    PROPERLY.                                                          ACRL0965
C                                                                       ACRL0966
      GO TO (1,2,3,4,9,21,22,23,24,25,26,28,30,31,34,35),N              ACRL0967
C                                                                       ACRL0968
C           THIS WRITES THE DATA HEADING ON A NEW PAGE                  ACRL0969
C                                                                       ACRL0970
   1  WRITE (L,39)                                                      ACRL0971
C                     THIS READS AND WRITES THE LABEL.                  ACRL0972
C                                                                       ACRL0973
      READ (K,40) (TMP1(I),I=1,8)                                       ACRL0974
      WRITE (L,41) (TMP1(I),I=1,8)                                      ACRL0975
      READ (K,42) LPR2,LSTOR                                            ACRL0976
      WRITE (L,43) LPR2,LSTOR                                           ACRL0977
C            THIS READS AND WRITES THE VALUES OF IAPMAX AND JAPMAX.     ACRL0978
C                                                                       ACRL0979
      READ (K,44) IA,IB                                                 ACRL0980
      WRITE (L,45) IA,IB                                                ACRL0981
      RETURN                                                            ACRL0982
C     THIS READS AND WRITES THE VALUES OF JPAR1, IEPS, LEVEL, JJ,       ACRL0983
C                     LEVELP, AND JJP.                                  ACRL0984
C                                                                       ACRL0985
   2  READ (K,44) IB,IC,((TMP1(I),TMP2(I),TMP3(I,1),TMP4(I,1)),I=1,IA)  ACRL0986
      WRITE (L,81) IB,IC,((TMP1(I),TMP2(I),TMP3(I,1),TMP4(I,1)),I=1,IA) ACRL0987
      RETURN                                                            ACRL0988
C           THIS WRITES DATA ON A NEW PAGE.                             ACRL0989
   3  IF (LZ) WRITE (L,39)                                              ACRL0990
C                     THIS READS AND WRITES THE LABEL.                  ACRL0991
C                                                                       ACRL0992
      READ (K,40) (TMP1(I),I=1,8)                                       ACRL0993
      IF (.NOT.LPR1) WRITE (L,78) (TMP1(I),I=1,8)                       ACRL0994
      IF (LPR1) WRITE (L,79) (TMP1(I),I=1,8)                            ACRL0995
C         THIS READS AND WRITES THE VALUES OF JTOT, JMIN, JDEG, TK,     ACRL0996
C                                 TMASS, AND WQ.                        ACRL0997
C                                                                       ACRL0998
      READ (K,46) IA,IB,IC,TMP2(1),TMP3(1,1),TMP4(1,1)                  ACRL0999
      IF (.NOT.LPR1) WRITE (L,47) IA,IB,IC,TMP2(1),TMP3(1,1),TMP4(1,1)  ACRL1000
      IF (LPR1) WRITE (L,84) IA,IB,IC,TMP2(1),TMP3(1,1),TMP4(1,1)       ACRL1001
      JMIN=IB                                                           ACRL1002
      JTOT=IA                                                           ACRL1003
      RETURN                                                            ACRL1004
C      THIS READS AND WRITES THE VALUES OF LPR3, LCHECK, NOREAD, N, J1, ACRL1005
C               NOPT, AND JINDEX FOR THE NEXT J-BLOCK TO BE READ.       ACRL1006
   4  IF (ID.LT.10) GO TO 5                                             ACRL1007
      READ (K,48) LPR3,LCHECK,LOREAD,IA,IB,IC,(TMP1(I),I=1,10)          ACRL1008
      IF (IA.LE.10) GO TO 6                                             ACRL1009
      WRITE (L,49) LPR3,LCHECK,LOREAD,IA,IB,IC,(TMP1(I),I=1,10)         ACRL1010
      READ (K,44) (TMP1(I),I=11,IA)                                     ACRL1011
      WRITE (L,82) (TMP1(I),I=11,IA)                                    ACRL1012
      GO TO 7                                                           ACRL1013
   5  READ (K,48) LPR3,LCHECK,LOREAD,IA,IB,IC,(TMP1(I),I=1,IA)          ACRL1014
   6  WRITE (L,49) LPR3,LCHECK,LOREAD,IA,IB,IC,(TMP1(I),I=1,IA)         ACRL1015
   7  J=IB-JMIN+1                                                       ACRL1016
      JY=0                                                              ACRL1017
      IF (LPR3) JY=-1                                                   ACRL1018
      J1=IB                                                             ACRL1019
      DO 8 I=1,IA                                                       ACRL1020
   8  TMP4(I,J)=TMP1(I)                                                 ACRL1021
      RETURN                                                            ACRL1022
   9  IK=IA                                                             ACRL1023
      JY=JY+1                                                           ACRL1024
      IF (LOREAD) GO TO 10                                              ACRL1025
C        THIS READS AND WRITES THE ROW REORDERING SEQUENCE IF REQUIRED. ACRL1026
C                                                                       ACRL1027
      READ (K,44) (TMP1(I),I=1,IK)                                      ACRL1028
      WRITE (L,83) (TMP1(I),I=1,IK)                                     ACRL1029
  10  IF (IB.GE.2) LZ=.TRUE.                                            ACRL1030
      IF ((IB.EQ.0).OR.(IB.EQ.3)) GO TO 13                              ACRL1031
C         THIS READS AND WRITES AN R J-BLOCK FOR NOPT=1 OR NOPT=2.      ACRL1032
C                                                                       ACRL1033
      DO 12 JK=1,IK                                                     ACRL1034
      IF (LZ) GO TO 11                                                  ACRL1035
      READ (K,51) (TMP4(I,JK),I=1,IK)                                   ACRL1036
      WRITE (L,52) (TMP4(I,JK),I=1,IK)                                  ACRL1037
      GO TO 12                                                          ACRL1038
  11  READ (K,51) (TMP3(I,JK),I=1,IK)                                   ACRL1039
      WRITE (L,52) (TMP3(I,JK),I=1,IK)                                  ACRL1040
  12  CONTINUE                                                          ACRL1041
      GO TO 15                                                          ACRL1042
C         THIS READS AND WRITES AN R J-BLOCK FOR NOPT=0 OR NOPT=3.      ACRL1043
C                                                                       ACRL1044
  13  IF (LZ) GO TO 14                                                  ACRL1045
      READ (K,51) ((TMP4(I,JK),I=1,IK),JK=1,IK)                         ACRL1046
      WRITE (L,52) ((TMP4(I,JK),I=1,IK),JK=1,IK)                        ACRL1047
      GO TO 15                                                          ACRL1048
  14  READ (K,51) ((TMP3(I,JK),I=1,IK),JK=1,IK)                         ACRL1049
      WRITE (L,52) ((TMP3(I,JK),I=1,IK),JK=1,IK)                        ACRL1050
  15  JINT=0                                                            ACRL1051
      IF (.NOT.LCHECK) GO TO 18                                         ACRL1052
C        THIS SECTION CHECKS THE SYMMETRY OF THE J-BLOCK IF REQUIRED.   ACRL1053
C                                                                       ACRL1054
      DO 17 I=1,IK                                                      ACRL1055
      DO 16 JK=I,IK                                                     ACRL1056
      IF (JK.EQ.I) GO TO 16                                             ACRL1057
      R=TMP4(I,JK)                                                      ACRL1058
      Q=TMP4(JK,I)                                                      ACRL1059
      IF (LZ) R=TMP3(I,JK)                                              ACRL1060
      IF (LZ) Q=TMP3(JK,I)                                              ACRL1061
      A=ABS(AMIN1(R,Q))                                                 ACRL1062
C        NORMALIZE THE SIZE OF THE ELEMENTS TO ABOUT ONE.               ACRL1063
      IF (A.LE.1.E-5) TAU=1.E+5                                         ACRL1064
      IF (A.GT.1.E-5) TAU=1.E+4                                         ACRL1065
      IF (A.GT.1.E-4) TAU=1.E+3                                         ACRL1066
      IF (A.GT.1.E-3) TAU=1.E+2                                         ACRL1067
      IF (A.GT.1.E-2) TAU=1.E+1                                         ACRL1068
      IF (A.GT.1.E-1) TAU=1.E+0                                         ACRL1069
      IF (A.GT.1.E+0) TAU=1.E-1                                         ACRL1070
      IF (A.GT.1.E+1) TAU=1.E-2                                         ACRL1071
      IF (A.GT.1.E+2) TAU=1.E-3                                         ACRL1072
      IF (ABS(R-Q)*TAU.LT.TMP2(1)) GO TO 16                             ACRL1073
C       THIS IS ENTERED WHEN AN ERROR IS FOUND                          ACRL1074
C                                                                       ACRL1075
C          THIS WRITES THE SYMMETRY CHECK ERROR MESSAGE.                ACRL1076
C                                                                       ACRL1077
      IF (JINT.EQ.0) WRITE (L,53)                                       ACRL1078
      JINT=JINT+1                                                       ACRL1079
C         PRINT THE ERROR AND SET LABORT .TRUE.                         ACRL1080
C              THIS WRITES THE SUBSCRIPT AND THE NONSYMMETRIC ELEMENTS. ACRL1081
C                                                                       ACRL1082
      WRITE (L,54) I,JK,Q,R                                             ACRL1083
      LABORT=.TRUE.                                                     ACRL1084
  16  CONTINUE                                                          ACRL1085
  17  CONTINUE                                                          ACRL1086
  18  IF (.NOT.LZ.OR.LABORT) GO TO 20                                   ACRL1087
      CALL TRNSFR (IK,TMP1,TMP3,TMP4,ID,IE,IG)                          ACRL1088
C       THIS WRITES END OF DATA, THE REORDERED J-BLOCK MESSAGE, AND     ACRL1089
C                      THE TOTAL ANGULAR MOMENTUM.                      ACRL1090
C                                                                       ACRL1091
      WRITE (L,56)                                                      ACRL1092
      IF (.NOT.LPR3) WRITE (L,57) IC                                    ACRL1093
      IF (LPR3) WRITE (L,58) IC                                         ACRL1094
C          THIS WRITES THE REORDERED J-BLOCK.                           ACRL1095
C                                                                       ACRL1096
      DO 19 JK=1,IK                                                     ACRL1097
      WRITE (L,55) (TMP4(I,JK),I=1,IK)                                  ACRL1098
  19  CONTINUE                                                          ACRL1099
C           THIS WRITES THE DATA HEADING ON THE NEXT LINE.              ACRL1100
C                                                                       ACRL1101
  20  IF ((LZ.AND.(.NOT.LABORT).OR.(JINT.NE.0)).AND.((J1.NE.JTOT).OR.(  ACRL1102
     1  JY.NE.1))) WRITE (L,59)                                         ACRL1103
C                                                                       ACRL1104
C            THIS WRITES END OF DATA.                                   ACRL1105
      IF ((.NOT.LZ.AND.(JINT.EQ.0)).AND.((J1.EQ.JTOT).AND.(JY.EQ.1)))   ACRL1106
     1  WRITE (L,56)                                                    ACRL1107
      RETURN                                                            ACRL1108
C       THIS WRITES N, J1, THE IMAGINARY PART OF THE T J-BLOCK, AND THE ACRL1109
C             REAL PART OF THE T J-BLOCK ON THE EXTERNAL DEVICE.        ACRL1110
C                                                                       ACRL1111
  21  WRITE (M) IA,IB,((TMP3(I,JK),I=1,IA),JK=1,IA),((TMP4(I,JK),I=1,IA)ACRL1112
     1,JK=1,IA)                                                         ACRL1113
      RETURN                                                            ACRL1114
C       THIS READS N, J1, THE IMAGINARY PART OF THE T J-BLOCK, AND THE  ACRL1115
C             REAL PART OF THE T J-BLOCK FROM THE EXTERNAL DEVICE.      ACRL1116
C                                                                       ACRL1117
  22  READ (M) IA,IB,((TMP3(I,JK),I=1,IA),JK=1,IA),((TMP4(I,JK),I=1,IA),ACRL1118
     1JK=1,IA)                                                          ACRL1119
      RETURN                                                            ACRL1120
C      THIS WRITES THE ROW ERROR MESSAGE AND THE TOTAL ANGULAR MOMENTUM.ACRL1121
C                                                                       ACRL1122
  23  WRITE (L,50) IA                                                   ACRL1123
      RETURN                                                            ACRL1124
  24  T1=TMP1(1)                                                        ACRL1125
      T2=TMP2(1)                                                        ACRL1126
      T3=TMP3(IA,1)                                                     ACRL1127
      T4=TMP4(IA,1)                                                     ACRL1128
C    THIS WRITES THE CURRENT VALUES OF S, SP, ALPHA(A), AND ALPHAP(AP). ACRL1129
C                                                                       ACRL1130
      WRITE (L,60) T3,T4,T1,T2                                          ACRL1131
      RETURN                                                            ACRL1132
C       THIS WRITES THE INITIAL KINETIC ENERGY IN VARIOUS SETS OF UNITS,ACRL1133
C               THE FACTOR ((2S+1)K2)**-1, AND THE RELATIVE VELOCITY.   ACRL1134
C                                                                       ACRL1135
  25  WRITE (L,61) TMP1(1),TMP2(1),TMP3(1,1),TMP4(1,1),IA               ACRL1136
      WRITE (L,62) IB,IC                                                ACRL1137
      RETURN                                                            ACRL1138
C        THIS WRITES THE TITLE LINE FOR THE EVEN AND ODD PARITY CONTRI- ACRL1139
C             BUTIONS TO THE TRANSITION PROBABILITY AND THE OPACITY.    ACRL1140
C                                                                       ACRL1141
  26  WRITE (L,63) T1,T3,T2,T4,T1,T3                                    ACRL1142
C       THIS WRITES THE RELATIVE ORBITAL ANGULAR MOMENTUM L, PROBL(L+1),ACRL1143
C               PROBLO(L+1), OPAC(L+1), AND OPACO(L+1).                 ACRL1144
C                                                                       ACRL1145
      DO 27 I=1,IA                                                      ACRL1146
      J=I-1                                                             ACRL1147
      WRITE (L,64) J,TMP1(I),TMP2(I),TMP3(I,1),TMP4(I,1)                ACRL1148
  27  CONTINUE                                                          ACRL1149
      RETURN                                                            ACRL1150
C       THIS WRITES THE TITLE LINE FOR THE TRANSITION PROBABILITY, THE  ACRL1151
C                 OPACITY, AND THE EVEN AND ODD PARITY CONTRIBUTIONS    ACRL1152
C           TO THE PARTIAL SUM OF THE PARTIAL INTEGRAL CROSS SECTIONS.  ACRL1153
C                                                                       ACRL1154
  28  IF ((IB.NE.0).AND.(IB.NE.4)) WRITE(L,65) T1,T3,T2,T4,T1,T3        ACRL1155
      IF ((IB.EQ.0).OR.(IB.EQ.4)) WRITE (L,66) T1,T3,T2,T4,T1,T3        ACRL1156
C         THIS WRITES THE RELATIVE ORBITAL ANGULAR MOMENTUM L, ITS TRAN-ACRL1157
C                SITION PROBABILITY AND OPACITY, AND THE EVEN AND ODD   ACRL1158
C       PARITY CONTRIBUTIONS TO THE PARTIAL INTEGRAL CROSS SECTION SUM. ACRL1159
C                                                                       ACRL1160
      DO 29 I=1,IA                                                      ACRL1161
      J=I-1                                                             ACRL1162
      R=TMP3(I,1)+TMP4(I,1)                                             ACRL1163
      WRITE (L,67) J,TMP1(I),TMP2(I),R,TMP3(I,1),TMP4(I,1)              ACRL1164
  29  CONTINUE                                                          ACRL1165
      RETURN                                                            ACRL1166
C   THIS WRITES THE BL COEFFICIENTS OF EQ.(4.5) AS CORRECTED BY HUBY.   ACRL1167
C                                                                       ACRL1168
  30  WRITE (L,68)                                                      ACRL1169
      WRITE (L,69) ((TMP1(I),TMP2(I)),I=1,IA)                           ACRL1170
      WQ=TMP3(1,1)                                                      ACRL1171
      RETURN                                                            ACRL1172
C     THIS WRITES THE INTEGRAL AND DIFFERENTIAL CROSS SECTIONS AND IF   ACRL1173
C                NECESSARY THE VALUE FOUND BY EQ.(3.16).                ACRL1174
C                                                                       ACRL1175
  31  IF (IB.EQ.9) WRITE (L,70)                                         ACRL1176
      IF ((IB.GT.0).AND.(IB.LT.4)) GO TO 32                             ACRL1177
      WRITE (L,71) TMP1(1),TMP3(1,1)                                    ACRL1178
      WRITE (L,72) IA,IA                                                ACRL1179
      WRITE (L,73)                                                      ACRL1180
      GO TO 33                                                          ACRL1181
  32  WRITE (L,74) TMP1(1),TMP3(1,1)                                    ACRL1182
      WRITE (L,72) IA,IA                                                ACRL1183
      WRITE (L,75) WQ                                                   ACRL1184
  33  WRITE (L,76) ((TMP2(I),TMP4(I,IC)),I=1,181)                       ACRL1185
      RETURN                                                            ACRL1186
C   THIS WRITES DATA ON A NEW PAGE AND READS AND WRITES THE NEW JAPMAX. ACRL1187
C                                                                       ACRL1188
  34  WRITE (L,39)                                                      ACRL1189
      READ (K,44) IA                                                    ACRL1190
      WRITE (L,80) IA                                                   ACRL1191
      GO TO 2                                                           ACRL1192
C            THIS WRITES THE ELAPSED TIME.                              ACRL1193
C                                                                       ACRL1194
  35  WRITE (L,77) TMP1(1)                                              ACRL1195
      RETURN                                                            ACRL1196
C                                                                       ACRL1197
  39  FORMAT (1H1,16H DATA RECORD NO.,3X,4HDATA,25X,69HSEE SECTION 2.4 OACRL1198
     1F THE TEXT FOR A DESCRIPTION OF THE DATA RECORD NOS.  )           ACRL1199
  40  FORMAT (8A10)                                                     ACRL1200
  41  FORMAT (9X,1H1,10X,8A10)                                          ACRL1201
  42  FORMAT (14L5)                                                     ACRL1202
  43  FORMAT (9X,1H2,10X,14L5)                                          ACRL1203
  44  FORMAT (14I5)                                                     ACRL1204
  45  FORMAT (9X,1H3,10X,14I5)                                          ACRL1205
  46  FORMAT (3I5,4E13.6)                                               ACRL1206
  47  FORMAT (9X,1H6,10X,3I5,4E13.6)                                    ACRL1207
  48  FORMAT (L1,2L2,13I5)                                              ACRL1208
  49  FORMAT (9X,1H7,10X,L1,2L2,13I5)                                   ACRL1209
  50  FORMAT (1H0,15(1H*),45H ROW LABEL ERROR FOR TOTAL ANGULAR MOMENTUMACRL1210
     1 =,I3,4X,37HCALCULATION STOPPED FOR THIS DATA SET)                ACRL1211
  51  FORMAT (5E13.6)                                                   ACRL1212
  52  FORMAT (9X,1H9,10X,5F13.6)                                        ACRL1213
  53  FORMAT (2X,11HEND OF DATA,15(1H*),100HTHE ABOVE MATRIX DATA CONTAIACRL1214
     1NS THE FOLLOWING ERRORS***THEREFORE THIS CALCULATION WAS ABORTED. ACRL1215
     2      /)                                                          ACRL1216
  54  FORMAT (15X,2I10,2E18.8)                                          ACRL1217
  55  FORMAT (1X,9F13.6)                                                ACRL1218
  56  FORMAT (1X,12H END OF DATA )                                      ACRL1219
  57  FORMAT (77H0  THE PROPERLY ORDERED J-BLOCK OF THE R MATRIX FOR TOTACRL1220
     1AL ANGULAR MOMENTUM = ,I3/)                                       ACRL1221
  58  FORMAT (77H0  THE PROPERLY ORDERED J-BLOCK OF THE T MATRIX FOR TOTACRL1222
     1AL ANGULAR MOMENTUM = ,I3/)                                       ACRL1223
  59  FORMAT (1H0,16H DATA RECORD NO.,3X,4HDATA)                        ACRL1224
  60  FORMAT (1H0,6X,12(1H*,1X),30H  INITIAL ANGULAR MOMENTUM S =,I3,5X,ACRL1225
     127HFINAL ANGULAR MOMENTUM SP =,I3,2X,12(1X,1H*)/58X,3HA =,I3,28X,4ACRL1226
     2HAP =,I3)                                                         ACRL1227
  61  FORMAT (15H0INITIAL K.E. =,E13.6,1X,4HA.U.,2X,1H=,E13.6,1X,3HEV.,2ACRL1228
     1X,1H=,E13.6,1X,12HERG/MOLECULE,2X,1H=,E13.6,1X,9HKCAL/MOLE,2X,1H=,ACRL1229
     2E13.6,1X,8H(CM)**-1/)                                             ACRL1230
  62  FORMAT (34X,1H2,1X,2H-1/26X,15H((2S+1)K )   = ,E13.6,1X,3HA02,20X,ACRL1231
     110HVELOCITY =,E13.6,1X,6HCM/SEC)                                  ACRL1232
  63  FORMAT (34H0RELATIVE ORBITAL ANGULAR MOMENTUM,3X,30HTRANSITION PROACRL1233
     1BABILITY FROM A=,I2,3H,S=,I2,7H TO AP=,I2,4H,SP=,I2,3X,30HTRANSITIACRL1234
     2ON PROBABILITY FROM A=,I2,3H,S=,I2/17X,1HL,7X,2(4X,24HEVEN PARITY ACRL1235
     3CONTRIBUTION,3X,23HODD PARITY CONTRIBUTION)/)                     ACRL1236
  64  FORMAT (14X,I4,14X,E16.6,11X,E15.6,12X,E16.6,11X,E15.6)           ACRL1237
  65  FORMAT (4H0  L,4X,20HTRANS. PROB. FROM A=,I2,3H,S=,I2,7H TO AP=,I2ACRL1238
     1,4H,SP=,I2,4X,20HTRANS. PROB. FROM A=,I2,3H,S=,I2,3X,49HPARTIAL INACRL1239
     2TER. INTEGRAL CR. SECTION SUMMED THRU L/23X,10HEVEN + ODD,31X,10HEACRL1240
     3VEN + ODD,16X,10HEVEN + ODD,7X,4HEVEN,10X,3HODD/)                 ACRL1241
  66  FORMAT (4H0  L,4X,20HTRANS. PROB. FROM A=,I2,3H,S=,I2,7H TO AP=,I2ACRL1242
     1,4H,SP=,I2,4X,20HTRANS. PROB. FROM A=,I2,3H,S=,I2,6X,44HPARTIAL INACRL1243
     2TEGRAL CROSS SECTION SUMMED THRU L /23X,10HEVEN + ODD,31X,10HEVEN ACRL1244
     3+ ODD,16X,10HEVEN + ODD,7X,4HEVEN,10X,3HODD/)                     ACRL1245
  67  FORMAT (1X,I3,14X,E16.6,26X,E15.6,13X,3(E13.6,1X))                ACRL1246
  68  FORMAT (/10X,111HTHE BL COEFFICIENTS OF THE LEGENDRE POLYNOMIALS IACRL1247
     1N THE DCS FORMULA OF BLATT AND BIEDENHARN AS CORRECTED BY HUBY/)  ACRL1248
  69  FORMAT (1X,7(I3,1X,E13.6,1X))                                     ACRL1249
  70  FORMAT (/3X,17(1H*),86H  SUM OF THE PROPERLY WEIGHTED INTERMEDIATEACRL1250
     1 INTEGRAL AND DIFFERENTIAL CROSS SECTIONS  ,17(1H*))              ACRL1251
  71  FORMAT (1H0,24X,35HINTEGRATED INTEGRAL CROSS SECTION =,E13.6,1X,3HACRL1252
     1A02,5X,11H4.0*PI*B0 =,E13.6,1X,3HA02)                             ACRL1253
  72  FORMAT (45X,2(14X,1H=,E13.6,1X,4HANG2))                           ACRL1254
  73  FORMAT (/45X,35HDIFFERENTIAL CROSS SECTION (A02/SR)/)             ACRL1255
  74  FORMAT (1H0,11X,48HINTEGRATED INTERMEDIATE INTEGRAL CROSS SECTION ACRL1256
     1=,E13.6,1X,3HA02,5X,11H4.0*PI*B0 =,E13.6,1X,3HA02)                ACRL1257
  75  FORMAT (/21X,76HINTERMEDIATE DIFFERENTIAL CROSS SECTION (A02/SR) IACRL1258
     1NCLUDING THE WEIGHT FACTOR,E14.6/)                                ACRL1259
  76  FORMAT (1X,8(1X,I3,E11.4,1X))                                     ACRL1260
  77  FORMAT (/1X,40HTOTAL ELAPSED TIME(SECONDS) IN CROSS =  ,E14.6/)   ACRL1261
  78  FORMAT (9X,1H5,10X,8A10)                                          ACRL1262
  79  FORMAT (8X,2H10,10X,8A10)                                         ACRL1263
  80  FORMAT (8X,2H12,10X,14I5)                                         ACRL1264
  81  FORMAT (9X,1H4,10X,14I5)                                          ACRL1265
  82  FORMAT (9X,1H7,10X,14I5)                                          ACRL1266
  83  FORMAT (9X,1H8,10X,14I5)                                          ACRL1267
  84  FORMAT (8X,2H11,10X,3I5,4E13.6)                                   ACRL1268
      END                                                               ACRL1269
      SUBROUTINE TRNSFR (NQ,KSTOR,H,R,NMAX,KMAX,MMAX)                   ACRL1270
C                                                                       ACRL1271
C         TRNSFR IS CALLED BY SUBPROGRAM INOUT WHICH MUST RETURN TO     ACRL1272
C    SUBPROGRAM CROSS A J-BLOCK WITH ITS ROWS IN THE CORRECT ORDER.     ACRL1273
C                                                                       ACRL1274
C          SUBPROGRAM CROSS REQUIRES THAT THE ROWS OF R OR TREAL BE     ACRL1275
C     ARRANGED IN ORDER OF ASCENDING RELATIVE ORBITAL ANGULAR MOMENTUM  ACRL1276
C     FOR EACH CHANNEL ANGULAR MOMENTUM AND TOTAL ANGULAR MOMENTUM.     ACRL1277
C                                                                       ACRL1278
      DIMENSION KSTOR(NMAX),H(NMAX,NMAX),R(KMAX,MMAX)                   ACRL1279
C                                                                       ACRL1280
C                         PARAMETER LIST                                ACRL1281
C                                                                       ACRL1282
C   NQ  -  AN INPUT VARIABLE  -  THE ROW SIZE OF THE CURRENT J-BLOCK    ACRL1283
C                                                                       ACRL1284
C   KSTOR  -  AN INPUT ARRAY  -  THE CORRECT ARRAY SUBSCRIPT SEQUENCE   ACRL1285
C                                                                       ACRL1286
C   H  -  AN INPUT ARRAY  -  THE J-BLOCK BEFORE ITS ROWS WERE REORDERED ACRL1287
C                                                                       ACRL1288
C   R  -  AN OUTPUT ARRAY  -  THE J-BLOCK WITH ITS ROWS REORDERED       ACRL1289
C                                                                       ACRL1290
C   NMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1291
C                                           STARTING ON CARD NO.191     ACRL1292
C   KMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1293
C   MMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1294
C                                                                       ACRL1295
      N=NQ                                                              ACRL1296
      DO 1 I=1,N                                                        ACRL1297
      II=KSTOR(I)                                                       ACRL1298
      DO 1 JK=1,N                                                       ACRL1299
      III=KSTOR(JK)                                                     ACRL1300
   1  R(II,III)=H(I,JK)                                                 ACRL1301
      RETURN                                                            ACRL1302
      END                                                               ACRL1303
      SUBROUTINE MXMPLY (A,B,C,IRA,JCA,JCB,IDA,IDB,IDC,TEMP)            ACRL1304
C                                                                       ACRL1305
C       MULTIPLIES A(IRA BY JCA) TIMES B(JCA BY JCB) AND RETURNS THE    ACRL1306
C               PRODUCT IN C(IRA BY JCB)                                ACRL1307
C       WHEN CALLING MXMPLY, C MAY BE THE SAME AS A BUT NOT THE SAME    ACRL1308
C         AS B.  B MAY BE THE SAME AS A IF C IS NOT THE SAME AS A.      ACRL1309
C                                                                       ACRL1310
C       TEMP IS A SCRATCH ARRAY (OF DIMENSION AT LEAST JCA) WHICH IS    ACRL1311
C         NECESSARY BECAUSE OF THE OPTIONS ON CARD NOS. 1308-1309.      ACRL1312
C                                                                       ACRL1313
      DIMENSION A(IDA,JCA),B(IDB,JCB),C(IDC,JCB),TEMP(JCA)              ACRL1314
      DO 2 K=1,IRA                                                      ACRL1315
      DO 1 I=1,JCB                                                      ACRL1316
      TEMP(I)=0.0                                                       ACRL1317
      DO 1 J=1,JCA                                                      ACRL1318
   1  TEMP(I)=A(K,J)*B(J,I)+TEMP(I)                                     ACRL1319
      DO 2 I=1,JCB                                                      ACRL1320
   2  C(K,I)=TEMP(I)                                                    ACRL1321
      RETURN                                                            ACRL1322
      END                                                               ACRL1323
      SUBROUTINE MXLNEQ (A,NN,IDA,DET,JRANK,EPS,IN)                     ACRL1324
C                                                                       ACRL1325
C        THIS IS A MODIFIED VERSION OF A SUBPROGRAM WRITTEN FOR THE     ACRL1326
C    UNIVERSITY COMPUTER CENTER OF THE UNIVERSITY OF MINNESOTA BY       ACRL1327
C    RICHARD L. HOTCHKISS AND MICHAEL J. FRISCH.                        ACRL1328
C        THIS SUBPROGRAM USES THE GAUSS-JORDAN METHOD TO FIND THE DETER-ACRL1329
C    MINANT AND INVERSE OF THE ARRAY A(NN BY NN).  NN MUST BE LESS THAN ACRL1330
C    OR EQUAL TO THE VARIABLE DIMENSION IDA.                            ACRL1331
C                                                                       ACRL1332
C       THE EXPLANATION OF THE /IOSEND/ VARIABLES BEGINS ON CARD NO.19. ACRL1333
      COMMON /IOSEND/ IREAD,IWRITE,ITAPE                                ACRL1334
      DIMENSION A(IDA,IDA),IN(IDA)                                      ACRL1335
C                                                                       ACRL1336
C                            PARAMETER LIST                             ACRL1337
C                                                                       ACRL1338
C   A  -  AN INPUT/OUTPUT ARRAY  -  THE MATRIX WHOSE INVERSE IS TO BE   ACRL1339
C                           TAKEN.  THE INVERSE OF A IS STORED IN A.    ACRL1340
C                                                                       ACRL1341
C   NN  -  AN INPUT VARIABLE  -  NUMBER OF ROWS OF A CURRENTLY IN USE   ACRL1342
C                                                                       ACRL1343
C   IDA  -  A VARIABLE DIMENSION  -  THE ROW DIMENSION OF THE ARRAY A   ACRL1344
C                               AND THE DIMENSION OF THE ARRAY IN       ACRL1345
C                                                                       ACRL1346
C   DET  -  AN OUTPUT VARIABLE  -  THE DETERMINANT OF THE ARRAY A       ACRL1347
C                                                                       ACRL1348
C   JRANK  -  AN OUTPUT VARIABLE  -  THE RANK OF THE ARRAY A            ACRL1349
C                                                                       ACRL1350
C   EPS  -  AN INPUT VARIABLE  -  1.E-N WHERE N IS THE NUMBER OF        ACRL1351
C                         SIGNIFICANT DIGITS FOR THE MACHINE IN USE     ACRL1352
C                                                                       ACRL1353
C   IN  -  A SCRATCH ARRAY                                              ACRL1354
C                                                                       ACRL1355
      DATA MAD,MMM/1,0/                                                 ACRL1356
      N=NN                                                              ACRL1357
      JRANK=N                                                           ACRL1358
      DET=1.0                                                           ACRL1359
      NFLAG=-1                                                          ACRL1360
      K1=1                                                              ACRL1361
      NM=NN                                                             ACRL1362
      DO 18 K=1,N                                                       ACRL1363
      PIV=0.0                                                           ACRL1364
      DO 2 I=K,N                                                        ACRL1365
      DO 2 J=K,N                                                        ACRL1366
      P=ABS(A(I,J))                                                     ACRL1367
      IF (P-PIV) 2,2,1                                                  ACRL1368
   1  PIV=P                                                             ACRL1369
      L=I                                                               ACRL1370
      M=J                                                               ACRL1371
   2  CONTINUE                                                          ACRL1372
      IF (NFLAG) 3,5,5                                                  ACRL1373
   3  IF (PIV-EPS) 4,4,7                                                ACRL1374
   4  JRANK=K-1                                                         ACRL1375
      WRITE (IWRITE,25) JRANK,N,EPS                                     ACRL1376
      NFLAG=0                                                           ACRL1377
   5  IF (PIV) 6,6,7                                                    ACRL1378
   6  JRANK=-K                                                          ACRL1379
      DET=0.                                                            ACRL1380
      WRITE (IWRITE,26)                                                 ACRL1381
      CALL EXIT                                                         ACRL1382
   7  PIVOT=A(L,M)                                                      ACRL1383
      DET=PIVOT*DET                                                     ACRL1384
      IN(K)=512*L+M                                                     ACRL1385
      IF (L-K) 10,10,8                                                  ACRL1386
   8  DET=-DET                                                          ACRL1387
      DO 9 J=K1,NM                                                      ACRL1388
      Z=A(L,J)                                                          ACRL1389
      A(L,J)=A(K,J)                                                     ACRL1390
   9  A(K,J)=Z                                                          ACRL1391
  10  IF (M-K) 13,13,11                                                 ACRL1392
  11  DET=-DET                                                          ACRL1393
      DO 12 I=1,N                                                       ACRL1394
      Z=A(I,M)                                                          ACRL1395
      A(I,M)=A(I,K)                                                     ACRL1396
  12  A(I,K)=Z                                                          ACRL1397
  13  PIVOT=1.0/PIVOT                                                   ACRL1398
      DO 14 J=K1,NM                                                     ACRL1399
  14  A(K,J)=PIVOT*A(K,J)                                               ACRL1400
      A(K,K)=0.0                                                        ACRL1401
      DO 17 I=K1,N                                                      ACRL1402
      Z=A(I,K)                                                          ACRL1403
      IF (Z) 15,17,15                                                   ACRL1404
  15  DO 16 J=K1,NM                                                     ACRL1405
  16  A(I,J)=-Z*A(K,J)+A(I,J)                                           ACRL1406
      A(I,K)=-Z*PIVOT                                                   ACRL1407
  17  CONTINUE                                                          ACRL1408
  18  A(K,K)=PIVOT                                                      ACRL1409
      IF (N.EQ.1) RETURN                                                ACRL1410
      K=N                                                               ACRL1411
      DO 24 J=2,N                                                       ACRL1412
      K=K-1                                                             ACRL1413
      L=MOD(IN(K),512)                                                  ACRL1414
      M=IN(K)/512                                                       ACRL1415
      IF (M-K) 21,21,19                                                 ACRL1416
  19  DO 20 I=1,N                                                       ACRL1417
      Z=A(I,K)                                                          ACRL1418
      A(I,K)=A(I,M)                                                     ACRL1419
  20  A(I,M)=Z                                                          ACRL1420
  21  IF (L-K) 24,24,22                                                 ACRL1421
  22  DO 23 I=K1,NM                                                     ACRL1422
      Z=A(K,I)                                                          ACRL1423
      A(K,I)=A(L,I)                                                     ACRL1424
  23  A(L,I)=Z                                                          ACRL1425
  24  CONTINUE                                                          ACRL1426
      RETURN                                                            ACRL1427
C                                                                       ACRL1428
  25  FORMAT (42H0****  MXLNEQ HAS FOUND ARGUMENT 5 JRANK =,I5,23H WHILEACRL1429
     1 ARGUMENT 2 NN = ,I5,22H FOR ARGUMENT 6 EPS = ,E20.10/30H  SOLUTIOACRL1430
     2N IS CONTINUED  **** )                                            ACRL1431
  26  FORMAT (68H0****  MXLNEQ HAS FOUND A ZERO DETERMINANT.  SOLUTION SACRL1432
     1TOPPED  **** )                                                    ACRL1433
      END                                                               ACRL1434
      SUBROUTINE OPACIT (JQ,JXQ,JPARQ,NVALUQ,NQ,KKQ,JSQ,JINDEX,G,R,OPAC,ACRL1435
     1  OPACO,NMAX,KMAX,MMAX,JMAXP)                                     ACRL1436
C                                                                       ACRL1437
C         THIS SUBPROGRAM CALCULATES THE OPACITY(TRANSITION PROBABILITY ACRL1438
C    OUT OF THE INITIAL STATE) FOR EACH L.                              ACRL1439
C                                                                       ACRL1440
      LOGICAL LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LT    ACRL1441
      COMMON/LCOM/LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LTACRL1442
C                                                                       ACRL1443
C         THE ONLY /LCOM/ VARIABLE USED BY SUBPROGRAM OPACIT IS LT WHICHACRL1444
C    IS EXPLAINED IN SUBPROGRAM PARITY. (SEE CARD NO.1514)              ACRL1445
C                                                                       ACRL1446
      DIMENSION JINDEX(NMAX),G(KMAX,MMAX),R(KMAX,MMAX),OPAC(JMAXP),     ACRL1447
     1          OPACO(JMAXP)                                            ACRL1448
C                                                                       ACRL1449
C                         PARAMETER LIST                                ACRL1450
C                                                                       ACRL1451
C   JQ  -  AN INPUT VARIABLE   -  THE INITIAL CHANNEL ANGULAR MOMENTUM  ACRL1452
C                                                                       ACRL1453
C   JXQ -  AN INPUT VARIABLE   -  THE PRESENT TOTAL ANGULAR MOMENTUM    ACRL1454
C                                                                       ACRL1455
C   JPARQ  -  AN INPUT VARIABLE  -  EXPLAINED IN CROSS (SEE CARD NO.337)ACRL1456
C                                                                       ACRL1457
C   NVALUQ  -  AN INPUT VARIABLE  -  100*(INITIAL STATE NON-ANGULAR     ACRL1458
C                                                  QUANTUM NUMBER)      ACRL1459
C                                                                       ACRL1460
C   NQ  -  AN INPUT VARIABLE  -  THE ROW SIZE OF THE CURRENT J-BLOCK    ACRL1461
C                                                                       ACRL1462
C   KKQ  -  AN INPUT VARIABLE  -  THE FIRST ROW OF THE INITIAL STATE    ACRL1463
C                                                                       ACRL1464
C   JSQ  -  AN INPUT VARIABLE  -  CHANNEL PARITY OF THE INITIAL STATE   ACRL1465
C                                                                       ACRL1466
C   JINDEX  -  AN INPUT ARRAY  -  THE J-BLOCK ROW LABEL                 ACRL1467
C                                                                       ACRL1468
C   G  -  AN INPUT ARRAY  -  THE REAL PART OF THE T MATRIX J-BLOCK      ACRL1469
C                                                                       ACRL1470
C   R  -  AN INPUT ARRAY  -  THE IMAGINARY PART OF THE T MATRIX J-BLOCK ACRL1471
C                                                                       ACRL1472
C   OPAC  -  AN OUTPUT ARRAY  -  EXPLAINED IN CROSS (SEE CARD NO.233)   ACRL1473
C                                                                       ACRL1474
C   OPACO  -  AN OUTPUT ARRAY  -  EXPLAINED IN CROSS (SEE CARD NO.237)  ACRL1475
C                                                                       ACRL1476
C   NMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1477
C   KMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1478
C   MMAX  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS     ACRL1479
C   JMAXP  -  A VARIABLE DIMENSION  -  EXPLAINED IN SUBPROGRAM CROSS    ACRL1480
C                                           STARTING AT CARD NO.187     ACRL1481
      J=JQ                                                              ACRL1482
      JPAR=JPARQ                                                        ACRL1483
      JX=JXQ                                                            ACRL1484
      NVALU=NVALUQ                                                      ACRL1485
      N=NQ                                                              ACRL1486
      KK=KKQ                                                            ACRL1487
      JS=JSQ                                                            ACRL1488
C                                                                       ACRL1489
C        FIND THE MAX AND MIN OF THE RELATIVE ANGULAR MOMENTUM + 1.     ACRL1490
      IF (JX-J) 1,2,1                                                   ACRL1491
   1  L1MIN=IABS(JX-J)+1                                                ACRL1492
      GO TO 3                                                           ACRL1493
   2  L1MIN=1                                                           ACRL1494
   3  L1MAX=JX+J+1                                                      ACRL1495
      DO 6 L1=L1MIN,L1MAX,JPAR                                          ACRL1496
      L2=L1-1                                                           ACRL1497
      CALL PARITY (JS,L2)                                               ACRL1498
      DO 5 LK=1,N                                                       ACRL1499
      IF ((JINDEX(LK)-NVALU).EQ.J) GO TO 5                              ACRL1500
      IF (LT) GO TO 4                                                   ACRL1501
      OPAC(L1)=OPAC(L1)+G(KK,LK)**2+R(KK,LK)**2                         ACRL1502
      GO TO 5                                                           ACRL1503
   4  OPACO(L1)=OPACO(L1)+G(KK,LK)**2+R(KK,LK)**2                       ACRL1504
   5  CONTINUE                                                          ACRL1505
   6  KK=KK+1                                                           ACRL1506
      RETURN                                                            ACRL1507
      END                                                               ACRL1508
      SUBROUTINE PARITY (IA,IB)                                         ACRL1509
C                                                                       ACRL1510
      LOGICAL LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LT    ACRL1511
      COMMON/LCOM/LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LTACRL1512
C                                                                       ACRL1513
C      LT IS THE ONLY /LCOM/ VARIABLE USED BY SUBPROGRAM PARITY.  THE   ACRL1514
C    SUBPROGRAM SETS LT .TRUE. IF THE PARITY IS ODD AND .FALSE. IF THE  ACRL1515
C    PARITY IS EVEN.                                                    ACRL1516
C                                                                       ACRL1517
C        IA AND IB ARE INTEGER INPUT VARIABLES FOR WHICH THE PARITY IS  ACRL1518
C    (-1.)**(IA+IB).                                                    ACRL1519
      LA=IA+IB                                                          ACRL1520
      LB=LA/2                                                           ACRL1521
      LC=LB*2                                                           ACRL1522
      IF (LA-LC) 1,2,1                                                  ACRL1523
   1  LT=.TRUE.                                                         ACRL1524
      RETURN                                                            ACRL1525
   2  LT=.FALSE.                                                        ACRL1526
      RETURN                                                            ACRL1527
      END                                                               ACRL1528
      FUNCTION P(L,X)                                                   ACRL1529
C                                                                       ACRL1530
C       THIS SUBPROGRAM CALCULATES THE LEGENDRE POLYNOMIAL P OF ORDER   ACRL1531
C   L AND ARGUMENT X.  IN THIS PROGRAM X IS COS(THETA).                 ACRL1532
C                                                                       ACRL1533
      LOGICAL LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LT    ACRL1534
      COMMON/LCOM/LPR1,LPR2,LPR3,LOREAD,LSTOR,LEGEND,LZ,LCHECK,LABORT,LTACRL1535
C                                                                       ACRL1536
C         THE ONLY VARIABLE IN /LCOM/ WHICH IS USED IN THIS SUBPROGRAM  ACRL1537
C    IS LEGEND WHICH IS AN INPUT VARIABLE.  LEGEND MAY BE SET .FALSE.   ACRL1538
C    FOR EACH CALL TO THIS SUBPROGRAM.  TO SAVE COMPUTER TIME, LEGEND   ACRL1539
C    SHOULD BE SET .TRUE. IF THE PREVIOUS CALL TO P WAS FOR A VALUE     ACRL1540
C    OF L ONE LESS THAN THE PRESENT VALUE AND FOR THE SAME VALUE OF X.  ACRL1541
C            P IS THE ONLY SUBPROGRAM WHICH USES LEGEND.                ACRL1542
C                                                                       ACRL1543
C                                                                       ACRL1544
      IF (LEGEND) GO TO 1                                               ACRL1545
      P1=X                                                              ACRL1546
      P2=0.5*(3.0*X*X-1.0)                                              ACRL1547
C                                                                       ACRL1548
C       CALCULATE P FOR L EQUALS 0, 1 OR 2.                             ACRL1549
   1  P=1.0                                                             ACRL1550
      IF (L.EQ.0) RETURN                                                ACRL1551
      IF (X.GT.0.99999999) RETURN                                       ACRL1552
      IF (L.GT.1) GO TO 2                                               ACRL1553
      P=X                                                               ACRL1554
      RETURN                                                            ACRL1555
   2  IF (L.GT.2) GO TO 3                                               ACRL1556
      P1=X                                                              ACRL1557
      P2=0.5*(3.0*X*X-1.0)                                              ACRL1558
      P=P2                                                              ACRL1559
      RETURN                                                            ACRL1560
C                                                                       ACRL1561
C      USE RECURSION RELATION FOR L.GT.2                                ACRL1562
   3  K=3                                                               ACRL1563
      IF (LEGEND) K=L                                                   ACRL1564
      DO 4 I=K,L                                                        ACRL1565
      XI=FLOAT(I)                                                       ACRL1566
      P3=((2.*XI-1.)*X*P2-(XI-1.0)*P1)/XI                               ACRL1567
      P1=P2                                                             ACRL1568
      P2=P3                                                             ACRL1569
   4  CONTINUE                                                          ACRL1570
      P=P2                                                              ACRL1571
      RETURN                                                            ACRL1572
      END                                                               ACRL1573
      FUNCTION Z(LA,LB,LC,LD,LL1,LEM)                                   ACRL1574
C                                                                       ACRL1575
C         THIS SUBPROGRAM CALCULATES THE Z COEFFICIENT OF BIEDENHARN,   ACRL1576
C    BLATT, AND ROSE, REV.MOD.PHYS. 24, 249(1952).  THE USER IS         ACRL1577
C    DIRECTED TO THAT PAPER FOR A DESCRIPTION OF THE PARAMETER LIST.    ACRL1578
C    THE PHASE FACTOR, I**(LEM-LA+LC), HAS BEEN REMOVED FROM THE Z COEF-ACRL1579
C    FICIENT DEFINITION AS RECOMMENDED BY R. HUBY. (SEE CARDS NOS. 55   ACRL1580
C    AND 1624)                                                          ACRL1581
C         THE Z COEFFICIENT, Z(LA,LB,LC,LD,LL1,LEM) IN THE REFERENCE,   ACRL1582
C    IS FOUND USING ITS DEFINITION IN TERMS OF WIGNER 6-J AND 3-J COEF- ACRL1583
C    FICIENTS.  THE 6-J COEFFICIENT IS FOUND IN SUBPROGRAM F6J AND THE  ACRL1584
C    3-J COEFFICIENT IS FOUND IN SUBPROGRAM F3J.  BOTH F6J AND F3J USE  ACRL1585
C    LOGARITHMS FOR THE FACTORIALS.  THEY RETURN THE LOGARITHM OF THE   ACRL1586
C    TRIANGLES IN PLOG AND THE 6-J COEFFICIENT SUM IN S.                ACRL1587
C                                                                       ACRL1588
C         SUBPROGRAM Z SHOULD BE CALLED BY SUBPROGRAM CROSS ONLY WHEN   ACRL1589
C    ALL THE TRIANGLE RELATIONS ARE SATISFIED.                          ACRL1590
C                                                                       ACRL1591
      COMMON/FACT/FL(322),NCALL,NFACT                                   ACRL1592
C       THE EXPLANATION OF THE /FACT/ VARIABLES BEGINS ON CARD NO.9.    ACRL1593
      Z=0.0                                                             ACRL1594
      L11=LA+LC+LEM                                                     ACRL1595
      L12=L11/2                                                         ACRL1596
      L1Z=L12*2                                                         ACRL1597
      IF (L11.NE.L1Z) GO TO 11                                          ACRL1598
C     DETERMINE WHETHER TO CALCULATE FL(N) S                            ACRL1599
      IF (NCALL+1867) 1,3,1                                             ACRL1600
   1  NCALL=-1867                                                       ACRL1601
      A=EXP(-64.)                                                       ACRL1602
C     CALCULATE FL(N) S                                                 ACRL1603
      FL(1)=0.0                                                         ACRL1604
      FL(2)=0.0                                                         ACRL1605
      DO 2 N=3,NFACT                                                    ACRL1606
      FN=FLOAT(N-1)                                                     ACRL1607
   2  FL(N)=FL(N-1)+ALOG(FN)                                            ACRL1608
   3  CALL F6J (LA,LB,LL1,LD,LC,LEM,S,PLOG,MIN)                         ACRL1609
      CALL F3J (LA,LC,LEM,PLOG)                                         ACRL1610
      IF (PLOG+64.0) 4,7,7                                              ACRL1611
   4  Q=PLOG+64.                                                        ACRL1612
      Q=EXP(Q)                                                          ACRL1613
      S6J=Q*S                                                           ACRL1614
      IF (ABS(S6J)-1.0) 11,11,6                                         ACRL1615
   6  S6J=S6J*A                                                         ACRL1616
      GO TO 8                                                           ACRL1617
   7  P=EXP(PLOG)                                                       ACRL1618
      S6J=P*S                                                           ACRL1619
   8  MIN2=MIN/2                                                        ACRL1620
      IF (MIN-2*MIN2) 9,10,9                                            ACRL1621
   9  S6J=-S6J                                                          ACRL1622
  10  CONTINUE                                                          ACRL1623
C     LSIGN=LA+LB+LD   --  WITHOUT THE R. HUBY CORRECTION               ACRL1624
      LSIGN=L12+LB+LD                                                   ACRL1625
      DT=(2.*FLOAT(LA)+1.)*(2.*FLOAT(LB)+1.)*(2.*FLOAT(LC)+1.)*         ACRL1626
     1(2.*FLOAT(LD)+1.)*(2.*FLOAT(LEM)+1.)                              ACRL1627
      Z=(-1.)**LSIGN*S6J*SQRT(DT)                                       ACRL1628
  11  RETURN                                                            ACRL1629
      END                                                               ACRL1630
      SUBROUTINE F6J (JD1,JD2,JD3,LD1,LD2,LD3,S,PLOG,MIN)               ACRL1631
C                                                                       ACRL1632
C         THIS SUBPROGRAM IS A MODIFIED FORM OF SUBPROGRAM S6J WRITTEN  ACRL1633
C    BY R.S.CRASWELL AND L.C.MAXIMON, U.S. DEPARTMENT OF COMMERCE       ACRL1634
C    NATIONAL BUREAU OF STANDARDS, NBS TECHNICAL NOTE 409 (1966).       ACRL1635
C                                                                       ACRL1636
      COMMON/FACT/FL(322),NCALL,NFACT                                   ACRL1637
C       THE EXPLANATION OF THE /FACT/ VARIABLES BEGINS ON CARD NO.9.    ACRL1638
C       THE EXPLANATION OF THE /IOSEND/ VARIABLES BEGINS ON CARD NO.19. ACRL1639
      COMMON /IOSEND/ IREAD,IWRITE,ITAPE                                ACRL1640
C                                                                       ACRL1641
C                            PARAMETER LIST                             ACRL1642
C                                                                       ACRL1643
C   JD1  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1644
C                                       (SEE CARD NO.1576)              ACRL1645
C   JD2  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1646
C   JD3  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1647
C   LD1  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1648
C   LD2  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1649
C   LD3  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1650
C                                                                       ACRL1651
C   S  -  AN OUTPUT VARIABLE  -  SUM IN THE 6-J COEFFICIENT DEFINITION  ACRL1652
C   (EQN. 6.3.7 OF EDMONDS, ANGULAR MOMENTUM IN QUANTUM MECHANICS)      ACRL1653
C                                                                       ACRL1654
C   PLOG  -  AN OUTPUT VARIABLE  -  LOGARITHM OF THE TRIANGLES IN THE   ACRL1655
C                    6-J COEFFICIENT DEFINITION(EQN.6.3.7 OF EDMONDS)   ACRL1656
C                                                                       ACRL1657
C   MIN  -  AN OUTPUT VARIABLE  -  USED BY SUBPROGRAM Z TO DETERMINE THEACRL1658
C                                     SIGN OF THE 6-J COEFFICIENT       ACRL1659
C                                                                       ACRL1660
      DIMENSION MA(4),MB(3),MED (12)                                    ACRL1661
      J1=JD1                                                            ACRL1662
      J2=JD2                                                            ACRL1663
      J3=JD3                                                            ACRL1664
      L1=LD1                                                            ACRL1665
      L2=LD2                                                            ACRL1666
      L3=LD3                                                            ACRL1667
      MED(1)=(-J1+J2+J3)                                                ACRL1668
      MED(2)=(+J1-J2+J3)                                                ACRL1669
      MED(3)=(+J1+J2-J3)                                                ACRL1670
      MED(4)=(-J1+L2+L3)                                                ACRL1671
      MED(5)=(+J1-L2+L3)                                                ACRL1672
      MED(6)=(+J1+L2-L3)                                                ACRL1673
      MED(7)=(-L1+J2+L3)                                                ACRL1674
      MED(8)=(+L1-J2+L3)                                                ACRL1675
      MED(9)=(+L1+J2-L3)                                                ACRL1676
      MED(10)=(-L1+L2+J3)                                               ACRL1677
      MED(11)=(+L1-L2+J3)                                               ACRL1678
      MED(12)=(+L1+L2-J3)                                               ACRL1679
C                                                                       ACRL1680
C         THE DO LOOP WHICH FOLLOWS CHECKS THE CONDITIONS NECESSARY FOR ACRL1681
C    A NONZERO RESULT.  THE CONDITIONS SHOULD ALL BE SATISFIED SINCE    ACRL1682
C    THEY HAVE ALL BEEN CHECKED IN SUBPROGRAM CROSS.                    ACRL1683
      DO 1 N=1,12                                                       ACRL1684
      IF (MED(N)) 10,1,1                                                ACRL1685
   1  CONTINUE                                                          ACRL1686
C                                                                       ACRL1687
      MA(1)=MED(1)+MED(2)+MED(3)                                        ACRL1688
      MA(2)=MED(4)+MED(5)+MED(6)                                        ACRL1689
      MA(3)=MED(7)+MED(8)+MED(9)                                        ACRL1690
      MA(4)=MED(10)+MED(11)+MED(12)                                     ACRL1691
      MB(1)=MA(1)+MED(12)                                               ACRL1692
      MB(2)=MA(1)+MED(4)                                                ACRL1693
      MB(3)=MA(1)+MED(8)                                                ACRL1694
C     DETERMINE MAXIMUM OF (J1+J2+J3),(J1+L2+L3),(L1+J2+L3),(L1+L2+J3)  ACRL1695
      MAX=MA(1)                                                         ACRL1696
      DO 3 N=2,4                                                        ACRL1697
      IF (MAX-MA(N)) 2,3,3                                              ACRL1698
   2  MAX=MA(N)                                                         ACRL1699
   3  CONTINUE                                                          ACRL1700
C     DETERMINE MINIMUM OF (J1+J2+L1+L2), (J2+J3+L2+L3),(J3+J1+L3+L1)   ACRL1701
      MIN=MB(1)                                                         ACRL1702
      DO 5 N=2,3                                                        ACRL1703
      IF (MIN-MB(N)) 5,5,4                                              ACRL1704
   4  MIN=MB(N)                                                         ACRL1705
   5  CONTINUE                                                          ACRL1706
      KMAX=MIN-MAX                                                      ACRL1707
      MINP1=MIN+1                                                       ACRL1708
      MINI=MINP1-MA(1)                                                  ACRL1709
      MIN2=MINP1-MA(2)                                                  ACRL1710
      MIN3=MINP1-MA(3)                                                  ACRL1711
      MIN4=MINP1-MA(4)                                                  ACRL1712
      MIN5=MINP1+1                                                      ACRL1713
      MIN6=MB(1)-MIN                                                    ACRL1714
      MIN7=MB(2)-MIN                                                    ACRL1715
      MIN8=MB(3)-MIN                                                    ACRL1716
C     SUM SERIES                                                        ACRL1717
      UK=1.E-15                                                         ACRL1718
      S=1.E-15                                                          ACRL1719
      IF (KMAX) 8,8,6                                                   ACRL1720
   6  DO 7 K=1,KMAX                                                     ACRL1721
      UK=-UK*FLOAT(MINI-K)*FLOAT(MIN2-K)*FLOAT(MIN3-K)*FLOAT(MIN4-K)/   ACRL1722
     1(FLOAT(MIN5-K)*FLOAT(MIN6+K)*FLOAT(MIN7+K)*FLOAT(MIN8+K))         ACRL1723
C     CUT OFF SERIES AT 1.0E-25                                         ACRL1724
      IF (ABS(UK)-1.E-25) 8,8,7                                         ACRL1725
   7  S=S+UK                                                            ACRL1726
   8  S=S*1.E+15                                                        ACRL1727
C     CALCULATE DELTA FUNCTIONS                                         ACRL1728
      DELOG=0.0                                                         ACRL1729
      DO 9 N=1,12                                                       ACRL1730
      NUM=MED(N)                                                        ACRL1731
   9  DELOG=DELOG+FL(NUM+1)                                             ACRL1732
      NUM1=MA(1)+2                                                      ACRL1733
      NUM2=MA(2)+2                                                      ACRL1734
      NUM3=MA(3)+2                                                      ACRL1735
      NUM4=MA(4)+2                                                      ACRL1736
      DELOG=DELOG-FL(NUM1)-FL(NUM2)-FL(NUM3)-FL(NUM4)                   ACRL1737
      DELOG=0.5*DELOG                                                   ACRL1738
      ULOG=FL(MIN5)-FL(MINI)-FL(MIN2)-FL(MIN3)-FL(MIN4)-FL(MIN6+1)-     ACRL1739
     1 FL(MIN7+1)-FL(MIN8+1)                                            ACRL1740
      PLOG=DELOG+ULOG                                                   ACRL1741
C                                                                       ACRL1742
C         TO SAVE COMPUTER TIME THE 6-J COEFFICIENT IS NOT EVALUATED BY ACRL1743
C    F6J BUT RATHER THE VALUES PLOG AND S ARE EVALUATED.  IF A 6-J COEF-ACRL1744
C    FICIENT IS WANTED THE FOLLOWING COMMENTED CARDS ARE THE STATEMENTS ACRL1745
C    NEEDED TO STORE A 6-J COEFFICIENT IN S6J.                          ACRL1746
C     IF(PLOG +64.0) 72,75,75                                           ACRL1747
C  72 Q=PLOG+64.                                                        ACRL1748
C     Q=EXP(Q)                                                          ACRL1749
C     S6J=Q*S                                                           ACRL1750
C     IF(ABS(S6J)-1.0) 73,73,74                                         ACRL1751
C 73  S6J=0.0                                                           ACRL1752
C     GO TO 90                                                          ACRL1753
C 74  S6J=S6J*EXP(-64.)                                                 ACRL1754
C     GO TO 78                                                          ACRL1755
C 75  P=EXP(PLOG)                                                       ACRL1756
C     S6J =P*S                                                          ACRL1757
C  78 MIN2=MIN/2                                                        ACRL1758
C     IF (MIN-2*MIN2) 80,90,80                                          ACRL1759
C  80  S6J = - S6J                                                      ACRL1760
C  90  CONTINUE                                                         ACRL1761
      RETURN                                                            ACRL1762
  10  WRITE (IWRITE,11) N                                               ACRL1763
      CALL EXIT                                                         ACRL1764
      RETURN                                                            ACRL1765
C                                                                       ACRL1766
  11  FORMAT (28H0ONE OF THE CONDITIONS (NO. ,I3,1X,50H) FOR A NON-ZERO ACRL1767
     16-J COEFFICIENT IS NOT SATISFIED )                                ACRL1768
      END                                                               ACRL1769
      SUBROUTINE F3J (JD1,JD2,JD3,PLOG)                                 ACRL1770
C                                                                       ACRL1771
C         THIS SUBPROGRAM CALCULATES THE LOGARITHM OF THE ABSOLUTE VALUEACRL1772
C    OF A 3-J COEFFICIENT FOR THE SPECIAL CASE GIVEN IN EQN. 3.7.17 OF  ACRL1773
C    EDMONDS, ANGULAR MOMENTUM IN QUANTUM MECHANICS, AND ADDS IT TO PLOGACRL1774
C                                                                       ACRL1775
      COMMON/FACT/FL(322),NCALL,NFACT                                   ACRL1776
C       THE EXPLANATION OF THE /FACT/ VARIABLES BEGINS ON CARD NO.9.    ACRL1777
C                                                                       ACRL1778
C                            PARAMETER LIST                             ACRL1779
C                                                                       ACRL1780
C   JD1  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1781
C                                       (SEE CARD NO.1576)              ACRL1782
C   JD2  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1783
C   JD3  -  AN INPUT VARIABLE  -  EXPLAINED IN SUBPROGRAM Z             ACRL1784
C                                                                       ACRL1785
C   PLOG  -  AN INPUT/OUTPUT VARIABLE  -  SEE F6J FOR THE VALUE OF PLOG ACRL1786
C                  UPON ENTERING F3J.  F3J ADDS TO PLOG THE LOGARITHM   ACRL1787
C                    OF THE ABSOLUTE VALUE OF THE 3-J COEFFICIENT       ACRL1788
C                          (EQN. 3.7.17 OF EDMUNDS, ANG.MOM. IN Q.M.)   ACRL1789
C                                                                       ACRL1790
      J1=JD1                                                            ACRL1791
      J2=JD2                                                            ACRL1792
      J3=JD3                                                            ACRL1793
      NUM1=J1+J2-J3                                                     ACRL1794
      NUM2=J1-J2+J3                                                     ACRL1795
      NUM3=J2-J1+J3                                                     ACRL1796
      NUM4=J1+J2+J3                                                     ACRL1797
      NUM5=NUM1/2                                                       ACRL1798
      NUM6=NUM2/2                                                       ACRL1799
      NUM7=NUM3/2                                                       ACRL1800
      NUM8=NUM4/2                                                       ACRL1801
      NUM1=NUM1+1                                                       ACRL1802
      NUM2=NUM2+1                                                       ACRL1803
      NUM3=NUM3+1                                                       ACRL1804
      NUM4=NUM4+2                                                       ACRL1805
      NUM5=NUM5+1                                                       ACRL1806
      NUM6=NUM6+1                                                       ACRL1807
      NUM7=NUM7+1                                                       ACRL1808
      NUM8=NUM8+1                                                       ACRL1809
      PLOG=PLOG+0.5*(FL(NUM1)+FL(NUM2)+FL(NUM3)-FL(NUM4))               ACRL1810
      ULOG=FL(NUM8)-FL(NUM5)-FL(NUM6)-FL(NUM7)                          ACRL1811
      PLOG=PLOG+ULOG                                                    ACRL1812
C         THE SIGN HAS BEEN ABSORBED INTO THE Z COEFFICIENT SUBPROGRAM. ACRL1813
      RETURN                                                            ACRL1814
      END                                                               ACRL1815
$                                                                       ACRL1816
DIFFERENTIAL CROSS SECTION TEST RUN DATA                                ACRL1817
 TRUE TRUE                                                              ACRL1818
    4    3                                                              ACRL1819
    6    3   -2    1    1    0   -2    1   -2    1   -2    1    3    2  ACRL1820
UNITARIZED BORN APPROX. ELECTRON-H 2P-1S CROSS SECTION AT 3.40125 EV.   ACRL1821
    4    1    0     2.50E-01     1.00E+00     1.00E+00                  ACRL1822
F F T    6    1    0  100  200  201  300  301  302                      ACRL1823
 5.965736E-01-2.212346E-01 3.533350E-02-9.196875E-02 8.711633E-03       ACRL1824
 1.616407E-03-2.212346E-01 2.568147E+00 6.250000E-01-5.069591E-01       ACRL1825
 2.608294E-01-1.873782E-01 3.533350E-02 6.250000E-01-1.818528E-01       ACRL1826
 3.243882E-01-4.517620E-01 3.426710E-02-9.196875E-02-5.069591E-01       ACRL1827
 3.243882E-01 5.039721E+00 6.123724E-01 6.187184E-01 8.711633E-03       ACRL1828
 2.608294E-01-4.517620E-01 6.123724E-01 1.039721E+00 1.190785E+00       ACRL1829
 1.616407E-03-1.873782E-01 3.426710E-02 6.187184E-01 1.190785E+00       ACRL1830
-7.102792E-01                                                           ACRL1831
F F T   12    2    0  100  200  201  201  201  300  301  301  301  302  ACRL1832
  302  302                                                              ACRL1833
 9.657359E-02-6.320988E-02-2.017947E-01 0.          -1.504254E-03       ACRL1834
-1.603125E-02-8.838552E-02 0.          -8.404572E-04-1.346999E-03       ACRL1835
 0.           3.241747E-04-6.320988E-02 8.181472E-01 3.608439E-01       ACRL1836
 0.           7.144345E-01-3.471368E-01-2.957003E-02 0.                 ACRL1837
 1.286503E-01 1.636035E-01 0.          -5.590082E-02-2.017947E-01       ACRL1838
 3.608439E-01 2.151481E+00 0.           5.892557E-02 7.863620E-02       ACRL1839
-5.927969E-01 0.          -4.825146E-02 2.135515E-01 0.                 ACRL1840
 2.593391E-02 0.           0.           0.           9.431472E-01       ACRL1841
 0.           0.           0.          -2.749139E-01 0.                 ACRL1842
 0.           3.735818E-02 0.          -1.504254E-03 7.144345E-01       ACRL1843
 5.892557E-02 0.          -1.401861E-01 1.346244E-01 1.068760E-02       ACRL1844
 0.          -1.356298E-01-3.431459E-01 0.          -3.205330E-02       ACRL1845
-1.603125E-02-3.471368E-01 7.863620E-02 0.           1.346244E-01       ACRL1846
 2.289721E+00 3.535534E-01 0.           8.750000E-01-5.031153E-01       ACRL1847
 0.           1.597524E-01-8.838552E-02-2.957003E-02-5.927969E-01       ACRL1848
 0.           1.068760E-02 3.535534E-01 4.789721E+00 0.                 ACRL1849
 4.419417E-01 3.557562E-01 0.          -2.743363E-01 0.                 ACRL1850
 0.           0.          -2.749139E-01 0.           0.                 ACRL1851
 0.           2.602221E+00 0.           0.           5.625000E-01       ACRL1852
 0.          -8.404572E-04 1.286503E-01-4.825146E-02 0.                 ACRL1853
-1.356298E-01 8.750000E-01 4.419417E-01 0.           3.522208E-01       ACRL1854
 3.633610E-01 0.           1.209554E+00-1.346999E-03 1.636035E-01       ACRL1855
 2.135515E-01 0.          -3.431459E-01-5.031153E-01 3.557562E-01       ACRL1856
 0.           3.633610E-01 9.772208E-01 0.          -5.103104E-02       ACRL1857
 0.           0.           0.           3.735818E-02 0.                 ACRL1858
 0.           0.           5.625000E-01 0.           0.                 ACRL1859
 4.772208E-01 0.           3.241747E-04-5.590082E-02 2.593391E-02       ACRL1860
 0.          -3.205330E-02 1.597524E-01-2.743363E-01 0.                 ACRL1861
 1.209554E+00-5.103104E-02 0.          -4.602792E-01                    ACRL1862
F F T   14    3    1  100  200  201  201  201  300  301  301  301  302  ACRL1863
  302  302  302  302                                                    ACRL1864
 1.713205E-02-1.462122E-02-1.334570E-01 0.          -3.450528E-03       ACRL1865
-2.001413E-03-3.800055E-02 0.          -5.655127E-04-2.031637E-02       ACRL1866
 0.          -1.046950E-03 0.           5.810262E-05                    ACRL1867
-1.462122E-02 2.503810E-01 5.533986E-01 0.           6.433680E-01       ACRL1868
-1.684218E-01-1.201232E-01 0.           4.121420E-02 2.442444E-01       ACRL1869
 0.           1.400358E-01 0.          -1.309572E-02                    ACRL1870
-1.334570E-01 5.533986E-01 4.931472E-01 0.          -6.123724E-02       ACRL1871
 8.047072E-02-3.456531E-01 0.          -3.523053E-02-5.942942E-01       ACRL1872
 0.           2.252175E-02 0.           1.207469E-02                    ACRL1873
 0.           0.           0.           4.431472E-01 0.                 ACRL1874
 0.           0.          -1.375311E-01 0.           0.                 ACRL1875
-2.160468E-01 0.           1.155349E-03 0.                              ACRL1876
-3.450528E-03 6.433680E-01-6.123724E-02 0.          -1.244067E-01       ACRL1877
 1.631661E-02 1.005309E-01 0.          -3.009078E-02-5.584610E-02       ACRL1878
 0.          -2.415717E-01 0.          -3.125936E-02                    ACRL1879
-2.001413E-03-1.684218E-01 8.047072E-02 0.           1.631661E-02       ACRL1880
 1.074747E+00 6.777721E-01 0.           1.059446E+00 2.766993E-01       ACRL1881
 0.          -5.197011E-01 0.          -5.201481E-02                    ACRL1882
-3.800055E-02-1.201232E-01-3.456531E-01 0.           1.005309E-01       ACRL1883
 6.777721E-01 1.977221E+00 0.           1.530931E-01 2.755676E-01       ACRL1884
 0.           5.672419E-01 0.          -1.309307E-01                    ACRL1885
 0.           0.           0.          -1.375311E-01 0.                 ACRL1886
 0.           0.           1.477221E+00 0.           0.                 ACRL1887
 2.515577E-01 0.           6.149187E-01 0.                              ACRL1888
-5.655127E-04 4.121420E-02-3.523053E-02 0.          -3.009078E-02       ACRL1889
 1.059446E+00 1.530931E-01 0.          -6.560256E-02-2.125000E-01       ACRL1890
 0.           4.631511E-01 0.           1.114602E+00                    ACRL1891
-2.031637E-02 2.442444E-01-5.942942E-01 0.          -5.584610E-02       ACRL1892
 2.766993E-01 2.755676E-01 0.          -2.125000E-01 4.339721E+00       ACRL1893
 0.           2.539861E-01 0.           8.017837E-02                    ACRL1894
 0.           0.           0.          -2.160468E-01 0.                 ACRL1895
 0.           0.           2.515577E-01 0.           0.                 ACRL1896
 2.152221E+00 0.           2.500000E-02 0.                              ACRL1897
-1.046950E-03 1.400358E-01 2.252175E-02 0.          -2.415717E-01       ACRL1898
-5.197011E-01 5.672419E-01 0.           4.631511E-01 2.539861E-01       ACRL1899
 0.           6.022208E-01 0.          -1.118034E-01                    ACRL1900
 0.           0.           0.           1.155349E-03 0.                 ACRL1901
 0.           0.           6.149187E-01 0.           0.                 ACRL1902
 2.500000E-02 0.           2.397208E-01 0.                              ACRL1903
 5.810262E-05-1.309572E-02 1.207469E-02 0.          -3.125936E-02       ACRL1904
-5.201481E-02-1.309307E-01 0.           1.114602E+00 8.017837E-02       ACRL1905
 0.          -1.118034E-01 0.          -3.516411E-01                    ACRL1906
F F T   14    4    1  100  200  201  201  201  300  301  301  301  302  ACRL1907
  302  302  302  302                                                    ACRL1908
 3.088937E-03-3.031913E-03-6.932256E-02 0.          -1.903316E-03       ACRL1909
-1.690661E-04-1.259689E-02 0.          -1.783280E-04-1.226990E-02       ACRL1910
 0.          -4.469763E-04 0.           8.590166E-06                    ACRL1911
-3.031913E-03 7.026996E-02 5.437453E-01 0.           5.338655E-01       ACRL1912
-6.417845E-02-1.628210E-01 0.           5.590961E-03 3.241552E-01       ACRL1913
 0.           9.719933E-02 0.          -1.075145E-03                    ACRL1914
-6.932256E-02 5.437453E-01 6.814718E-02 0.          -7.692976E-02       ACRL1915
 3.768810E-02-1.363088E-01 0.          -1.369596E-02-6.063916E-01       ACRL1916
 0.          -1.713321E-02 0.           3.777207E-03                    ACRL1917
 0.           0.           0.           2.422164E-01 0.                 ACRL1918
 0.           0.          -6.614666E-02 0.           0.                 ACRL1919
-1.916587E-01 0.          -5.774255E-03 0.                              ACRL1920
-1.903316E-03 5.338655E-01-7.692976E-02 0.          -9.024609E-02       ACRL1921
-2.388023E-02 1.211259E-01 0.          -6.045809E-04-1.193303E-01       ACRL1922
 0.          -1.575185E-01 0.          -2.045355E-02                    ACRL1923
-1.690661E-04-6.417845E-02 3.768810E-02 0.          -2.388023E-02       ACRL1924
 4.547272E-01 8.953950E-01 0.           1.075714E+00 1.045825E-01       ACRL1925
 0.          -4.918830E-01 0.          -1.190795E-01                    ACRL1926
-1.259689E-02-1.628210E-01-1.363088E-01 0.           1.211259E-01       ACRL1927
 8.953950E-01 7.540065E-01 0.          -5.949533E-02 6.827993E-01       ACRL1928
 0.           6.465338E-01 0.          -2.886531E-02                    ACRL1929
 0.           0.           0.          -6.614666E-02 0.                 ACRL1930
 0.           0.           8.966362E-01 0.           0.                 ACRL1931
 3.674842E-01 0.           5.939127E-01 0.                              ACRL1932
-1.783280E-04 5.590961E-03-1.369596E-02 0.          -6.045809E-04       ACRL1933
 1.075714E+00-5.949533E-02 0.          -1.924052E-01-1.106567E-01       ACRL1934
 0.           4.533113E-01 0.           9.601793E-01                    ACRL1935
-1.226990E-02 3.241552E-01-6.063916E-01 0.          -1.193303E-01       ACRL1936
 1.045825E-01 6.827993E-01 0.          -1.106567E-01 1.396864E+00       ACRL1937
 0.          -6.492384E-12 0.           4.374089E-02                    ACRL1938
 0.           0.           0.          -1.916587E-01 0.                 ACRL1939
 0.           0.           3.674842E-01 0.           0.                 ACRL1940
 1.057578E+00 0.          -2.823462E-02 0.                              ACRL1941
-4.469763E-04 9.719933E-02-1.713321E-02 0.          -1.575185E-01       ACRL1942
-4.918830E-01 6.465338E-01 0.           4.533113E-01-6.492384E-12       ACRL1943
 0.           3.105541E-01 0.          -1.088058E-01                    ACRL1944
 0.           0.           0.          -5.774255E-03 0.                 ACRL1945
 0.           0.           5.939127E-01 0.           0.                 ACRL1946
-2.823462E-02 0.           1.276360E-01 0.                              ACRL1947
 8.590166E-06-1.075145E-03 3.777207E-03 0.          -2.045355E-02       ACRL1948
-1.190795E-01-2.886531E-02 0.           9.601793E-01 4.374089E-02       ACRL1949
 0.          -1.088058E-01 0.          -2.561448E-01                    ACRL1950
UNITARIZED BORN APPROX. ELECTRON-H 2P-2P CROSS SECTION AT 3.40125 EV.   ACRL1951
    4    1    0     2.50E-01     1.00E+00     1.00E+00                  ACRL1952
UNITARIZED BORN APPROX. ELECTRON-H 2P-3D CROSS SECTION AT 3.40125 EV.   ACRL1953
    4    1    0     2.50E-01     1.00E+00     1.00E+00                  ACRL1954
    1                                                                   ACRL1955
    5    2    1    0   -2    1                                          ACRL1956
C.C. APPROX. ELECTRON-H SINGLET STATE 1S-2P CROSS SECTION AT 13.605 EV. ACRL1957
    4    1    3     1.00E+00     1.00E+00     2.50E-01                  ACRL1958
F T T    6    1    0  100  200  201  300  301  302                      ACRL1959
 8.107000E-01 5.133900E-01-2.648800E-01-4.442700E-02-3.117200E-01       ACRL1960
 8.160700E-01 5.134900E-01 1.016100E+00 6.667200E-02 8.096500E-01       ACRL1961
 1.797700E-01 3.475900E-01-2.643000E-01 6.672500E-02-8.491700E-01       ACRL1962
-3.283000E-01-4.711500E-01 1.266800E+00-4.444800E-02 8.096300E-01       ACRL1963
-3.284100E-01 4.586500E-01 2.995100E-01 1.547400E+00-3.117800E-01       ACRL1964
 1.797600E-01-4.711500E-01 2.995000E-01 1.106200E+00-1.279000E+00       ACRL1965
 8.162000E-01 3.475500E-01 1.266600E+00 1.547400E+00-1.278900E+00       ACRL1966
 2.010700E+00                                                           ACRL1967
F T F    9    2    3  100  200  201  201  300  301  301  302  302       ACRL1968
    1    2    4    3    5    7    6    9    8                           ACRL1969
-3.498700E-01 1.253800E+00-4.143200E-01 1.372000E+00-1.110900E-01       ACRL1970
-6.812300E-01-3.908400E-02-5.724800E-01-9.600900E-01 1.254000E+00       ACRL1971
-5.200900E+00 2.667800E+00-5.206400E+00 4.462700E-01 2.912900E+00       ACRL1972
 7.954800E-02 2.015700E+00 3.556400E+00-4.143400E-01 2.667600E+00       ACRL1973
-2.921300E+00 1.713200E+00-1.719700E-01-4.635600E-01 8.383700E-01       ACRL1974
-8.934100E-01-1.871400E+00 1.372200E+00-5.206600E+00 1.713500E+00       ACRL1975
-5.418400E+00 4.543800E-01 3.033100E+00 5.703600E-01 2.405800E+00       ACRL1976
 4.024900E+00-1.111100E-01 4.463100E-01-1.720300E-01 4.543500E-01       ACRL1977
 9.431100E-02-9.386100E-01 2.940400E-01 4.345800E-01-2.786700E-01       ACRL1978
-6.814400E-01 2.913200E+00-4.638400E-01 3.033400E+00-9.386700E-01       ACRL1979
-1.321900E+00-1.126900E+00-1.128800E+00-2.334000E+00-3.913300E-02       ACRL1980
 7.972200E-02 8.382800E-01 5.705000E-01 2.940000E-01-1.127000E+00       ACRL1981
-5.206900E-01-3.430700E-01 1.701900E-01-5.726100E-01 2.015800E+00       ACRL1982
-8.934600E-01 2.405900E+00 4.345900E-01-1.128700E+00-3.430700E-01       ACRL1983
-3.185300E-02-1.067300E+00-9.601800E-01 3.558200E+00-1.871400E+00       ACRL1984
 4.024500E+00-2.786200E-01-2.333600E+00 1.703500E-01-1.067100E+00       ACRL1985
-2.468500E+00                                                           ACRL1986
F T F   10    3    2  100  200  201  201  300  301  301  302  302  302  ACRL1987
    1    2    4    3    5    7    6   10    9    8                      ACRL1988
 1.890700E-01-2.205700E-01 1.363900E-01-3.867000E-01 2.573700E-01       ACRL1989
-1.872100E-01 3.602700E-01-9.727200E-02-1.375800E-01-5.530500E-01       ACRL1990
-2.206200E-01 2.190200E+00-2.669800E-01 1.557900E+00-3.428900E-01       ACRL1991
-1.543700E+00-5.021900E-01-3.241800E-01 1.505900E+00 2.618700E+00       ACRL1992
 1.363700E-01-2.669700E-01-1.569900E-02-7.265100E-01 7.844100E-02       ACRL1993
 4.110400E-01-5.761100E-01 4.356400E-01-1.312200E+00-8.534300E-01       ACRL1994
-3.867800E-01 1.557800E+00-7.266100E-01 1.082000E+00-9.350400E-01       ACRL1995
 8.725800E-01-1.560500E+00 6.220600E-01 2.686100E-01 2.291000E+00       ACRL1996
 2.574400E-01-3.432400E-01 7.849300E-02-9.352100E-01 9.675800E-01       ACRL1997
 8.611800E-01 2.447700E-01 2.502200E-01-7.313500E-01-1.713500E+00       ACRL1998
-1.871200E-01-1.544400E+00 4.109300E-01 8.720000E-01 8.612600E-01       ACRL1999
 2.049200E+00 2.613900E+00 5.797600E-01-8.052100E-01-1.494900E+00       ACRL2000
 3.603200E-01-5.024200E-01-5.761900E-01-1.560700E+00 2.446200E-01       ACRL2001
 2.614000E+00-1.298300E+00 1.046400E+00-1.239900E+00-1.157000E+00       ACRL2002
-9.722300E-02-3.247100E-01 4.355800E-01 6.217500E-01 2.503300E-01       ACRL2003
 5.799000E-01 1.046500E+00 1.306700E+00-1.662100E-01 2.168200E-01       ACRL2004
-1.376100E-01 1.505900E+00-1.312200E+00 2.687100E-01-7.311700E-01       ACRL2005
-8.049200E-01-1.239700E+00-1.660700E-01 1.388100E+00 2.432800E+00       ACRL2006
-5.531300E-01 2.618700E+00-8.534400E-01 2.291100E+00-1.713100E+00       ACRL2007
-1.494200E+00-1.156700E+00 2.173100E-01 2.432800E+00 3.820400E+00       ACRL2008
F T T   10    4    2  100  200  201  201  300  301  301  302  302  302  ACRL2009
 4.091400E-02-2.814800E-02-6.045700E-03-6.657600E-02 1.236300E-01       ACRL2010
 7.890200E-03 5.281000E-02-5.095700E-02-6.550500E-02-1.884100E-02       ACRL2011
-2.807100E-02-1.672200E+00 7.577300E-01 3.544600E+00 1.144500E+01       ACRL2012
 8.421100E-01 5.688300E+00-4.230400E+00-4.736600E+00 9.312800E-02       ACRL2013
-6.058900E-03 7.581700E-01-3.772200E-01-5.937400E-01-2.228000E+00       ACRL2014
-1.728700E-01-1.069300E+00 8.634000E-01 4.568500E-01-1.128500E-01       ACRL2015
-6.667300E-02 3.547000E+00-5.936000E-01-4.054500E+00-1.658700E+01       ACRL2016
-1.191200E+00-7.879500E+00 6.517000E+00 7.351700E+00 1.010400E+00       ACRL2017
 1.232800E-01 1.145200E+01-2.227100E+00-1.658400E+01-6.344700E+01       ACRL2018
-3.959400E+00-3.120900E+01 2.386900E+01 2.750400E+01 1.273000E+00       ACRL2019
 7.886800E-03 8.424500E-01-1.728500E-01-1.190900E+00-3.959200E+00       ACRL2020
-5.881600E-01-2.430100E+00 1.119200E+00 2.553800E+00-7.766200E-02       ACRL2021
 5.263600E-02 5.692200E+00-1.068900E+00-7.878400E+00-3.121000E+01       ACRL2022
-2.430300E+00-1.564500E+01 1.170000E+01 1.413600E+01 2.941900E-01       ACRL2023
-5.082700E-02-4.232900E+00 8.629400E-01 6.515500E+00 2.386700E+01       ACRL2024
 1.119200E+00 1.169800E+01-8.240000E+00-9.998800E+00-1.187500E+00       ACRL2025
-6.534200E-02-4.739300E+00 4.563000E-01 7.349400E+00 2.750100E+01       ACRL2026
 2.553700E+00 1.413400E+01-9.998200E+00-1.200800E+01-8.027600E-01       ACRL2027
-1.883800E-02 9.341800E-02-1.130100E-01 1.009200E+00 1.268400E+00       ACRL2028
-7.793700E-02 2.919100E-01-1.185800E+00-8.008100E-01-5.999900E-01       ACRL2029
    1                                                                   ACRL2030
    5    3    1    0   -2    1                                          ACRL2031
C.C. APPROX. ELECTRON-H TRIPLET STATE 1S-2P CROSS SECTION AT 13.605 EV. ACRL2032
    4    1    1     1.00E+00     1.00E+00     7.50E-01                  ACRL2033
F F T    6    1    0  100  200  201  300  301  302                      ACRL2034
 5.291100E+00-3.405900E+00-2.457000E+00 1.099700E-01-5.055700E-02       ACRL2035
 1.474100E+00-3.406400E+00-1.479700E+01-1.002300E+01 9.910100E-01       ACRL2036
-1.410100E-02 6.140100E+00-2.457400E+00-1.002300E+01-6.768800E+00       ACRL2037
 3.158700E-01-2.027200E-01 4.217800E+00 1.102000E-01 9.919600E-01       ACRL2038
 3.164700E-01-7.443200E-02 9.821900E-01-1.900800E+00-5.030600E-02       ACRL2039
-1.298900E-02-2.020600E-01 9.821200E-01 9.619000E-01-1.138200E+00       ACRL2040
 1.474200E+00 6.139600E+00 4.217500E+00-1.900400E+00-1.137700E+00       ACRL2041
-3.391200E+00                                                           ACRL2042
F F T    9    2    0  100  200  201  201  300  301  301  302  302       ACRL2043
 2.741100E-01 9.060500E-01 6.546000E-01-1.295300E+00 1.829000E+00       ACRL2044
 1.157600E+00-5.474300E-01-1.720900E+00-3.389700E-01 9.062600E-01       ACRL2045
-5.032200E+00-3.232900E+00 6.069700E+00-6.942600E+00-4.751900E+00       ACRL2046
 3.267400E+00 5.665100E+00 9.010900E-01 6.546500E-01-3.232400E+00       ACRL2047
-3.399700E+00 6.470600E+00-8.886900E+00-4.826400E+00 1.490300E+00       ACRL2048
 9.442100E+00 1.750300E+00-1.295400E+00 6.068800E+00 6.470700E+00       ACRL2049
-1.191100E+01 1.799400E+01 1.035300E+01-2.683200E+00-1.905100E+01       ACRL2050
-4.070000E+00 1.829200E+00-6.941100E+00-8.887100E+00 1.799400E+01       ACRL2051
-2.704100E+01-1.515000E+01 2.727900E+00 2.995500E+01 6.887000E+00       ACRL2052
 1.157800E+00-4.751500E+00-4.827000E+00 1.035400E+01-1.515100E+01       ACRL2053
-9.535400E+00 2.519100E+00 1.683800E+01 4.228600E+00-5.475700E-01       ACRL2054
 3.267600E+00 1.490600E+00-2.683700E+00 2.728400E+00 2.519300E+00       ACRL2055
-1.437600E+00-2.884200E+00 1.949700E-01-1.721200E+00 5.664000E+00       ACRL2056
 9.443100E+00-1.905300E+01 2.995700E+01 1.683800E+01-2.883900E+00       ACRL2057
-3.225200E+01-8.498500E+00-3.390200E-01 9.006600E-01 1.750500E+00       ACRL2058
-4.070200E+00 6.887400E+00 4.228500E+00 1.950800E-01-8.498400E+00       ACRL2059
-2.794300E+00                                                           ACRL2060
F F T   10    3    1  100  200  201  201  300  301  301  302  302  302  ACRL2061
 8.945300E-02-8.082200E-02 4.196700E-02 2.805800E-02-5.845900E-02       ACRL2062
 5.876200E-02-1.748700E-02-6.211600E-02-1.901000E-02-6.731300E-02       ACRL2063
-8.084000E-02 1.130300E+00 8.686800E-01-2.010500E-01 6.040800E-01       ACRL2064
 1.332000E-01 3.808100E-01 1.325200E-01-1.126500E-01 3.883100E-01       ACRL2065
 4.197200E-02 8.685200E-01-2.923700E+00-1.818900E-01 7.416200E-01       ACRL2066
 5.440600E-02-2.443800E-01 6.292200E-01 6.892200E-01 5.798200E-01       ACRL2067
 2.806200E-02-2.010900E-01-1.818800E-01-6.854900E-02 2.189300E-01       ACRL2068
-3.122600E-01 2.196800E-01 5.208300E-01 1.915700E-01 5.192200E-01       ACRL2069
-5.847700E-02 6.040300E-01 7.418200E-01 2.189500E-01 8.072900E-01       ACRL2070
-3.217600E-01 7.749100E-02 5.271900E-01 8.821100E-02 6.622900E-02       ACRL2071
 5.878600E-02 1.333300E-01 5.398300E-02-3.123000E-01-3.217400E-01       ACRL2072
-1.177800E+00 9.717400E-01-5.868400E-01-6.250900E-01 3.137600E-01       ACRL2073
-1.749200E-02 3.807100E-01-2.443200E-01 2.196700E-01 7.741300E-02       ACRL2074
 9.717900E-01-6.016500E-01-4.705600E-01-5.430200E-01 7.947300E-01       ACRL2075
-6.213600E-02 1.326600E-01 6.292400E-01 5.208400E-01 5.272400E-01       ACRL2076
-5.868400E-01-4.704900E-01 8.655300E-01-6.758700E-01 2.686900E-01       ACRL2077
-1.901900E-02-1.124400E-01 6.891300E-01 1.915800E-01 8.837900E-02       ACRL2078
-6.251700E-01-5.429700E-01-6.757800E-01-1.275400E+00-5.681800E-01       ACRL2079
-6.731300E-02 3.880400E-01 5.798400E-01 5.191600E-01 6.598200E-02       ACRL2080
 3.137700E-01 7.946800E-01 2.685100E-01-5.683100E-01 2.159300E+00       ACRL2081
T F F   10    4    2  100  200  201  201  300  301  301  302  302  302  ACRL2082
    1    2    4    3    5    7    6   10    9    8                      ACRL2083
-2.115300E-02-2.607100E-03-1.811000E-03 1.110600E-01-6.435700E-02       ACRL2084
 1.863300E-02 7.389100E-02-9.802100E-03-3.066400E-02 3.236100E-02       ACRL2085
-2.609600E-03-5.498900E-01-9.708200E-02-1.515200E-01 6.893900E-02       ACRL2086
-4.343200E-02 1.104100E-01-3.508100E-02 5.303400E-02-3.976800E-02       ACRL2087
-1.814100E-03-9.710200E-02-3.484100E-01-1.581100E-01-1.847300E-01       ACRL2088
 6.162500E-02 1.322500E-01 8.577300E-02 2.765300E-01 3.843200E-02       ACRL2089
 1.110800E-01-1.515400E-01-1.581100E-01-1.428300E+00-6.295300E-02       ACRL2090
-1.281500E-02-2.094900E-02-2.552200E-02-1.965300E-01-1.205400E-01       ACRL2091
-6.437000E-02 6.914700E-02-1.846700E-01-6.292900E-02-1.537000E+00       ACRL2092
 1.557700E-01-6.123100E-02 1.931300E-02-1.156300E-01 1.986900E-01       ACRL2093
 1.863500E-02-4.330500E-02 6.165700E-02-1.276000E-02 1.557600E-01       ACRL2094
-2.646300E-01-3.170700E-01-3.369800E-01-1.125000E-01-1.208900E-01       ACRL2095
 7.390300E-02 1.104300E-01 1.322700E-01-2.092000E-02-6.125400E-02       ACRL2096
-3.170900E-01-1.375400E+00-1.543400E-01 7.351200E-02-4.073600E-01       ACRL2097
-9.803300E-03-3.507300E-02 8.582500E-02-2.548200E-02 1.933300E-02       ACRL2098
-3.369700E-01-1.543200E-01-1.695000E+00 1.276600E-01 5.512300E-02       ACRL2099
-3.067400E-02 5.309200E-02 2.765700E-01-1.965100E-01-1.156400E-01       ACRL2100
-1.125000E-01 7.354000E-02 1.276700E-01-1.508500E+00 6.023000E-02       ACRL2101
 3.237100E-02-3.982700E-02 3.839500E-02-1.205900E-01 1.986700E-01       ACRL2102
-1.209400E-01-4.073700E-01 5.509600E-02 6.020400E-02-2.991900E-01       ACRL2103
 7.242600E-02 3.215100E-02-1.044600E-02-7.246600E-02 2.849100E-02       ACRL2104
-3.269800E-02 4.500700E-02 1.141100E-02 9.603200E-03-4.132600E-02       ACRL2105
 3.215100E-02 7.499100E-01 2.782600E-01-6.902700E-02 5.149000E-02       ACRL2106
 1.567600E-01 1.076300E-01-1.233600E-01 8.753900E-03 2.020900E-01       ACRL2107
-1.045000E-02 2.782500E-01-3.543900E-01 1.633800E-01 1.536900E-01       ACRL2108
-2.371300E-01-3.033300E-02 2.870300E-02-2.901000E-01 6.846100E-02       ACRL2109
-7.247600E-02-6.901000E-02 1.634100E-01 7.707600E-01-3.523700E-03       ACRL2110
 3.360100E-02 1.668700E-01-8.650900E-02-1.063100E-01 1.390200E-01       ACRL2111
 2.849400E-02 5.152500E-02 1.536600E-01-3.460600E-03 5.655800E-01       ACRL2112
 2.664200E-02-1.062700E-01 4.132800E-01-5.675500E-04-2.316500E-01       ACRL2113
-3.270700E-02 1.566200E-01-2.371700E-01 3.364500E-02 2.663900E-02       ACRL2114
 2.439900E-01-2.580900E-02-1.170900E-01 7.907100E-02-1.534000E-01       ACRL2115
 4.500700E-02 1.077200E-01-3.030200E-02 1.669700E-01-1.062700E-01       ACRL2116
-2.577200E-02-5.404900E-01 2.678200E-01 3.251100E-01 1.100900E-02       ACRL2117
 1.141000E-02-1.233600E-01 2.871800E-02-8.645700E-02 4.133000E-01       ACRL2118
-1.171300E-01 2.677900E-01-2.583100E-01 5.367500E-02 3.185000E-02       ACRL2119
 9.604800E-03 8.781100E-03-2.901200E-01-1.063200E-01-5.735400E-04       ACRL2120
 7.910200E-02 3.250900E-01 5.372400E-02 5.901300E-01 1.115400E-01       ACRL2121
-4.133500E-02 2.021400E-01 6.845600E-02 1.390800E-01-2.316700E-01       ACRL2122
-1.534000E-01 1.100600E-02 3.184600E-02 1.115200E-01-3.258200E-01       ACRL2123
    2                                                                   ACRL2124
    2    3    1    0    1    2    1    2    1    2                      ACRL2125
C.C. APPROX. HE-H2 J=0 TO J=2 CROSS SECTION AT 2.087E-3 EV.             ACRL2126
    6    1    0   3.7445E-01    2.441E+03     1.00E+00                  ACRL2127
T F T    2    1    1  100  102                                          ACRL2128
-.1693495E+01 .5323425E+00                                              ACRL2129
 .5323425E+00-.1695150E+00                                              ACRL2130
-.4570551E+00 .1636344E+00                                              ACRL2131
 .1636344E+00 .1139452E-01                                              ACRL2132
T F T    3    2    1  100  102  102                                     ACRL2133
-.7097577E+00-.2183917E-01 .5630469E-01                                 ACRL2134
-.2183917E-01-.7036191E+00-.4817463E-01                                 ACRL2135
 .5630469E-01-.4817463E-01-.1120148E-01                                 ACRL2136
 .2021023E+00 .9317109E+00 .5640690E-01                                 ACRL2137
 .9317109E+00-.1797029E+00-.9477518E-01                                 ACRL2138
 .5640690E-01-.9477518E-01 .6798824E-01                                 ACRL2139
T F T    4    3    1  100  102  102  102                                ACRL2140
-.1631003E+01 .4239192E+00-.1747616E+00 .1175974E-01                    ACRL2141
 .4239192E+00-.6267108E+00-.1847001E+00-.1549880E-01                    ACRL2142
-.1747616E+00-.1847001E+00-.1570284E+00-.5728969E-02                    ACRL2143
 .1175974E-01-.1549880E-01-.5728969E-02-.5914525E-03                    ACRL2144
-.4848438E+00 .3764711E+00 .1203968E+00 .1216714E-01                    ACRL2145
 .3764711E+00 .5426095E+00 .4583819E+00 .1799443E-01                    ACRL2146
 .1203968E+00 .4583819E+00 .9966603E-02 .1258322E-02                    ACRL2147
 .1216714E-01 .1799443E-01 .1258322E-02 .1725710E-01                    ACRL2148
T F T    4    4    1  100  102  102  102                                ACRL2149
-.5161516E+00 .3875456E+00-.1022936E+00 .4026469E-03                    ACRL2150
 .3875456E+00-.4148307E+00 .1222162E+00-.2198170E-03                    ACRL2151
-.1022936E+00 .1222162E+00-.5531648E-01-.2518338E-04                    ACRL2152
 .4026469E-03-.2198170E-03-.2518338E-04-.1649628E-04                    ACRL2153
 .3595631E+00-.6412597E+00 .2544315E+00 .1635454E-03                    ACRL2154
-.6412597E+00 .2823707E+00 .3873358E-01 .1005711E-02                    ACRL2155
 .2544315E+00 .3873358E-01 .1262398E+00 .2764965E-03                    ACRL2156
 .1635454E-03 .1005711E-02 .2764965E-03 .5627305E-02                    ACRL2157
T F T    4    5    1  100  102  102  102                                ACRL2158
-.3244886E+00 .1947003E+00-.1292030E-02-.1222172E-04                    ACRL2159
 .1947003E+00-.1243915E+00 .1022102E-02 .9586214E-05                    ACRL2160
-.1292030E-02 .1022102E-02-.1443689E-03-.8785323E-06                    ACRL2161
-.1222172E-04 .9586214E-05-.8785323E-06-.2653204E-05                    ACRL2162
 .5766874E+00-.4161516E+00 .4892093E-02 .4452926E-04                    ACRL2163
-.4161516E+00 .1490413E+00 .2244424E-02 .1817823E-04                    ACRL2164
 .4892093E-02 .2244424E-02 .1603204E-01 .8032801E-04                    ACRL2165
 .4452926E-04 .1817823E-04 .8032801E-04 .2301607E-02                    ACRL2166
T F T    4    6    1  100  102  102  102                                ACRL2167
-.1019423E+01-.1311549E+00-.1547393E-03 .2059899E-05                    ACRL2168
-.1311549E+00-.2209323E-01-.4854122E-04 .2869992E-06                    ACRL2169
-.1547393E-03-.4854122E-04-.1364904E-04 .4496582E-08                    ACRL2170
 .2059899E-05 .2869992E-06 .4496582E-08-.5880764E-06                    ACRL2171
-.9812391E+00-.1399672E+00-.2215424E-03 .2041282E-05                    ACRL2172
-.1399672E+00 .8310020E-01 .5063346E-03-.1302553E-06                    ACRL2173
-.2215424E-03 .5063346E-03 .5192894E-02-.1297100E-05                    ACRL2174
 .2041286E-05-.1302559E-06-.1297073E-05 .1084501E-02                    ACRL2175
C.C. APPROX. HE-H2 J=2 TO J=2 CROSS SECTION AT 1.905E-5 EV.             ACRL2176
    6    1    0   3.4174E-03    2.441E+03     1.00E+00                  ACRL2177
                                                                        ACRL****
ACRL0001 ACRL ADAPTED FOR IBM360/370. PROGRAM ACRL TO CALCULATE         ACRLA000
1   DIFFERENTIAL AND INTEGRAL CROSS SECTIONS ADAPTED TO RUN ON IBM      ACRLA000
2   COMPUTERS.  BRANDT, M.A., TRUHLAR, D.G., SMITH, R.L.                ACRLA000
REF. IN COMP. PHYS. COMMUN. 7 (1974) 172                                ACRLA000
C         REPLACE CARD 41 WITH NEXT 5 CARDS                             ACRLA001
      NCALL = 6                                                         ACRLA002
      NFACT = 70                                                        ACRLA003
      IREAD = 5                                                         ACRLA004
      IWRITE = 6                                                        ACRLA005
      ITAPE = 9                                                         ACRLA006
C         REPLACE CARD 169 WITH NEXT CARD                               ACRLA007
     6         ,LABEL(80)                                               ACRLA008
C         REPLACE CARD 295 WITH NEXT CARD                               ACRLA009
      CALL INOUT (1,LABEL,TMP2,TMP3,TMP4,IAPMAX,JAPMAX,IC,80,IE,IG)     ACRLA010
C         REPLACE CARD 394 WITH NEXT CARD                               ACRLA011
      CALL INOUT (3,LABEL,TK,TMASS,WQ,JTOT,JMIN,JDEG,80,IE,IG)          ACRLA012
C         REPLACE CARDS 974 TO 975 WITH NEXT 2 CARDS                    ACRLA013
      READ (K,40) (TMP1(I),I=1,80)                                      ACRLA014
      WRITE (L,41) (TMP1(I),I=1,80)                                     ACRLA015
C         REPLACE CARDS 993 TO 995 WITH NEXT 3 CARDS                    ACRLA016
      READ (K,40) (TMP1(I),I=1,80)                                      ACRLA017
      IF (.NOT.LPR1) WRITE (L,78) (TMP1(I),I=1,80)                      ACRLA018
      IF (LPR1) WRITE (L,79) (TMP1(I),I=1,80)                           ACRLA019
C         INSERT NEXT CARD AFTER CARD 1001                              ACRLA020
      IF (LPR1) WRITE (L,56)                                            ACRLA021
C         REPLACE CARDS 1200 TO 1201 WITH NEXT 2 CARDS                  ACRLA022
  40  FORMAT (80A1)                                                     ACRLA023
  41  FORMAT (9X,1H1,10X,80A1)                                          ACRLA024
C         REPLACE CARDS 1262 TO 1263 WITH NEXT 2 CARDS                  ACRLA025
  78  FORMAT (9X,1H5,10X,80A1)                                          ACRLA026
  79  FORMAT (8X,2H10,10X,80A1)                                         ACRLA027
                                                                        ACRL****
