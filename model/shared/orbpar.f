      Subroutine ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
C****
!@sum ORBPAR calculates the three orbital parameters as a function of
!@+   YEAR.  The source of these calculations is: Andre L. Berger,
!@+   1978, "Long-Term Variations of Daily Insolation and Quaternary
!@+   Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!@+   Berger, May 1978, "A Simple Algorithm to Compute Long Term
!@+   Variations of Daily Insolation", published by Institut
!@+   D'Astronomie de Geophysique, Universite Catholique de Louvain,
!@+   Louvain-la Neuve, No. 18.
!@+
!@+   Tables and equations refer to the first reference (JAS).  The
!@+   corresponding table or equation in the second reference is
!@+   enclosed in parentheses.  The coefficients used in this
!@+   subroutine are slightly more precise than those used in either
!@+   of the references.  The generated orbital parameters are precise
!@+   within plus or minus 1000000 years from present.
C****
!@auth Gary L. Russell (with extra terms from D. Thresher)
C****
      USE CONSTANT, only : twopi,PI180=>radian
      IMPLICIT NONE
C**** Input:
!@var YEAR   = years C.E. are positive, B.C.E are -ve (i.e 4BCE = -3)
      REAL*8, INTENT(IN) :: YEAR
C**** Output:
!@var ECCEN  = eccentricity of orbital ellipse
!@var OBLIQ  = latitude of Tropic of Cancer in degrees
!@var OMEGVP = longitude of perihelion =
!@+          = spatial angle from vernal equinox to perihelion
!@+            in degrees with sun as angle vertex
      REAL*8, INTENT(OUT) :: ECCEN,OBLIQ,OMEGVP
C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      REAL*8, PARAMETER, DIMENSION(3,47) :: TABL1 = RESHAPE( (/
     1            -2462.2214466d0, 31.609974d0, 251.9025d0,
     2             -857.3232075d0, 32.620504d0, 280.8325d0,
     3             -629.3231835d0, 24.172203d0, 128.3057d0,
     4             -414.2804924d0, 31.983787d0, 292.7252d0,
     5             -311.7632587d0, 44.828336d0,  15.3747d0,
     6              308.9408604d0, 30.973257d0, 263.7951d0,
     7             -162.5533601d0, 43.668246d0, 308.4258d0,
     8             -116.1077911d0, 32.246691d0, 240.0099d0,
     9              101.1189923d0, 30.599444d0, 222.9725d0,
     O              -67.6856209d0, 42.681324d0, 268.7809d0,
     1               24.9079067d0, 43.836462d0, 316.7998d0,
     2               22.5811241d0, 47.439436d0, 319.6024d0,
     3              -21.1648355d0, 63.219948d0, 143.8050d0,
     4              -15.6549876d0, 64.230478d0, 172.7351d0,
     5               15.3936813d0,  1.010530d0,  28.9300d0,
     6               14.6660938d0,  7.437771d0, 123.5968d0,
     7              -11.7273029d0, 55.782177d0,  20.2082d0,
     8               10.2742696d0,   .373813d0,  40.8226d0,
     9                6.4914588d0, 13.218362d0, 123.4722d0,
     O                5.8539148d0, 62.583231d0, 155.6977d0,
     1               -5.4872205d0, 63.593761d0, 184.6277d0,
     2               -5.4290191d0, 76.438310d0, 267.2772d0,
     3                5.1609570d0, 45.815258d0,  55.0196d0,
     4                5.0786314d0,  8.448301d0, 152.5268d0,
     5               -4.0735782d0, 56.792707d0,  49.1382d0,
     6                3.7227167d0, 49.747842d0, 204.6609d0,
     7                3.3971932d0, 12.058272d0,  56.5233d0,
     8               -2.8347004d0, 75.278220d0, 200.3284d0,
     9               -2.6550721d0, 65.241008d0, 201.6651d0,
     O               -2.5717867d0, 64.604291d0, 213.5577d0,
     1               -2.4712188d0,  1.647247d0,  17.0374d0,
     2                2.4625410d0,  7.811584d0, 164.4194d0,
     3                2.2464112d0, 12.207832d0,  94.5422d0,
     4               -2.0755511d0, 63.856665d0, 131.9124d0,
     5               -1.9713669d0, 56.155990d0,  61.0309d0,
     6               -1.8813061d0, 77.448840d0, 296.2073d0,
     7               -1.8468785d0,  6.801054d0, 135.4894d0,
     8                1.8186742d0, 62.209418d0, 114.8750d0,
     9                1.7601888d0, 20.656133d0, 247.0691d0,
     O               -1.5428851d0, 48.344406d0, 256.6114d0,
     1                1.4738838d0, 55.145460d0,  32.1008d0,
     2               -1.4593669d0, 69.000539d0, 143.6804d0,
     3                1.4192259d0, 11.071350d0,  16.8784d0,
     4               -1.1818980d0, 74.291298d0, 160.6835d0,
     5                1.1756474d0, 11.047742d0,  27.5932d0,
     6               -1.1316126d0,  0.636717d0, 348.1074d0,
     7                1.0896928d0, 12.844549d0,  82.6496d0/),(/3,47/) )
C**** Table 2 (4).  Precessional parameter: ECCEN sin(omega) (unused)
      REAL*8, PARAMETER, DIMENSION(3,46) :: TABL2 = RESHAPE(  (/
     1     .0186080D0,  54.646484D0,   32.012589D0,
     2     .0162752D0,  57.785370D0,  197.181274D0,
     3    -.0130066D0,  68.296539D0,  311.699463D0,
     4     .0098883D0,  67.659821D0,  323.592041D0,
     5    -.0033670D0,  67.286011D0,  282.769531D0,
     6     .0033308D0,  55.638351D0,   90.587509D0,
     7    -.0023540D0,  68.670349D0,  352.522217D0,
     8     .0014002D0,  76.656036D0,  131.835892D0,
     9     .0010070D0,  56.798447D0,  157.536392D0,
     O     .0008570D0,  66.649292D0,  294.662109D0,
     1     .0006499D0,  53.504456D0,  118.253082D0,
     2     .0005990D0,  67.023102D0,  335.484863D0,
     3     .0003780D0,  68.933258D0,  299.806885D0,
     4    -.0003370D0,  56.630219D0,  149.162415D0,
     5     .0003334D0,  86.256454D0,  283.915039D0,
     6     .0003334D0,  23.036499D0,  320.110107D0,
     7     .0002916D0,  89.395340D0,   89.083817D0,
     8     .0002916D0,  26.175385D0,  125.278732D0,
     9     .0002760D0,  69.307068D0,  340.629639D0,
     O    -.0002330D0,  99.906509D0,  203.602081D0,
     1    -.0002330D0,  36.686569D0,  239.796982D0,
     2     .0001820D0,  67.864838D0,  155.484787D0,
     3     .0001772D0,  99.269791D0,  215.494690D0,
     4     .0001772D0,  36.049850D0,  251.689606D0,
     5    -.0001740D0,  56.625275D0,  130.232391D0,
     6    -.0001240D0,  68.856720D0,  214.059708D0,
     7     .0001153D0,  87.266983D0,  312.845215D0,
     8     .0001153D0,  22.025970D0,  291.179932D0,
     9     .0001008D0,  90.405869D0,  118.013870D0,
     O     .0001008D0,  25.164856D0,   96.348694D0,
     1     .0000912D0,  78.818680D0,  160.318298D0,
     2     .0000912D0,  30.474274D0,   83.706894D0,
     3    -.0000806D0, 100.917038D0,  232.532120D0,
     4    -.0000806D0,  35.676025D0,  210.866943D0,
     5     .0000798D0,  81.957565D0,  325.487061D0,
     6     .0000798D0,  33.613159D0,  248.875565D0,
     7    -.0000638D0,  92.468735D0,   80.005234D0,
     8    -.0000638D0,  44.124329D0,    3.393823D0,
     9     .0000612D0, 100.280319D0,  244.424728D0,
     O     .0000612D0,  35.039322D0,  222.759552D0,
     1    -.0000603D0,  98.895981D0,  174.672028D0,
     2    -.0000603D0,  35.676025D0,  210.866943D0,
     3     .0000597D0,  87.248322D0,  342.489990D0,
     4     .0000597D0,  24.028381D0,   18.684967D0,
     5     .0000559D0,  86.630264D0,  324.737793D0,
     6     .0000559D0,  22.662689D0,  279.287354D0/), (/3,46/) )
C**** Table 3 (5).  Eccentricity: ECCEN (unused)
      REAL*8, PARAMETER, DIMENSION(3,42) :: TABL3 = RESHAPE( (/
     1     .01102940D0,   3.138886D0,  165.168686D0,
     2    -.00873296D0,  13.650058D0,  279.687012D0,
     3    -.00749255D0,  10.511172D0,  114.518250D0,
     4     .00672394D0,  13.013341D0,  291.579590D0,
     5     .00581229D0,   9.874455D0,  126.410858D0,
     6    -.00470066D0,   0.636717D0,  348.107422D0,
     7    -.00254464D0,  12.639528D0,  250.756897D0,
     8     .00231485D0,   0.991874D0,   58.574905D0,
     9    -.00221955D0,   9.500642D0,   85.588211D0,
     O     .00201868D0,   2.147012D0,  106.593765D0,
     1    -.00172371D0,   0.373813D0,   40.822647D0,
     2    -.00166112D0,  12.658154D0,  221.112030D0,
     3     .00145096D0,   1.010530D0,   28.930038D0,
     4     .00131342D0,  12.021467D0,  233.004639D0,
     5     .00101442D0,   0.373813D0,   40.822647D0,
     6    -.00088343D0,  14.023871D0,  320.509521D0,
     7    -.00083395D0,   6.277772D0,  330.337402D0,
     8     .00079475D0,   6.277772D0,  330.337402D0,
     9     .00067546D0,  27.300110D0,  199.373871D0,
     O    -.00066447D0,  10.884985D0,  155.340912D0,
     1     .00062591D0,  21.022339D0,  229.036499D0,
     2     .00059751D0,  22.009552D0,   99.823303D0,
     3    -.00053262D0,  27.300110D0,  199.373871D0,
     4    -.00052983D0,   5.641055D0,  342.229980D0,
     5    -.00052983D0,   6.914489D0,  318.444824D0,
     6     .00052836D0,  12.002811D0,  262.649414D0,
     7     .00051457D0,  16.788940D0,   84.855621D0,
     8    -.00050748D0,  11.647654D0,  192.181992D0,
     9    -.00049048D0,  24.535049D0,   75.027847D0,
     O     .00048888D0,  18.870667D0,  294.654541D0,
     1     .00046278D0,  26.026688D0,  223.159103D0,
     2     .00046212D0,   8.863925D0,   97.480820D0,
     3     .00046046D0,  17.162750D0,  125.678268D0,
     4     .00042941D0,   2.151964D0,  125.523788D0,
     5     .00042342D0,  37.174576D0,  325.784668D0,
     6     .00041713D0,  19.748917D0,  252.821732D0,
     7    -.00040745D0,  21.022339D0,  229.036499D0,
     8    -.00040569D0,   3.512699D0,  205.991333D0,
     9    -.00040569D0,   1.765073D0,  124.346024D0,
     O    -.00040385D0,  29.802292D0,   16.435165D0,
     1     .00040274D0,   7.746099D0,  350.172119D0,
     2     .00040068D0,   1.142024D0,  273.759521D0/), (/3,42/) )
C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      REAL*8, PARAMETER, DIMENSION(3,19) :: TABL4 = RESHAPE( (/
     1     .01860798D0,   4.207205D0,   28.620089D0,
     2     .01627522D0,   7.346091D0,  193.788772D0,
     3    -.01300660D0,  17.857263D0,  308.307024D0,
     4     .00988829D0,  17.220546D0,  320.199637D0,
     5    -.00336700D0,  16.846733D0,  279.376984D0,
     6     .00333077D0,   5.199079D0,   87.195000D0,
     7    -.00235400D0,  18.231076D0,  349.129677D0,
     8     .00140015D0,  26.216758D0,  128.443387D0,
     9     .00100700D0,   6.359169D0,  154.143880D0,
     O     .00085700D0,  16.210016D0,  291.269597D0,
     1     .00064990D0,   3.065181D0,  114.860583D0,
     2     .00059900D0,  16.583829D0,  332.092251D0,
     3     .00037800D0,  18.493980D0,  296.414411D0,
     4    -.00033700D0,   6.190953D0,  145.769910D0,
     5     .00027600D0,  18.867793D0,  337.237063D0,
     6     .00018200D0,  17.425567D0,  152.092288D0,
     7    -.00017400D0,   6.186001D0,  126.839891D0,
     8    -.00012400D0,  18.417441D0,  210.667199D0,
     9     .00001250D0,   0.667863D0,   72.108838D0/), (/3,19/) )
C**** Table 5 (3).  General precession in longitude: psi
      REAL*8, PARAMETER, DIMENSION(3,78) :: TABL5 = RESHAPE( (/
     1    7391.0225890d0,  31.609974d0,   251.9025d0,
     2    2555.1526947d0,  32.620504d0,   280.8325d0,
     3    2022.7629188d0,  24.172203d0,   128.3057d0,
     4   -1973.6517951d0,   0.636717d0,   348.1074d0,
     5    1240.2321818d0,  31.983787d0,   292.7252d0,
     6     953.8679112d0,   3.138886d0,   165.1686d0,
     7    -931.7537108d0,  30.973257d0,   263.7951d0,
     8     872.3795383d0,  44.828336d0,    15.3747d0,
     9     606.3544732d0,   0.991874d0,    58.5749d0,
     O    -496.0274038d0,   0.373813d0,    40.8226d0,
     1     456.9608039d0,  43.668246d0,   308.4258d0,
     2     346.9462320d0,  32.246691d0,   240.0099d0,
     3    -305.8412902d0,  30.599444d0,   222.9725d0,
     4     249.6173246d0,   2.147012d0,   106.5937d0,
     5    -199.1027200d0,  10.511172d0,   114.5182d0,
     6     191.0560889d0,  42.681324d0,   268.7809d0,
     7    -175.2936572d0,  13.650058d0,   279.6869d0,
     8     165.9068833d0,   0.986922d0,    39.6448d0,
     9     161.1285917d0,   9.874455d0,   126.4108d0,
     O     139.7878093d0,  13.013341d0,   291.5795d0,
     1    -133.5228399d0,   0.262904d0,   307.2848d0,
     2     117.0673811d0,   0.004952d0,    18.9300d0,
     3     104.6907281d0,   1.142024d0,   273.7596d0,
     4      95.3227476d0,  63.219948d0,   143.8050d0,
     5      86.7824524d0,   0.205021d0,   191.8927d0,
     6      86.0857729d0,   2.151964d0,   125.5237d0,
     7      70.5893698d0,  64.230478d0,   172.7351d0,
     8     -69.9719343d0,  43.836462d0,   316.7998d0,
     9     -62.5817473d0,  47.439436d0,   319.6024d0,
     O      61.5450059d0,   1.384343d0,    69.7526d0,
     1     -57.9364011d0,   7.437771d0,   123.5968d0,
     2      57.1899832d0,  18.829299d0,   217.6432d0,
     3     -57.0236109d0,   9.500642d0,    85.5882d0,
     4     -54.2119253d0,   0.431696d0,   156.2147d0,
     5      53.2834147d0,   1.160090d0,    66.9489d0,
     6      52.1223575d0,  55.782177d0,    20.2082d0,
     7     -49.0059908d0,  12.639528d0,   250.7568d0,
     8     -48.3118757d0,   1.155138d0,    48.0188d0,
     9     -45.4191685d0,   0.168216d0,     8.3739d0,
     O     -42.2357920d0,   1.647247d0,    17.0374d0,
     1     -34.7971099d0,  10.884985d0,   155.3409d0,
     2      34.4623613d0,   5.610937d0,    94.1709d0,
     3     -33.8356643d0,  12.658184d0,   221.1120d0,
     4      33.6689362d0,   1.010530d0,    28.9300d0,
     5     -31.2521586d0,   1.983748d0,   117.1498d0,
     6     -30.8798701d0,  14.023871d0,   320.5095d0,
     7      28.4640769d0,   0.560178d0,   262.3602d0,
     8     -27.1960802d0,   1.273434d0,   336.2148d0,
     9      27.0860736d0,  12.021467d0,   233.0046d0,
     O     -26.3437456d0,  62.583231d0,   155.6977d0,
     1      24.7253740d0,  63.593761d0,   184.6277d0,
     2      24.6732126d0,  76.438310d0,   267.2772d0,
     3      24.4272733d0,   4.280910d0,    78.9281d0,
     4      24.0127327d0,  13.218362d0,   123.4722d0,
     5      21.7150294d0,  17.818769d0,   188.7132d0,
     6     -21.5375347d0,   8.359495d0,   180.1364d0,
     7      18.1148363d0,  56.792707d0,    49.1382d0,
     8     -16.9603104d0,   8.448301d0,   152.5268d0,
     9     -16.1765215d0,   1.978796d0,    98.2198d0,
     O      15.5567653d0,   8.863925d0,    97.4808d0,
     1      15.4846529d0,   0.186365d0,   221.5376d0,
     2      15.2150632d0,   8.996212d0,   168.2438d0,
     3      14.5047426d0,   6.771027d0,   161.1199d0,
     4     -14.3873316d0,  45.815258d0,    55.0196d0,
     5      13.1351419d0,  12.002811d0,   262.6495d0,
     6      12.8776311d0,  75.278220d0,   200.3284d0,
     7      11.9867234d0,  65.241008d0,   201.6651d0,
     8      11.9385578d0,  18.870667d0,   294.6547d0,
     9      11.7030822d0,  22.009553d0,    99.8233d0,
     O      11.6018181d0,  64.604291d0,   213.5577d0,
     1     -11.2617293d0,  11.498094d0,   154.1631d0,
     2     -10.4664199d0,   0.578834d0,   232.7153d0,
     3      10.4333970d0,   9.237738d0,   138.3034d0,
     4     -10.2377466d0,  49.747842d0,   204.6609d0,
     5      10.1934446d0,   2.147012d0,   106.5938d0,
     6     -10.1280191d0,   1.196895d0,   250.4676d0,
     7      10.0289441d0,   2.133898d0,   332.3345d0,
     8     -10.0034259d0,   0.173168d0,    27.3039d0/), (/3,78/) )
C****
      REAL*8 :: YM1950,SUMC,ARG,ESINPI,ECOSPI,PIE,PSI,FSINFD
      INTEGER :: I
C****
      YM1950 = YEAR-1950.
C****
C**** Obliquity from Table 1 (2):
C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
C****   OBLIQ  = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
C****
      SUMC = 0.
      DO I=1,47
        ARG  = PI180*(YM1950*TABL1(2,I)/3600.+TABL1(3,I))
        SUMC = SUMC + TABL1(1,I)*COS(ARG)
      END DO
      OBLIQ = 23.320556D0 + SUMC/3600.
!      OBLIQ  = OBLIQ*PI180 ! not needed for output in degrees
C****
C**** Eccentricity from Table 4 (1):
C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
C****
      ESINPI = 0.
      ECOSPI = 0.
      DO I=1,19
        ARG    = PI180*(YM1950*TABL4(2,I)/3600.+TABL4(3,I))
        ESINPI = ESINPI + TABL4(1,I)*SIN(ARG)
        ECOSPI = ECOSPI + TABL4(1,I)*COS(ARG)
      END DO
      ECCEN  = SQRT(ESINPI*ESINPI+ECOSPI*ECOSPI)
C****
C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
C****
      PIE = ATAN2(ESINPI,ECOSPI)
      FSINFD = 0.
      DO I=1,78
        ARG    = PI180*(YM1950*TABL5(2,I)/3600.+TABL5(3,I))
        FSINFD = FSINFD + TABL5(1,I)*SIN(ARG)
      END DO
      PSI    = PI180*(3.392506D0+(YM1950*50.439273D0+FSINFD)/3600.)
      OMEGVP = MOD(PIE+PSI+.5*TWOPI,TWOPI)
      IF(OMEGVP.lt.0.)  OMEGVP = OMEGVP + TWOPI
      OMEGVP = OMEGVP/PI180  ! for output in degrees
C****
      RETURN
      END SUBROUTINE ORBPAR
