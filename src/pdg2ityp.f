cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ityp(pdgid)
c
c
c     input pdgid  : Particle-ID according to Particle Data Group  
c 
c     converts PDG-ID to UrQMD-ID 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer pdgid
      integer ityp
      integer charge

      integer tab_size
      parameter (TAB_SIZE = 183)

      logical anti
      integer abs_ityp
      integer norm_iso3
      integer idtab(3,TAB_SIZE)
      integer first
      integer last
      integer next
      integer isoit

      integer i

c----------------------------------------------------------------------- 

      data idtab/
c Neutron
     .       1,  0,  2112,  
c Proton
     .       1,  1,  2212,
c N*
     .       2,  0,  12112,       2,  1,  12212,
     .       3,  0,   1214,       3,  1,   2124, 
     .       4,  0,  22112,       4,  1,  22212,
     .       5,  0,  32112,       5,  1,  32212,
     .       6,  0, 102116,       6,  1, 102216,
     .       7,  0,   2116,       7,  1,   2216,
     .       8,  0,  12116,       8,  1,  12216,
     .       9,  0,  21214,       9,  1,  22124, 
     .      10,  0,  42112,      10,  1,  42212, 
     .      11,  0,  31214,      11,  1,  32124, 
     .      12,  0, 9902114,     12,  1, 9902214, 
     .      13,  0, 9912114,     13,  1, 9912214,   
     .      14,  0,   1218,      14,  1,   2128, 
c Delta
     .      17, -1,   1114,  17, 0,  2114,   17, 1,  2214,   17, 2,   2224,
     .      18, -1,  31114,  18, 0, 32114,   18, 1, 32214,   18, 2,  32224,
     .      19, -1,   1112,  19, 0,  1212,   19, 1,  2212,   19, 2,   2222,
     .      20, -1,  11114,  20, 0, 12114,   20, 1, 12214,   20, 2,  12224,
     .      21, -1,   1116,  21, 0,  1216,   21, 1,  2126,   21, 2,   2226,
     .      22, -1,  21112,  22, 0, 21212,   22, 1, 22122,   22, 2,  22222,
     .      23, -1,  21114,  23, 0, 22114,   23, 1, 22214,   23, 2,  22224,
     .      24, -1,  11116,  24, 0, 11216,   24, 1, 12126,   24, 2,  12226,
     .      25, -1,   1118,  25, 0,  2118,   25, 1,  2218,   25, 2,   2228,
     .      26, -1, 201118,  26, 0, 202118,  26, 1, 202218,  26, 2, 202228,
c Lambda
     .      27,  0,  3122,
     .      28,  0, 13122,   
     .      29,  0,  3124,   
     .      30,  0, 23122,   
     .      31,  0, 33122,
     .      32,  0, 13124,
     .      33,  0, 43122,   
     .      34,  0, 53122,   
     .      35,  0,  3126,   
     .      36,  0, 13126,   
     .      37,  0, 23124,   
     .      38,  0,  3128,   
     .      39,  0, 23126,   
c Sigma
     .      40, -1,  3112,    40,  0,  3212,    40,  1,  3222,
     .      41, -1,  3114,    41,  0,  3214,    41,  1,  3224,
     .      42, -1, 13112,    42,  0, 13212,    42,  1, 13222,
     .      43, -1, 13114,    43,  0, 13214,    43,  1, 13224,
     .      44, -1, 23112,    44,  0, 23212,    44,  1, 23222,
     .      45, -1,  3116,    45,  0,  3216,    45,  1,  3226,
     .      46, -1, 13116,    46,  0, 13216,    46,  1, 13226,
     .      47, -1, 23114,    47,  0, 23214,    47,  1, 23224,
     .      48, -1,  3118,    48,  0,  3218,    48,  1,  3228,
c Xi
     .      49, -1,  3312,     49,  1,  3322,
     .      50, -1,  3314,     50,  1,  3324,
     .      52, -1, 13314,     52,  1, 13324,
c Omega
     .      55,  0,  3334,
c gamma
     .     100,  0,    22,    100,  7, 9000221,  100,  8, 11,   100,  9, 13,
c pion
     .     101, -1,  -211,    101,  0,   111,    101,  1,   211,
c eta
     .     102,  0,   221,
c omega
     .     103,  0,   223,
c rho
     .     104, -1,  -213,    104,  0,   113,    104,  1,   213,
c f0(980)
     .     105,  0, 9010221,
c kaon
     .     106,  0,   311,    106,  1,   321,
c eta'
     .     107,  0,   331,
c k*(892)
     .     108,  0,   313,    108,  1,   323,
c phi
     .     109,  0,   333,
c k0*(1430)
     .     110,  0, 10313,    110,  1, 10323,
c a0(980)
     .     111, -2, -9000211,    111,  0, 900111,    111,  2, 900211,
c f0(1370)
     .     112,  0, 10221,
c k1(1270)
     .     113, -1, 10313,    113,  1, 10323,
c a1(1260)
     .     114, -2,-20213,    114,  0, 20113,    114,  2, 20213,
c f1(1285)
     .     115,  0, 20223,
c f1'(1510)
     .     116,  0, 40223,
c k2*(1430)
     .     117, -1,   315,    117,  1,   325,
c a2(1329)
     .     118, -2,  -215,    118,  0,   115,    118,  2,   215,
c f2(1270)
     .     119,  0,  225,
c f2'(1525)
     .     120,  0,   335,
c k1(1400)
     .     121, -1, 20313,    121,  1, 20323,
c b1
     .     122, -2,-10213,    122,  0, 10113,    122,  2, 10213, 
c h1
     .     123,  0, 10223,
c K* (1410)
     .     125, -1, 100313,   125,  1, 100323,
c rho (1450)
     .     126, -2,-40213,    126,  0, 40113,    126,  2, 40213,
c omega (1420)
     .     127,  0, 50223,
c phi(1680)
     .     128,  0, 10333,
c k*(1680)
     .     129, -1, 40313,    129,  1, 40323,
c rho(1700)
     .     130, -2,-30213,    130,  0, 30113,    130,  2, 30213,
c omega(1600)
     .     131,  0, 60223,
c phi(1850)     
     .     132,  0,   337,
c D
     .     133,  -1,   421,    133,   1,   411,  
c D*
     .     134,  -1, 10421,    134,   1, 10411,  
c J/Psi
     .     135,  0, 443,
c Chi_c
     .     136,  0, 100443,
c Psi'
     .     137,  0, 10441,
c Ds
     .     138,   0,   431,  
c Ds*
     .     139,   0,   433 /

c----------------------------------------------------------------------- 

      ityp=0
      charge=0

      do 101 i=1,TAB_SIZE
        if (idtab(3,i).eq.pdgid) then
          ityp=idtab(1,i)
          charge=idtab(2,i)
c          write(0,*)'old, new ityp =',pdgid,ityp !Debug only
          goto 200
        endif
 101  continue


      do 102 i=1,TAB_SIZE
        if (idtab(3,i).eq.abs(pdgid)) then
          ityp=-idtab(1,i)
          charge=-idtab(2,i)
c          write(0,*)'old, new ityp =',pdgid,ityp !Debug only
          goto 200
        endif
 102  continue

      if(pdgid.gt.9900000.OR.
     &  (pdgid.gt.1000.AND.pdgid.lt.9999).OR.
     &  (pdgid.gt.11000.AND.pdgid.lt.19999).OR.
     &  (pdgid.gt.21000.AND.pdgid.lt.29999).OR.
     &  (pdgid.gt.31000.AND.pdgid.lt.39999).OR.
     &  (pdgid.gt.41000.AND.pdgid.lt.49999).OR.
     &  (pdgid.gt.51000.AND.pdgid.lt.59999).OR.
     &  (pdgid.gt.101000.AND.pdgid.lt.119999).OR.
     &  (pdgid.gt.201000.AND.pdgid.lt.219999)) then
          ityp=60
          charge=0 
          goto 200
      endif

      if(abs(pdgid).gt.9900000.OR.
     &  (abs(pdgid).gt.1000.AND.abs(pdgid).lt.9999).OR.
     &  (abs(pdgid).gt.11000.AND.abs(pdgid).lt.19999).OR.
     &  (abs(pdgid).gt.21000.AND.abs(pdgid).lt.29999).OR.
     &  (abs(pdgid).gt.31000.AND.abs(pdgid).lt.39999).OR.
     &  (abs(pdgid).gt.41000.AND.abs(pdgid).lt.49999).OR.
     &  (abs(pdgid).gt.51000.AND.abs(pdgid).lt.59999).OR.
     &  (abs(pdgid).gt.101000.AND.abs(pdgid).lt.119999).OR.
     &  (abs(pdgid).gt.201000.AND.abs(pdgid).lt.219999)) then
          ityp=-60
          charge=0 
          goto 200
      endif

c IF NOT OTHERWISE DEFINED, ASSUME PARTICLE IS MESON AND UNCHARGED
      if(abs(pdgid).lt.10000000) then
       ityp=140
       charge=0
       goto 200
      endif

      write(0,*)'In function ITYP: Error in tablelookup'
      write(0,*)'PDG-Ityp =',pdgid
      ityp = 0
      stop
    
 200  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function charge(pdgid)
c
c
c     input pdgid  : Particle-ID according to Particle Data Group  
c 
c     converts PDG-ID to UrQMD-ID 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer pdgid
      integer ityp
      integer charge

      integer tab_size
      parameter (TAB_SIZE = 183)

      logical anti
      integer abs_ityp
      integer norm_iso3
      integer idtab(3,TAB_SIZE)
      integer first
      integer last
      integer next
      integer isoit

      integer i

c----------------------------------------------------------------------- 


      data idtab/
c Neutron
     .       1,  0,  2112,  
c Proton
     .       1,  1,  2212,
c N*
     .       2,  0,  12112,       2,  1,  12212,
     .       3,  0,   1214,       3,  1,   2124, 
     .       4,  0,  22112,       4,  1,  22212,
     .       5,  0,  32112,       5,  1,  32212,
     .       6,  0, 102116,       6,  1, 102216,
     .       7,  0,   2116,       7,  1,   2216,
     .       8,  0,  12116,       8,  1,  12216,
     .       9,  0,  21214,       9,  1,  22124, 
     .      10,  0,  42112,      10,  1,  42212, 
     .      11,  0,  31214,      11,  1,  32124, 
     .      12,  0, 9902114,     12,  1, 9902214, 
     .      13,  0, 9912114,     13,  1, 9912214,   
     .      14,  0,   1218,      14,  1,   2128, 
c Delta
     .      17, -1,   1114,  17, 0,  2114,   17, 1,  2214,   17, 2,   2224,
     .      18, -1,  31114,  18, 0, 32114,   18, 1, 32214,   18, 2,  32224,
     .      19, -1,   1112,  19, 0,  1212,   19, 1,  2212,   19, 2,   2222,
     .      20, -1,  11114,  20, 0, 12114,   20, 1, 12214,   20, 2,  12224,
     .      21, -1,   1116,  21, 0,  1216,   21, 1,  2126,   21, 2,   2226,
     .      22, -1,  21112,  22, 0, 21212,   22, 1, 22122,   22, 2,  22222,
     .      23, -1,  21114,  23, 0, 22114,   23, 1, 22214,   23, 2,  22224,
     .      24, -1,  11116,  24, 0, 11216,   24, 1, 12126,   24, 2,  12226,
     .      25, -1,   1118,  25, 0,  2118,   25, 1,  2218,   25, 2,   2228,
     .      26, -1, 201118,  26, 0, 202118,  26, 1, 202218,  26, 2, 202228,
c Lambda
     .      27,  0,  3122,
     .      28,  0, 13122,   
     .      29,  0,  3124,   
     .      30,  0, 23122,   
     .      31,  0, 33122,
     .      32,  0, 13124,
     .      33,  0, 43122,   
     .      34,  0, 53122,   
     .      35,  0,  3126,   
     .      36,  0, 13126,   
     .      37,  0, 23124,   
     .      38,  0,  3128,   
     .      39,  0, 23126,   
c Sigma
     .      40, -1,  3112,    40,  0,  3212,    40,  1,  3222,
     .      41, -1,  3114,    41,  0,  3214,    41,  1,  3224,
     .      42, -1, 13112,    42,  0, 13212,    42,  1, 13222,
     .      43, -1, 13114,    43,  0, 13214,    43,  1, 13224,
     .      44, -1, 23112,    44,  0, 23212,    44,  1, 23222,
     .      45, -1,  3116,    45,  0,  3216,    45,  1,  3226,
     .      46, -1, 13116,    46,  0, 13216,    46,  1, 13226,
     .      47, -1, 23114,    47,  0, 23214,    47,  1, 23224,
     .      48, -1,  3118,    48,  0,  3218,    48,  1,  3228,
c Xi
     .      49, -1,  3312,     49,  1,  3322,
     .      50, -1,  3314,     50,  1,  3324,
     .      52, -1, 13314,     52,  1, 13324,
c Omega
     .      55,  0,  3334,
c gamma
     .     100,  0,    22,    100,  7, 9000221,  100,  8, 11,   100,  9, 13,
c pion
     .     101, -1,  -211,    101,  0,   111,    101,  1,   211,
c eta
     .     102,  0,   221,
c omega
     .     103,  0,   223,
c rho
     .     104, -1,  -213,    104,  0,   113,    104,  1,   213,
c f0(980)
     .     105,  0, 9010221,
c kaon
     .     106,  0,   311,    106,  1,   321,
c eta'
     .     107,  0,   331,
c k*(892)
     .     108,  0,   313,    108,  1,   323,
c phi
     .     109,  0,   333,
c k0*(1430)
     .     110,  0, 10313,    110,  1, 10323,
c a0(980)
     .     111, -2, -9000211,    111,  0, 900111,    111,  2, 900211,
c f0(1370)
     .     112,  0, 10221,
c k1(1270)
     .     113, -1, 10313,    113,  1, 10323,
c a1(1260)
     .     114, -2,-20213,    114,  0, 20113,    114,  2, 20213,
c f1(1285)
     .     115,  0, 20223,
c f1'(1510)
     .     116,  0, 40223,
c k2*(1430)
     .     117, -1,   315,    117,  1,   325,
c a2(1329)
     .     118, -2,  -215,    118,  0,   115,    118,  2,   215,
c f2(1270)
     .     119,  0,  225,
c f2'(1525)
     .     120,  0,   335,
c k1(1400)
     .     121, -1, 20313,    121,  1, 20323,
c b1
     .     122, -2,-10213,    122,  0, 10113,    122,  2, 10213, 
c h1
     .     123,  0, 10223,
c K* (1410)
     .     125, -1, 100313,   125,  1, 100323,
c rho (1450)
     .     126, -2,-40213,    126,  0, 40113,    126,  2, 40213,
c omega (1420)
     .     127,  0, 50223,
c phi(1680)
     .     128,  0, 10333,
c k*(1680)
     .     129, -1, 40313,    129,  1, 40323,
c rho(1700)
     .     130, -2,-30213,    130,  0, 30113,    130,  2, 30213,
c omega(1600)
     .     131,  0, 60223,
c phi(1850)     
     .     132,  0,   337,
c D
     .     133,  -1,   421,    133,   1,   411,  
c D*
     .     134,  -1, 10421,    134,   1, 10411,  
c J/Psi
     .     135,  0, 443,
c Chi_c
     .     136,  0, 100443,
c Psi'
     .     137,  0, 10441,
c Ds
     .     138,   0,   431,  
c Ds*
     .     139,   0,   433 /

c----------------------------------------------------------------------- 
      ityp=0
      charge=0

      do 101 i=1,TAB_SIZE
        if (idtab(3,i).eq.pdgid) then
          ityp=idtab(1,i)
          charge=idtab(2,i)
c          write(0,*)'charge =',ityp,charge !Debug only
          goto 200
        endif
 101  continue

c..If not in list, check for ANTIPARTICLES
      do 102 i=1,TAB_SIZE
        if (idtab(3,i).eq.abs(pdgid)) then
          ityp=-idtab(1,i)
          charge=-idtab(2,i)
c          write(0,*)'charge =',ityp,charge !Debug only
          goto 200
        endif
 102  continue


      if(pdgid.gt.9900000.OR.
     &  (pdgid.gt.1000.AND.pdgid.lt.9999).OR.
     &  (pdgid.gt.11000.AND.pdgid.lt.19999).OR.
     &  (pdgid.gt.21000.AND.pdgid.lt.29999).OR.
     &  (pdgid.gt.31000.AND.pdgid.lt.39999).OR.
     &  (pdgid.gt.41000.AND.pdgid.lt.49999).OR.
     &  (pdgid.gt.51000.AND.pdgid.lt.59999).OR.
     &  (pdgid.gt.101000.AND.pdgid.lt.119999).OR.
     &  (pdgid.gt.201000.AND.pdgid.lt.219999)) then
          ityp=60
          charge=0 
          goto 200
      endif

      if(abs(pdgid).gt.9900000.OR.
     &  (abs(pdgid).gt.1000.AND.abs(pdgid).lt.9999).OR.
     &  (abs(pdgid).gt.11000.AND.abs(pdgid).lt.19999).OR.
     &  (abs(pdgid).gt.21000.AND.abs(pdgid).lt.29999).OR.
     &  (abs(pdgid).gt.31000.AND.abs(pdgid).lt.39999).OR.
     &  (abs(pdgid).gt.41000.AND.abs(pdgid).lt.49999).OR.
     &  (abs(pdgid).gt.51000.AND.abs(pdgid).lt.59999).OR.
     &  (abs(pdgid).gt.101000.AND.abs(pdgid).lt.119999).OR.
     &  (abs(pdgid).gt.201000.AND.abs(pdgid).lt.219999)) then
          ityp=-60
          charge=0 
          goto 200
      endif

c IF NOT OTHERWISE DEFINED, ASSUME PARTICLE IS MESON AND UNCHARGED
      if(abs(pdgid).lt.10000000) then
       ityp=140
       charge=0
       goto 200
      endif

      write(0,*)'In function ITYP: Error in tablelookup'
      write(0,*)'PDG-Ityp =',pdgid
      ityp = 0
      stop

      
 200  continue

      return
      end

