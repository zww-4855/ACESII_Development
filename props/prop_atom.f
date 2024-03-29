










      subroutine  atom_gfac(znuc,gfactor,isotope,abun,niso)

c     Function will extract nuclear g-factor for a given atom

      implicit none
c  this need to modify
c     This include file contains atomic information.
c
c     For each isotope of each atom we have:
c     - mass
c     - abundance
c     - nuclear spin
c     - nuclear g-factor
c     - nuclear quadrupole moment
c
      integer          iii,jjj
      integer          atomlocator(2,118)
      double precision atomprop(5,323)
c
c     For each atom number in the periodic system: 
c     (# of isotopes, loc. in atomprop)
c
      data atomlocator/
c
c       1:Hydrogen        2:Helium          3:Lithium         4:Beryllium 
     &  2,1,              2,3,              2,5,              1,7,
c
c       5:Boron           6:Carbon          7:Nitrogen        8:Oxygen
     &  2,8,              2,10,             2,12,             3,14,
c
c       9:Fluorine       10:Neon           11:Sodium         12:Magnesium
     &  1,17,            3,18,             1,21,             3,22,
c
c      13:Aluminium      14:Silicon        15:Phosphorous    16:Sulphur 
     & 1,25,             3,26,             1,29,             4,30,
c
c      17:Chlorine       18:Argon          19:Potassium      20:Calcium 
     & 2,34,             3,36,             3,39,             6,42,
c
c      21:Scandium       22:Titanium       23:Vanadium       24:Chromium
     & 1,48,             5,49,             2,54,             4,56,
c
c      25:Manganese      26:Iron           27:Cobalt         28:Nickel
     & 1,60,             4,61,             1,65,             5,66,
c
c      29:Copper         30:Zinc           31:Gallium        32:Germanium
     & 2,71,             5,73,             2,78,             5,80,
c
c      33:Arsenic        34:Selenium       35:Bromine        36:Krypton
     & 1,85,             6,86,             2,92,             6,94,
c
c      37:Rubidium       38:Strontium      39:Yttrium        40:Zirconium
     & 2,100,            4,102,            1,106,            5,107,
c
c      41:Niobium        42:Molybdenum     43:Technetium     44:Ruthenium
     & 1,112,            7,113,            1,120,            7,121,
c
c      45:Rhodium        46:Palladium      47:Silver         48:Cadmium
     & 1,128,            6,129,            2,135,            8,137,
c
c      49:Indium         50:Tin            51:Antinomy       52:Tellurium 
     & 2,145,            10,147,           2,157,            8,159,
c
c      53:Iodine         54:Xenon          55:Caesium        56:Barium         
     & 1,167,            9,168,            1,177,            7,178,
c
c      57:Lanthanum      58:Cerium
     & 2,185,            4,187,
c
c      59:Praseodymium   60:Neodymiumium   61:Promethium     62:Samarium
     & 1,191,            7,192,            1,199,            7,200,
c
c      63:Europium       64:Gadolinium     65:Terbium        66:Dysprosium
     & 2,207,            7,209,            1,216,            7,217,
c
c      67:Holmium        68:Erbium         69:Thulium        70:Ytterbium 
     & 1,224,            6,225,            1,231,            7,232,
c
c      71:Lutetium       72:Hafnium        73:Tantalum       74:Tungsten
     & 2,239,            6,241,            2,247,            5,249,
c
c      75:Rhenium        76:Osmium         77:Iridium        78:Platinum
     & 2,254,            7,256,            2,263,            6,265,
c
c      79:Gold           80:Mercury        81:Thallium       82:Lead
     & 1,271,            7,272,            2,279,            4,281,
c
c      83:Bismuth        84:Polonium       85:Astatine       86:Radon
     & 1,285,            1,286,            1,287,            2,288,
c
c      87:Francium       88:Radium         89:Actinium       90:Thorium
     & 1,290,            1,291,            1,292,            1,293,
c
c      91:Protoactinium  92:Uranium        93:Neptunium      94:Plutonium
     & 1,294,            3,295,            1,298,            1,299,
c
c      95:Americium      96:Curium         97:Berkelium      98:Californium
     & 1,300,            1,301,            1,302,            1,303,
c
c      99:Einsteinium   100:Fermium       101:Mendelevium   102:Nobelium
     & 1,304,           1,305,            1,306,            1,307,
c
c     103:Lawrencium    104:Rutherfordium 105:Dubnium       106:Seaborgium
     &1,308,            1,309,            1,310,            1,311,
c
c     107:Bohrium       108:Hassium       109:Meitnerium    110:Darmstadtium
     &1,312,            1,313,            1,314,            1,315,
c
c     111:Roentgenium   112:Ununbium      113:Ununtrium     114:Ununquadium
     &1,316,            1,317,            1,318,            1,319,
c
c     115:Ununpentium   116:Ununhexium    117:Ununseptium   118:Ununoctium
     &1,320,            1,321,            1,322,            1,323/
c
c
c     Atomic properties table (for naturaly occuring isotopes), with:
c     - Atomic mass
c     - Abundance
c     - Nuclear spin
c     - Nuclear magnetic moment (divide by spin to get nuclear g-factor)
c     - Nuclear Quadrupole Moment (Barn)
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=1,20) /
c
c     Hydrogen
c
     &  1.00782504d0,   99.9885d0, 0.5d0,  2.7928446d0,      0.0d0,  !1
     &  2.01410178d0,    0.0115d0, 1.0d0,  0.8574376d0, 0.002860d0,    
c
c     Helium
c
     &  3.01602931d0,  0.000137d0, 0.5d0, -2.1276240d0,      0.0d0,  !3
     &  4.00260324d0, 99.999863d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Lithium
c
     &  6.01512140d0,      7.59d0, 1.0d0,  0.8220467d0,-0.000808d0,  !5
     &  7.01600300d0,     92.41d0, 1.5d0,  3.2564240d0,-0.040100d0,
c
c     Beryllium 
c
     &  9.01218220d0,     100.0d0, 1.5d0, -1.1779000d0, 0.052880d0,  !7
c
c     Boron
c 
     & 10.0129369d0,      19.9d0, 3.0d0,  1.8006500d0, 0.084590d0,   !8
     & 11.0093054d0,      80.1d0, 1.5d0,  2.6886370d0, 0.040590d0,
c     
c     Carbon
c
     & 12.0000000d0,     98.93d0, 0.0d0,        0.0d0,      0.0d0,  !10
     & 13.0033548d0,      1.07d0, 0.5d0,  0.7024110d0,      0.0d0,
c 
c     Nitrogen
c
     & 14.0030740d0,    99.632d0, 1.0d0,  0.4037607d0,  0.02044d0,  !12
     & 15.0001090d0,     0.368d0, 0.5d0, -0.2831892d0,      0.0d0,
c
c     Oxygen
c
     & 15.9949146d0,    99.757d0, 0.0d0,        0.0d0,      0.0d0,  !14
     & 16.9991312d0,     0.038d0, 2.5d0, -1.8938000d0, -0.02558d0,
     & 17.9991603d0,     0.205d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Fluorine       
c
     & 18.9984932d0,     100.0d0, 0.5d0,  2.6288670d0, -0.09420d0,  !17
c
c     Neon           
c
     & 19.9924356d0,     90.48d0, 0.0d0,        0.0d0,      0.0d0,  !18
     & 20.9938428d0,      0.27d0, 1.5d0, -0.6617960d0,  0.10155d0,
     & 21.9913831d0,      9.25d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=21,38) /
c
c     Sodium         
c
     & 22.9897677d0,     100.0d0, 1.5d0,  2.2175200d0,  0.10400d0,  !21
c
c     Magnesium
c
     & 23.9859423d0,     78.99d0, 0.0d0,        0.0d0,      0.0d0,  !22
     & 24.9858374d0,     10.00d0, 2.5d0, -0.8554600d0,  0.19940d0,
     & 25.9825937d0,     11.01d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Aluminium     
c
     & 26.9815386d0,     100.0d0, 2.5d0,  3.6415040d0,  0.13660d0,  !25
c
c     Silicon       
c
     & 27.9769271d0,   92.2297d0, 0.0d0,        0.0d0,      0.0d0,  !26
     & 28.9764949d0,    4.6832d0, 0.5d0, -0.5552900d0,      0.0d0,
     & 29.9737707d0,    3.0872d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Phosphorous   
c
     & 30.9737620d0,     100.0d0, 0.5d0,  1.1316000d0,      0.0d0,  !29
c
c     Sulphur 
c
     & 31.9720707d0,     94.93d0, 0.0d0,        0.0d0,      0.0d0,  !30
     & 32.9714584d0,      0.76d0, 1.5d0,  0.6438210d0, -0.06780d0,
     & 33.9678667d0,      4.29d0, 0.0d0,        0.0d0,      0.0d0,
     & 35.9670806d0,      0.02d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Chlorine      
c
     & 34.9688527d0,     75.78d0, 1.5d0,  0.8218736d0, -0.08165d0,  !34
     & 36.9659026d0,     24.22d0, 1.5d0,  0.6841230d0, -0.06435d0,
c
c     Argon         
c
     & 35.9675455d0,    0.3365d0, 0.0d0,        0.0d0,      0.0d0,  !36
     & 37.9627325d0,    0.0632d0, 0.0d0,        0.0d0,      0.0d0,
     & 39.9623837d0,   99.6003d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=39,99) /
c
c     Potassium     
c
     & 38.9637074d0,   93.2581d0, 1.5d0,  0.3914658d0,  0.05850d0,  !39
     & 39.9639992d0,    0.0117d0, 4.0d0, -1.2980990d0, -0.07300d0,
     & 40.9618254d0,    6.7302d0, 1.5d0,  0.2148699d0,  0.07110d0,
c
c     Calcium 
c
     & 39.9625906d0,    96.941d0, 0.0d0,        0.0d0,      0.0d0,  !42
     & 41.9586176d0,     0.647d0, 0.0d0,        0.0d0,      0.0d0,
     & 42.9587662d0,     0.135d0, 3.5d0, -1.3172700d0, -0.04080d0, 
     & 43.9554806d0,     2.086d0, 0.0d0,        0.0d0,      0.0d0,
     & 45.9536890d0,     0.004d0, 0.0d0,        0.0d0,      0.0d0,
     & 47.9525330d0,     0.187d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Scandium      
c
     & 44.9559100d0,     100.0d0, 3.5d0,  4.7564830d0, -0.22000d0,  !48
c
c     Titanium      
c
     & 45.9526294d0,      8.25d0, 0.0d0,        0.0d0,      0.0d0,  !49
     & 46.9517640d0,      7.44d0, 2.5d0, -0.7884800d0,  0.30200d0,
     & 47.9479473d0,     73.72d0, 0.0d0,        0.0d0,      0.0d0,
     & 48.9478711d0,      5.41d0, 3.5d0, -1.1041700d0,  0.24700d0,
     & 49.9447921d0,      5.18d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Vanadium      
c
     & 49.9471609d0,     0.250d0, 3.0d0,  3.3474500d0,  0.21000d0,  !54
     & 50.9439617d0,    99.750d0, 3.5d0,  5.1514000d0, -0.05200d0,
c
c     Chromium
c
     & 49.9460464d0,     4.345d0, 0.0d0,        0.0d0,      0.0d0,  !56
     & 51.9405098d0,    83.789d0, 0.0d0,        0.0d0,      0.0d0,
     & 52.9406513d0,     9.501d0, 1.5d0, -0.4745400d0, -0.15000d0,
     & 53.9388825d0,     2.365d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Manganese     
c
     & 54.9380471d0,     100.0d0, 2.5d0,  3.4532000d0,  0.33000d0,  !60
c
c     Iron          
c
     & 53.9396127d0,     5.845d0, 0.0d0,        0.0d0,      0.0d0,  !61
     & 55.9349393d0,    91.754d0, 0.0d0,        0.0d0,      0.0d0,
     & 56.9353958d0,     2.119d0, 0.5d0,  0.0906229d0,  0.16000d0,
     & 57.9332773d0,     0.282d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Cobalt        
c
     & 58.9331976d0,     100.0d0, 3.5d0,  4.6270000d0,  0.42000d0,  !65
c
c     Nickel
c
     & 57.9353462d0,   68.0769d0, 0.0d0,        0.0d0,      0.0d0,  !66
     & 59.9307884d0,   26.2231d0, 0.0d0,        0.0d0,      0.0d0,
     & 60.9310579d0,    1.1399d0, 1.5d0, -0.7500200d0,  0.16200d0,
     & 61.9283461d0,    3.6345d0, 0.0d0,        0.0d0,      0.0d0,
     & 63.9279679d0,    0.9256d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Copper        
c
     & 62.9295989d0,     69.17d0, 1.5d0,  2.2233000d0, -0.22000d0,  !71
     & 64.9277929d0,     30.83d0, 1.5d0,  2.3817000d0, -0.20400d0,
c
c     Zinc          
c
     & 63.9291448d0,     48.63d0, 0.0d0,        0.0d0,      0.0d0,  !73
     & 65.9260347d0,     28.90d0, 0.0d0,        0.0d0,      0.0d0,
     & 66.9271291d0,      4.10d0, 2.5d0,  0.8754790d0,  0.15000d0,
     & 67.9248459d0,     18.75d0, 0.0d0,        0.0d0,      0.0d0,
     & 69.9253250d0,      0.62d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Gallium       
c
     & 68.9255800d0,    60.108d0, 1.5d0,  2.0165900d0,  0.17100d0,  !78
     & 70.9247005d0,    39.892d0, 1.5d0,  2.5622700d0,  0.10700d0,
c
c     Germanium
c
     & 69.9242497d0,     20.84d0, 0.0d0,        0.0d0,      0.0d0,  !80
     & 71.9220789d0,     27.54d0, 0.0d0,        0.0d0,      0.0d0,
     & 72.9234626d0,      7.73d0, 4.5d0, -0.8794669d0, -0.19600d0,
     & 73.9211774d0,     36.28d0, 0.0d0,        0.0d0,      0.0d0,
     & 75.9214016d0,      7.61d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Arsenic       
c
     & 74.9215942d0,     100.0d0, 1.5d0,  1.4394700d0,  0.31400d0,  !85
c
c     Selenium      
c
     & 73.9224746d0,      0.89d0, 0.0d0,        0.0d0,      0.0d0,  !86
     & 75.9192120d0,      9.37d0, 0.0d0,        0.0d0,      0.0d0,
     & 76.9199125d0,      7.63d0, 0.5d0,  0.5350600d0,  0.76000d0,
     & 77.9173076d0,     23.77d0, 0.0d0,        0.0d0,      0.0d0,
     & 79.9165196d0,     49.61d0, 0.0d0,        0.0d0,      0.0d0,
     & 81.9166978d0,      8.73d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Bromine       
c
     & 78.9183361d0,     50.69d0, 1.5d0,  2.1063990d0,  0.31300d0,  !92
     & 80.9162890d0,     49.31d0, 1.5d0,  2.2705600d0,  0.26200d0,
c
c     Krypton
c
     & 77.9203960d0,      0.35d0, 0.0d0,        0.0d0,      0.0d0,  !94
     & 79.9163800d0,      2.28d0, 0.0d0,        0.0d0,      0.0d0,
     & 81.9134820d0,     11.58d0, 0.0d0,        0.0d0,      0.0d0,
     & 82.9141350d0,     11.49d0, 4.5d0, -0.9706690d0,  0.25900d0,
     & 83.9115070d0,     57.00d0, 0.0d0,        0.0d0,      0.0d0,
     & 85.9106160d0,     17.30d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=100,176) /
c
c     Rubidium      
c
     & 84.9117940d0,     71.17d0, 2.5d0,  1.3530300d0,  0.27600d0,  !100
     & 86.9091870d0,     27.83d0, 1.5d0,  2.7512400d0,  0.13350d0,
c
c     Strontium     
c
     & 84.9134300d0,      0.56d0, 0.0d0,        0.0d0,      0.0d0,  !102
     & 85.9092672d0,      9.86d0, 0.0d0,        0.0d0,      0.0d0,
     & 86.9088841d0,      7.00d0, 4.5d0, -1.0928300d0,  0.33500d0,
     & 87.9056188d0,     82.58d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Yttrium       
c
     & 88.9058490d0,     100.0d0, 0.5d0, -0.1374153d0,      0.0d0,  !106
c
c     Zirconium
c
     & 89.9047026d0,     51.45d0, 0.0d0,        0.0d0,      0.0d0,  !107
     & 90.9056439d0,     11.22d0, 2.5d0, -1.3036200d0, -0.17600d0,
     & 91.9050386d0,     17.15d0, 0.0d0,        0.0d0,      0.0d0,
     & 93.9063148d0,     17.38d0, 0.0d0,        0.0d0,      0.0d0,
     & 95.9082750d0,      2.80d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Niobium       
c
     & 92.9063772d0,     100.0d0, 4.5d0,  6.1705000d0, -0.32000d0,  !112
c
c     Molybdenum    
c
     & 91.9068090d0,     14.84d0, 0.0d0,        0.0d0,      0.0d0,  !113
     & 93.9050853d0,      9.25d0, 0.0d0,        0.0d0,      0.0d0,
     & 94.9058411d0,     15.92d0, 2.5d0, -0.9142000d0, -0.02200d0,
     & 95.9046785d0,     16.68d0, 0.0d0,        0.0d0,      0.0d0,
     & 96.9060205d0,      9.55d0, 2.5d0, -0.9335000d0,  0.25500d0,
     & 97.9054073d0,     24.13d0, 0.0d0,        0.0d0,      0.0d0,
     & 99.9074770d0,      9.63d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Technetium    
c
     & 98.9070000d0,     100.0d0, 4.5d0,  5.6847000d0, -0.12900d0,  !120
c
c     Ruthenium
c
     & 95.9075990d0,      5.54d0, 0.0d0,        0.0d0,      0.0d0,  !121
     & 97.9052870d0,      1.87d0, 0.0d0,        0.0d0,      0.0d0,
     & 98.9059389d0,     12.76d0, 2.5d0, -0.6413000d0,  0.07900d0,
     & 99.9042192d0,     12.60d0, 0.0d0,        0.0d0,      0.0d0,
     &100.9042192d0,     17.06d0, 2.5d0, -0.7189000d0,  0.45700d0,
     &101.9043485d0,     31.55d0, 0.0d0,        0.0d0,      0.0d0,
     &103.9054240d0,     18.62d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Rhodium       
c
     &102.9055000d0,     100.0d0, 0.5d0, -0.0884000d0,      0.0d0,  !128
c
c     Palladium     
c
     &101.9056340d0,      1.02d0, 0.0d0,        0.0d0,      0.0d0,  !129
     &103.9040290d0,     11.14d0, 0.0d0,        0.0d0,      0.0d0,
     &104.9050790d0,     22.33d0, 2.5d0, -0.6420000d0,  0.66000d0,
     &105.9034780d0,     27.33d0, 0.0d0,        0.0d0,      0.0d0,
     &107.9038950d0,     26.46d0, 0.0d0,        0.0d0,      0.0d0,
     &109.9051670d0,     11.72d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Silver        
c
     &106.9050920d0,    51.839d0, 0.5d0, -0.1135700d0,  0.98000d0,  !135
     &108.9047560d0,    48.161d0, 0.5d0, -0.1306905d0,      0.0d0,
c
c     Cadmium
c
     &105.9064610d0,      1.25d0, 0.0d0,        0.0d0,      0.0d0,  !137
     &107.9041760d0,      0.89d0, 0.0d0,        0.0d0,      0.0d0,
     &109.9030050d0,     12.49d0, 0.0d0,        0.0d0,      0.0d0,
     &110.9041820d0,     12.80d0, 0.5d0, -0.5948857d0, -0.85000d0,
     &111.9027570d0,     24.13d0, 0.0d0,        0.0d0,      0.0d0,
     &112.9044000d0,     12.22d0, 0.5d0, -0.6223005d0,      0.0d0,
     &113.9033570d0,     28.73d0, 0.0d0,        0.0d0,      0.0d0,
     &115.9047550d0,      7.49d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Indium        
c
     &112.9040610d0,      4.29d0, 4.5d0,  5.5289000d0,   0.7990d0,  !145
     &114.9038820d0,     95.71d0, 4.5d0,  5.5408000d0,   0.8100d0,
c
c     Tin           
c
     &111.9048260d0,      0.97d0, 0.0d0,        0.0d0,      0.0d0,  !147
     &113.9027840d0,      0.66d0, 0.0d0,        0.0d0,      0.0d0,
     &114.9033480d0,      0.34d0, 0.5d0, -0.9188400d0,      0.0d0,
     &115.9017470d0,     14.54d0, 0.0d0,        0.0d0,      0.0d0,
     &116.9029560d0,      7.68d0, 0.5d0, -1.0010500d0,      0.0d0,
     &117.9016090d0,     24.22d0, 0.0d0,        0.0d0,      0.0d0,
     &118.9033110d0,      8.59d0, 0.5d0, -1.0472900d0, -0.12800d0,
     &119.9021991d0,     32.58d0, 0.0d0,        0.0d0,      0.0d0,
     &121.9034404d0,      4.63d0, 0.0d0,        0.0d0,      0.0d0,
     &123.9052743d0,      5.79d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Antinomy      
c
     &120.9038212d0,     57.36d0, 2.5d0,  3.3634000d0, -0.36000d0,  !157
     &122.9042160d0,     42.64d0, 3.5d0,  2.5498000d0, -0.49000d0,
c
c     Tellurium 
c
     &119.9040480d0,      0.09d0, 0.0d0,        0.0d0,      0.0d0,  !159
     &121.9030500d0,      2.55d0, 0.0d0,        0.0d0,      0.0d0,
     &122.9042710d0,      0.89d0, 0.5d0, -0.7367900d0,      0.0d0,
     &123.9028180d0,      4.74d0, 0.0d0,        0.0d0,      0.0d0,
     &124.9044285d0,      7.07d0, 0.5d0, -0.8882800d0, -0.31000d0,
     &125.9033095d0,     18.84d0, 0.0d0,        0.0d0,      0.0d0,
     &127.9044630d0,     31.74d0, 0.0d0,        0.0d0,      0.0d0,
     &129.9062290d0,     34.08d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Iodine        
c
     &126.9044730d0,     100.0d0, 2.5d0,  2.8132800d0, -0.71000d0,  !167
c
c     Xenon         
c
     &123.9058942d0,      0.09d0, 0.0d0,        0.0d0,      0.0d0,  !168
     &125.9042810d0,      0.09d0, 0.0d0,        0.0d0,      0.0d0,
     &127.9035312d0,      1.92d0, 0.0d0,        0.0d0,      0.0d0,
     &128.9047801d0,     26.44d0, 0.5d0, -0.7779770d0, -0.39300d0,
     &129.9035094d0,      4.08d0, 0.0d0,        0.0d0,      0.0d0,
     &130.9050720d0,     21.18d0, 1.5d0,  0.6918610d0, -0.11400d0,
     &131.9041440d0,     26.89d0, 0.0d0,        0.0d0,      0.0d0,
     &133.9053950d0,     10.44d0, 0.0d0,        0.0d0,      0.0d0,
     &135.9072140d0,      8.87d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=177,238) /
c
c     Caesium       
c
     &132.9054290d0,     100.0d0, 3.5d0,  2.5829240d0, -0.00343d0,  !177
c
c     Barium         
c
     &129.9062820d0,     0.106d0, 0.0d0,        0.0d0,      0.0d0,  !178
     &131.9050420d0,     0.101d0, 0.0d0,        0.0d0,      0.0d0,
     &133.9044860d0,     2.417d0, 0.0d0,        0.0d0,      0.0d0,
     &134.9056650d0,     6.592d0, 1.5d0,  0.8379430d0,  0.16000d0,
     &135.9045530d0,     7.854d0, 0.0d0,        0.0d0,      0.0d0,
     &136.9058120d0,    11.232d0, 1.5d0,  0.9373650d0,  0.24500d0,
     &137.9052320d0,    71.698d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Lanthanum     
c
     &137.9071050d0,     0.090d0, 5.0d0,  3.7139000d0,  0.45000d0,  !185
     &138.9063470d0,    99.910d0, 3.5d0,  2.7832000d0,  0.20000d0,
c
c     Cerium
c
     &135.9071400d0,     0.185d0, 0.0d0,        0.0d0,      0.0d0,  !187
     &137.9059850d0,     0.251d0, 0.0d0,        0.0d0,      0.0d0,
     &139.9054330d0,    88.450d0, 0.0d0,        0.0d0,      0.0d0,
     &141.9092410d0,    11.114d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Praseodymium  
c
     &140.9076470d0,     100.0d0, 2.5d0,  4.1360000d0, -0.05890d0,  !191
c
c     Neodymiumium  
c
     &141.9077190d0,      27.2d0, 0.0d0,        0.0d0,      0.0d0,  !192
     &142.9098100d0,      12.2d0, 3.5d0, -1.0650000d0, -0.63000d0,
     &143.9100830d0,      23.8d0, 0.0d0,        0.0d0,      0.0d0,
     &144.9125700d0,       8.3d0, 3.5d0, -0.6560000d0, -0.33000d0,
     &145.9131130d0,      17.2d0, 0.0d0,        0.0d0,      0.0d0,
     &147.9168890d0,       5.7d0, 0.0d0,        0.0d0,      0.0d0,
     &149.9208870d0,       5.6d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Promethium    
c
     &146.9151340d0,     100.0d0, 3.5d0,  2.6000000d0,  0.74000d0,  !199
c
c     Samarium
c
     &143.9119980d0,      3.07d0, 0.0d0,        0.0d0,      0.0d0,  !200
     &146.9148940d0,     14.99d0, 3.5d0, -0.8149000d0, -0.25900d0,
     &147.9148190d0,     11.24d0, 0.0d0,        0.0d0,      0.0d0,
     &148.9171800d0,     13.82d0, 3.5d0, -0.6718000d0,  0.07400d0,
     &149.9172730d0,      7.38d0, 0.0d0,        0.0d0,      0.0d0,
     &151.9197280d0,     26.75d0, 0.0d0,        0.0d0,      0.0d0,
     &153.9222050d0,     22.75d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Europium      
c
     &150.9197020d0,     47.81d0, 2.5d0,  3.4718000d0,  0.90300d0,  !207
     &152.9212250d0,     52.19d0, 2.5d0,  1.5331000d0,  2.41200d0,
c
c     Gadolinium    
c
     &151.9197860d0,      0.20d0, 0.0d0,        0.0d0,      0.0d0,  !209
     &153.9208610d0,      2.18d0, 0.0d0,        0.0d0,      0.0d0,
     &154.9226180d0,     14.80d0, 1.5d0, -0.2591000d0,  1.27000d0,
     &155.9221180d0,     20.47d0, 0.0d0,        0.0d0,      0.0d0,
     &156.9239560d0,     15.65d0, 1.5d0, -0.3399000d0,  1.35000d0,
     &157.9240190d0,     24.84d0, 0.0d0,        0.0d0,      0.0d0,
     &159.9270490d0,     21.86d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Terbium       
c
     &158.9253420d0,     100.0d0, 1.5d0,  2.0140000d0,  1.43200d0,  !216
c
c     Dysprosium
c
     &155.9242770d0,      0.06d0, 0.0d0,        0.0d0,      0.0d0,  !217
     &157.9244030d0,      0.10d0, 0.0d0,        0.0d0,      0.0d0,
     &159.9251930d0,      2.34d0, 0.0d0,        0.0d0,      0.0d0,
     &160.9269300d0,     18.91d0, 2.5d0, -0.4806000d0,  2.50700d0,
     &161.9267950d0,     25.51d0, 0.0d0,        0.0d0,      0.0d0,
     &162.9287280d0,     24.90d0, 2.5d0,  0.6726000d0,  2.64800d0,
     &163.9291710d0,     28.18d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Holmium       
c
     &164.9303190d0,     100.0d0, 3.5d0,  4.1730000d0,  3.58000d0,  !224
c
c     Erbium        
c
     &161.9287750d0,      0.14d0, 0.0d0,        0.0d0,      0.0d0,  !225
     &163.9291980d0,      1.61d0, 0.0d0,        0.0d0,      0.0d0,
     &165.9302900d0,     33.61d0, 0.0d0,        0.0d0,      0.0d0,
     &166.9320460d0,     22.93d0, 3.5d0, -0.5665000d0,  3.56500d0,
     &167.9323680d0,     26.78d0, 0.0d0,        0.0d0,      0.0d0,
     &169.9354610d0,     14.93d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Thulium       
c
     &168.9342120d0,     100.0d0, 0.5d0, -0.2316000d0, -1.20000d0,  !231
c
c     Ytterbium 
c
     &167.9338940d0,      0.13d0, 0.0d0,        0.0d0,      0.0d0,  !232
     &169.9347590d0,      3.04d0, 0.0d0,        0.0d0,      0.0d0,
     &170.9363230d0,     14.28d0, 0.5d0,  0.4919000d0,      0.0d0,
     &171.9363780d0,     21.83d0, 0.0d0,        0.0d0,      0.0d0,
     &172.9382080d0,     16.13d0, 2.5d0, -0.6776000d0,  2.80000d0,
     &173.9388590d0,     31.83d0, 0.0d0,        0.0d0,      0.0d0,
     &175.9425640d0,     12.76d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=239,289) /
c
c     Lutetium      
c
     &174.9407700d0,     97.41d0, 3.5d0,  2.2327000d0,  3.49000d0,  !239
     &175.9426790d0,      2.59d0, 7.0d0,  3.1900000d0,  4.97000d0,
c
c     Hafnium       
c
     &173.9400440d0,      0.16d0, 0.0d0,        0.0d0,      0.0d0,  !241
     &175.9414060d0,      5.26d0, 0.0d0,        0.0d0,      0.0d0,
     &176.9432170d0,     18.60d0, 3.5d0,  0.7936000d0,  3.36500d0,
     &177.9436960d0,     27.28d0, 0.0d0,        0.0d0,      0.0d0,
     &178.9458122d0,     13.62d0, 4.5d0, -0.6409000d0,  3.79300d0,
     &179.9465457d0,     35.08d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Tantalum      
c
     &179.9474620d0,     0.012d0, 8.0d0,        0.0d0,      0.0d0,  !247
     &180.9479920d0,    99.988d0, 3.5d0,  2.3710000d0,  3.17000d0,
c
c     Tungsten
c
     &179.9467010d0,      0.12d0, 0.0d0,        0.0d0,      0.0d0,  !249
     &181.9482020d0,     26.50d0, 0.0d0,        0.0d0,      0.0d0,
     &182.9502200d0,     14.31d0, 0.5d0,  0.1177847d0,      0.0d0,
     &183.9509280d0,     30.64d0, 0.0d0,        0.0d0,      0.0d0,
     &185.9543570d0,     28.43d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Rhenium       
c
     &184.9529510d0,     37.40d0, 2.5d0,  3.1871000d0,  2.18000d0,  !254
     &186.9557440d0,     62.60d0, 2.5d0,  3.2197000d0,  2.07000d0,
c
c     Osmium        
c
     &183.9524880d0,      0.02d0, 0.0d0,        0.0d0,      0.0d0,  !256
     &185.9538300d0,      1.59d0, 0.0d0,        0.0d0,      0.0d0,
     &186.9557410d0,      1.96d0, 0.5d0, 0.06465185d0,      0.0d0,
     &187.9558300d0,     13.23d0, 0.0d0,        0.0d0,      0.0d0,
     &188.9581370d0,     16.15d0, 1.5d0,  0.6599330d0,  0.85600d0,
     &189.9584360d0,     26.26d0, 0.0d0,        0.0d0,      0.0d0,
     &191.9614670d0,     40.78d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Iridium       
c
     &190.9605840d0,      37.3d0, 1.5d0,  0.1462000d0,  0.81600d0,  !263
     &192.9629170d0,      62.7d0, 1.5d0,  0.1592000d0,  0.75100d0,
c
c     Platinum
c
     &189.9599170d0,     0.014d0, 0.0d0,        0.0d0,      0.0d0,  !265
     &191.9610190d0,     0.782d0, 0.0d0,        0.0d0,      0.0d0,
     &193.9626550d0,    32.967d0, 0.0d0,        0.0d0,      0.0d0,
     &194.9647660d0,    33.832d0, 0.5d0,  0.6095000d0,      0.0d0,
     &195.9649260d0,    25.242d0, 0.0d0,        0.0d0,      0.0d0,
     &197.9678690d0,     7.163d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Gold          
c
     &196.9665430d0,     100.0d0, 1.5d0,  0.1481590d0,  0.54700d0,  !271
c
c     Mercury       
c
     &195.9658070d0,      0.15d0, 0.0d0,        0.0d0,      0.0d0,  !272
     &197.9667430d0,      9.97d0, 0.0d0,        0.0d0,      0.0d0,
     &198.9682540d0,     16.87d0, 0.5d0,  0.5058852d0,  0.67400d0,
     &199.9683000d0,     23.10d0, 0.0d0,        0.0d0,      0.0d0,
     &200.9702770d0,     13.18d0, 1.5d0, -0.5602250d0,  0.38600d0,
     &201.9706170d0,     29.86d0, 0.0d0,        0.0d0,      0.0d0,
     &203.9734670d0,      6.87d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Thallium      
c
     &202.9723200d0,    29.524d0, 0.5d0,  1.6222570d0,      0.0d0,  !279
     &204.9744010d0,    70.476d0, 0.5d0,  1.6382135d0,      0.0d0,  
c
c     Lead
c
     &203.9730200d0,       1.4d0, 0.0d0,        0.0d0,      0.0d0,  !281
     &205.9744400d0,      24.1d0, 0.0d0,        0.0d0,      0.0d0,
     &206.9758720d0,      22.1d0, 0.5d0,  0.5821900d0, -0.26900d0,
     &207.9766270d0,      52.4d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Bismuth       
c
     &208.9803740d0,     100.0d0, 4.5d0,  4.1106000d0, -0.51600d0,  !285
c
c     Polonium      
c
     &208.9824040d0,       0.0d0, 0.5d0,        0.0d0,      0.0d0,  !286
c
c     Astatine
c
     &209.9871260d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !287
c
c     Radon
c
     &208.9903800d0,       0.0d0, 2.5d0,  0.8388000d0,  0.31100d0,  !288
     &222.0175710d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0/
c
      data ((atomprop(iii,jjj),iii=1,5),jjj=290,323) /
c
c     Francium      
c
     &223.0197330d0,       0.0d0, 1.5d0,  1.1700000d0,  1.17000d0,  !290
c
c     Radium        
c
     &226.0254030d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !291
c
c     Actinium      
c
     &227.0277500d0,       0.0d0, 1.5d0,  1.1000000d0,  1.70000d0,  !292
c
c     Thorium
c
     &232.0380508d0,     100.0d0, 0.0d0,        0.0d0,      0.0d0,  !293
c
c     Protoactinium 
c
     &231.0358800d0,     100.0d0, 1.5d0,  2.0100000d0, -1.72000d0,  !294
c
c     Uranium       
c
     &234.0409468d0,    0.0055d0, 0.0d0,        0.0d0,      0.0d0,  !295
     &235.0439242d0,    0.7200d0, 3.5d0, -0.3500000d0,  4.93600d0,
     &238.0507847d0,   99.2745d0, 0.0d0,        0.0d0,      0.0d0,
c
c     Neptunium     
c
     &237.0481678d0,       0.0d0, 2.5d0,  3.1400000d0,  3.88600d0,  !298
c
c     Plutonium
c
     &244.0641990d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !299
c
c     Americium     
c
     &243.0613750d0,       0.0d0, 2.5d0,  1.6100000d0,  4.21000d0,  !300
c
c     Curium        
c
     &247.0703470d0,       0.0d0, 4.5d0,  0.3700000d0,      0.0d0,  !301
c
c     Berkelium     
c
     &247.0703000d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !302
c
c     Californium
c
     &251.0795800d0,       0.0d0, 0.5d0,        0.0d0,      0.0d0,  !303
c
c     Einsteinium   
c
     &252.0829440d0,       0.0d0, 5.0d0,        0.0d0,  6.70000d0,  !304
c
c     Fermium       
c
     &257.0950990d0,       0.0d0, 4.5d0,        0.0d0,      0.0d0,  !305
c
c     Mendelevium   
c
     &258.0985700d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !306
c
c     Nobelium
c
     &259.1009310d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !307
c
c     Lawrencium    
c
     &260.1053200d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !308
c
c     Rutherfordium 
c
     &261.1086900d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !309
c
c     Dubnium       
c
     &262.1137600d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !310
c
c     Seaborgium
c
     &263.1182200d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !311
c
c     Bohrium       
c
     &262.1229300d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !312
c
c     Hassium       
c
     &265.1301600d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !313
c
c     Meitnerium    
c
     &268.1435000d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !314
c
c     Darmstadtium
c
     &269.1451000d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !315
c
c     Roentgenium   
c
     &272.1535000d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !316
c
c     Ununbium      
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !317
c
c     Ununtrium     
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !318
c
c     Ununquadium
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !319
c
c     Ununpentium   
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !320
c
c     Ununhexium    
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !321
c
c     Ununseptium   
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0,  !322
c
c     Ununoctium
c
     &        0.0d0,       0.0d0, 0.0d0,        0.0d0,      0.0d0/  !323
c
c
      double precision znuc       ! [in]  nuclear charge
      double precision gfactor    ! [out] nuclear g-factor
c
      integer iatom, niso, i, iloc, isotope
      double precision abun, spin, gfac
c
      iatom = atomlocator(2,nint(znuc))
      niso = atomlocator(1,nint(znuc))
c
c     Determine isotope with spin and highest abundance
c
      abun = 0.0d0
      spin = 0.0d0
      gfac = 0.0d0
      isotope = 0
      do i = 1, niso
         iloc = iatom+i-1
         if (atomprop(3,iloc).ne.0.0d0) then
            if (atomprop(2,iloc).ge.abun) then
               isotope = nint(atomprop(1,iloc))
               abun = atomprop(2,iloc)
               spin = atomprop(3,iloc)
               gfac = atomprop(4,iloc)
            endif
         endif
      enddo
      gfactor = 0.0d0
      if (spin.gt.0.0d0) gfactor = gfac/spin
c
c
      return
      end      
