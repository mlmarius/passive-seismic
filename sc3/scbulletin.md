After playback
==============
$ scbulletin -i ga2017qpupjg_reconstructed.xml -3

ERROR: ga2017qpupjg_reconstructed.xml: no origin and no event in eventparameters found


Before playback
===============
Event `ga2017qpupjg` was dumped using `scxmldump --debug -f -E ga2017qpupjg -P -A -M -F -d $DBFLAG > ga2017qpupjg.xml`.
 
$ scbulletin -i ga2017qpupjg.xml -3

Event:
    Public ID              ga2017qpupjg
    Description
      region name: Gascoyne, Gascoyne, Western Australia
Origin:
    Date                   2017-08-24
    Time                   05:30:11.2  +/-  1.8 s
    Latitude               -25.26 deg  +/-   11 km
    Longitude              116.76 deg  +/-   17 km
    Depth                     282 km   +/-   26 km
    Agency                 NA
    Mode                   automatic
    Status                 NOT SET
    Residual RMS              1.07 s
    Azimuthal gap             172 deg

5 Network magnitudes:
    mB        4.63            1
    Mw(mB)    3.84            1
    mb        3.81 +/- 0.59   4
    M         4.32            7 preferred
    MLv       4.57 +/- 0.46   7

11 Phase arrivals:
    sta  net   dist azi  phase   time         res     wt  sta
    MEEK  AU    2.2 130  P       05:30:58.0   0.1 A  1.0  MEEK
    MEEK  AU    2.2 130  S       05:31:47.5  12.2 AX 0.0  MEEK
    GIRL  AU    3.5 318  P       05:31:10.8  -0.1 A  1.0  GIRL
    MORW  AU    3.8 189  P       05:31:15.7   0.9 A  1.0  MORW
    BLDU  AU    5.3 180  P       05:31:29.6  -2.5 A  1.0  BLDU
    BLDU  AU    5.3 180  S       05:31:48.2 -49.1 AX 0.0  BLDU
    KLBR  AU    6.4 172  P       05:31:44.9   0.4 A  1.0  KLBR
    KLBR  AU    6.4 172  S       05:31:48.8 -71.3 AX 0.0  KLBR
    MUN   AU    6.7 184  P       05:31:49.6   0.9 A  1.0  MUN
    MUN   AU    6.7 184  S       05:32:24.0 -43.7 AX 0.0  MUN
    NWAO  AU    7.7 177  P       05:32:00.6   0.2 A  1.0  NWAO

12 Station magnitudes:
    sta  net   dist azi  type   value   res        amp per
    MEEK  AU    2.2 130  MLv     2.75 -1.82   0.112566
    GIRL  AU    3.5 318  MLv     4.05 -0.52   0.415849
    MORW  AU    3.8 189  MLv     5.05  0.47    3.06219
    BLDU  AU    5.3 180  MLv     4.72  0.14   0.603897
    BLDU  AU    5.3 180  mb      4.84  1.03    202.069 0.93
    KLBR  AU    6.4 172  MLv     4.87  0.30   0.477181
    KLBR  AU    6.4 172  mB      4.63  0.00    1148.89
    KLBR  AU    6.4 172  mb      3.75 -0.06      16.77 0.96
    MUN   AU    6.7 184  MLv     4.59  0.02   0.205489
    MUN   AU    6.7 184  mb      3.57 -0.23    14.1566 1.29
    NWAO  AU    7.7 177  MLv     4.81  0.24   0.197757
    NWAO  AU    7.7 177  mb      3.35 -0.46    14.4165 2.87

