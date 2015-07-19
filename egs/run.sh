#!/bin/bash
# Copyright 2015 SpeechLab at FEE CTU Prague (Author: Petr Mizera)
# Apache 2.0
#
# Usage:
#   bash example.sh 1
#------------------------------------------
#echo "copy ../Default/ctucopy4 ."
[ ! -f ../src/bin/ctucopy4 ] && echo "The ctucopy4 tool was not compiled. See README" && exit 1
[ ! -f ctucopy4 ] && cp ../src/bin/ctucopy4 .

case "$1" in
    0)
      echo "00) clean all - data/*";
      rm -r data/
      rm ctucopy4
      echo "========================="
      ;;
    1)
      echo "01) create param. MFCC ";
      rm -r data/01
      mkdir -p data/01
      ./ctucopy4 -C conf/01_mfcc_0_16k_2510_30.ctuconf  -S conf/01_list.scp
      echo "========================="
      echo "It was created:"
      tree data/01
      echo "========================="
      ;;
    2)
      echo "02) only create stat.cvmn from raw/wav signal => data/02/cmvn.stat ";
      echo "  -example of list: ";
      echo "         path/to/signal speaker_ID";
      echo "         --------------------------";
      echo "         /data/SPEECON/ADULT1CS/BLOCK00/SES000/SA000CB1.CS0 SA000";
      echo "   or: ";
      echo "         /data/SPEECON/ADULT1CS/BLOCK00/SES000/SA000CB1.CS0 data/02/SA000CB1.CS0.mfcc_0_16k_2510_30 SA000";
      rm -r data/02
      mkdir -p data/02
      ./ctucopy4 -C conf/02_mfcc_0_16k_2510_30.ctuconf  -S conf/02_list.scp
      echo "========================="
      echo "It was created:"
      tree data/02
      echo "========================="
    ;;
    3)
      echo "03) only create stat.cvmn from features (HTK) => data/03/cmvn.stat ";
      echo "  -example of list: ";
      echo "         path/to/features speaker_ID";
      echo "         --------------------------";
      echo "         data/01/SA000CB1.CS0.mfcc_0_16k_2510_30 SA000";
      echo "   or: ";
      echo "         data/01/SA000CB1.CS0.mfcc_0_16k_2510_30 data/03/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/03
      mkdir -p data/03
      ./ctucopy4 -C conf/03_mfcc_0_16k_2510_30.ctuconf  -S conf/03_list.scp
      echo "========================="
      echo "It was created:"
      tree data/03
      echo "========================="
    ;;
    4)
      echo "04) Apply cmvn, format_in is raw/wav. Stat. cmvn open from => data/02/cmvn.stat ";
      echo "  -example of list: ";
      echo "         path/to/signal path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         /data/SPEECON/ADULT1CS/BLOCK00/SES000/SA000CB1.CS0 data/04/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/04
      mkdir -p data/04
      ./ctucopy4 -C conf/04_mfcc_0_16k_2510_30.ctuconf  -S conf/04_list.scp
      echo "========================="
      echo "It was created:"
      tree data/04
      echo "========================="
    ;;
    5)
      echo "05) Apply cmvn, format_in is htk. Stat. cmvn open from => data/03/cmvn.stat ";
      echo "  -example of list: ";
      echo "         path/to/features path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         data/01/SA000CB1.CS0.mfcc_0_16k_2510_30 data/05/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/05
      mkdir -p data/05
      ./ctucopy4 -C conf/05_mfcc_0_16k_2510_30.ctuconf  -S conf/05_list.scp
      echo "========================="
      echo "It was created:"
      tree data/05
      echo "========================="
    ;;
    6)
      echo "06) Apply cmvn, format_in is raw/wav. Stat. cmvn is not found.";
      echo "  -example of list: ";
      echo "         path/to/signal path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         /data/SPEECON/ADULT1CS/BLOCK00/SES000/SA000CB1.CS0 data/06/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/06
      mkdir -p data/06
      ./ctucopy4 -C conf/06_mfcc_0_16k_2510_30.ctuconf  -S conf/06_list.scp
      echo "========================="
      echo "It was created:"
      tree data/06
      echo "========================="
    ;;
    7)
      echo "07) Apply cmvn, format_in is htk (features). Stat. cmvn is not found.";
      echo "  -example of list: ";
      echo "         path/to/features path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         data/01/SA000CB1.CS0.mfcc_0_16k_2510_30 data/07/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/07
      mkdir -p data/07
      ./ctucopy4 -C conf/07_mfcc_0_16k_2510_30.ctuconf  -S conf/07_list.scp
      #kdyz neni uveden soubor se stat. cmvn, tak se ozve : OPTS: Syntax error in option "-apply_cmvn".
      echo "========================="
      echo "It was created:"
      tree data/07
      echo "========================="
    ;;
    8)
      echo "08) Apply cmvn, format_in is raw/wav. and add operator -stat_cmvn .";
      echo "  -example of list: ";
      echo "         path/to/signal path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         /data/SPEECON/ADULT1CS/BLOCK00/SES000/SA000CB1.CS0 data/08/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/08
      mkdir -p data/08
      ./ctucopy4 -C conf/08_mfcc_0_16k_2510_30.ctuconf  -S conf/08_list.scp
      echo "========================="
      echo "It was created:"
      tree data/08
      echo "========================="
    ;;
    9)
      echo "09) Apply cmvn, format_in is htk. and add operator -stat_cmvn .";
      echo "  -example of list: ";
      echo "         path/to/features path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         path/to/features data/09/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/09
      mkdir -p data/09
      ./ctucopy4 -C conf/09_mfcc_0_16k_2510_30.ctuconf  -S conf/09_list.scp
      echo "========================="
      echo "It was created:"
      tree data/09
      echo "========================="
    ;;
    10)
      echo "10) Apply cmvn, format_in is htk, format_out=ark . and add operator -stat_cmvn .";
      echo "  -example of list: ";
      echo "         path/to/features path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         path/to/features data/10/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/10
      mkdir -p data/10
      ./ctucopy4 -C conf/10_mfcc_0_16k_2510_30.ctuconf  -S conf/10_list.scp
      echo "========================="
      echo "It was created:"
      tree data/10
      echo "========================="
    ;;
    11)
      echo "11) create param. MFCC and save to ark (kaldi format) ";
      rm -r data/11
      mkdir -p data/11
      ./ctucopy4 -C conf/11_mfcc_0_16k_2510_30.ctuconf  -S conf/11_list.scp
      echo "========================="
      echo "It was created:"
      tree data/11
      echo "========================="
      ;;
    12)
      echo "12) Apply cmvn, format_in is raw/wav. and add operator -stat_cmvn. large list of speakers";
      echo "  -example of list: ";
      echo "         path/to/features path/to/norm_features speaker_ID";
      echo "         --------------------------";
      echo "         path/to/features data/12/SA000CB1.CS0.mfcc_0_16k_2510_30.norm SA000";
      rm -r data/12
      mkdir -p data/12/MFCC
      ./ctucopy4 -C conf/12_mfcc_0_16k_2510_30.ctuconf  -S conf/12_list.scp
      echo "========================="
      echo "It was created:"
      tree data/12
      echo "========================="
    ;;
    13)
      echo "13) delta MFCC ";
      rm -r data/15
      mkdir -p data/15
      ./ctucopy4 -C conf/15_mfcc_0_16k_2510_30.ctuconf  -S conf/15_list.scp
      echo "========================="
      echo "It was created:"
      tree data/15
      echo "========================="
      ;;
    14)
      echo "14) delta-delta MFCC to ark + CMS exp.";
      rm -r data/16
      mkdir -p data/16
      ./ctucopy4 -C conf/16_mfcc_0_16k_2510_30.ctuconf  -S conf/16_list.scp
      echo "========================="
      echo "It was created:"
      tree data/16
      echo "========================="
      ;;
    15)
      echo "15) Create context information from features defined -fea_kind. save to ark";
      rm -r data/17
      mkdir -p data/17
      ./ctucopy4 -C conf/17_mfcc_0_16k_2510_30.ctuconf  -S conf/17_list.scp
      echo "========================="
      echo "It was created:"
      tree data/17
      echo "========================="
      ;;
    16)
      echo "16) Create context information from input htk features. save to ark";
      rm -r data/18
      mkdir -p data/18
      ./ctucopy4 -C conf/18_mfcc_0_16k_2510_30.ctuconf  -S conf/18_list.scp
      echo "========================="
      echo "It was created:"
      tree data/18
      echo "========================="
      ;;
    17)
      echo "17) Create ark from input htk features.";
      rm -r data/19
      mkdir -p data/19
      ./ctucopy4 -C conf/19_mfcc_0_16k_2510_30.ctuconf  -S conf/19_list.scp
      echo "========================="
      echo "It was created:"
      tree data/19
      echo "========================="
      ;;
    18)
      echo "18) Make td-iir-mfcc based features from raw input to HTK format.";
      rm -r data/20
      mkdir -p data/20
      ./ctucopy4 -C conf/20_td-iir-mfcc.ctuconf  -S conf/20_list.scp
      echo "========================="
      echo "It was created:"
      tree data/20
      echo "========================="
      ;;
    19)
      echo "19) SS exten - raw input";
      rm -r data/21
      mkdir -p data/21
      ./ctucopy4 -C conf/21_exten.ctuconf -S conf/21_list.scp
      echo "========================="
      echo "It was created:"
      tree data/21
      echo "========================="
      ;;
esac
