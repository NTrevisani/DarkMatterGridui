
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=100

TEST="test"

NJETS=0

CHANNELS="OF"

PROOFMODE="Cluster"

MULTIPLICATEXS=$1

#"All SF OF EE MuE EMu MuMu "

SAMPLES="
Dark1              \
Dark10             \
Dark100            \
Dark500            \
Dark1000           \
ZH                 \
HWW                \
WW                 \
"

#QCD                \ 
#Top                \
#WJets              \
#TTJets             \
#VBF                \
#"

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFiles/AllJet/
	mkdir rootFiles/AllJet/${CHANNEL}	
	mkdir rootFiles/AllJet/${CHANNEL}/${MULTIPLICATEXS}	
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\",\"$PROOFMODE\","$MULTIPLICATEXS")"
	mv ${SAMPLE}.root rootFiles/AllJet/${CHANNEL}/${MULTIPLICATEXS}
  
    done

done


