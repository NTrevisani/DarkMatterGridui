
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=19.365

TEST="test"

NJETS=$1

CHANNELS="OF"

PROOFMODE="Cluster"
#"All SF OF EE MuE EMu MuMu "


SAMPLES="
Dark1              \
Dark10             \
Dark100            \
Dark500            \
Dark1000           \
ZH                 \
QCD                \ 
Top                \
WJets              \
TTJets             \
HWW                \
WW                 \
VBF                \
"

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFiles/AllJet/
	mkdir rootFiles/AllJet/${CHANNEL}	
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\",\"$PROOFMODE\")"
	mv ${SAMPLE}.root rootFiles/AllJet/${CHANNEL}
  
    done

done


