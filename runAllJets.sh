
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
#"All SF OF EE MuE EMu MuMu "

PROOFMODE="Cluster"

SAMPLES="
Dark100            \
Dark10             \
Dark1              \
ZH                 \
HWW                \
WW                 \
"
#TTJets             \
#DY                 \
#"
#Dark100            \
#Dark500            \
#Dark1000           \
#Dark1              \
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
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\",\"$PROOFMODE\")"
	mv ${SAMPLE}.root rootFiles/AllJet/${CHANNEL}/${MULTIPLICATEXS}
  
    done

done


