
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=0.004008

NJETS=$1

CHANNELS="OF"
#"All SF OF EE MuE EMu MuMu "

PROOFMODE="Sequential"

SAMESIGN="OS" 

MUONIDS="MediumIDTighterIP"
#"MediumIDTighterIP MediumID TightID TightIDTighterIP"

SAMPLES="
Data201550
WW50
HWW25
monoH_500GeV25
monoH_1GeV25
monoH_10GeV25
monoH_100GeV25
monoH_1000GeV25
"

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	for MUONID in $MUONIDS; do 

	    mkdir rootFiles/
	    mkdir rootFiles/${CHANNEL}	
	    mkdir rootFiles/${CHANNEL}/${MUONID}	
	    root -l -b -q "RunMuonAnalyzer.C(\"$SAMPLE\",\"$CHANNEL\",\"$SAMESIGN\",\"$PROOFMODE\",$LUMINOSITY,\"$MUONID\")"  

	done

    done

done


