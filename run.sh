#!/bin/bash
VERSION="gaussianfrankenspline-test-new"
OUTPUT_ROOT="/home/krzysiek/output_store"
CHAINS=5

OUTPUT_DIR=$OUTPUT_ROOT/packages/gaussianfrankenspline/$VERSION
mkdir -p $OUTPUT_DIR
echo -ne "Output to: $OUTPUT_DIR\n\n\n"
FILE_INFO_FILE=$OUTPUT_DIR/$VERSION-file-data.txt
rm -f $FILE_INFO_FILE

for i in $(seq "$CHAINS")
do
	INIT_INPUT_FILE=circular-rw.inits
	OPTIM_OUTPUT_FILE=$OUTPUT_DIR/$VERSION-optimize-try-$i-output.csv
	OPTIM_DIAGNOSTIC_FILE=$OUTPUT_DIR/$VERSION-optimize-try-$i-diagnostics.csv
	OPTIM_TERM_FILE=$OUTPUT_DIR/$VERSION-optimize-try-$i-term.csv
	SAMPLER_OUTPUT_FILE=$OUTPUT_DIR/$VERSION-sampler-chain-$i-output.csv
	SAMPLER_DIAGNOSTIC_FILE=$OUTPUT_DIR/$VERSION-sampler-chain-$i-diagnostics.csv
	SAMPLER_TERM_FILE=$OUTPUT_DIR/$VERSION-sampler-chain-$i-term.csv
	
	echo -ne "output-dir:$OUTPUT_DIR\n" >> $FILE_INFO_FILE
	echo -ne "optim-try-$i-output-file:$OPTIM_OUTPUT_FILE\n" >> $FILE_INFO_FILE
	echo -ne "optim-try-$i-diagnostic-file:$OPTIM_DIAGNOSTIC_FILE\n" >> $FILE_INFO_FILE
	echo -ne "sampler-chain-$i-output-file:$SAMPLER_OUTPUT_FILE\n" >> $FILE_INFO_FILE
	echo -ne "sampler-chain-$i-diagnostic-file:$SAMPLER_DIAGNOSTIC_FILE\n" >> $FILE_INFO_FILE
	
	./gaussfrankenspline optimize id=$i \
		data file=circular-rw.data init=circular-rw.inits \
		output \
			file=$OPTIM_OUTPUT_FILE diagnostic_file=$OPTIM_DIAGNOSTIC_FILE   \
			refresh=10 > $OPTIM_TERM_FILE &
	
	./gaussfrankenspline sample id=$i \
		data file=circular-rw.data init=circular-rw.inits \
		output \
			file=$SAMPLER_OUTPUT_FILE diagnostic_file=$SAMPLER_DIAGNOSTIC_FILE   \
			refresh=10 > $SAMPLER_TERM_FILE &
done

