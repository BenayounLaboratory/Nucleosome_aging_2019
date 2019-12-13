#!/bin/bash

if [[ "$#" -lt 2 ]]
then
    echo "$(basename $0) [Ass] [BedDir] "  1>&2
    echo "   [Ass]: assembly" 1>&2 
    echo "   [BedDir]: folder with bam files" 1>&2

    exit 1
fi

Ass=$(echo $1 | sed 's:/$::g')
BedDir=$(echo $2 | sed 's:/$::g')

# make output directory if it doesnt exist
[[ ! -d "${BedDir}/GO_RAND" ]] && mkdir "${oDir}/GO_RAND"

##########################################
# 1. convert SAM into BAM, clean up alignemnts
for f in $(find "${BedDir}" -name '*.bed')
do
	filePath="${BedDir}/GO_RAND"
	
	# filter to retain mapped reads
	DirName=$(basename "${f}" | sed 's/\.bed//g');
    oDname="${filePath}/${DirName}"
    GenomeOntBas="${oDname}/basic.genomeOntology.txt"
    GenomeOntRep="${oDname}/repeats.genomeOntology.txt"
	
	# run genome ontology HOMER
    annotatePeaks.pl $f $Ass -genomeOntology $oDname > temp.xls
	
	fileName1=$(basename "${f}" | sed 's/\.bed/_basic_genomeOntology\.txt/g');
    oFname1="${filePath}/${fileName1}"
    
	fileName2=$(basename "${f}" | sed 's/\.bed/_repeat_genomeOntology\.txt/g');
    oFname2="${filePath}/${fileName2}"

	mv $GenomeOntBas $oFname1;
	mv $GenomeOntRep $oFname2;

	rm -r $oDname

    
done
echo "Finished Genome Ontology annotation"
