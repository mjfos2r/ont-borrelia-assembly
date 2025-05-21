samtools view -H <file1.bam> | grep '^@RG'
samtools stats --split RG <file1.bam> 
docker images | awk 'NR != 1 {print "IMG:  "$1"\nID:   "$3"\nSize: "$7"\n"}'
