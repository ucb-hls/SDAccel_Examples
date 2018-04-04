#!/bin/bash 

#TARGET=15

for TARGET in `seq 1 14`;do
    aws s3 ls s3://adamacc/inputs
    aws s3 cp s3://adamacc/inputs/IR_INPUTS/ch$TARGET-ir-targets.tar.bz2 .
    aws s3 cp s3://adamacc/inputs/IR_INPUTS/schedules/${TARGET}valid.txt .
    aws s3 cp s3://adamacc/inputs/IR_INPUTS/target-dirs/${TARGET}.ir.targets.tsv .
    tar xvf ch$TARGET-ir-targets.tar.bz2
    #rm ${TARGET}.ir.targets.tsv
    #aws s3 cp s3://adamacc/inputs/IR_INPUTS/target-dirs/target-parse.pl .
    perl target-parse.pl ${TARGET}.ir.targets.tsv ./${TARGET} 
done
