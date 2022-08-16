#!/bin/sh

[ ! -d "bin" ] && mkdir bin
[ ! -d "tmp" ] && mkdir tmp

printf "IMPACTA version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (Alec made this)" >> version.txt

export PETSC_DIR=~/petsc
cd setup
make 
