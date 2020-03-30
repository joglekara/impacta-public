#!/bin/sh

mkdir bin
mkdir tmp

printf "IMPACTA version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (Archis made this)" >> version.txt

cd setup
make debug