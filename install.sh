#!/bin/sh

[ ! -d "bin" ] && mkdir bin
[ ! -d "tmp" ] && mkdir tmp

printf "IMPACTA version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (Alec made this)" >> version.txt

export PETSC_DIR=/opt/homebrew/Cellar/petsc/3.18.1
export BOOST_INCLUDE=/opt/homebrew/Cellar/boost/1.80.0/include
export BOOST_LIB=/opt/homebrew/Cellar/boost/1.80.0/lib

cd setup
make
