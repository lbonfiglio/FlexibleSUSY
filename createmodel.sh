#!/bin/sh

# configure script for FlexibleSUSY
# Author: Alexander Voigt

if test $# -lt 1 ; then
    echo "Error: Too few arguments"
    echo "Usage: ./`basename $0` <model-name>"
    exit 1
fi

# directory of this script
BASEDIR=$(dirname $0)

# model name
model=$1
if test -d "models/$model" ; then
    echo "Error: directory models/$model already exists!"
    echo "Please chose another model name."
    exit 1
fi

# target directory
targetdir="models/$model"
echo -n "creating model directory $targetdir ... "
mkdir $targetdir
echo "done"

# target makefile module
MAKEFILE=module.mk
MAKEFILE_TMPL=config/module_template.mk

sed -e "s|@DIR@|$targetdir|g"         \
    -e "s|@MODEL@|$model|g"           \
    < $MAKEFILE_TMPL > $targetdir/$MAKEFILE
