#!/bin/bash

rm -r build/doctrees
rm -r build/html
rm source/plinder.*

cp ../examples/*.ipynb source/
sphinx-apidoc -o source -d 10 -f ../src/plinder/
cp plinder.rst source/
make html

if [[ -n "$1" ]]; then
  open build/html/index.html
else
  echo "View docs at build/html/index.html"
fi
rm -f source/*.ipynb
