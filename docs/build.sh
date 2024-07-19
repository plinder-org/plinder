#!/bin/bash

rm source/plinder.*
cp ../examples/*.ipynb .
sphinx-apidoc -o source -d 10 -f ../src/plinder/
cp plinder.rst source/
rm -r build/doctrees
rm -r build/html
make html

if [[ -n "$1" ]]; then
  open build/html/index.html
else
  echo "View docs at build/html/index.html"
fi
rm *.ipynb
