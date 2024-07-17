#!/bin/bash

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <PLINDER_RELEASE> <PLINDER_ITERATION> <PLINDER_VERSION>"
  exit 1
fi

PLINDER_RELEASE="$1"
PLINDER_ITERATION="$2"
PLINDER_VERSION="$3"

mkdir -p releases/
rm -f releases/plinder.tar.gz
cp ~/Downloads/plinder-${PLINDER_VERSION}.tar.gz releases/plinder.tar.gz
# curl -L https://github.com/aivant/plinder/archive/refs/tags/v${PLINDER_VERSION}.tar.gz -o releases/plinder.tar.gz
pushd releases
tar -xzf plinder.tar.gz
rm -r plinder-${PLINDER_VERSION}/examples/
rm -r plinder-${PLINDER_VERSION}/analysis/
rm -r plinder-${PLINDER_VERSION}/dockerfiles/
rm -r plinder-${PLINDER_VERSION}/scripts/
rm -r plinder.tar.gz
tar -czvf plinder.tar.gz plinder-${PLINDER_VERSION}/
popd
gsutil cp releases/plinder.tar.gz gs://plinder/${PLINDER_RELEASE}/${PLINDER_ITERATION}/plinder.tar.gz
#
