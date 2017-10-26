#!/bin/bash

git rm .gitmodules
git rm --cached src/SeqLib
rm -rf .git/modules/src/SeqLib
rm -rf src/SeqLib
