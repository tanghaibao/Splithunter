#!/bin/bash

git rm .gitmodules
git rm --cached SeqLib
rm -rf .git/modules/SeqLib
rm -rf SeqLib
