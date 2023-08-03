#!/bin/bash
#author:siayouyang

grep -v "^@" $1 |cut -f9 |sort -n |uniq -c
