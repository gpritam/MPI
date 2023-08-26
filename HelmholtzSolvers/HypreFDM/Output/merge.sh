#!/bin/bash

cat $(printf "Solution-*") > $(printf "Solution.tec")
rm $(printf "Solution-*")
