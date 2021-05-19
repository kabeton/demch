#!/bin/bash

steps=(20 40 80 160 320 640 1280)
echo > log.txt

for i in "${steps[@]}"; do
  echo computing with $i
  $(./a.out $i >> log.txt)
  echo done
done
