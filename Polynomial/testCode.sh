#!/bin/bash

g++ Test.cpp -o TestedCode
g++ Testing/Gen.cpp -o Gen
g++ Testing/Checker.cpp -o Checker
ver="Accepted"
limit=10000

for((i = 0; i < limit; i++)); do
	./Gen > int
	./TestedCode < int > out2
	./Checker int out2 out1 ver
	ver=$(cat ver)
	if [ "$ver" != "Accepted" ]; then
		tc=$i
		break
	fi
done

if [ "$ver" == "Accepted" ]; then
	echo "Accepted. $limit test cases passed"
else
	echo "$ver on test $tc"
fi
