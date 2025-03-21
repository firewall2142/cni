#!/bin/bash

for (( c=1; c<=5; c++ ))
do

echo "###################################"
echo "########## TEST - $c ###############"
echo "###################################"
./fast_opt.py test_files/lab_test_data_00$c.csv test_files/district_test_data_00$c.csv test_files/solution_00$c.csv


echo ""
echo ""
echo "Checking constraints for $c"
echo "--------------------------"

./constraint_checks_v2.py test_files/lab_test_data_00$c.csv test_files/district_test_data_00$c.csv test_files/solution_00$c.csv

echo ""
echo ""
echo ""

done



echo "###############################"
echo "########## DONE ###############"
echo "###############################"
