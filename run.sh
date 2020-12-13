#!/bin/bash
touch analysis/trash && \
rm analysis/* && \
./parse_output.py fast_output.csv && \
./constraint_checks_v2.py fast_output.csv
