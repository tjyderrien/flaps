#!/bin/bash

tail -n1 TimeMax.dat | awk '{ print $10 }'
