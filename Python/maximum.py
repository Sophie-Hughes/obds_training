#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:32:31 2020

Find maximum number from list of numbers 

@author: Sophie
"""


number_list = [26, 54, 93, 17, 77, 31, 44, 44, 20] 
number_max = number_list[0] 

for number in number_list:
    if number > number_max:
        number_max = number
       
print('The maximum number is', number_max)


number_max.selection_sort 