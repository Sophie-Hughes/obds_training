#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:04:02 2020
Sorted 
@author: Sophie
"""

number_list = [26, 54, 93, 17, 77, 31, 44, 44, 20] 


def max_number(number_list):
    number_max = number_list[0]
    max_index = 0
    loop_index = 0 
    for number in number_list:
        if number > number_max:
            number_max = number
            max_index = loop_index
        #print(number, loop_index)
        loop_index += 1 
        
    return [number_max, max_index]

#print(max_number(number_list))
       
#sort of number_list 

for i in range(len(number_list), 1, -1):
    biggest_number, max_index = max_number(number_list[:i])
    print(number_list[:i])
    temp = number_list[i-1] 
    number_list[i-1] = biggest_number
    number_list[max_index] = temp 
    print(i, biggest_number, max_index)
    
print(number_list)
    
    

    