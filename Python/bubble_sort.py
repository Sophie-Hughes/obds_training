#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:41:43 2020

@author: Sophie
"""

number_list = [26, 54, 93, 17, 77, 31, 44, 55, 20] 

       
# bubble sort of number_list 
# Outer loop - stepping loop - j defines length of list, going down to 0 starts length of list and goes down to 1 by increments of 1
# Stepping loop 
for j in range(len(number_list), 1, -1): 
    for i in range(0, j-1, 1): # inner loop - swapping loop - inner loop pushes largest number to end  
        if number_list[i]>number_list[i+1]: # pairwise comparison 
            temp = number_list[i] # creating temporary variable 
            number_list[i] = number_list[i+1] # swapping varibles part 1: create duplicate
            number_list[i+1] = temp     
            print(number_list)
    #print(number_list[i], number_list[i+1])
    
# tip - instead of using i and j - give it a meaningful variable name - makes following code easier 
# better for debugging
# j = step iteration 
# i = current index of position in list 


   