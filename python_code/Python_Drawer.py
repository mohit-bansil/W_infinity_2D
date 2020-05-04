# -*- coding: utf-8 -*-
"""
Created on Sun May  3 23:40:15 2020

@author: Mohit
"""

import numpy as np
import cv2 as cv
import random as rand

mode = "all seperate"
#all seperate means that the program will display each cell one at a time
#combined means that the program will display all cells at once
W = 700

inputsmin = -6
inputsmax = 8


window = "Cell Diagram"
# Create black empty images
size = W, W, 3
image = np.zeros(size, dtype=np.uint8)

def image_coordinates(x):
    return int( (x - inputsmin)*W / (inputsmax - inputsmin) )

f = open("Cell_Data.txt", 'r')
N = int(f.readline())

for j in range(2**N - 1):
    count = int(f.readline())
    
    color = (rand.randint(0, 255), rand.randint(0, 255), rand.randint(0, 255))
    
    for i in range(count):
        x0, x1, y0, y1 = f.readline().split(" ")
        x0 = image_coordinates(float(x0))
        x1 = image_coordinates(float(x1))
        y0 = W - image_coordinates(float(y0))
        y1 = W - image_coordinates(float(y1))
        
        cv.rectangle(image, (x0,y0), (x1,y1), color, -1, 8)
        
    if(mode == "all seperate"):
        cv.imshow(window, image)
        cv.moveWindow(window, W, 200)
        cv.waitKey(0)
        cv.destroyAllWindows()

if(mode == "combined"):
    cv.imshow(window, image)
    
    cv.moveWindow(window, W, 200)
    cv.waitKey(0)
    cv.imwrite("Cell_Diagram.png", image)
    cv.destroyAllWindows()
