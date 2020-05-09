# -*- coding: utf-8 -*-
"""
Created on Sun May  3 23:40:15 2020

@author: Mohit Bansil

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import cv2 as cv
import random as rand

cellDataFileLocation = "Cell_Data.txt"

#all seperate means that the program will display each cell one at a time
#combined means that the program will display all cells at once
mode = "all seperate"

#This is the resolution of the image
W = 700

#These should be set to the smallest and largest possible 
#x,y coordinates in the cells
inputsmin = -6
inputsmax = 8

rand.seed(2)

#Rescale coordinates to correct scale for the diagram
def image_coordinates(x):
    return int( (x - inputsmin)*W / (inputsmax - inputsmin) )

def read_line(file):
    for line in file:
        if(not line.startswith("#")):
           return line
    print("Warning: End of File Reached")
    return
                             
def print_cell_diagram():
    
    # Create black empty images
    size = W, W, 3
    cellWindow = "Cell Diagram"
    cellImage = np.zeros(size, dtype=np.uint8)
    
    cellDataFile = open(cellDataFileLocation, 'r')
    N = int(read_line(cellDataFile))
    
    '''
    used just to visualize support of mu'''
    x0 = image_coordinates(float(0))
    x1 = image_coordinates(float(4))
    y0 = W - image_coordinates(float(0))
    y1 = W - image_coordinates(float(4))

    cv.rectangle(cellImage, (x0,y0), (x1,y1), (0,255,255), -1, 8)
    
    
    for j in range(2**N - 1):
        count = int(read_line(cellDataFile))
        
        color = (rand.randint(0, 255), rand.randint(0, 255), rand.randint(0, 255))
        
        for i in range(count):
            x0, x1, y0, y1 = read_line(cellDataFile).split(" ")
            x0 = image_coordinates(float(x0))
            x1 = image_coordinates(float(x1))
            y0 = W - image_coordinates(float(y0))
            y1 = W - image_coordinates(float(y1))
            
            cv.rectangle(cellImage, (x0,y0), (x1,y1), color, -1, 8)
            
        if(mode == "all seperate"):
            cv.imshow(cellWindow, cellImage)
            cv.moveWindow(cellWindow, W, 200)
            cv.waitKey(0)
            cv.destroyAllWindows()
    
    if(mode == "combined"):
        cv.imshow(cellWindow, cellImage)
        
        cv.moveWindow(cellWindow, W, 200)
        cv.waitKey(0)
        cv.imwrite("Cell_Diagram.png", cellImage)
        cv.destroyAllWindows()
    
    cellDataFile.close()


#start main program

print_cell_diagram()













