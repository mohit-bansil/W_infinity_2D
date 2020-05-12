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
import os

cellDataFileLocation = "Cell_Data.txt"
cellDataWithMuFileLocation = "Cell_Data_With_Mu.txt"
transportPlanDataFileLocation = "Optimal_Transport_Plan.txt"

#This is the folder where the images should be stored

#This is the resolution of the image
W = 700

#These should be set to the smallest and largest possible 
#x,y coordinates in the cells
inputsmin = -.5
inputsmax = 4.5

tol = 0.00001
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
                             
def read_rectangle(file):
    x0, x1, y0, y1 = read_line(file).split(" ")
    x0 = image_coordinates(float(x0))
    x1 = image_coordinates(float(x1))
    y0 = W - image_coordinates(float(y0))
    y1 = W - image_coordinates(float(y1)) 
    return x0, x1, y0, y1

def show_image(window, image):
    cv.imshow(window, image)
    cv.moveWindow(window, W, 200)
    cv.waitKey(0)
    cv.destroyAllWindows()
  

#for mode, all seperate means that the program will display each cell one at a time
#combined means that the program will display all cells at once
#showImages tells whether or not the images should be displayed on the screen
def print_cell_diagram(mode, showImages):
    
    if mode == "all seperate" and not os.path.exists('Cell_Diagrams'):
        os.makedirs('Cell_Diagrams')
    
    # Create black empty images
    size = W, W, 3
    cellWindow = "Cell Diagram"
    cellImage = np.full(size, 255, dtype=np.uint8)
    
    cellDataFile = open(cellDataFileLocation, 'r')
    N = int(read_line(cellDataFile))
    
    '''
    used just to visualize support of mu
    x0 = image_coordinates(float(0))
    x1 = image_coordinates(float(4))
    y0 = W - image_coordinates(float(0))
    y1 = W - image_coordinates(float(4)) 
    cv.rectangle(cellImage, (x0,y0), (x1,y1), (0,255,255), -1, 8)
    '''
    
    for j in range(1, 2**N):
        count = int(read_line(cellDataFile))
        
        color = (rand.randint(0, 255), rand.randint(0, 255), rand.randint(0, 255))
        
        for i in range(count):
            x0, x1, y0, y1 = read_rectangle(cellDataFile)
            cv.rectangle(cellImage, (x0,y0), (x1,y1), color, -1, 8)
            
        if(mode == "all seperate"):
            if(showImages):
                show_image(cellWindow, cellImage)
            cv.imwrite("Cell_Diagrams/" + "Cell_Diagram"+ str(j) +".png", cellImage)
            cellImage = np.full(size, 255, dtype=np.uint8)
    
    if(mode == "combined"):
        if(showImages):
            show_image(cellWindow, cellImage)     
        cv.imwrite("Cell_Diagram.png", cellImage)
    
    cellDataFile.close()

def print_transport_diagram():
    
    if not os.path.exists('Warehouse_Diagrams'):
        os.makedirs('Warehouse_Diagrams')

    cellDataWithMuFile = open(cellDataWithMuFileLocation, 'r')
    transportPlanDataFile = open(transportPlanDataFileLocation, 'r')
    N = int(read_line(cellDataWithMuFile))
    read_line(transportPlanDataFile) #kills the N in that file
    read_line(transportPlanDataFile) #kills empty cell in that file
    size = W, W, 3
    transportWindow = "Transport Diagram"
    warehouseImages = [np.full(size, 255, dtype=np.uint8) for i in range(N)]
    
    for j in range(2**N - 1):
        size = float(read_line(cellDataWithMuFile))
        count = int(read_line(cellDataWithMuFile))
        
        transportAmounts = read_line(transportPlanDataFile).split(" ")
        
        cell_rectangles = []
        
        for i in range(count):
            cell_rectangles.append((read_rectangle(cellDataWithMuFile)))
        
        for k in range(N):
            if(float(transportAmounts[k]) < tol):
                shade = 255
            else:
                shade = 255 - 255* float(transportAmounts[k]) / size
        
            for i in range(count):
                cv.rectangle(warehouseImages[k], (cell_rectangles[i][0],cell_rectangles[i][2]), (cell_rectangles[i][1],cell_rectangles[i][3]), (shade, shade, shade), -1, 8)
    
    for k in range(N):
        show_image(transportWindow, warehouseImages[k])     
        cv.imwrite("Warehouse_Diagrams/" + "Warehouse_Diagram" + str(k + 1) + ".png", warehouseImages[k])
    
    cellDataWithMuFile.close()
    transportPlanDataFile.close()

    

#start main program
print_transport_diagram()
print_cell_diagram("all seperate", False)













