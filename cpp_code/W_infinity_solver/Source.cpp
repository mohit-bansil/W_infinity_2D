/*
W_infinity_solver

Author: Mohit Bansil

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
*/

#include<iostream>
#include<fstream>
#include <algorithm>
#include <bitset> 

using namespace std;

struct rectangle;
rectangle intersect(rectangle A, rectangle B);

const int maxN = 8;
const double tolerance = 0.0000001;
//const double infinity = 9999999;

int N;
double y[maxN][2];
double lambda[maxN];
rectangle* pointerMu;


ifstream inputFile("Input_Data.txt");
ofstream cellOutputFile("Cell_Data.txt");
ofstream cellOutputFileWithMu("Cell_Data_With_Mu.txt");
ofstream optimalTransportPlanFile("Optimal_Transport_Plan.txt");

/*
Rectangle object is charaterized by its left and right x coordinates and its upper and lower y coordinates
Currently methods are
1. Determine if rectangle is non-empty (non-degenerate)
2. Find area
3. Find mu area (area of rectangle intersected with mu)
4. Determine if a point is inside the rectangle
5. Print the coordinates of the rectangle to cellOutputFile
6. Print the coordinates of the rectangle to cellOutputFileWithMu
7. Intersect two rectangles*/

struct rectangle
{
	//data
	double x0, x1, y0, y1;

	//methods
	bool is_nondegenerate()
	{
		if (x1 < x0 || y1 < y0)
			return false;
		else
			return true;
	}
	double area()
	{
		if (is_nondegenerate())
			return (x1 - x0)*(y1 - y0);
		else
			return 0;
	}
	double mu_area()
	{
		return intersect(*this, *pointerMu).area() / ((*pointerMu).area());
	}
	bool is_inside(double x, double y_point)
	{
		if (x < x0 || x > x1 || y_point < y0 || y_point > y1)
			return false;
		else
			return true;
	}
	void print()
	{
		cellOutputFile << x0 << " " << x1 << " " << y0 << " " << y1 << endl;
	}
	void print_mu()
	{
		cellOutputFileWithMu << x0 << " " << x1 << " " << y0 << " " << y1 << endl;
	}

	//constructors
	rectangle() {}
	rectangle(double x0_, double x1_, double y0_, double y1_)
	{
		x0 = x0_;
		x1 = x1_;
		y0 = y0_;
		y1 = y1_;
	}
};
rectangle intersect(rectangle A, rectangle B)
{
	return rectangle(max(A.x0, B.x0), min(A.x1, B.x1), max(A.y0, B.y0), min(A.y1, B.y1));
}

/*This function returns the rectangle given by the points of distance 
 at most omega from the warehouse y_i */ 
rectangle generate_warehouse_rectangle(const int i, const double omega)
{
	return rectangle(y[i][0] - omega, y[i][0] + omega, y[i][1] - omega, y[i][1] + omega);
}

/*
A cell is modeled as an array of (disjoint) rectangles. There is also a count variable that tells us how many rectangles there are.

Currently the methods are:
1. intersect the cell with a given rectangle.
2. intersect the cell with the complement of a given rectangle (Outersect).
3. Compute the mu area of the cell. 
4. Print the coordinates of each rectangle in the cell to cellOutputFile
5. Print the coordinates of each rectangle in the cell and the mu area to cellOutputFileWithMu
6. Delete all degenerate rectangles from cell.

WARNING!!! The way Outersect is written, it assumes that the outersecting rectangle 
is at least as big as each of the rectangles in the cell! 
This is true for our purposes but may not be true in general. Be careful!

A cell can also be contrustructed from a rectangle. 
*/
struct cell
{
	//data
	int count;
	//Note that an application of outersect_with_rectangle can increase the number of pieces by at most 8
	rectangle pieces[10 * maxN + 4];
	

	//methods
	void intersect_with_rectangle(rectangle A)
	{

		for (int i = 0; i < count; i++)
		{
			pieces[i] = intersect(pieces[i], A);
		}

		return;
	}

	void outersect_with_rectangle(rectangle A)
	{
		int warningFlag;

		int newPiecesAdded = 0;

		for (int i = 0; i < count; i++)
		{
			warningFlag = 0;

			if (intersect(pieces[i], A).is_nondegenerate() != true)
				continue;

			if (pieces[i].is_inside(A.x0, A.y0))
			{
				pieces[count + newPiecesAdded] = rectangle(pieces[i].x0, A.x0, A.y0, pieces[i].y1);
				newPiecesAdded++;
				pieces[count + newPiecesAdded] = rectangle(A.x0, pieces[i].x1, pieces[i].y0, A.y0);
				newPiecesAdded++;

				pieces[i].x1 = A.x0;
				pieces[i].y1 = A.y0;

				continue;
			}
			if (pieces[i].is_inside(A.x1, A.y0))
			{
				pieces[count + newPiecesAdded] = rectangle(A.x1, pieces[i].x1, A.y0, pieces[i].y1);
				newPiecesAdded++;
				pieces[count + newPiecesAdded] = rectangle(pieces[i].x0, A.x1, pieces[i].y0, A.y0);
				newPiecesAdded++;

				pieces[i].x0 = A.x1;
				pieces[i].y1 = A.y0;

				continue;
			}
			if (pieces[i].is_inside(A.x0, A.y1))
			{
				pieces[count + newPiecesAdded] = rectangle(pieces[i].x0, A.x0, pieces[i].y0, A.y1);
				newPiecesAdded++;
				pieces[count + newPiecesAdded] = rectangle(A.x0, pieces[i].x1, A.y1, pieces[i].y1);
				newPiecesAdded++;

				pieces[i].x1 = A.x0;
				pieces[i].y0 = A.y1;

				continue;
			}
			if (pieces[i].is_inside(A.x1, A.y1))
			{
				pieces[count + newPiecesAdded] = rectangle(pieces[i].x0, A.x1, A.y1, pieces[i].y1);
				newPiecesAdded++;
				pieces[count + newPiecesAdded] = rectangle(A.x1, pieces[i].x1, pieces[i].y0, A.y1);
				newPiecesAdded++;

				pieces[i].x0 = A.x1;
				pieces[i].y0 = A.y1;

				continue;
			}

			if(A.y0 < pieces[i].y0 && A.y1 > pieces[i].y1)
			{
				if (pieces[i].x1 > A.x0 && A.x1 > pieces[i].x1)
				{
					pieces[i].x1 = A.x0;
					warningFlag++;
				}
				if (pieces[i].x0 < A.x1 && A.x0 < pieces[i].x0)
				{
					pieces[i].x0 = A.x1;
					warningFlag++;
				}
			}

			if (A.x0 < pieces[i].x0 && A.x1 > pieces[i].x1)
			{
				if (pieces[i].y1 > A.y0 && A.y1 > pieces[i].y1)
				{
					pieces[i].y1 = A.y0;
					warningFlag++;
				}
				if (pieces[i].y0 < A.y1 && A.y0 < pieces[i].y0)
				{
					pieces[i].y0 = A.y1;
					warningFlag++;
				}
			}

/*
			if (warningFlag > 1)
				cerr << "Warning something might be wrong. Multiple sides were cut from a cell piece." << endl;
			if (warningFlag == 0)
				cerr << "Warning something might be wrong. No sides were cut from a cell piece." << endl;
*/
		}

		count += newPiecesAdded;

		return;
	}
	double mu_area()
	{
		double answer = 0;
		for (int i = 0; i < count; i++)
		{
			answer += pieces[i].mu_area();
		}

		return answer;
	}
	void print()
	{
		cellOutputFile << count << endl;
		for (int i = 0; i < count; i++)
			pieces[i].print();

	}
	void print_mu()
	{
		cellOutputFileWithMu << (*this).mu_area() << endl;
		cellOutputFileWithMu << count << endl;
		for (int i = 0; i < count; i++)
			pieces[i].print_mu();
	}
	void remove_degenerate_rectangles()
	{
		for (int i = 0; i < count; i++)
			if (pieces[i].is_nondegenerate() == false)
			{
				count--;
				for (int j = i; j < count; j++)
					pieces[j] = pieces[j + 1];
				i--;
			}
	}


	//constructors
	cell() {
		count = 0;
	}
	cell(rectangle A_) {
		count = 1;
		rectangle temp = A_;
		pieces[0] = temp;
	}
};

/*This function takes in an omega stores the left verticies of the omega-transport graph in the array cellsizes.
It also prints the cells to the Output File */
void compute_cell_graph(const double omega, double cellsizes[], bool printCells, bool printCellsMu)
{
	rectangle warehouseRectangles[maxN];

	for (int i = 0; i < N; i++)
	{
		warehouseRectangles[i] = generate_warehouse_rectangle(i, omega);
	}

	cell currentCell;
	double currentSumOfCellMuSizes = 0;

	for (int i = 1; i < (1 << N); i++)
	{
		bitset<32>b(i);

		int firstRectangle;

		for (int j = 0; j < N; j++)
		{
			if (b[j] == 1)
			{
				firstRectangle = j;
				break;
			}
		}

		currentCell = cell(warehouseRectangles[firstRectangle]);

		for (int j = 0; j < N; j++)
		{
			if (j == firstRectangle)
				continue;
			else if (b[j] == 1)
				currentCell.intersect_with_rectangle(warehouseRectangles[j]);
			else
				currentCell.outersect_with_rectangle(warehouseRectangles[j]);
		}

		cellsizes[i] = currentCell.mu_area();

		currentSumOfCellMuSizes += cellsizes[i];

		if (printCells)
		{
			currentCell.remove_degenerate_rectangles();
			currentCell.print();
		}

		if (printCellsMu)
		{
			currentCell.intersect_with_rectangle(*pointerMu);
			currentCell.remove_degenerate_rectangles();
			currentCell.print_mu();
		}
	}

	cellsizes[0] = 1 - currentSumOfCellMuSizes;

	return;
}

/* This function takes a transport graph and determines if there is a perfect matching*/
bool is_there_perfect_matching(double cellsizes[])
{
	//while this isn't needed it should make the code much faster
	if (cellsizes[0] > tolerance)
		return false;

	double leftsum, rightsum;

	for (int i = 1; i < (1 << N); i++)
	{
		leftsum = 0;
		rightsum = 0;
		bitset<32>currentSubsetBeingChecked(i);
		
		//Compute sum of weights of the right verticies in the current subset being checked
		for (int j = 0; j < N; j++)
		{
			if (currentSubsetBeingChecked[j] == 1)
				rightsum += lambda[j];
		}

		//Compute sum of weights of neighbor (left) verticies in the current subset being checked
		for (int j = 1; j < (1 << N); j++)
		{
			bitset<32>currentCellBeingChecked(j);

			for (int k = 0; k < N; k++)
			{
				if (currentCellBeingChecked[k] == 1 && currentSubsetBeingChecked[k] == 1)
				{
					leftsum += cellsizes[j];
					break;
				}
			}
		}

		if (leftsum < rightsum - tolerance)
			return false;

	}

	return true;
}

/* This function reads the values of N, lambda, y from the input file*/
void input_data()
{
	inputFile >> N;

	for (int i = 0; i < N; i++)
		inputFile >> lambda[i];

	for (int i = 0; i < N; i++)
		inputFile >> y[i][0] >> y[i][1];
}

/* This function solves the W_infinity transport problem to within an error of desiredError.
The array cellSizes is written with the resulting transport graph vertex weights
The resulting cell decomposition is written to the cellOutputFile and cellOutputFileWithMu*/
void solve(double lowerOmegaBound, double upperOmegaBound, double desiredError, double cellSizes[])
{
	double currentOmegaGuess;

	while (upperOmegaBound - lowerOmegaBound > desiredError)
	{
		currentOmegaGuess = (lowerOmegaBound + upperOmegaBound) / 2;
		compute_cell_graph(currentOmegaGuess, cellSizes, false, false);
		if (is_there_perfect_matching(cellSizes))
			upperOmegaBound = currentOmegaGuess;
		else
			lowerOmegaBound = currentOmegaGuess;
	}

	compute_cell_graph(upperOmegaBound, cellSizes, true, true);

//	cellOutputFile << endl << lowerOmegaBound << " " << upperOmegaBound << endl;

	return;
}

// bool for left or right and a number
struct vertex {
	int num;
	bool isLeft;

	vertex() {	}
	vertex(int num_, bool isLeft_)
	{
		num = num_;
		isLeft = isLeft_;
	}
};

/* This function finds the shortest path between the source and sink in the residual flow graph.
The return value is the vertex number of the last right vertex of the shortest path. -1 is returned if no path is found.
The path data found is stored in leftVerticiesPreviousVertex and rightVerticiesPreviousVertex, 
For example if k is the number of a left vertex in the shortest path leftVerticiesPreviousVertex[k] is the number of the right vertex that preceeds it
-2 is used for the source vertex
*/
int shortestPath(const double leftVertexWeights[], const double rightVertexWeights[], const double partialMatching[][maxN], const double partialLeftSums[], const double partialRightSums[], int leftVerticiesPreviousVertex[], int rightVerticiesPreviousVertex[])
{
	vertex queue[1 << maxN];
	int queueCurrentPosition = 0;
	int queueSize = 0;

	//Initialize all previous verticies to -1
	for (int i = 0; i < (1 << N); i++)
		leftVerticiesPreviousVertex[i] = -1;
	for (int i = 0; i < N; i++)
		rightVerticiesPreviousVertex[i] = -1;

	//Set all left verticies connected to the source to have previous vertex -2 (for source) and put them in the queue
	for (int i = 1; i < (1 << N); i++)
	{
		if (partialLeftSums[i] < leftVertexWeights[i] - tolerance)
		{
			leftVerticiesPreviousVertex[i] = -2;
			queue[queueSize] = vertex(i, true);
			queueSize++;
		}
	}

	for (; queueCurrentPosition < queueSize; queueCurrentPosition++)
	{
		int currentVertexNum = queue[queueCurrentPosition].num;

		//current vertex in queue is a left vertex
		if (queue[queueCurrentPosition].isLeft)
		{
			bitset<32>currentLeftVertex(currentVertexNum);
			for (int j = 0; j < N; j++)
			{
				if (currentLeftVertex[j] == 1 && rightVerticiesPreviousVertex[j] == -1)
				{
					rightVerticiesPreviousVertex[j] = currentVertexNum;
					queue[queueSize] = vertex(j, false);
					queueSize++;

					//Check if we can get to the sink!!!
					if (partialRightSums[j] < rightVertexWeights[j] - tolerance)
					{
						return j;
					}


				}
			}
		}
		//current vertex in queue is a right vertex
		else
		{
			for (int j = 1; j < (1 << N); j++)
			{
				bitset<32>leftVertex(j);

				if (partialMatching[j][currentVertexNum] > tolerance && leftVertex[currentVertexNum] == 1 && leftVerticiesPreviousVertex[j] == -1)
				{
					leftVerticiesPreviousVertex[j] = currentVertexNum;
					queue[queueSize] = vertex(j, true);
					queueSize++;
				}
			}
		}
	}

	return -1;
}

/* This function finds the a maximal matching in the transport graph given by leftVertexWeights and rightVertexWeights
It is stored in partialMatching
*/
void edmondKarp(const double leftVertexWeights[], const double rightVertexWeights[], double partialMatching[][maxN])
{
	double partialLeftSums[1 << maxN];
	double partialRightSums[maxN];

	int shortestPathLeftVerticiesPreviousVertex[1 << maxN];
	int shortestPathRightVerticiesPreviousVertex[maxN];

	//Initialize partialMatching and its left and right sums
	for (int i = 0; i < (1 << N); i++)
	{
		for (int j = 0; j < N; j++)
			partialMatching[i][j] = 0;
	}
	for (int i = 0; i < (1 << N); i++)
		partialLeftSums[i] = 0;
	for (int j = 0; j < N; j++)
		partialRightSums[j] = 0;


	int shortestPathLastVextex = shortestPath(leftVertexWeights, rightVertexWeights, partialMatching, partialLeftSums, partialRightSums, shortestPathLeftVerticiesPreviousVertex, shortestPathRightVerticiesPreviousVertex);
	
	double capacityMatched;
	int currentRightVertexNum;
	int currentLeftVertexNum;
	while (shortestPathLastVextex != -1)
	{
		//Update partial matching

		//First we find the capacity of our path
		capacityMatched = rightVertexWeights[shortestPathLastVextex] - partialRightSums[shortestPathLastVextex];
		currentRightVertexNum = shortestPathLastVextex;
		while (true)
		{
			currentLeftVertexNum = shortestPathRightVerticiesPreviousVertex[currentRightVertexNum];
			if (shortestPathLeftVerticiesPreviousVertex[currentLeftVertexNum] == -2)
				break;
			capacityMatched = min(capacityMatched, partialMatching[currentLeftVertexNum][shortestPathLeftVerticiesPreviousVertex[currentLeftVertexNum] ]);
			currentRightVertexNum = shortestPathLeftVerticiesPreviousVertex[currentLeftVertexNum] ;
		}

		capacityMatched = min(capacityMatched, leftVertexWeights[currentLeftVertexNum] - partialLeftSums[currentLeftVertexNum]);

		//Now we update our partial matching with the new path and its capacity computed above
		partialLeftSums[currentLeftVertexNum] += capacityMatched;
		partialRightSums[shortestPathLastVextex] += capacityMatched;

		currentRightVertexNum = shortestPathLastVextex;
		while (currentRightVertexNum != -2)
		{
			currentLeftVertexNum = shortestPathRightVerticiesPreviousVertex[currentRightVertexNum];
			
			partialMatching[currentLeftVertexNum][currentRightVertexNum] += capacityMatched;
			
			currentRightVertexNum = shortestPathLeftVerticiesPreviousVertex[currentLeftVertexNum];

			if(currentRightVertexNum != -2)
				partialMatching[currentLeftVertexNum][shortestPathLeftVerticiesPreviousVertex[currentLeftVertexNum]] -= capacityMatched;

		}

		//Grab new shortest path
		shortestPathLastVextex = shortestPath(leftVertexWeights, rightVertexWeights, partialMatching, partialLeftSums, partialRightSums, shortestPathLeftVerticiesPreviousVertex, shortestPathRightVerticiesPreviousVertex);
	}

	return;
}

/* This function prints the Optimal Transport Plan to the optimalTransportPlan output file*/
void writeOptimalTransportPlan(double partialMatching[][maxN])
{
	optimalTransportPlanFile << N << endl;

	for (int i = 0; i < (1 << N); i++)
	{
		for (int j = 0; j < N; j++)
			optimalTransportPlanFile << partialMatching[i][j] << " ";
		optimalTransportPlanFile << endl;
	}

	optimalTransportPlanFile << endl;
	return;
}

void prepareCellOutputFiles()
{
	cellOutputFile << "#This file contains the data on the cell decomposition." << endl << "#" << endl;	
	cellOutputFile << "#The first line is N, which is the number of warehouses." << endl;
	cellOutputFile << "#Then there are 2^N -1 blocks each representing one cell." << endl;
	cellOutputFile << "#The first line of each block gives b_i the number of disjoint rectangles in the cell." << endl;
	cellOutputFile << "#Each of the next b_i lines in a block give the coordinates of one rectangle in the cell." << endl;
	cellOutputFile << "#The rectangles are given as x0, x1, y0, y1." << endl;
	cellOutputFile << "#Note that the cell corresponding to the empty set isn't printed!" << endl;
	cellOutputFile << "#" << endl;
	cellOutputFile << N << endl;

	cellOutputFileWithMu << "#This file contains the data on the cell decomposition after intersection with the support of mu." << endl << "#" << endl;
	cellOutputFileWithMu << "#The first line is N, which is the number of warehouses." << endl;
	cellOutputFileWithMu << "#Then there are 2^N -1 blocks each representing one cell." << endl;
	cellOutputFileWithMu << "#The first line of each block gives the mu measure of the cell." << endl;
	cellOutputFileWithMu << "#The second line of each block gives b_i the number of disjoint rectangles in the cell." << endl;
	cellOutputFileWithMu << "#Each of the next b_i lines in a block give the coordinates of one rectangle in the cell." << endl;
	cellOutputFileWithMu << "#The rectangles are given as x0, x1, y0, y1." << endl;
	cellOutputFileWithMu << "#Note that the cell corresponding to the empty set isn't printed!" << endl;
	cellOutputFileWithMu << "#" << endl;
	cellOutputFileWithMu << N << endl;
}

int main()
{
	double cellSizes[1 << maxN];

	input_data();
	prepareCellOutputFiles();

	rectangle mu = rectangle(0, 4, 0, 4);
	pointerMu = &mu;

	solve(0, 10, 0.00001, cellSizes);

	cellOutputFile.close();

	double partialMatching[1 << maxN][maxN];
	edmondKarp(cellSizes, lambda, partialMatching);

	writeOptimalTransportPlan(partialMatching);

	return 0;
}