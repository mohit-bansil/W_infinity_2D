#include<iostream>
#include<fstream>
#include <algorithm>
#include <bitset> 

using namespace std;

const int maxN = 10;
const double tolerance = 0.0000001;

int N;
double y[maxN][2];
double lambda[maxN];

ifstream inputFile("Input_Data.txt");
ofstream cellOutputFile("Cell_Data.txt");

/*
Rectangle object is charaterized by its left and right x coordinates and its upper and lower y coordinates
Currently methods are
1. Determine if rectangle is non-empty (non-degenerate)
2. Find area
3. Find mu area (area of rectangle intersected with mu)
4. Determine if a point is inside the rectangle
5. Print the coordinates of the rectangle to cellOutputFile
6. intersect two rectangles

Also note that the definition of mu is inside of the mu_area function
*/
struct rectangle;
rectangle intersect(rectangle A, rectangle B);
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
		return intersect(*this, rectangle(0, 4, 0, 4)).area() / 16;
	}
	bool is_inside(double x, double y)
	{
		if (x < x0 || x > x1 || y < y0 || y > y1)
			return false;
		else
			return true;
	}
	void print()
	{
		cellOutputFile << x0 << " " << x1 << " " << y0 << " " << y1 << endl;
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
5. Delete all degenerate rectangles from cell.

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

			if (pieces[i].x1 > A.x0)
			{
				pieces[i].x1 = A.x0;
				warningFlag++;
			}
			if (pieces[i].x0 < A.x1)
			{
				pieces[i].x0 = A.x1;
				warningFlag++;
			}
			if (pieces[i].y1 > A.y0)
			{
				pieces[i].y1 = A.y0;
				warningFlag++;
			}
			if (pieces[i].y0 < A.y1)
			{
				pieces[i].y0 = A.y1;
				warningFlag++;
			}

			if (warningFlag > 1)
				cerr << "Warning something might be wrong. Multiple sides were cut from a cell piece." << endl;
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
void compute_cell_graph(const double omega, double cellsizes[], bool printCells)
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

		if (leftsum > rightsum)
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
The resulting cell decomposition and the upper/lower omega bounds are written to the cellOutputFile*/
void solve(double lowerOmegaBound, double upperOmegaBound, double desiredError, double cellSizes[])
{
	double currentOmegaGuess;

	while (upperOmegaBound - lowerOmegaBound > desiredError)
	{
		currentOmegaGuess = (lowerOmegaBound + upperOmegaBound) / 2;
		compute_cell_graph(currentOmegaGuess, cellSizes, false);
		if (is_there_perfect_matching(cellSizes))
			upperOmegaBound = currentOmegaGuess;
		else
			lowerOmegaBound = currentOmegaGuess;
	}

	compute_cell_graph(upperOmegaBound, cellSizes, true);

	cellOutputFile << endl << lowerOmegaBound << " " << upperOmegaBound << endl;

	return;
}


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

void edmondKarp(double leftVertexWeights[], double rightVertexWeights[])
{


	return;
}

bool shortestPath(const double leftVertexWeights[], const double rightVertexWeights[], double partialMatching[][maxN], double partialLeftSums[], double partialRightSums[])
{
	vertex queue[1 << maxN];
	int queueCurrentPosition = 0;
	int queueSize = 0;


	int leftVerticiesDistanceFromSource[1 << maxN];
	int rightVerticiesDistanceFromSource[maxN];

	//Initialize all distances to +infinity (represented by -1)
	for (int i = 0; i < (1 << N); i++)
		leftVerticiesDistanceFromSource[i] = -1;
	for (int i = 0; i < N; i++)
		rightVerticiesDistanceFromSource[i] = -1;

	//Set all left verticies connected to the source to have distance one and put them in the queue
	for (int i = 1; i < (1 << N); i++)
	{
		if (leftVertexWeights[i] < partialLeftSums[i] - tolerance)
		{
			leftVerticiesDistanceFromSource[i] = 1;
			queue[queueSize] = vertex(i, true);
			queueSize++;
		}
	}

	for (queueCurrentPosition; queueCurrentPosition < queueSize; queueCurrentPosition++)
	{
		int currentVertexNum = queue[queueCurrentPosition].num;

		//current vertex in queue is a left vertex
		if (queue[queueCurrentPosition].isLeft)
		{
			bitset<32>currentLeftVertex(currentVertexNum);
			for (int j = 0; j < N; j++)
			{
				if (currentLeftVertex[j] == 1 && rightVerticiesDistanceFromSource[j] == -1)
				{
					//Check if we can get to the sink!!!
					if (partialRightSums[j] < rightVertexWeights[j] - tolerance)
					{
						cout << "I found the shorted path!" << endl;
						return true;
					}

					rightVerticiesDistanceFromSource[j] = leftVerticiesDistanceFromSource[currentVertexNum] + 1;
					queue[queueSize] = vertex(j, false);
					queueSize++;
				}
			}
		}
		//current vertex in queue is a right vertex
		else
		{
			for (int j = 1; j < (1 << N); j++)
			{
				bitset<32>leftVertex(j);

				if (partialMatching[j][currentVertexNum] > tolerance && leftVertex[currentVertexNum] == 1 && leftVerticiesDistanceFromSource[j] == -1)
				{
					leftVerticiesDistanceFromSource[j] = rightVerticiesDistanceFromSource[currentVertexNum] + 1;
					queue[queueSize] = vertex(j, true);
					queueSize++;
				}
			}
		}
	}

	cout << "There is no path" << endl;

	return false;
}


int main()
{
	double cellSizes[1 << maxN];

	input_data();

	cellOutputFile << N << endl;

	solve(0, 100, 0.1, cellSizes);

	cellOutputFile.close();

	//system("pause");
	return 0;
}