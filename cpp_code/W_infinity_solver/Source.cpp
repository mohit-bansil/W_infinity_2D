#include<iostream>
#include<fstream>
#include <algorithm>
#include <bitset> 

using namespace std;

const int N = 3;
double y[N][2] = { { 0,0 },{ 1,1 },{ 2,2 } };
double lambda[N];

ofstream OutputFile("Cell_Data.txt");

/*
Rectangle object is charaterized by its left and right x coordinates and its upper and lower y coordinates
Currently methods are
1. Determine if rectangle is non-empty (non-degenerate)
2. Find area
3. Find mu area (area of rectangle intersected with mu)
4. Determine if a point is inside the rectangle
5. Print the coordinates of the rectangle to OutputFile
6. Intersect two rectangles

Also note that the definition of mu is inside of the mu_area function
*/
struct rectangle;
rectangle Intersect(rectangle A, rectangle B);
struct rectangle
{
	//data
	double x0, x1, y0, y1;

	//methods
	bool isNondegenerate()
	{
		if (x1 < x0 || y1 < y0)
			return false;
		else
			return true;
	}
	double area()
	{
		if (isNondegenerate())
			return (x1 - x0)*(y1 - y0);
		else
			return 0;
	}
	double mu_area()
	{
		return Intersect(*this, rectangle(0, 4, 0, 4)).area() / 16;
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
		OutputFile << x0 << " " << x1 << " " << y0 << " " << y1 << endl;
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
rectangle Intersect(rectangle A, rectangle B)
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
1. Intersect the cell with a given rectangle.
2. Intersect the cell with the complement of a given rectangle (Outersect).
3. Compute the mu area of the cell. 
4. Print the coordinates of each rectangle in the cell to OutputFile
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
	//Note that an application of Outsect_with_Rectangle can increase the number of pieces by at most 8
	rectangle pieces[10 * N + 4];
	

	//methods
	void Intersect_with_Rectangle(rectangle A)
	{

		for (int i = 0; i < count; i++)
		{
			pieces[i] = Intersect(pieces[i], A);
		}

		return;
	}
	void Outersect_with_Rectangle(rectangle A)
	{
		int warning_flag;

		int new_pieces_added = 0;

		for (int i = 0; i < count; i++)
		{
			warning_flag = 0;

			if (Intersect(pieces[i], A).isNondegenerate() != true)
				continue;

			if (pieces[i].is_inside(A.x0, A.y0))
			{
				pieces[count + new_pieces_added] = rectangle(pieces[i].x0, A.x0, A.y0, pieces[i].y1);
				new_pieces_added++;
				pieces[count + new_pieces_added] = rectangle(A.x0, pieces[i].x1, pieces[i].y0, A.y0);
				new_pieces_added++;

				pieces[i].x1 = A.x0;
				pieces[i].y1 = A.y0;

				continue;
			}
			if (pieces[i].is_inside(A.x1, A.y0))
			{
				pieces[count + new_pieces_added] = rectangle(A.x1, pieces[i].x1, A.y0, pieces[i].y1);
				new_pieces_added++;
				pieces[count + new_pieces_added] = rectangle(pieces[i].x0, A.x1, pieces[i].y0, A.y0);
				new_pieces_added++;

				pieces[i].x0 = A.x1;
				pieces[i].y1 = A.y0;

				continue;
			}
			if (pieces[i].is_inside(A.x0, A.y1))
			{
				pieces[count + new_pieces_added] = rectangle(pieces[i].x0, A.x0, pieces[i].y0, A.y1);
				new_pieces_added++;
				pieces[count + new_pieces_added] = rectangle(A.x0, pieces[i].x1, A.y1, pieces[i].y1);
				new_pieces_added++;

				pieces[i].x1 = A.x0;
				pieces[i].y0 = A.y1;

				continue;
			}
			if (pieces[i].is_inside(A.x1, A.y1))
			{
				pieces[count + new_pieces_added] = rectangle(pieces[i].x0, A.x1, A.y1, pieces[i].y1);
				new_pieces_added++;
				pieces[count + new_pieces_added] = rectangle(A.x1, pieces[i].x1, pieces[i].y0, A.y1);
				new_pieces_added++;

				pieces[i].x0 = A.x1;
				pieces[i].y0 = A.y1;

				continue;
			}

			if (pieces[i].x1 > A.x0)
			{
				pieces[i].x1 = A.x0;
				warning_flag++;
			}
			if (pieces[i].x0 < A.x1)
			{
				pieces[i].x0 = A.x1;
				warning_flag++;
			}
			if (pieces[i].y1 > A.y0)
			{
				pieces[i].y1 = A.y0;
				warning_flag++;
			}
			if (pieces[i].y0 < A.y1)
			{
				pieces[i].y0 = A.y1;
				warning_flag++;
			}

			if (warning_flag > 1)
				cerr << "Warning something is wrong. Multiple sides were cut from a cell piece." << endl;
		}

		count += new_pieces_added;

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
		OutputFile << count << endl;
		for (int i = 0; i < count; i++)
			pieces[i].print();

	}
	void remove_degenerate_rectangles()
	{
		for (int i = 0; i < count; i++)
			if (pieces[i].isNondegenerate() == false)
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

/*This function takes in an omega stores the omega-transport graph in the array cellsizes.
It also prints the cells to the Output File */
void compute_cell_graph(const double omega, double cellsizes[])
{
	rectangle warehouse_rectangles[N];

	for (int i = 0; i < N; i++)
	{
		warehouse_rectangles[i] = generate_warehouse_rectangle(i, omega);
	}

	cellsizes[0] = 0;

	cell currentCell;

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

		currentCell = cell(warehouse_rectangles[firstRectangle]);

		for (int j = 0; j < N; j++)
		{
			if (j == firstRectangle)
				continue;
			else if (b[j] == 1)
				currentCell.Intersect_with_Rectangle(warehouse_rectangles[j]);
			else
				currentCell.Outersect_with_Rectangle(warehouse_rectangles[j]);
		}

		cellsizes[i] = currentCell.mu_area();

		currentCell.remove_degenerate_rectangles();
		currentCell.print();

	}

	return;
}

int main()
{
	OutputFile << N << endl;

	double cellsizes[1 << N];

	compute_cell_graph(0.6, cellsizes);

	OutputFile.close();

	//system("pause");
	return 0;
}