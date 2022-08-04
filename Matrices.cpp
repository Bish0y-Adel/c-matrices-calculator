#include <iostream>
#include <limits>
#include <cmath>
using namespace std;
void matrix_cin(float** matrix, int num_row, int num_col);// 1
void sum_matrix(float** mat1, float** mat2, float** result, int& num_row, int& num_col);// 2
void sub_matrix(float** mat1, float** mat2, float** result, int& num_row, int& num_col);// 3
void print(float** result, int& num_row, int& num_col);// 4
void product(float** mat1, float** mat2, float** result, int& mat1_row, int& mat1_col, int& mat2_row, int& mat2_col);// 5
float determinant(float** matrix, int n);// 6
void getCofactor(float** matrix, float** temp, int p, int q, int n); // 7
void adjoint(float** matrix, float** adj, int N);// 8
bool inverse(float** matrix, float** inverse, int N);// 9
void transpose(float** matrix, int rowsize, int colsize);// 10
void conssum(float** matrix, float** result, int rowsize, int colsize, int cons);// 11
void consmult(float** matrix, float** result, int rowsize, int colsize, int cons);// 12
void consdiv(float** matrix, float** result, int rowsize, int colsize, int cons);// 13

int main()
{
	int  operation, row1size = 0, row2size = 0, col1size = 0, col2size = 0, cons=0;
	cout << "\n************************************************************************************\n\n";
	cout << "\t\t\t\t============================";
	cout << "\n\t\t\t\t|| Built by ' RTCC ' Team || \n";
	cout << "\t\t\t\t============================";
	cout << "\n\n************************************************************************************";
	while (1)
	{
		cout << "\n\n****************************************MENUE***************************************\n";
		
		cout << " press ' 1 ' to add tow matrices \n press ' 2 ' to subtract tow matrices \n press ' 3 ' to multiply tow matrices\n press ' 4 ' to find the determinant of matrix\n"
		<< " press ' 5 ' to find tarnspose of matrix \n press ' 6 ' to find Adjoint of matrix \n press ' 7 ' to find Inverse of matrix\n"
		<< " press ' 8 ' to add matrix to constant value \n press ' 9 ' to multiply matrix by const value\n press '10 ' to devide matrix by const value\n press '11 ' to exit the program\n\n";

		cout << "************************************************************************************\n";
		//cout << "Select which operations you want to perform from the previous sentences ?  ";
		//cin >> operation;
		//-----------------------------------------------------------------
		do {
			cout << "\nSelect which operations you want to perform from the previous sentences ?  ";
			if (cin >> operation) break;
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
		} while (true);
		//-----------------------------------------------------------------
		//accept order of a matrix 1 from user; 
		if (operation == 4 || operation == 6 || operation == 7)
		{
			//------------------------------------------------------------------------------------------
			do {
				cout << "\nEnter the order of matrix: ";
				if (cin >> row1size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//------------------------------------------------------------------------------------------
			//cout << "Enter the order of matrix: ";
			//cin >> row1size;
			col1size = row1size;
			col2size = row1size;
		}
		else if (operation == 5 || operation == 9 || operation == 8 || operation == 10)
		{
			//-------------------------------------------------------------------------------
			do {
				cout << "\nEnter number of rows of this matrix: ";
				if (cin >> row1size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//-------------------------------------------------------------------------------
			do {
				cout << "\nEnter number of columns of this matrix: ";
				if (cin >> col1size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//-------------------------------------------------------------------------------
			//cout << "Enter number of rows and columns of this matrix: \n";
			//cin >> row1size;
			//cin >> col1size;
			col2size = col1size;
			if (operation == 8 || operation == 9 || operation == 10)
			{
				//-----------------------------------------------------------------
				do {
					cout << "\nEnter the const value: ";
					if (cin >> cons) break;
					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
				} while (true);
				//-----------------------------------------------------------------
				//cout << "Enter the const value: ";
				//cin >> cons;
			}
		}
		else if (operation == 1 || operation == 2 || operation == 3)
		{
			
			//---------------------------------------------------------------------------
			do {
				cout << "\nEnter number of rows of the first matrix: ";
				if (cin >> row1size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//---------------------------------------------------------------------------
			do {
			cout << "\nEnter number of columns of the first matrix: ";
			if (cin >> col1size) break;
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
		    } while (true);
			//---------------------------------------------------------------------------
			do {
				cout << "\nEnter number of rows of the socend matrix: ";
				if (cin >> row2size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//---------------------------------------------------------------------------
			do {
				cout << "\nEnter number of columns of the socend matrix: ";
				if (cin >> col2size) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//---------------------------------------------------------------------------
			//cout << "Enter number of rows and columns of the first matrix: \n";
			//cin >> row1size;
			//cin >> col1size;
			//cout << "Enter number of rows and columns of the socend matrix: \n";
			//cin >> row2size;
			//cin >> col2size;
		}
		else if (operation==11)
			break;
		else
			cout << "\n\t\t\t\t\tInvalid Operayion\n\n\t\t\t\t     Enter a valid operation\n\n";
		//********************************************************************
					// matrix one decleratin
		float** matrix1 = new float* [row1size];
		for (int x = 0; x < row1size; x++)
			matrix1[x] = new float[col1size];
		// matrix tow decleratin
		float** matrix2 = new float* [row2size];
		for (int x = 0; x < row2size; x++)
			matrix2[x] = new float[col2size];
		// resultant matrix decleratin
		float** result = new float* [row1size];
		for (int x = 0; x < row1size; x++)
			result[x] = new float[col2size];
		// temporary matrix function
		float** temp = new float* [row1size];
		for (int x = 0; x < row1size; x++)
			temp[x] = new float[col1size];
		// adj matrix declaration
		float** adj = new float* [row1size];
		for (int x = 0; x < row1size; x++)
			adj[x] = new float[row1size];
		//  inverse matrix declaration
		float** inversee = new float* [row1size];
		for (int x = 0; x < row1size; x++)
			inversee[x] = new float[row1size];
//*********************************************************************
//*********************************************************************
		switch (operation)
		{
		case 1:  // sum of tow matrix
			if (row1size == row2size && col1size == col2size) {
				//********************************************************************
				//matrix 1 cin
				matrix_cin(matrix1, row1size, col1size);
				//matrix 2 cin
				matrix_cin(matrix2, row2size, col2size);

				sum_matrix(matrix1, matrix2, result, row1size, col1size);
				cout << "The sum of the Tow matric =  "; print(result, row1size, col1size);

				break;
			}
			else
				cout << "\n\t\t\t\t\t\tE_R_R_O_R\n\n\t\t\tTow matrices should have the same order\n\n";
			break;
//*******************************************************************************************************************************
//*******************************************************************************************************************************
		case 2:
			//subtertion of tow matrices
			if (row1size == row2size && col1size == col2size) {
				//********************************************************************
				//matrix 1 cin
				matrix_cin(matrix1, row1size, col1size);
				//matrix 2 cin
				matrix_cin(matrix2, row2size, col2size);
				//*********************************************************************
				sub_matrix(matrix1, matrix2, result, row1size, col1size);
				cout << "The subtraction of the Tow matrices = "; print(result, row1size, col1size);
				break;
			}
			else
				cout << "\n\t\t\t\t\t\tE_R_R_O_R\n\n\t\t\tTow matrices should have the same order\n\n";
			break;
//*******************************************************************************************************************************
		case 3:
			//product of tow matrices
			if (col1size == row2size) {
				//********************************************************************
				//matrix 1 cin
				matrix_cin(matrix1, row1size, col1size);
				//matrix 2 cin
				matrix_cin(matrix2, row2size, col2size);
				//*********************************************************************
				cout << "The product of the Tow matrices = ";
				product(matrix1, matrix2, result, row1size, col1size, row2size, col2size);
				print(result, row1size, col2size);

				break;
			}

			else
				cout << "\n\t\t\t\t\t\tE_R_R_O_R\n\n\t\t Num of columns of the first matrix should equal to the number of rows of the socend matrix\n\n";
			break;
//*****************************************************************************************************************************
		case 4:
			// determinant of matrix
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			cout << "\n Determinant of the matrix = " << determinant(matrix1, row1size) << "\n\n";
			break;
//****************************************************************************************************************************
		case 5:
			// transpose of matrix
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			transpose(matrix1, row1size, col1size);
			break;
//****************************************************************************************************************************
		case 6:
			//adjoint of matrix
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			adjoint(matrix1, adj, row1size);
			cout << "\nthe Adjoint matrix = \n";
			print(adj, row1size, row1size);
			break;
//****************************************************************************************************************************
		case 7:
			//Inverse of matrix
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			inverse(matrix1, inversee, row1size);
			cout << "\nthe Inverse matrix = \n";
			print(inversee, row1size, row1size);
			break;
//****************************************************************************************************************************
		case 8:
			//add by constant
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			cout << "\nthe rusultant matrix = \n";
			conssum(matrix1, result, row1size, col1size, cons);
			print(result, row1size, col1size);
			break;
//****************************************************************************************************************************
		case 9:
			//add by constant
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			cout << "\nthe rusultant matrix = \n";
			consmult(matrix1, result, row1size, col1size, cons);
			print(result, row1size, col1size);
			break;
//*****************************************************************************************************************************
		case 10:
			//add by constant
			matrix_cin(matrix1, row1size, col1size);
			cout << "The entered matrix is:" << endl;
			print(matrix1, row1size, col1size);
			cout << "\nthe rusultant matrix = \n";
			consdiv(matrix1, result, row1size, col1size, cons);
			print(result, row1size, col1size);
			break;

		}
//*********************************************************************
		/*
		// deleting all matrices
		for (int x = 0; x < row1size; x++)
			delete[] matrix1[x];
		delete[] matrix1;
		for (int x = 0; x < row2size; x++)
			delete[] matrix2[x];
		delete[] matrix2;
		for (int x = 0; x < row1size; x++)
			delete[] result[x];
		delete[] result;
		for (int x = 0; x < row1size; x++)
			delete[] temp[x];
		delete[] temp;
		for (int x = 0; x < row1size; x++)
			delete[] adj[x];
		delete[] adj;
		for (int x = 0; x < row1size; x++)
			delete[] inversee[x];
		delete[] inversee;
//*********************************************************************
*/
	}
	
	

	cout << "\n\n\t\t\t\t=============================================";
	cout << "\n\t\t\t\t|| 2022 All Rights Reserved, ' RTCC ' Team || \n";
	cout << "\t\t\t\t=============================================\n\n\n\n";
	system("pause");
}

/*
            // deleting all matrices
			for (int x = 0; x < row1size; x++)
				delete[] matrix1[x];
			delete[] matrix1;
			for (int x = 0; x < row2size; x++)
				delete[] matrix2[x];
			delete[] matrix2;
			for (int x = 0; x < row1size; x++)
				delete[] result[x];
			delete[] result;
			for (int x = 0; x < row1size; x++)
				delete[] temp[x];
			delete[] temp;
			for (int x = 0; x < row1size; x++)
				delete[] adj[x];
			delete[] adj;
			for (int x = 0; x < row1size; x++)
				delete[] inversee[x];
			delete[] inversee;

*/

//input function
void matrix_cin(float** matrix, int num_row, int num_col)
{
	cout << "Enter matrix elements: \n";
	for (int row = 0; row < num_row; row++)
	{
		for (int col = 0; col < num_col; col++)
		{
			//-------------------------------------------------------------------
			do {
				cout << endl;
				if (cin >> matrix[row][col]) break;
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\t\t\t\t    ERORR\n\t\t\tEnter a valid integer number";
			} while (true);
			//-------------------------------------------------------------------
			//cin >> matrix[row][col];
		}
			
		cout << endl;
	}
}
//-----------------------------------------------------------------------------------------------------------
// Adding functin
void sum_matrix(float** mat1, float** mat2, float** result, int& num_row, int& num_col)
{
	for (int row = 0; row < num_row; row++)
	{
		for (int col = 0; col < num_col; col++)
		{
			result[row][col] = mat1[row][col] + mat2[row][col];
		}
	}
}
//-----------------------------------------------------------------------------------------------------------
//Subtraction function
void sub_matrix(float** mat1, float** mat2, float** result, int& num_row, int& num_col)
{
	for (int row = 0; row < num_row; row++)
	{
		for (int col = 0; col < num_col; col++)
		{
			result[row][col] = mat1[row][col] - mat2[row][col];
		}
	}
}
//-----------------------------------------------------------------------------------------------------------
//output function
void print(float** result, int& num_row, int& num_col)
{
	cout << endl;
	for (int row = 0; row < num_row; row++)
	{
		for (int col = 0; col < num_col; col++)
		{
			cout << result[row][col];
			cout << "\t";
		}
		cout << endl;
	}

}
//-----------------------------------------------------------------------------------------------------------
// multiplication function
void product(float** mat1, float** mat2, float** result, int& mat1_row, int& mat1_col, int& mat2_row, int& mat2_col)
{

	for (int row = 0; row < mat1_row; row++) {
		for (int col = 0; col < mat2_col; col++) {
			float sum = 0.0;
			for (int x = 0; x < mat1_col; x++)
				sum += (mat1[row][x] * mat2[x][col]);
			result[row][col] = sum;
		}
	}
}
//------------------------------------------------------------------------------------------------------------
// Function to get cofactor of A[p][q] in temp[][].
void getCofactor(float** matrix, float** temp, int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = matrix[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}

}
// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(float** matrix, float** adj, int N)
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	float** temp = new float* [N];
	for (int x = 0; x < N; x++)
		temp[x] = new float[N];
	int sign = 1;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// Get cofactor of A[i][j]
			getCofactor(matrix, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(temp, N - 1));
		}
	}

}
// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(float** matrix, float** inverse, int N)
{
	// Find determinant of A[][]
	float det = determinant(matrix, N);
	if (det == 0)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	float** adj = new float* [N];
	for (int x = 0; x < N; x++)
		adj[x] = new float[N];
	adjoint(matrix, adj, N);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			inverse[i][j] = adj[i][j] / float(det);

	return true;
}
// function to get transpose
void transpose(float** matrix, int rowsize , int colsize)
{
	float** temp = new float* [colsize];
	for (int x = 0; x < colsize; x++)
		temp[x] = new float[rowsize];

	for (int row = 0; row < rowsize; row++)
		for (int col = 0; col < colsize; col++)
		{
			temp[col][row] = matrix[row][col];
		}
	cout << "the transpose of this matrix is: \n\n";
	for (int i = 0; i < colsize; i++) {
		for (int j = 0; j < rowsize; j++)
			cout << temp[i][j] << " ";
		cout << endl;
	}
} 
//************************************************************************************************************
// Recursive function for finding determinant of matrix.
// n is current dimension of A[][].
float determinant(float** matrix, int n)
{
	float D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1)
		return matrix[0][0];
	// To store cofactors
	float** temp = new float* [n];
	for (int x = 0; x < n; x++)
		temp[x] = new float[n];

	int sign = 1;  // To store sign multiplier

	 // Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f]
		getCofactor(matrix, temp, 0, f, n);
		D += sign * matrix[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}
//************************************************************************************************************
//function matrix multiplication by constant
void consmult (float** matrix,float** result ,int rowsize,int colsize,int cons)
{
	for (int row = 0; row < rowsize; row++) {
		for (int col = 0; col < colsize; col++)
		{
			result[row][col] = cons * matrix[row][col];
		}
	}
}
//function matrix add by constant
void conssum(float** matrix, float** result, int rowsize, int colsize, int cons)
{
	for (int row = 0; row < rowsize; row++) {
		for (int col = 0; col < colsize; col++)
		{
			result[row][col] = cons + matrix[row][col];
		}
	}
}
//************************************************************************************************************
void consdiv(float** matrix, float** result, int rowsize, int colsize, int cons)
{
	for (int row = 0; row < rowsize; row++) {
		for (int col = 0; col < colsize; col++)
		{
			result[row][col] = matrix[row][col] / cons;
		}
	}
}
//RTCC means: Runtime Terror Chasers Clan ;