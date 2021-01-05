#include "iostream"
#include <math.h>

#define MATRIX_SIZE 2048
#define NUM_COLUMNS MATRIX_SIZE
#define NUM_ROWS MATRIX_SIZE

class Matrix
{
public:
	int num_columns;
	int num_rows;
	//number of elements between the beginnings of adjacent
	// rows in the memory layout (useful for representing sub-matrices)
	int pitch;
	//Pointer to the first element of the matrix represented
	float* elements;
public:
	Matrix()
	{
		
	}

	/* УМНОЖЕНИЕ МАТРИЦ */
	Matrix matrix_multiply(const Matrix A, const Matrix B)
	{
		Matrix C;
		C.num_columns = C.pitch = A.num_columns;
		C.num_rows = A.num_rows;
		unsigned int size = C.num_rows * C.num_columns;
		C.elements = (float*)malloc(size * sizeof(float));

		for (unsigned int i = 0; i < A.num_columns; i++)
			for (unsigned int j = 0; j < B.num_rows; j++) {
				double sum = 0.0f;
				for (unsigned int k = 0; k < A.num_columns; k++) {
					double a = A.elements[i * A.num_columns + k];
					double b = B.elements[k * B.num_rows + j];
					sum += a * b;
				}
				C.elements[i * B.num_rows + j] = (float)sum;
			}
		return C;
	}

	/* ВЫВОД МАТРИЦЫ */
	void print_matrix(const Matrix M) {
		for (unsigned int i = 0; i < M.num_rows; i++) {
			for (unsigned int j = 0; j < M.num_columns; j++)
				printf("%f ", M.elements[i * M.num_rows + j]);
			printf("\n");
		}
		printf("\n");
	}

	/* This function implements a row-oriented Cholesky decomposition on the input matrix A to generate an upper triangular matrix U
		such that A = U^TU.
	 */
	int chol_gold(const Matrix A, Matrix U)
	{
		int i, j, k;
		int size = A.num_rows * A.num_columns;

		// Копируем содержимое матрицы A в рабочую матрицу U
		for (i = 0; i < size; i++)
			U.elements[i] = A.elements[i];

		// Выполняем разложение Холецкого на месте матрицы U
		for (k = 0; k < U.num_rows; k++)
		{
			// Возьмите квадратный корень из диагонального элемента
			U.elements[k * U.num_rows + k] = sqrt(U.elements[k * U.num_rows + k]);
			if (U.elements[k * U.num_rows + k] <= 0) {
				printf("Cholesky decomposition failed. \n");
				return 0;
			}

			// Шаг деления
			for (j = (k + 1); j < U.num_rows; j++)
				U.elements[k * U.num_rows + j] /= U.elements[k * U.num_rows + k]; // Шаг деления


			// Шаг устранения
			for (i = (k + 1); i < U.num_rows; i++)
				for (j = i; j < U.num_rows; j++)
					U.elements[i * U.num_rows + j] -= U.elements[k * U.num_rows + i] * U.elements[k * U.num_rows + j];

		}

		// В качестве последнего шага обнуляем нижнюю треугольную часть U
		for (i = 0; i < U.num_rows; i++)
			for (j = 0; j < i; j++)
				U.elements[i * U.num_rows + j] = 0.0;

		printf("The Upper triangular matrix is: \n");
		print_matrix(U);

		return 1;
	}

	/* Helper function which checks the correctness of Cholesky decomposition by attempting to recover the original matrix A using U^TU. */
	int check_chol(const Matrix A, const Matrix U)
	{
		Matrix U_transpose;
		U_transpose.num_columns = U_transpose.pitch = U.num_columns;
		U_transpose.num_rows = U.num_rows;
		unsigned int size = U_transpose.num_rows * U_transpose.num_columns;
		U_transpose.elements = (float*)malloc(size * sizeof(float));

		// Determine the transpose of U
		unsigned int i, j;
		for (i = 0; i < U.num_rows; i++)
			for (j = 0; j < U.num_columns; j++)
				U_transpose.elements[i * U.num_rows + j] = U.elements[j * U.num_columns + i];

		// Multiply U and U_transpose
		Matrix A_recovered = matrix_multiply(U_transpose, U);
		// print_matrix(A_recovered);

		// Compare the two matrices A and A_recovered
		for (i = 0; i < size; i++)
			if (fabs(A.elements[i] - A_recovered.elements[i]) > 0.01)
				return 0;

		return 1;
	}

};

/* Функция проверяет, является ли матрица симметричной. */
int check_if_symmetric(const Matrix M)
{
	for (unsigned int i = 0; i < M.num_rows; i++)
		for (unsigned int j = 0; j < M.num_columns; j++)
			if (M.elements[i * M.num_rows + j] != M.elements[j * M.num_columns + i])
				return 0;
	return 1;
}

/* Функция проверяет, является ли матрица доминирующей по диагонали. */
int check_if_diagonal_dominant(const Matrix M)
{
	float diag_element;
	float sum;
	for (unsigned int i = 0; i < M.num_rows; i++) {
		sum = 0.0;
		diag_element = M.elements[i * M.num_rows + i];
		for (unsigned int j = 0; j < M.num_columns; j++) {
			if (i != j)
				sum += abs(M.elements[i * M.num_rows + j]);
		}
		if (diag_element <= sum)
			return 0;
	}

	return 1;
}

Matrix create_positive_definite_matrix(unsigned int num_rows, unsigned int num_columns)
{
	Matrix M;
	M.num_columns = M.pitch = num_columns;
	M.num_rows = num_rows;
	unsigned int size = M.num_rows * M.num_columns;
	M.elements = (float*)malloc(size * sizeof(float));

	// Step 1: Create a matrix with random numbers between [-.5 and .5]
	printf("Creating a %d x %d matrix with random numbers between [-.5, .5]...", num_rows, num_columns);
	unsigned int i;
	unsigned int j;
	for (i = 0; i < size; i++)
		M.elements[i] = ((float)rand() / (float)RAND_MAX) - 0.5;
	printf("done. \n");
	// print_matrix(M);
	// getchar();

	// Step 2: Make the matrix symmetric by adding its transpose to itself
	printf("Generating the symmetric matrix...");
	Matrix transpose;
	transpose.num_columns = transpose.pitch = num_columns;
	transpose.num_rows = num_rows;
	size = transpose.num_rows * transpose.num_columns;
	transpose.elements = (float*)malloc(size * sizeof(float));

	for (i = 0; i < M.num_rows; i++)
		for (j = 0; j < M.num_columns; j++)
			transpose.elements[i * M.num_rows + j] = M.elements[j * M.num_columns + i];
	// print_matrix(transpose);

	for (i = 0; i < size; i++)
		M.elements[i] += transpose.elements[i];
	if (check_if_symmetric(M))
		printf("done. \n");
	else {
		printf("error. \n");
		free(M.elements);
		M.elements = NULL;
	}
	// print_matrix(M);
	// getchar();

	// Step 3: Make the diagonal entries large with respect to the row and column entries
	printf("Generating the positive definite matrix...");
	for (i = 0; i < num_rows; i++)
		for (j = 0; j < num_columns; j++) {
			if (i == j)
				M.elements[i * M.num_rows + j] += 0.5 * M.num_rows;
		}
	if (check_if_diagonal_dominant(M))
		printf("done. \n");
	else {
		printf("error. \n");
		free(M.elements);
		M.elements = NULL;
	}
	// print_matrix(M);
	// getchar();

	free(transpose.elements);

	// M is diagonally dominant and symmetric!
	return M;
}



int main()
{
//    A[0][0] = 6.25; A[1][0] = -1; A[2][0] = 0.5;
//    A[0][1] = -1; A[1][1] = 5;  A[2][1] = 2.12;
//    A[0][2] = 0.5; A[1][2] = 2.12; A[2][2] = 3.6;
//    
//    B[0] = 7.5;
//    B[1] = -8.68;
//    B[2] = -0.24;
//
	Matrix M;
	M.num_columns = M.pitch = 4;
	M.num_rows = 3;
	int size = M.num_rows * M.num_columns;
	M.elements = (float*)malloc(size * sizeof(float));

	M.elements[0] = 6.25;
	M.elements[1] = -1;
	M.elements[2] = 0.5;
	M.elements[3] = 7.5;

	M.elements[4] = -1;
	M.elements[5] = 5;
	M.elements[6] = 2.12;
	M.elements[7] = -8.68;

	M.elements[8] = 0.5;
	M.elements[9] = 2.12;
	M.elements[10] = 3.6;
	M.elements[11] = -0.24;

	Matrix U;
	U.num_columns = U.pitch = 4;
	U.num_rows = 3;
	size = U.num_rows * U.num_columns;
	U.elements = (float*)malloc(size * sizeof(float));

	M.chol_gold(M, U);

}


//
//
//#include <iostream>
//
//using namespace std;
//
//double FindDet(double** a, int n) 
//{
//    double p, kst = 0;
//    int t = 0;
//    p = 0;
//    for (int i = 0; i < n - 1; i++)
//    {
//        t = 1;
//        while (a[i][i] == 0)
//        {
//            for (int j = 0; j < n; j++)
//            {
//                a[i][j] = kst;
//                a[i][j] = a[i + t][j];
//                a[i + t][j] = kst;
//            }
//            p++;
//            t++;
//        }
//
//        for (int k = i + 1; k < n; k++)
//        {
//            kst = a[k][i] / a[i][i];
//            for (int j = 0; j < n; j++)
//                a[k][j] -= a[i][j] * kst;
//        }
//    }
//
//    kst = pow(-1, p);
//    for (int i = 0; i < n; i++)
//        kst *= a[i][i];
//   
//    return kst;
//}
//
//void PRINT(double** A, int N) 
//{
//    cout << endl << "Matrix:" << endl;
//    
//    for (int i = 0; i < N; i++) 
//    {
//        for (int j = 0; j < N; j++) 
//        {
//            cout << A[i][j] << " ";
//        }
//        cout << endl;
//    }
//}
//
//int main()
//{
//    int N = 3;
//   
//    double** A = new double* [N];
//    double** L = new double* [N];
//    double** Lt = new double* [N];
//    
//    for (int i = 0; i < N; i++) 
//    {
//        A[i] = new double[N];
//        L[i] = new double[N];
//        Lt[i] = new double[N];
//    }
//
//    for (int i = 0; i < N; i++) 
//    {
//        for (int j = 0; j < N; j++) 
//        {
//            A[i][j] = 0.0;
//            L[i][j] = 0.0;
//            Lt[i][j] = 0.0;
//        }
//    }
//
//    double* Y = new double[N];
//    double* B = new double[N];
//    double* Det = new double[N];
//    double* X = new double[N];
//    
//    for (int i = 0; i < N; i++) 
//    {
//        Y[i] = 0;
//        B[i] = 0;
//        Det[i] = 0;
//        X[i] = 0;
//    }
//        
//
//    A[0][0] = 6.25; A[1][0] = -1; A[2][0] = 0.5;
//    A[0][1] = -1; A[1][1] = 5;  A[2][1] = 2.12;
//    A[0][2] = 0.5; A[1][2] = 2.12; A[2][2] = 3.6;
//    
//    B[0] = 7.5;
//    B[1] = -8.68;
//    B[2] = -0.24;
//
//    for (int i = 0; i < N; i++) 
//    {
//        for (int j = 0; j < N; j++)
//        {
//            cout << A[i][j] << " ";
//        }
//        cout << B[i] << endl;
//    }
//
//    for (int i = 0; i < N; ++i) 
//    {
//        const int n = i + 1;
//        
//        double** temp = new double* [N];
//        for (int i = 0; i < N; i++)
//        {
//            temp[i] = new double[N];
//        }
//
//        for (int i = 0; i < N; i++)
//        {
//            for (int j = 0; j < N; j++)
//            {
//                temp[i][j] = 0.0;
//            }
//        }
//
//        for (int j = 0; j < n; j++)
//        {
//            for (int k = 0; k < n; k++)
//            {
//                temp[j][k] = A[j][k];
//            }
//        }
//        PRINT(temp, n);
//        Det[i] = FindDet(temp, n);
//    }
//
//    for (int i = 0; i < N; i++)
//    {
//        if (Det[i] < 0)
//            return -1;
//    }
//
//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < N; j++)
//        {
//            if ((i == 0) && (j == 0))
//                L[i][j] = sqrt(A[i][j]);
//            else if (j == 0) 
//            {
//                if (L[0][0] != 0)
//                    L[i][j] = A[i][0] / L[0][0];
//            }
//            else if (i == j) 
//            {
//                double temp = A[i][j];
//                for (int i1 = 0; i1 < i; ++i1) 
//                {
//                    temp -= pow(L[i][i1], 2);
//                }
//                L[i][j] = sqrt(temp);
//            }
//            else 
//            {
//                double temp = A[i][j];
//                for (int k = 0; k < i; ++k) 
//                {
//                    temp -= L[i][k] * L[j][k]; //check
//                }
//                if (L[j][j] != 0)
//                    L[i][j] = temp / L[j][j];
//            }
//        }
//    }
//
//    PRINT(L, N);
//    double _temp = 0.0;
//    int j;
//    
//    for (int i = N - 1; i >= 0; --i) 
//    {
//        _temp = 0.0;
//        for (j = N - 1; j > i; --j) 
//        {
//            _temp += L[N - i - 1][N - j - 1] * Y[N - j - 1];
//        }
//        if ((N - i - 1 >= 0) && (N - i - 1 < 3))
//            Y[N - i - 1] = (B[N - i - 1] - _temp) / L[N - i - 1][N - j - 1];
//        //    cout << Y[N - i - 1] << " ";
//    }
//    
//    for (int i = 0; i < N; ++i) 
//    {
//        cout << Y[i] << " ";
//    }
//
//    for (int i = 0; i < N; ++i) 
//    {
//        for (int j = 0; j < N; ++j)
//        {
//            Lt[i][j] = L[j][i];
//        }
//    }
//
//    PRINT(Lt, N);
//
//    for (int i = N - 1; i >= 0; --i) 
//    {
//        _temp = 0.0;
//        for (j = N - 1; j > i; --j) 
//        {
//            _temp += Lt[i][j] * X[j];
//        }
//        if ((N - i - 1 >= 0) && (N - i - 1 < 3))
//            X[i] = (Y[i] - _temp) / Lt[i][j];
//        //    cout << Y[N - i - 1] << " ";
//    }
//
//    for (int i = 0; i < N; ++i)
//    {
//        cout << X[i] << " ";
//    }
//}