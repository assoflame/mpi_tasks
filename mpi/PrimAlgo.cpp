#include <iostream>
#include <mpi.h>

#define n 5

//void createAndShareAdjacencyMatrix(int rank)
//{
//	
//}

int main(int argc, char** argv)
{
	int rank, size;
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int matrix[n][n];

	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
			for (int j = i; j < n; ++j)
			{
				matrix[i][j] = i == j ? 0 : rand() % 21;
				matrix[j][i] = matrix[i][j];
			}
	}
	MPI_Bcast(&matrix[0][0], n * n, MPI_INT, 0, MPI_COMM_WORLD);




	MPI_Finalize();

	return EXIT_SUCCESS;
}