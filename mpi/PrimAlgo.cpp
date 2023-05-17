#include <iostream>
#include <mpi.h>
#include <set>

#define n 4    // |V|
#define maxWeight 5


int main(int argc, char** argv)
{
	int rank, size;
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int matrix[n][n];
	bool visited[n];
	std::pair<int, int> tree[n - 1];  // save edges that we are include in the 'minimum spanning tree' (MST), (our result)

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int localCount = rest / size;
	displs[0] = 0;
	counts[0] = localCount;
	for (int i = 1; i < size; ++i)
	{
		rest -= localCount;
		localCount = rest / (size - i);
		counts[i] = localCount;
		displs[i] = counts[i - 1] + displs[i - 1];
	}

	if (rank == 0)
	{
		printf("counts:\n");
		for (int i = 0; i < size; ++i)
			printf("%d ", counts[i]);
		printf("\ndispls:\n");
		for (int i = 0; i < size; ++i)
			printf("%d ", displs[i]);
		printf("\n\n");
	}

	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
		{
			visited[i] = false;
			for (int j = i; j < n; ++j)
			{
				matrix[i][j] = i == j ? INT_MAX : rand() % (maxWeight + 1);
				matrix[j][i] = matrix[i][j];
			}
		}
		visited[0] = true;  // 0 is the number of the first vertex we are considering (can be any from 0 to n - 1)
	}
	MPI_Bcast(&matrix[0][0], n * n, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&visited[0], n, MPI_C_BOOL, 0, MPI_COMM_WORLD);

	for (int i = 0; i < n - 1; ++i)
	{
		std::pair<int, int>* currentEdges = new std::pair<int, int>[size];
		std::pair<int, int> currentEdge = std::make_pair(-1, -1);

		int minWeight = INT_MAX;
		for (int v = 0; v < counts[rank]; ++v)
		{
			if (!visited[v + displs[rank]])
			{
				for (int j = 0; j < n; ++j)
				{
					if (visited[j] && minWeight > matrix[v + displs[rank]][j])
					{
						minWeight = matrix[v + displs[rank]][j];
						currentEdge = std::make_pair(v + displs[rank], j);
					}
				}
			}
		}
		/*
		else {
			for (int j = 0; j < size; ++j)
				currentEdges[j] = std::make_pair(-1, -1);
		}*/
		/*MPI_Barrier(MPI_COMM_WORLD);*/
		MPI_Gather(&currentEdge, sizeof(std::pair<int, int>), MPI_BYTE, currentEdges, sizeof(std::pair<int, int>), MPI_BYTE, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			std::pair<int, int> edgeWithMinWeight;
			int minWeight = INT_MAX;
			printf("current edges from procs:\n");
			for (int j = 0; j < size; ++j)
			{
				printf("(%d, %d) ", currentEdges[j].first, currentEdges[j].second);
				if (currentEdges[j].first != -1 && minWeight > matrix[currentEdges[j].first][currentEdges[j].second])
				{
					minWeight = matrix[currentEdges[j].first][currentEdges[j].second];
					edgeWithMinWeight = currentEdges[j];
				}
			}
			printf("\n");

			visited[edgeWithMinWeight.first] = true;
			tree[i] = std::make_pair(edgeWithMinWeight.first, edgeWithMinWeight.second);
		}
		delete[] currentEdges;
		MPI_Bcast(&visited[0], n, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		printf("\n");
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				printf("%d ", matrix[i][j] != INT_MAX ? matrix[i][j] : -1);
			printf("\n");
		}
		printf("\n");

		for (int i = 0; i < n - 1; ++i)
			printf("(%d, %d)\n", tree[i].first, tree[i].second);
	}

	delete[] counts;
	delete[] displs;

	MPI_Finalize();

	return EXIT_SUCCESS;
}

// delete memory !!!