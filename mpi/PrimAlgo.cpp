#include <iostream>
#include <mpi.h>

#define n 4    // |V|
#define maxWeight 10



void fillMatrixAndVisitedVector(int matrix[n][n], bool visited[n]);

unsigned long getTreeWeight(std::pair<int, int> edges[n - 1], int matrix[n][n]);

void printMatrix(int matrix[n][n]);

void distributeVertices(int** counts, int** displs, int size);

std::pair<int, int> getLocalMinWeightEdge(int* counts, int* displs, int rank,
	int matrix[n][n], bool* visited);

std::pair<int, int> getGlobalMinWeightEdge(int size, std::pair<int, int>* currentEdges,
	int matrix[n][n]);

std::pair<int, int>* getMinimumSpanningTree(int size, int rank, int* counts, int* displs,
	int matrix[n][n], bool* visited);


int main(int argc, char** argv)
{
	int rank, size;
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int matrix[n][n];
	bool visited[n];
	//std::pair<int, int> tree[n - 1];  // save edges that we are include in the 'minimum spanning tree' (MST), (our result)

	int* counts = NULL, *displs = NULL;
	distributeVertices(&counts, &displs, size);

	/*if (rank == 0)
	{
		printf("counts:\n");
		for (int i = 0; i < size; ++i)
			printf("%d ", counts[i]);
		printf("\ndispls:\n");
		for (int i = 0; i < size; ++i)
			printf("%d ", displs[i]);
		printf("\n\n");
	}*/

	if (rank == 0)
	{
		fillMatrixAndVisitedVector(matrix, visited);
	}
	MPI_Bcast(&matrix[0][0], n * n, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&visited[0], n, MPI_C_BOOL, 0, MPI_COMM_WORLD);

	std::pair<int, int>* tree = getMinimumSpanningTree(size, rank, counts, displs, matrix, visited);

	if (rank == 0)
	{
		printMatrix(matrix);
		printf("\ntree weight = %lu", getTreeWeight(tree, matrix));
	}

	delete[] counts;
	delete[] displs;

	MPI_Finalize();

	return EXIT_SUCCESS;
}



void fillMatrixAndVisitedVector(int matrix[n][n], bool visited[n])
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

void printMatrix(int matrix[n][n])
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			printf("%d ", matrix[i][j] != INT_MAX ? matrix[i][j] : -1);
		printf("\n");
	}
	printf("\n");
}

unsigned long getTreeWeight(std::pair<int, int> tree[n - 1], int matrix[n][n])
{
	unsigned long treeWeight = 0;
	for (int i = 0; i < n - 1; ++i)
	{
		treeWeight += matrix[tree[i].first][tree[i].second];
		printf("(%d, %d)\n", tree[i].first, tree[i].second);
	}

	return treeWeight;
}

void distributeVertices(int** counts, int** displs, int size)
{
	*counts = new int[size];
	*displs = new int[size];
	int rest = n;
	int localCount = rest / size;
	(*displs)[0] = 0;
	(*counts)[0] = localCount;
	for (int i = 1; i < size; ++i)
	{
		rest -= localCount;
		localCount = rest / (size - i);
		(*counts)[i] = localCount;
		(*displs)[i] = (*counts)[i - 1] + (*displs)[i - 1];
	}
}

std::pair<int, int> getLocalMinWeightEdge(int* counts, int* displs, int rank, int matrix[n][n], bool* visited)
{
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

	return currentEdge;
}

std::pair<int, int> getGlobalMinWeightEdge(int size, std::pair<int, int>* currentEdges, int matrix[n][n])
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

	return edgeWithMinWeight;
}

std::pair<int, int>* getMinimumSpanningTree(int size, int rank, int* counts, int* displs, int matrix[n][n], bool* visited)
{
	std::pair<int, int>* tree = new std::pair<int, int>[n - 1];
	for (int i = 0; i < n - 1; ++i)
	{
		if(rank == 0)
			printf("iteration %d\n", i);
		std::pair<int, int>* currentEdges = new std::pair<int, int>[size];
		std::pair<int, int> currentEdge = getLocalMinWeightEdge(counts, displs, rank, matrix, visited);

		MPI_Gather(&currentEdge, sizeof(std::pair<int, int>), MPI_BYTE, currentEdges, sizeof(std::pair<int, int>), MPI_BYTE, 0, MPI_COMM_WORLD);

		if (rank == 0)
		{
			std::pair<int, int> edgeWithMinWeight = getGlobalMinWeightEdge(size, currentEdges, matrix);
			visited[edgeWithMinWeight.first] = true;
			tree[i] = std::make_pair(edgeWithMinWeight.first, edgeWithMinWeight.second);
		}

		delete[] currentEdges;
		MPI_Bcast(&visited[0], n, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	}

	return tree;
}