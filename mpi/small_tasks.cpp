#include <iostream>
#include <mpi.h>

void task1(int* argc, char*** argv);
void task2(int* argc, char*** argv);
void task3(int* argc, char*** argv);
void task4(int* argc, char*** argv);
void task5(int* argc, char*** argv);
void task6(int* argc, char*** argv);
void task7(int* argc, char*** argv);
void task8(int* argc, char*** argv);
void task9(int* argc, char*** argv);
void task10(int* argc, char*** argv);
void task11and12(int* argc, char*** argv);
void task11and12_loop(int rank, int size, MPI_Comm);
void task13(int* argc, char*** argv);

double rnd(double min, double max);

int main(int argc, char* argv[])
{
	task1(&argc, &argv);
	/*task2(&argc, &argv);
	task3(&argc, &argv);
	task4(&argc, &argv);
	task5(&argc, &argv);
	task6(&argc, &argv);
	task7(&argc, &argv);
	task8(&argc, &argv);
	task9(&argc, &argv);
	task10(&argc, &argv);
	task11and12(&argc, &argv);
	task13(&argc, &argv);*/

	return EXIT_SUCCESS;
}

void task1(int* argc, char*** argv) {
	int rank, size;
	MPI_Init(argc, argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("Hello from %d process out of %d!", rank, size);

	MPI_Finalize();
}

void task2(int* argc, char*** argv)
{
	srand(time(NULL));
	int rank, size;
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int arrayLength = 16;
	int* counts = new int[size];
	int* displs = new int[size];
	int rest = arrayLength;
	int proccesN = rest / size;
	displs[0] = 0;
	counts[0] = proccesN;
	for (int i = 1; i < size; ++i)
	{
		rest -= proccesN;
		proccesN = rest / (size - i);
		counts[i] = proccesN;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int* singleProcessArray = new int[counts[rank]];
	int singleProcessCount = counts[rank];

	int mainArray[arrayLength];
	int result = 0;

	if (rank == 0)
	{
		for (int i = 0; i < arrayLength; ++i)
			mainArray[i] = rand() % 15;
	}
	MPI_Scatterv(mainArray, counts, displs, MPI_INT, singleProcessArray, singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);

	int max = singleProcessArray[0];
	for (int i = 1; i < singleProcessCount; ++i)
		if (max < singleProcessArray[i])
			max = singleProcessArray[i];

	
	printf("process %d\n", rank);
	for (int i = 0; i < singleProcessCount; ++i)
		printf("%d ", singleProcessArray[i]);
	printf("\n");

	MPI_Reduce(&max, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf("max = %d\n", result);
	}

	delete[] singleProcessArray;
	delete[] counts;
	delete[] displs;

	MPI_Finalize();
}

void task3(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double width = 10; // square side size
	int pointsCount = (int)1e8;
	int singleProcessCount = pointsCount / size;
	int circleCount = 0;

	for (int i = 0; i < singleProcessCount; ++i)
	{
		double x = rnd(-width / 2, width / 2);
		double y = rnd(-width / 2, width / 2);
		if (x * x + y * y <= (width / 2) * (width / 2))
			++circleCount;
	}

	int resultCircleCount = 0;

	MPI_Reduce(&circleCount, &resultCircleCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf("Pi = %lf\n", 4 * (resultCircleCount / (double)(singleProcessCount * size)));
	}

	MPI_Finalize();
}

void task4(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 16;
	int array[n];
	for (int i = 0; i < n; ++i)
		array[i] = rand() % 15 - 8;

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int proccesN = n / size;
	displs[0] = 0;
	counts[0] = proccesN;
	for (int i = 1; i < size; ++i)
	{
		rest -= proccesN;
		proccesN = rest / (size - i);
		counts[i] = proccesN;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int* singleProcessArray = new int[counts[rank]];
	int singleProcessCount = counts[rank];

	MPI_Scatterv(&array[0], counts, displs, MPI_INT, singleProcessArray, singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);

	int singleProcessResult[2];
	singleProcessResult[0] = 0;
	singleProcessResult[1] = 0;

	for (int i = 0; i < singleProcessCount; ++i)
	{
		if (singleProcessArray[i] > 0)
		{
			singleProcessResult[0] += singleProcessArray[i];
			++singleProcessResult[1];
		}
	}

	int result[2];
	result[0] = 0;
	result[1] = 0;
	MPI_Reduce(singleProcessResult, result, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{	
		for (int i = 0; i < n; ++i)
			printf("%d ", array[i]);
		printf("\n%d\n", (double)result[0] / result[1]);
	}

	MPI_Finalize();

	delete[] counts;
	delete[] displs;
	delete[] singleProcessArray;
}

void task5(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 6;
	int firstVector[n];
	int secondVector[n];
	int result = 0;
	for (int i = 0; i < n; ++i)
	{
		firstVector[i] = rand() % 10;
		secondVector[i] = rand() % 10;
	}

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int proccesN = n / size;
	displs[0] = 0;
	counts[0] = proccesN;
	for (int i = 1; i < size; ++i)
	{
		rest -= proccesN;
		proccesN = rest / (size - i);
		counts[i] = proccesN;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int* singleProcessFirstVector = new int[counts[rank]];
	int* singleProcessSecondVector = new int[counts[rank]];
	int singleProcessCount = counts[rank];

	MPI_Scatterv(&firstVector[0], counts, displs, MPI_INT, singleProcessFirstVector, singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(&secondVector[0], counts, displs, MPI_INT, singleProcessSecondVector, singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);

	int singleProcessResult = 0;
	for (int i = 0; i < singleProcessCount; ++i)
		singleProcessResult += singleProcessFirstVector[i] * singleProcessSecondVector[i];

	MPI_Reduce(&singleProcessResult, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	printf("process %d: count = %d, displ = %d\n", rank, counts[rank], displs[rank]);

	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
			printf("%d ", firstVector[i]);
		printf("\n");
		for (int i = 0; i < n; ++i)
			printf("%d ", secondVector[i]);
		printf("\n");

		printf("result = %d\n", result);
	}

	MPI_Finalize();

	delete[] counts;
	delete[] displs;
	delete[] singleProcessFirstVector;
	delete[] singleProcessSecondVector;
}

void task6(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 9, m = 9;
	int matrix[n][m];
	int globalMaxMin = INT_MIN;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			matrix[i][j] = rand() % 10;

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int proccesN = n / size;
	displs[0] = 0;
	counts[0] = proccesN * m;
	for (int i = 1; i < size; ++i)
	{
		rest -= proccesN;
		proccesN = rest / (size - i);
		counts[i] = proccesN * m;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int* localMatrix = new int[counts[rank]];

	MPI_Scatterv(matrix[0], counts, displs, MPI_INT, localMatrix, counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	printf("process %d:\n", rank);
	for (int i = 0; i < counts[rank] / m; ++i)
	{
		for (int j = 0; j < m; ++j)
			printf("%d ", localMatrix[i * m + j]);
		printf("\n");
	}
	printf("\n");

	int* mins = new int[counts[rank] / m];
	for (int i = 0; i < counts[rank] / m; ++i)
	{
		int min = localMatrix[m * i];
		for (int j = 1; j < m; ++j)
		{
			if (min > localMatrix[m * i + j])
				min = localMatrix[m * i + j];
		}
		mins[i] = min;
	}
	int maxMin = mins[0];
	for (int i = 1; i < counts[rank] / m; ++i)
		if (maxMin < mins[i])
			maxMin = mins[i];

	MPI_Reduce(&maxMin, &globalMaxMin, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				printf("%d ", matrix[i][j]);
			printf("\n");
		}

		printf("maxmin = %d", globalMaxMin);
	}

	delete[] counts;
	delete[] displs;
	delete[] localMatrix;

	MPI_Finalize();
}

void task7(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 10;
	const int m = 10;

	MPI_Datatype column, columnType;
	MPI_Type_vector(n, 1, m, MPI_INT, &column);
	MPI_Type_commit(&column);
	MPI_Type_create_resized(column, 0, sizeof(int), &columnType);
	MPI_Type_commit(&columnType);

	int result = 0;
	int matrix[n][m];
	int vector[m];
	int resultVector[n];
	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				matrix[i][j] = rand() % 10;

		for (int i = 0; i < m; ++i)
		{
			vector[i] = rand() % 10;
		}

		printf("matrix:\n");
		for (int i = 0; i < n; ++i)
		{
			resultVector[i] = 0;
			for (int j = 0; j < m; ++j)
				printf("%d ", matrix[i][j]);
			printf("\n");
		}
		printf("\nvector:\n");
		for (int i = 0; i < m; ++i)
			printf("%d\n", vector[i]);
		printf("\n");
	}

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = m;
	int localCount = rest / size;
	displs[0] = 0;
	counts[0] = localCount;
	for (int i = 1; i < size; ++i)
	{
		rest -= localCount;
		localCount = rest / (size - i);
		counts[i] = localCount;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int localColumns[n][m];
	int localVector[m];
	localCount = counts[rank];

	MPI_Scatterv(&matrix[0][0], counts, displs, columnType, &localColumns[0][0], localCount, columnType, 0, MPI_COMM_WORLD);
	MPI_Scatterv(&vector[0], counts, displs, MPI_INT, &localVector[0], localCount, MPI_INT, 0, MPI_COMM_WORLD);

	int localResult[n];
	for (int i = 0; i < n; ++i)
	{
		localResult[i] = 0;
		for (int j = 0; j < localCount; ++j)
		{
			localResult[i] += localColumns[i][j] * localVector[j];
		}
	}

	MPI_Reduce(&localResult[0], &resultVector[0], n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf("result:\n");
		for (int i = 0; i < n; ++i)
			printf("%d\n", resultVector[i]);
	}


	MPI_Finalize();

	delete[] counts;
	delete[] displs;
}

void task8(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status status;
	const int n = 24;
	if (rank == 0)
	{
		int vector[n];
		int* localVector = new int[n / size];
		for (int i = 0; i < n; ++i)
		{
			vector[i] = rand() % 10;
			printf("%d ", vector[i]);
		}
		printf("\n");

		for (int i = 0; i < n / size; ++i)
			localVector[i] = vector[i];

		for (int i = 1; i < size; ++i)
		{
			MPI_Send(&vector[0] + i * (n / size), n / size, MPI_INT, i, i, MPI_COMM_WORLD);
		}

		printf("process %d:\n", rank);
		for (int i = 0; i < n / size; ++i)
			printf("%d ", localVector[i]);
		printf("\n");

		int* lastVector = new int[n];

		for (int i = 0; i < n / size; ++i)
			lastVector[i] = localVector[i];

		for (int i = 1; i < size; ++i)
			MPI_Recv(lastVector + i * (n / size), n / size, MPI_INT, i, i, MPI_COMM_WORLD, &status);

		printf("{");
		for (int i = 0; i < n; ++i)
			printf("%d ", lastVector[i]);

		printf("}\n");

		delete[] localVector;
		delete[] lastVector;
	}
	else {
		int* localVector = new int[n / size];
		MPI_Recv(&localVector[0], n / size, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);

		printf("process %d:\n", rank);
		for (int i = 0; i < n / size; ++i)
			printf("%d ", localVector[i]);
		printf("\n");

		MPI_Send(localVector, n / size, MPI_INT, 0, rank, MPI_COMM_WORLD);

		delete[] localVector;
	}

	MPI_Finalize();
}

void task9(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 14;
	int vector[n];
	for (int i = 0; i < n; ++i)
		vector[i] = rand() % 10;

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int proccesN = n / size;
	displs[0] = 0;
	counts[0] = proccesN;
	for (int i = 1; i < size; ++i)
	{
		rest -= proccesN;
		proccesN = rest / (size - i);
		counts[i] = proccesN;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int* localVector = new int[counts[rank]];
	int localCount = counts[rank];

	MPI_Scatterv(&vector[0], counts, displs, MPI_INT, localVector, localCount, MPI_INT, 0, MPI_COMM_WORLD);

	printf("proess %d\n", rank);
	for (int i = 0; i < counts[rank]; ++i)
		printf("%d ", localVector[i]);
	printf("\n");

	for (int i = 0; i < localCount / 2; ++i)
	{
		int tmp = localVector[i];
		localVector[i] = localVector[localCount - 1 - i];
		localVector[localCount - 1 - i] = tmp;
	}

	printf("proess %d\n", rank);
	for (int i = 0; i < counts[rank]; ++i)
		printf("%d ", localVector[i]);
	printf("\n");

	if (rank != 0)
	{
		MPI_Send(&localVector[0], counts[rank], MPI_INT, 0, rank, MPI_COMM_WORLD);
	}
	else {
		for (int i = 0; i < n; ++i)
			printf("%d ", vector[i]);
		printf("\n");

		for (int i = 0; i < counts[0]; ++i)
			vector[n - counts[0] + i] = localVector[i];

		MPI_Status status;
		int step = counts[0];
		for (int i = 1; i < size; ++i)
		{
			step += counts[i];
			MPI_Recv(&vector[0] + n - step, counts[i], MPI_INT, i, i, MPI_COMM_WORLD, &status);
		}

		for (int i = 0; i < n; ++i)
			printf("%d ", vector[i]);
		printf("\n");
	}

	delete[] localVector;
	delete[] displs;
	delete[] counts;

	MPI_Finalize();
}

void task10(int* argc, char*** argv)
{
	int size, rank;
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n = (int)1e5;
	MPI_Status status;

	int* buffer = new int[n + MPI_BSEND_OVERHEAD];
	int bsize = sizeof(int) * (n + MPI_BSEND_OVERHEAD);
	MPI_Buffer_attach(buffer, bsize);

	if (rank == 0)
	{
		int* array = new int[n];
		for (int i = 0; i < n; ++i)
			array[i] = i;

		double time = MPI_Wtime();

		/*MPI_Send(&array[0], n, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&array[0], n, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		printf("Send time = %lf\n", MPI_Wtime() - time);*/

		/*MPI_Rsend(&array[0], n, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&array[0], n, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		printf("Rsend time = %lf\n", MPI_Wtime() - time);*/

		
		MPI_Ssend(&array[0], n, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&array[0], n, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		printf("Ssend time = %lf\n", MPI_Wtime() - time);

		
		/*MPI_Bsend(&array[0], n, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&array[0], n, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
		printf("Bsend time = %lf\n", MPI_Wtime() - time);*/
	}
	else {
		int* recv = new int[n];

		/*MPI_Recv(&recv[0], n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&recv[0], n, MPI_INT, 0, 1, MPI_COMM_WORLD);*/

		/*MPI_Recv(&recv[0], n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Rsend(&recv[0], n, MPI_INT, 0, 1, MPI_COMM_WORLD);*/

		MPI_Recv(&recv[0], n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Ssend(&recv[0], n, MPI_INT, 0, 1, MPI_COMM_WORLD);

	/*	MPI_Recv(&recv[0], n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Bsend(&recv[0], n, MPI_INT, 0, 1, MPI_COMM_WORLD);*/
	}
	MPI_Finalize();
}

void task11and12(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	task11and12_loop(rank, size, MPI_COMM_WORLD);

	MPI_Group group, newGroup;
	MPI_Comm newComm;

	MPI_Comm_group(MPI_COMM_WORLD, &group);
	MPI_Group_incl(group, 4, new int[4] {0, 1, 2, 3}, &newGroup);
	MPI_Comm_create(MPI_COMM_WORLD, newGroup, &newComm);

	if (newComm != MPI_COMM_NULL)
	{
		MPI_Comm_size(newComm, &size);
		MPI_Comm_rank(newComm, &rank);
		task11and12_loop(rank, size, newComm);
	}

	MPI_Finalize();
}

void task11and12_loop(int rank, int size, MPI_Comm comm)
{
	MPI_Status status;
	if (rank == 0)
	{
		int message = 0;
		MPI_Send(&message, 1, MPI_INT, rank + 1, rank, comm);
		MPI_Recv(&message, 1, MPI_INT, size - 1, size - 1, comm, &status);
		++message;
		printf("for size %d message = %d\n", size, message);
	}
	else
	{
		int recv;
		MPI_Recv(&recv, 1, MPI_INT, rank - 1, rank - 1, comm, &status);
		++recv;
		MPI_Send(&recv, 1, MPI_INT, (rank + 1) % size, rank, comm);
	}
}

void task13(int* argc, char*** argv)
{
	int size, rank;
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int n = 5;

	MPI_Datatype column, columnType;
	MPI_Type_vector(n, 1, n, MPI_INT, &column);
	MPI_Type_commit(&column);
	MPI_Type_create_resized(column, 0, sizeof(int), &columnType);
	MPI_Type_commit(&columnType);

	int result = 0;
	int matrix[n][n];
	
	if (rank == 0)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				matrix[i][j] = i + j;

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				printf("%d ", matrix[i][j]);
			printf("\n");
		}
	}

	int* counts = new int[size];
	int* displs = new int[size];
	int rest = n;
	int localCount = n / size;
	displs[0] = 0;
	counts[0] = localCount;
	for (int i = 1; i < size; ++i)
	{
		rest -= localCount;
		localCount = rest / (size - i);
		counts[i] = localCount;
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	int localRows[n][n];
	int localColumns[n][n];
	localCount = counts[rank];

	MPI_Scatterv(&matrix[0][0], counts, displs, columnType, &localColumns[0][0], localCount, columnType, 0, MPI_COMM_WORLD);

	for (int i = 0; i < size; ++i)
	{
		displs[i] *= n;
		counts[i] *= n;
	}

	localCount *= n;
	MPI_Scatterv(&matrix[0][0], counts, displs, MPI_INT, &localRows[0][0], localCount, MPI_INT, 0, MPI_COMM_WORLD);

	bool symmetric = true;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < localCount / n; ++j)
			if (localColumns[i][j] != localRows[j][i])
			{
				symmetric = false;
				break;
			}
	int localResult = symmetric;

	MPI_Reduce(&localResult, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf(result == size ? "matrix is symmetric\n" : "matrix is not symmetric\n");
	}

	MPI_Finalize();
}

double rnd(double min, double max)
{
	return (double)(rand()) / RAND_MAX * (max - min) + min;
}