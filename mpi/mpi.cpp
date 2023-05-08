#include <iostream>
#include <mpi.h>

void task1(int* argc, char*** argv);
void task2(int* argc, char*** argv);
void task2_with_reduction(int* argc, char*** argv);
void task3(int* argc, char*** argv);
void task4(int* argc, char*** argv);
void task5(int* argc, char*** argv);
void task6(int* argc, char*** argv);
void task7(int* argc, char*** argv);
void task8(int* argc, char*** argv);
void task9(int* argc, char*** argv);
void task10(int* argc, char*** argv);
void task11(int* argc, char*** argv);
void task11_loop(int rank, int size, MPI_Comm);

int sign(int number);
double rnd(double min, double max);

int main(int argc, char* argv[])
{
	/*task1(&argc, &argv);*/
	//task2(&argc, &argv);
	/*task2_with_reduction(&argc, &argv);*/
	/*task3(&argc, &argv);*/
	/*task4(&argc, &argv);*/
	/*task5(&argc, &argv);*/
	/*task6(&argc, &argv);*/
	/*task7(&argc, &argv);*/
	/*task8(&argc, &argv);*/
	/*task9(&argc, &argv);*/
	/*task10(&argc, &argv);*/
	task11(&argc, &argv);

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

	const int arrayLength = 25;
	int singleProcessCount = arrayLength / size + sign(arrayLength % size);
	int mainArray[arrayLength];
	int* array = new int[singleProcessCount];
	int* result = new int[size];

	for (int i = 0; i < arrayLength; ++i)
		mainArray[i] = rand() % 15;

	MPI_Scatter(&mainArray[0], singleProcessCount, MPI_INT, &array[0], singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);
	
	int max = array[0];
	for (int i = 1; i < singleProcessCount; ++i)
		if (max < array[i])
			max = array[i];

	std::cout << "process " << rank << std::endl;
	for (int i = 0; i < singleProcessCount; ++i)
		std::cout << array[i] << " ";
	std::cout << std::endl;

	MPI_Gather(&max, 1, MPI_INT, &result[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (int i = 0; i < arrayLength; ++i)
			std::cout << mainArray[i] << " ";
		std::cout << std::endl;
		
		for (int i = 0; i < size; ++i)
			std::cout << result[i] << " ";
		std::cout << std::endl;

		int mainMax = result[0];
		for (int i = 1; i < size; ++i)
			if (mainMax < result[i])
				mainMax = result[i];

		std::cout << "max = " << mainMax;
	}

	delete[] result;
	delete[] array;

	MPI_Finalize();
}

void task2_with_reduction(int* argc, char*** argv)
{
	srand(time(NULL));
	int rank, size;
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int arrayLength = 25;
	int singleProcessCount = arrayLength / size + sign(arrayLength % size);
	int mainArray[arrayLength];
	int* array = new int[singleProcessCount];
	int result = 0;

	for (int i = 0; i < arrayLength; ++i)
		mainArray[i] = rand() % 15;

	MPI_Scatter(&mainArray[0], singleProcessCount, MPI_INT, &array[0], singleProcessCount, MPI_INT, 0, MPI_COMM_WORLD);

	int max = array[0];
	for (int i = 1; i < singleProcessCount; ++i)
		if (max < array[i])
			max = array[i];

	std::cout << "process " << rank << std::endl;
	for (int i = 0; i < singleProcessCount; ++i)
		std::cout << array[i] << " ";
	std::cout << std::endl;

	MPI_Reduce(&max, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		std::cout << "max = " << result << std::endl;
	}

	delete[] array;

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
		std::cout << "Pi = " << 4 * (resultCircleCount / (double)(singleProcessCount * size)) << std::endl;
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
			std::cout << array[i] << " ";
		std::cout << std::endl;

		std::cout << (double)result[0] / result[1];
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
			std::cout << firstVector[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < n; ++i)
			std::cout << secondVector[i] << " ";
		std::cout << std::endl;

		std::cout << "result = " << result << std::endl;
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

	const int n = 18, m = 5;
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
			printf("%d ", localMatrix[i + j]);
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

	MPI_Status status;

	const int n = 10;
	const int m = 10;

	MPI_Datatype ColumnType;
	MPI_Type_vector(n, 1, n, MPI_INT, &ColumnType);
	MPI_Type_commit(&ColumnType);

	if (rank == 0)
	{
		int matrix[n][m];
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				matrix[i][j] = rand() % 10;

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				printf("%d ", matrix[i][j]);
			printf("\n");
		}

		MPI_Send(&matrix[0][0], 2, ColumnType, 3, 99, MPI_COMM_WORLD);
	}
	else if (rank == 3)
	{
		int localMatrix[n][m];
		MPI_Recv(&localMatrix[0][0], 2, ColumnType, 0, 99, MPI_COMM_WORLD, &status);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				printf("%d ", localMatrix[i][j]);
			printf("\n");
		}
	}

	MPI_Finalize();
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
	}
	else {
		int* localVector = new int[n / size];
		MPI_Recv(&localVector[0], n / size, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);

		printf("process %d:\n", rank);
		for (int i = 0; i < n / size; ++i)
			printf("%d ", localVector[i]);
		printf("\n");

		MPI_Send(localVector, n / size, MPI_INT, 0, rank, MPI_COMM_WORLD);
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

	MPI_Finalize();
}

void task10(int* argc, char*** argv)
{

}

void task11(int* argc, char*** argv)
{
	int size, rank;
	srand(time(NULL));
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	task11_loop(rank, size, MPI_COMM_WORLD);

	MPI_Group group, newGroup;
	MPI_Comm newComm;

	MPI_Comm_group(MPI_COMM_WORLD, &group);
	MPI_Group_incl(group, 4, new int[4] {0, 1, 2, 3}, &newGroup);
	MPI_Comm_create(MPI_COMM_WORLD, newGroup, &newComm);

	if (newComm != MPI_COMM_NULL)
	{
		MPI_Comm_size(newComm, &size);
		MPI_Comm_rank(newComm, &rank);
		task11_loop(rank, size, newComm);
	}

	MPI_Finalize();
}

void task11_loop(int rank, int size, MPI_Comm comm)
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

int sign(int number)
{
	return number > 0
		? 1
		: number == 0 
			? 0
			: -1;
}

double rnd(double min, double max)
{
	return (double)(rand()) / RAND_MAX * (max - min) + min;
}