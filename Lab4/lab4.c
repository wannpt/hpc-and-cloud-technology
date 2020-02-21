#include<mpi.h>
#include<stdlib.h>
#include<stdio.h>

int main(int argc,char* argv[])
{
	int p;
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int row_mat1, col_mat1, row_mat2, col_mat2;
	int i, j, k;
	double *mat1, *mat2, *resMat, *mat1Buf, *resMatBuf;
    char line[100];
	
	if (rank == 0)
	{
		FILE* inputFile1 = fopen(argv[1], "r");
		fscanf(inputFile1, "%d %d", &row_mat1, &col_mat1);
		mat1 = (double*)malloc(row_mat1 * col_mat1 * sizeof(double));
		while (fgets(line,sizeof(line),inputFile1)!=NULL)
		{
			for (i = 0; i < row_mat1; i++)
			{
				for (j = 0; j < col_mat1; j++)
				{
					fscanf(inputFile1, "%lf", &mat1[(i * col_mat1) + j]);
				}
			}
		}
		fclose(inputFile1);
		FILE* inputFile2 = fopen(argv[2], "r");
		fscanf(inputFile2, "%d %d", &row_mat2, &col_mat2);
		mat2 = (double*)malloc(row_mat2 * col_mat2 * sizeof(double));

		while (fgets(line,sizeof(line),inputFile2)!=NULL)
		{
			for (i = 0; i < row_mat2; i++)
			{
				for (j = 0; j < col_mat2; j++)
				{
					fscanf(inputFile2, "%lf", &mat2[(j * row_mat2) + i]);
				}
			}
		}

        int temp = row_mat2;
        row_mat2 = col_mat2;
        col_mat2 = temp;

		fclose(inputFile2);
		resMat = (double*)malloc(row_mat1 * row_mat2 * sizeof(double));
	}

	MPI_Bcast(&row_mat1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&col_mat1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&row_mat2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank != 0)
	{
		mat2 = (double*)malloc(col_mat1 * row_mat2 * sizeof(double));
	}

	mat1Buf = (double*)malloc(row_mat1 / p * col_mat1 * sizeof(double));
	resMatBuf = (double*)malloc(row_mat1 / p * row_mat2 * sizeof(double));
	
	MPI_Bcast(mat2, col_mat1 * row_mat2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(mat1, row_mat1 / p * col_mat1, MPI_DOUBLE, mat1Buf, row_mat1 / p * col_mat1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < row_mat1 / p; i++)
	{
		for (j = 0; j < row_mat2; j++)
		{
			resMatBuf[i * row_mat2 + j] = 0;

			for (k = 0; k < col_mat1; k++)
			{
				resMatBuf[i * row_mat2 + j] += mat1Buf[i * col_mat1 + k] * mat2[j * col_mat1 + k];
			}
		}
	}
	MPI_Gather(resMatBuf, row_mat1 / p * row_mat2, MPI_DOUBLE, resMat, row_mat1 / p * row_mat2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        if(row_mat1%p!=0)
        {
            for (i = row_mat1/p * p; i < (row_mat1/p *p)+row_mat1%p; i++)
	        {
                for (j = 0; j < row_mat2; j++)
                {
                    resMat[i * row_mat2 + j] = 0;

                    for (k = 0; k < col_mat1; k++)
                    {
                        resMat[i * row_mat2 + j] += mat1[i * col_mat1 + k] * mat2[j * col_mat1 + k];
                    }
                }
            }
        }

        FILE *outputFile = fopen(argv[3],"w+");
        fprintf(outputFile,"%d %d\n",row_mat1,row_mat2);
        for (i = 0; i < row_mat1; i++)
        {
            for (j = 0; j < row_mat2; j++)
            {
                fprintf(outputFile,"%.10lf ", resMat[i * row_mat2 + j]);
            }
            fprintf(outputFile,"\n");
        }

        free(mat1);
        free(resMat);
    }
        free(mat2);
        free(mat1Buf);
        free(resMatBuf);

    MPI_Finalize();
}