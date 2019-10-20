#include<stdio.h>
#include<math.h>
#define SIZE 5
double A[SIZE][SIZE] = {1.9940  ,  1.7463  ,  1.3490  ,  1.1347  ,  1.0048,
1.7463  ,  2.5237  ,  1.3500  ,  1.8107  ,  1.3668,
1.3490  ,  1.3500  ,  1.1141  ,  0.8219  ,  0.6519,
1.1347  ,  1.8107  ,  0.8219  ,  1.7730  ,  0.9674,
1.0048  ,  1.3668  ,  0.6519  ,  0.9674  ,  0.8289};
double AnswerL[SIZE][SIZE]={1.4121 ,  1.2367 ,  0.9553 ,  0.8035 ,  0.7116,
	 0 ,  0.9972 ,  0.1691 ,  0.8193 ,  0.4882,
	 0 ,  0 ,  0.4158,   -0.2026,   -0.2657,
	 0 ,  0 ,  0 ,  0.6442,   -0.0905,
0 ,  0 ,  0 ,  0,  0.0734};

void PrintMatrix(double A[SIZE][SIZE])
{
	for(int i=0;i<SIZE;i++)
	{
		for(int j=0;j<SIZE;j++)
		{
			printf("%2.6lf\t",A[i][j]);
		}
		printf("\n");
	}
}
void PrintVector(double x[SIZE])
{
	for(int i=0;i<SIZE;i++)
		printf("%2.6lf\n",x[i]);
}

int MatrixEqual(const double A[SIZE][SIZE],const double B[SIZE][SIZE])
{
	int t=0;
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			if(A[i][j]!=B[i][j])
				t++;
	return t;
}

int VectorEqual(const double A[SIZE],const double B[SIZE])
{
	int t=0;
	for(int i=0;i<SIZE;i++)
			if(A[i]!=B[i])
				t++;
	return t;
}

double VectorNorm(const double A[SIZE])
{
	double sum=0.0;
	for(int i=0;i<SIZE;i++)
		sum+=A[i]*A[i];
	
	return sqrt(sum);
}

void SubtractVector(const double A[SIZE],const double B[SIZE],double C[SIZE])
{
	for(int i=0;i<SIZE;i++)
		C[i]=A[i]-B[i];
}

void AbsoluteMatrix(const double A[SIZE][SIZE], double B[SIZE][SIZE])
{
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			B[i][j]=fabs(A[i][j]);
}
int IsZeroMatrix(const double A[SIZE][SIZE])
{
	int t=0;
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			if(fabs(A[i][j])!=0)
				t++;
	//printf("t = %d\n",t);
	return t;
}

void MatrixMult(const double S[SIZE][SIZE], const double T[SIZE][SIZE],double C[SIZE][SIZE])
{	
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			for(int k=0;k<SIZE;k++)
				C[i][j]+=S[i][k]*T[k][j];
}


void MatrixSubtraction(const double A[SIZE][SIZE],const double B[SIZE][SIZE], double C[SIZE][SIZE])
{
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			C[i][j]=A[i][j]-B[i][j];
}


void MatrixVector(const double S[SIZE][SIZE], const double X[SIZE],double Y[SIZE])
{	
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
				Y[i]+=S[i][j]*X[j];
	PrintVector(Y);
}


double MatrixNorm(const double A[SIZE][SIZE])
{
	double sum=0.0;
	for(int i=0;i<SIZE;i++)
	{
		for(int j=0;j<SIZE;j++)
		{
			if(i==j)
				sum=sum + fabs(A[i][j]*A[i][j]);
			else
				sum+=A[i][j]*A[i][j];
		}
	}
	
	return sqrt(sum);
}
void Transpose(const double A[SIZE][SIZE],double AT[SIZE][SIZE])
{
	for(int i=0;i<SIZE;i++)
		for(int j=0;j<SIZE;j++)
			AT[j][i]=A[i][j];
}


void ForwardSubst(double A[SIZE][SIZE],const double b[SIZE],double x[SIZE])
{
	x[0]=b[0]/A[0][0];
	for(int i=1;i<SIZE;i++)
	{
		double sum=b[i];
		for(int j=0;j<=i-1;j++)
		{
			sum-=A[i][j]*x[j];
		}
		x[i]=sum/A[i][i];
	}
}




void BackwardSubst(const double A[SIZE][SIZE],const double b[SIZE],double x[SIZE])
{
	x[SIZE-1]=b[SIZE-1]/A[SIZE-1][SIZE-1];
	for(int i=SIZE-2;i>=0;i--)
	{
		double sum=b[i];
		for(int j=i+1;j<SIZE;j++)
		{
			sum-=A[i][j]*x[j];
		}
		x[i]=sum/A[i][i];
	}
}
void Cholesky(const double A[SIZE][SIZE],double L[SIZE][SIZE])
{
	for(int k=0;k<SIZE;k++)
	{
		double sum=0;
		for(int s=0;s<=k-1;s++)
		{
			sum=sum+L[k][s]*L[k][s];
		}
		L[k][k]=sqrt(A[k][k]-sum);
		for(int i=k+1;i<SIZE;i++)
		{
			double sum=0;
			for(int s=0;s<=k-1;s++)
			{
				sum+=L[i][s]*L[k][s];
			}
			L[i][k]=(A[i][k]-sum)/L[k][k];
		}	
	}	
}


int main()
{
	PrintMatrix(A);
	double L[SIZE][SIZE]={0};	
	Cholesky(A,L);
	double LT[SIZE][SIZE]={0};	
	Transpose(L,LT);
	PrintMatrix(LT);
	double Mult[SIZE][SIZE]={0};
	MatrixMult(L,LT,Mult);
	double Final[SIZE][SIZE]={0};
	MatrixSubtraction(A,Mult,Final);
	double matnorm=MatrixNorm(Final);
	printf("%2.3lf\n",matnorm);
	PrintMatrix(Final);
	int Mark=0;
	if(matnorm<0.001)
		Mark=5;
	else
		Mark=1;
	printf("You got %d Marks for this Problem\n",Mark);
	//double x[SIZE]={0};	
	return 0;
}


