#include <iostream>
#include <cmath>

int main()
{
	double a = -1.;
	double b = 1.;
	int np = 4;
	int ne = 4;

	//Element boundary points:
	double dx = (b-a) / ne;
	double x[ne+1];
	for( int i=0; i<ne+1; i++ )
	{
		x[i] = a + i*dx;
	}

	//Total number of nodes (element boundary nodes are repeated):
	int N = np*ne;

	//GLL nodes on standard interval from -1 to 1:
	double xGLL[np];
	if( np == 2 )
	{
		xGLL[0] = -1.;
		xGLL[1] = 1.;
	}
	else if( np == 3 )
	{
		xGLL[0] = -1.;
		xGLL[1] = 0.;
		xGLL[2] = 1.;
	}
	else if( np == 4 ) 
	{
		xGLL[0] = -1.;
		xGLL[1] = -sqrt(5.)/5.;
		xGLL[2] = sqrt(5.)/5;
		xGLL[3] = 1.;
	}

	//Location of the center of each element:
	double xc[ne];
	for( int i=0; i<ne; i++ )
	{
		xc[i] = ( x[i] + x[i+1] ) / 2;
	}

	//All of the x-coordinates in one long array:
	double X[N];
	for( int i=0; i<ne; i++ )
	{
		for( int j=0; j<np; j++ )
		{
			X[np*i+j] = xc[i] + dx/2 * xGLL[j];
		}
	}
	//Print X to make sure it is correct:
	for( int i=0; i<N; i++ )
	{
		std::cout << "X[" << i << "] = " << X[i] << std::endl;
	}
	

}
