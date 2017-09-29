#ifndef RBF_HPP
#define RBF_HPP

namespace interpolation
{
	class rbf
	{
		private:
			int r;
			int param;

		public:
			rbf();
			rbf( char );
			rbf( char, double );
			rbf( char, double, double )
			rbf d(rbf,int);
	};
}

#endif