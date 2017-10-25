void addVecs( const int& n, double X[], double Y[], double Z[] ) {
    for( int i=0; i<n; i++ ) {
        Z[i] = X[i] + Y[i];
    }
}
