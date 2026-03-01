#define KALMAN_ROLL_PITCH_H
#include <math.h>

typedef struct{

    double roll; // variables that wer are trying to filter for 
    double pitch;

    double P[4]; // 2x2 covariance variance matrix
    double Q[2]; // 2x2, but its completely diagonal
    double R[3]; // 3x3, but completely diagonal


} KRP;

void KalmanRPI(KRP *krp, double Pinit, double*Q, double *R);
void KalmanRPP(KRP *krp, double );
void KalmanRP(KRP *krp);