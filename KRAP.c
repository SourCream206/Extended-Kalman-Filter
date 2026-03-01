#include "kalmanrollpitch.h"

// In the actual code, mpu.get_gyro is never used in the actual thing so that needs to be done, Gstruct has double XAxis, YAxis, and ZAxis
// Gstruct imu_gyr 
// mpu.get_gyro(0, &imu_gyr); // comes from double GyroRange[4]={131.0,65.5,32.8,16.4};, picking 0,1,2,3 chooses sensitivity of the gyro
// 
void KalmanRPI(KRP *krp, double Pinit, double *Q, double *R){
    krp->roll = 0.0; // phi
    krp->pitch=0.0; // theta

    krp->P[0] = Pinit; krp->P[1] = 0.0;
    krp ->P[3] =  0.0; krp -> P[4] = Pinit;
    krp->Q[0] = Q[0]; krp->Q[1] = Q[1]; 
    
    krp->R[0] = R[0]; krp->R[1] = R[1]; krp->R[2] = R[2]; 

}

void KalmanRPP(KRP *krp, struct Gstruct *gyro, double t){

    double roll  = krp->roll;
    double pitch = krp->pitch;

    double p = gyro->XAxis;
    double q = gyro->YAxis;
    double r = gyro->ZAxis;

    krp->roll  = roll  + t*(p + tan(pitch)*(q*sin(roll) + r*cos(roll)));
    krp->pitch = pitch + t*(q*cos(roll) - r*sin(roll));

    double A[4];
    // yakubian matrix
    A[0] = tan(pitch)*(q*cos(roll) - r*sin(roll));
    A[1] = (r*cos(roll) + q*sin(roll)) / (cos(pitch)*cos(pitch));
    A[2] = -(r*cos(roll) + q*sin(roll));
    A[3] = 0.0;


}