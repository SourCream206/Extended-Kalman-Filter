#include "kalmanrollpitch.h"

// In the actual code, mpu.get_gyro is never used in the actual thing so that needs to be done, Gstruct has double XAxis, YAxis, and ZAxis
// Gstruct imu_gyr 
// mpu.get_gyro(0, &imu_gyr); // comes from double GyroRange[4]={131.0,65.5,32.8,16.4};, picking 0,1,2,3 chooses sensitivity of the gyro

void KalmanRPI(KRP *krp, double Pinit, double *Q, double *R){
    krp->roll = 0.0; // phi
    krp->pitch=0.0; // theta

    krp->P[0] = Pinit; krp->P[1] = 0.0;
    krp ->P[3] =  0.0; krp -> P[4] = Pinit;
    krp->Q[0] = Q[0]; krp->Q[3] = Q[3]; krp -> Q[2] = 0.0; krp -> Q[1] = 0.0;
    
    krp->R[0] = R[0]; krp->R[1] = R[1]; krp->R[2] = R[2]; 

}

void mat2Multily(double *A, double *B, double *C){
    C[0] = A[0]*B[0] + A[1]*B[2]; 
    C[1] = A[0]*B[1] + A[1]*B[3]; 
    C[2] = A[2]*B[0] + A[3]*B[2]; 
    C[3] = A[2]*B[1] + A[3]*B[3]; 
}

void mat2Transpose(double *A, double *AT){
    AT[0] = A[0]; AT[1] = A[2];
    AT[2] = A[1]; AT[3] = A[3];
}

void mat2Add(double *A, double *B, double *C){
    C[0] = A[0] + B[0]; C[1] = A[1] + B[1];
    C[2] = A[2] + B[2]; C[3] = A[3] + B[3];
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
        
    double F[4];
    // F = I + A*t, approximation for F without using heavy matrix operations 
    F[0] = 1.0 + A[0]*t;
    F[1] = A[1]*t;
    F[2] = A[2]*t;
    F[3] = 1.0 + A[3]*t;
    
    double FP[4];
    double Ft[4];
    double FPFt[4];
    
    // Calculate P = FPFt + Q 
    mat2Multily(F,krp->P,FP);
    mat2Transpose(F,Ft);
    mat2Multily(FP,Ft,FPFt);
    mat2Add(FPFt, krp->Q, krp->P); 
}

void KalmanRP(KRP *krp, struct Astruct *acc){
    
    double z0 = acc->XAxis;
    double z1 = acc->YAxis;
    
    double x0 = krp->roll;
    double x1 = krp->pitch;


    krp->P[0] += krp->Q[0];
    krp->P[1] += krp->Q[1];
    krp->P[2] += krp->Q[2];
    krp->P[3] += krp->Q[3];


    
}