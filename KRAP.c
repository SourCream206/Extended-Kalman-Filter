#include <math.h>

typedef struct{
    double roll; // variables that wer are trying to filter for 
    double pitch;
    double P[4]; // 2x2 covariance variance matrix
    double Q[4]; // 2x2, 
    double R[9]; // 3x3,     
} KRP;

// In the actual code, mpu.get_gyro is never used in the actual thing so that needs to be done, Gstruct has double XAxis, YAxis, and ZAxis
// Gstruct imu_gyr 
// mpu.get_gyro(0, &imu_gyr); // comes from double GyroRange[4]={131.0,65.5,32.8,16.4};, picking 0,1,2,3 chooses sensitivity of the gyro

void KalmanRPI(KRP *krp, double Pinit, double *Q, double *R){
    krp->roll = 0.0; // phi
    krp->pitch=0.0; // theta

    krp->P[0] = Pinit; 
    krp->P[1] = 0.0;
    krp ->P[2] =  0.0; 
    krp -> P[3] = Pinit;
    
    krp->Q[0] = Q[0]; krp->Q[3] = Q[3]; krp -> Q[2] = 0.0; krp -> Q[1] = 0.0;
    
    krp->R[0] = R[0]; krp->R[4] = R[4]; krp->R[8] = R[8]; 
    krp->R[1] = 0.0;  krp->R[2] = 0.0;  krp->R[3] = 0.0;  krp->R[5] = 0.0; krp->R[6] = 0.0;   krp->R[7] = 0.0;  
}

void matMult(double *A, double *B, double *C, int rA, int cA, int cB){
    for(int i = 0; i<rA; i++){
        for(int j = 0; j<cB; j++){
            C[i*cB + j ] = 0.0;
            for(int k = 0; k<cA; k++){
                C[i*cB + j] += A[i*cA + k] * B[k*cB + j];
            }
        }
    }
}

void matTransgender(double *A, double *AT, int N, int M){
    for (int i = 0; i< N; i++){
        for (int j = 0; j<M; j++){
            AT[j*N + i] = A[i*M + j];
        }
    }
}

void matAdd(double *A, double *B, double *C, int N, int M){
    for(int i = 0; i< N*M; i++){
        C[i] = A[i] + B[i];
    }
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
    matMult(F,krp->P,FP, 2,2,2);
    matTransgender(F,Ft, 2,2);
    matMult(FP,Ft,FPFt, 2,2,2);
    matAdd(FPFt, krp->Q, krp->P, 2, 2); 
}

void matInv(double *A, double *Ainv){
    double a = A[0]; // for apples
    double b = A[1]; // for balls
    double c = A[2]; // for cats
    double d = A[3]; // for dogs
    double e = A[4]; // for elephants
    double f = A[5]; // for fish 
    double g = A[6]; // for goats
    double h = A[7]; // for hats
    double i = A[8]; // for igloos

    double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);

    double invdet = 1.0 /det;

    Ainv[0] = (e*i - f*h) * invdet;
    Ainv[1] = (c*h - b*i) * invdet;
    Ainv[2] = (b*f - c*e) * invdet;
    Ainv[3] = (f*g - d*i) * invdet;
    Ainv[4] = (a*i - c*g) * invdet;
    Ainv[5] = (c*d - a*f) * invdet;
    Ainv[6] = (d*h - e*g) * invdet;
    Ainv[7] = (b*g - a*h) * invdet;
    Ainv[8] = (a*e - b*d) * invdet;

}

void KalmanRP(KRP *krp, struct Astruct *acc){
    
    double accX = acc->XAxis;
    double accY = acc->YAxis;
    double accZ = acc->ZAxis;
    
    double roll  = krp->roll;
    double pitch = krp->pitch;

    double rollAcc = atan2(accY, accZ);
    double pitchAcc = atan2(-accX, sqrt(accY*accY + accZ*accZ));

    double g = 9.81;

    //Measurement Function 
    double h[3];
    h[0] = g*sin(pitch);
    h[1] = -g*cos(pitch)*sin(roll);
    h[2] = -g*cos(pitch)*cos(roll);
    // Measurement Jacobian 
    double H[6];
    double Ht[6];
    H[0]= 0.0;
    H[1] = -g * cos(pitch);

    H[2] = g * cos(pitch) * cos(roll);
    H[3] = -g * sin(pitch) * sin(roll);
    
    H[4] = g * cos(pitch) * sin(roll); 
    H[5] = -g * sin(pitch) * cos(roll);

    
    
    matTransgender(H, Ht, 3, 2);
    // Kalman gain = P * H' / ( HPH' + R)
    double HP[6];   // 3x2 * 2x2 = 3x2
    matMult(H, krp->P, HP, 3, 2, 2);

    double HPHt[9]; // 3x2 * 2x3 = 3x3
    matMult(HP, Ht, HPHt, 3, 2, 3);  // this is O(n^3) 🤤🤤🤤()

    double S[9];
    double Sinv[9];
    matAdd(HPHt, krp->R, S, 3, 3);
    matInv(S, Sinv);
        
    double K[6];  // 2x3
    matMult(krp->P, Ht, K, 2,2,3); // P*H'  🤤🤤🤤()
    matMult(K, Sinv, K, 2,3,3);  // PH' * S-1

    // P = (I - KH)P
        
    double KH[4];  // 2x3 * 3x2 = 2x2
    matMult(K,H,KH,2,3,2);
    KH[0] = 1- KH[0];
    KH[1] = -KH[1];
    KH[2] = -KH[2];
    KH[3] = 1 - KH[3];
    double Ptemp[4];
    matMult(KH, krp->P, Ptemp,2,2,2);
    
    for(int i=0;i<4;i++) 
        KH[i] = Ptemp[i];

    krp->roll = krp->roll + K[0] * (accX - h[0]) + K[1]*(accY - h[1]) + K[2]*(accZ - h[2]);
    krp->pitch = krp->pitch + K[3]*(accX - h[0]) + K[4]*(accY - h[1]) + K[5]*(accZ - h[2]);

}