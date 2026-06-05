#include <math.h>
typedef struct {
    double x[3];  
    double P[9];  
    double Q[9];  
    double R;     
} KVel;

void KVelInit(KVel *kv, double Pinit,
              double q_vx, double q_vy, double q_vz,
              double r_vz) {
    kv->x[0] = 0.0;
    kv->x[1] = 0.0;
    kv->x[2] = 0.0;

    // P = Pinit * I
    for (int i = 0; i < 9; i++) kv->P[i] = 0.0;
    kv->P[0] = Pinit;
    kv->P[4] = Pinit;
    kv->P[8] = Pinit;

    // Q diagonal
    for (int i = 0; i < 9; i++) kv->Q[i] = 0.0;
    kv->Q[0] = q_vx;
    kv->Q[4] = q_vy;
    kv->Q[8] = q_vz;

    kv->R = r_vz;
}
void KVelPredict(KVel *kv, KRP *krp, struct Astruct *acc, double dt) {
    double cr = cos(krp->roll),  sr = sin(krp->roll);
    double cp = cos(krp->pitch), sp = sin(krp->pitch);

    double ax_w =  cp          * acc->XAxis
                 + sp*sr       * acc->YAxis
                 + sp*cr       * acc->ZAxis;

    double ay_w =  cr          * acc->YAxis
                 - sr          * acc->ZAxis;

    double az_w = -sp          * acc->XAxis
                 + cp*sr       * acc->YAxis
                 + cp*cr       * acc->ZAxis
                 - 9.81;        // subtract gravity

    kv->x[0] += ax_w * dt;
    kv->x[1] += ay_w * dt;
    kv->x[2] += az_w * dt;

    // F = I  =>  P = P + Q
    for (int i = 0; i < 9; i++)
        kv->P[i] += kv->Q[i];
}

void KVelUpdate(KVel *kv, double ka_vz) {

    double S    = kv->P[8] + kv->R;
    double Sinv = 1.0 / S;

    double K0 = kv->P[2] * Sinv;   // gain for vx
    double K1 = kv->P[5] * Sinv;   // gain for vy
    double K2 = kv->P[8] * Sinv;   // gain for vz

    double y = ka_vz - kv->x[2];   // innovation

    kv->x[0] += K0 * y;
    kv->x[1] += K1 * y;
    kv->x[2] += K2 * y;

    double p00 = kv->P[0], p01 = kv->P[1], p02 = kv->P[2];
    double p10 = kv->P[3], p11 = kv->P[4], p12 = kv->P[5];
    double p20 = kv->P[6], p21 = kv->P[7], p22 = kv->P[8];

    kv->P[0] = p00 - K0 * p20;
    kv->P[1] = p01 - K0 * p21;
    kv->P[2] = p02 - K0 * p22;

    kv->P[3] = p10 - K1 * p20;
    kv->P[4] = p11 - K1 * p21;
    kv->P[5] = p12 - K1 * p22;

    kv->P[6] = p20 - K2 * p20;
    kv->P[7] = p21 - K2 * p21;
    kv->P[8] = p22 - K2 * p22;
}

// ----- Usage sketch (call ts from main loop btw) --------------------------------
//
// KVel kvel;
// KVelInit(&kvel, 1.0, 0.1, 0.1, 0.05, 0.2);
//
// // each loop iteration:
// KVelPredict(&kvel, &krp, &imu_acc, dt);
// KAltPredict(&kalt, &krp, &imu_acc, dt);
// KAltUpdate (&kalt, pressure_pa, temp_k, p0_pa);
// KVelUpdate (&kvel, kalt.x[1]);          // wire KAlt vz into velocity filter
//
// // read out:
// double vx = kvel.x[0];
// double vy = kvel.x[1];
// double vz = kvel.x[2];
//
//
// ------------------------------------------------------------------------------

