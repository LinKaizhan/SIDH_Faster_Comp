/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
* Copyright (c) Microsoft Corporation
*
* Website: https://github.com/microsoft/PQCrypto-SIDH
* Released under MIT license
*
* Abstract: pairing computation for compression
*********************************************************************************************/


#define t_points  2


static void final_exponentiation_2_torsion(f2elm_t f, const f2elm_t finv, f2elm_t fout)
{ // The final exponentiation for pairings in the 2^eA-torsion group. Raising the value f to the power (p^2-1)/2^eA.
    felm_t one = {0};
    f2elm_t temp = {0};
    unsigned int i; 

    fpcopy((digit_t*)&Montgomery_one, one);
    
    // f = f^p
    fp2_conj(f, temp);
    fp2mul_mont(temp, finv, temp);              // temp = f^(p-1)

    for (i = 0; i < OBOB_EXPON; i++) {
        cube_Fp2_cycl(temp, one);
    }
    fp2copy(temp, fout);
}


static void final_exponentiation_3_torsion(f2elm_t f, const f2elm_t finv, f2elm_t fout)
{ // The final exponentiation for pairings in the 3-torsion group. Raising the value f to the power (p^2-1)/3^eB.
    felm_t one = {0};
    f2elm_t temp;
    unsigned int i; 

    fpcopy((digit_t*)&Montgomery_one, one);
    
    // f = f^p
    fp2_conj(f, temp); 
    fp2mul_mont(temp, finv, temp);              // temp = f^(p-1)

    for (i = 0; i < OALICE_BITS; i++) {
        sqr_Fp2_cycl(temp, one);
    }
    fp2copy(temp, fout);
}


void Tate3_pairings(point_full_proj_t* Qj, f2elm_t* f)
{
    felm_t* x3, * x2plam12, * lam1p2, * x1x1px2plam12plam1p2y1, * l1, * c;
    f2elm_t XQ2S[t_points], finv[2 * t_points], one = { 0 };
    f2elm_t g, h, tf, t0, t1;

    fpcopy((digit_t*)&Montgomery_one, one[0]);

    for (int j = 0; j < t_points; j++)
    {
        fp2copy(one, f[j]);
        fp2copy(one, f[j + t_points]);
        fp2sqr_mont(Qj[j]->X, XQ2S[j]);
    }

    for (int k = 0; k < OBOB_EXPON - 1; k++) {
        x3 = (felm_t*)T3_0 + k;
        x2plam12 = (felm_t*)T3_1 + k;
        lam1p2 = (felm_t*)T3_2 + k;
        x1x1px2plam12plam1p2y1 = (felm_t*)T3_3 + k;
        for (int j = 0; j < t_points; j++) {
            fpmul_mont(Qj[j]->X[0], *x2plam12, t0[0]);
            fpmul_mont(Qj[j]->X[1], *x2plam12, t0[1]);
            fp2add(XQ2S[j], t0, g);
            fpmul_mont(Qj[j]->Y[0], *lam1p2, t1[0]);
            fpmul_mont(Qj[j]->Y[1], *lam1p2, t1[1]);
            fp2sub(g, t1, g);
            fpadd(g[0], *x1x1px2plam12plam1p2y1, g[0]);
            fpsub(Qj[j]->X[0], *x3, h[0]);
            fpcopy(Qj[j]->X[1], h[1]);
            fp2_conj(h, h);
            fp2mul_mont(g, h, g);
            fp2sqr_mont(f[j], tf);
            fp2mul_mont(f[j], tf, f[j]);
            fp2mul_mont(f[j], g, f[j]);

            fp2sub(XQ2S[j], t0, g);
            fpadd(g[1], t1[0], g[1]);
            fpsub(g[0], t1[1], g[0]);
            fpadd(g[0], *x1x1px2plam12plam1p2y1, g[0]);
            fpadd(Qj[j]->X[0], *x3, h[0]);
            fp2mul_mont(g, h, g);
            fp2sqr_mont(f[j + t_points], tf);
            fp2mul_mont(f[j + t_points], tf, f[j + t_points]);
            fp2mul_mont(f[j + t_points], g, f[j + t_points]);
        }
    }
    for (int j = 0; j < t_points; j++) {
        l1 = (felm_t*)T3_4 + 0;
        c = (felm_t*)T3_4 + 1;

        fpmul_mont(Qj[j]->X[0], *l1, t0[0]);
        fpmul_mont(Qj[j]->X[1], *l1, t0[1]);
        fp2sub(t0, Qj[j]->Y, g);
        fpsub(g[0], *c, g[0]);

        fp2sqr_mont(f[j], tf);
        fp2mul_mont(f[j], tf, f[j]);
        fp2mul_mont(f[j], g, f[j]);

        fpsub(t0[1], Qj[j]->Y[0], g[0]);
        fpadd(t0[0], Qj[j]->Y[1], g[1]);
        fpneg(g[1]);
        fpsub(g[1], *c, g[1]);

        fp2sqr_mont(f[j + t_points], tf);
        fp2mul_mont(f[j + t_points], tf, f[j + t_points]);
        fp2mul_mont(f[j + t_points], g, f[j + t_points]);
    }

    // Final exponentiation:
    mont_n_way_inv(f, 2 * t_points, finv);
    for (int j = 0; j < 2 * t_points; j++) {
        final_exponentiation_3_torsion(f[j], finv[j], f[j]);
    }
}


void Tate2_pairings(point_full_proj_t* Qj, f2elm_t* f)
{
    felm_t* lam1, * lam2, * lam1x1y1, * lam2x2y2;
    f2elm_t finv[2 * t_points], one = { 0 };
    f2elm_t l1_first, g, h;
    felm_t t0, t1;

    fpcopy((digit_t*)&Montgomery_one, one[0]);

    for (int j = 0; j < t_points; j++) {
        fp2copy(one, f[j]);
        fp2copy(one, f[j + t_points]);
    }
    // Pairings with P
    fpcopy(((felm_t*)T4_6)[0], t0);
    fpcopy(((felm_t*)T4_6)[1], t1);
    fpcopy(((felm_t*)T4_0)[2], l1_first[0]);
    fpcopy(((felm_t*)T4_0)[3], l1_first[1]);

    for (int j = 0; j < t_points; j++) {
        fpsub(Qj[j]->X[0], t0, h[0]);
        fpcopy(Qj[j]->X[1], h[1]);
        fp2copy(h, g);
        fp2mul_mont(g, l1_first, g);
        fp2sub(g, Qj[j]->Y, g);
        fpsub(g[1], t1, g[1]);
        fp2_conj(h, h);
        fp2mul_mont(g, h, g);
        fp2copy(g, f[j]);
    }

    for (int k = 0; k < OALICE_BITS / 2 - 1; k++) {
        lam1 = (felm_t*)T4_7 + k;
        lam1x1y1 = (felm_t*)T4_8 + k;
        lam2 = (felm_t*)T4_9 + k;
        lam2x2y2 = (felm_t*)T4_10 + k;

        for (int j = 0; j < t_points; j++) {
            fpmul_mont(*lam1, Qj[j]->X[0], g[1]);
            fpmul_mont(*lam1, Qj[j]->X[1], g[0]);
            fpneg(g[0]);
            fp2sub(g, Qj[j]->Y, g);
            fpadd(g[1], *lam1x1y1, g[1]);
            fpmul_mont(*lam2, Qj[j]->X[0], h[1]);
            fpmul_mont(*lam2, Qj[j]->X[1], h[0]);
            fpneg(h[0]);
            fp2sub(h, Qj[j]->Y, h);
            fpadd(h[1], *lam2x2y2, h[1]);
            fp2_conj(h, h);

            fp2sqr_mont(f[j], f[j]);
            fp2mul_mont(g, f[j], f[j]);
            fp2sqr_mont(f[j], f[j]);
            fp2mul_mont(f[j], h, f[j]);
        }
    }
    // Last iteration
#if (OALICE_BITS % 2 == 0)
    for (int j = 0; j < t_points; j++) {
        fpcopy(((felm_t*)T4_11)[1], t0);
        fpsub(Qj[j]->X[0], t0, g[0]);
        fpcopy(Qj[j]->X[1], g[1]);
        fp2sqr_mont(f[j], f[j]);
        fp2mul_mont(f[j], g, f[j]);
    }
#else
    for (int j = 0; j < t_points; j++) {
        fpcopy(Qj[j]->X[1], h[1]);
        fpcopy(((felm_t*)T4_11)[2], t0);
        fpsub(Qj[j]->X[0], t0, h[0]);
        fpcopy(((felm_t*)T4_11)[3], t1);
        fpmul_mont(h[0], t1, g[1]);
        fpmul_mont(h[1], t1, g[0]);
        fpneg(g[0]);
        fp2sub(g, Qj[j]->Y, g);
        fp2_conj(h, h);

        fp2sqr_mont(f[j], f[j]);
        fp2mul_mont(g, f[j], f[j]);
        fp2sqr_mont(f[j], f[j]);
        fp2mul_mont(f[j], h, f[j]);
    }
#endif  

    // Pairings with Q
    fpcopy(((felm_t*)T4_1)[0], t0);
    fpcopy(((felm_t*)T4_1)[1], t1);
    fpcopy(((felm_t*)T4_0)[0], l1_first[0]);
    fpcopy(((felm_t*)T4_0)[1], l1_first[1]);

    for (int j = 0; j < t_points; j++) {
        fpsub(Qj[j]->X[0], t0, h[0]);
        fpcopy(Qj[j]->X[1], h[1]);
        fp2copy(h, g);
        fp2mul_mont(g, l1_first, g);
        fp2sub(g, Qj[j]->Y, g);
        fpsub(g[0], t1, g[0]);
        fp2_conj(h, h);
        fp2mul_mont(g, h, g);
        fp2copy(g, f[j + t_points]);
    }
    for (int k = 0; k < OALICE_BITS / 2 - 1; k++) {
        lam1 = (felm_t*)T4_2 + k;
        lam1x1y1 = (felm_t*)T4_3 + k;
        lam2 = (felm_t*)T4_4 + k;
        lam2x2y2 = (felm_t*)T4_5 + k;

        for (int j = 0; j < t_points; j++) {
            fpmul_mont(*lam1, Qj[j]->X[0], g[0]);
            fpmul_mont(*lam1, Qj[j]->X[1], g[1]);
            fp2sub(g, Qj[j]->Y, g);
            fpadd(g[0], *lam1x1y1, g[0]);
            fpmul_mont(*lam2, Qj[j]->X[0], h[0]);
            fpmul_mont(*lam2, Qj[j]->X[1], h[1]);
            fp2sub(h, Qj[j]->Y, h);
            fpadd(h[0], *lam2x2y2, h[0]);
            fp2_conj(h, h);

            fp2sqr_mont(f[j + t_points], f[j + t_points]);
            fp2mul_mont(g, f[j + t_points], f[j + t_points]);
            fp2sqr_mont(f[j + t_points], f[j + t_points]);
            fp2mul_mont(f[j + t_points], h, f[j + t_points]);
        }
    }
    // Last iteration
#if (OALICE_BITS % 2 == 0)
    for (int j = 0; j < t_points; j++) {
        fpcopy(((felm_t*)T4_11)[0], t0);
        fpsub(Qj[j]->X[0], t0, g[0]);
        fpcopy(Qj[j]->X[1], g[1]);
        fp2sqr_mont(f[j + t_points], f[j + t_points]);
        fp2mul_mont(f[j + t_points], g, f[j + t_points]);
    }
#else
    for (int j = 0; j < t_points; j++) {
        fpcopy(Qj[j]->X[1], h[1]);
        fpcopy(((felm_t*)T4_11)[0], t0);
        fpsub(Qj[j]->X[0], t0, h[0]);
        fpcopy(((felm_t*)T4_11)[1], t1);
        fpmul_mont(h[0], t1, g[0]);
        fpmul_mont(h[1], t1, g[1]);
        fp2sub(g, Qj[j]->Y, g);
        fp2_conj(h, h);
        fp2sqr_mont(f[j + t_points], f[j + t_points]);
        fp2mul_mont(g, f[j + t_points], f[j + t_points]);
        fp2sqr_mont(f[j + t_points], f[j + t_points]);
        fp2mul_mont(f[j + t_points], h, f[j + t_points]);;
    }
#endif  

    // Final exponentiation:
    mont_n_way_inv(f, 2 * t_points, finv);
    for (int j = 0; j < 2 * t_points; j++) {
        final_exponentiation_2_torsion(f[j], finv[j], f[j]);
    }
    }

