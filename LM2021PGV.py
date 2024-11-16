# this model is developd by Chenying Liu
# using BCHydro GMM functional form for PGV for NGAsub database
# last verified Jan 26 2022
import numpy as np
from numpy import log, exp

print("The estimation of PGV for large M intraslab EQ is higher than BChydro CGMM")


def encode_reg(region):
    if region == "Alaska":
        MB_slab = 7.2
        MB_face = 8.4
        reg = 0
    elif region == "Cascadia":
        MB_slab = 6.4
        MB_face = 7.7
        reg = 1
    elif region == "CAM":
        MB_slab = 6.7
        MB_face = 7.3
        reg = 2
    elif region == "Japan":
        MB_slab = 7.6
        MB_face = 8.1
        reg = 3
    elif region == "NewZealand":
        MB_slab = 7.5
        MB_face = 8.3
        reg = 4
    elif region == "SouthAmerica":
        MB_slab = 7
        MB_face = 8.5
        reg = 5
    elif region == "Taiwan":
        MB_slab = 7.1
        MB_face = 7.1
        reg = 6
    elif region == "global":
        MB_slab = 7.6
        MB_face = 7.9
        reg = 7
    else:
        print("invalid string for region")

    return MB_slab, MB_face, reg


def LM2021PGV(M, Rrup, Vs30, ZTOR, mechanism, region='global'):
    # PGV in cm/s
    Vlin = 400
    n = 1.18
    b = -1.186
    c = 1.88
    [MB_slab, MB_face, regid] = encode_reg(region)
    if mechanism == "interface":
        Fs = 0
    if mechanism == "intraslab":
        Fs = 1
    vssmall = Vs30 <= Vlin
    vsbig = Vs30 > Vlin
    Vs30_new = min(Vs30, 1000)
    C1 = Fs * MB_slab + (1 - Fs) * MB_face
    C1slab_C1face = MB_slab - MB_face
    Msmall = M <= C1
    Mbig = 1 - Msmall
    Zsmall = ZTOR <= 100
    Zbig = 1 - Zsmall

    a1 = np.zeros(8)
    a6 = np.zeros(8)
    a12 = np.zeros(8)

    a1[0] = 9.306274418239529
    a1[1] = 9.431237996123231
    a1[2] = 9.888323824222548
    a1[3] = 10.113212187851387
    a1[4] = 9.801095367543883
    a1[5] = 9.224550504063641
    a1[6] = 10.574316427946735
    a1[7] = 9.750898391003714  # mean of a1

    a6[0] = -0.0003293625672081765
    a6[1] = -0.00033879598141665625
    a6[2] = -0.0001764763651728594
    a6[3] = -0.001490966079195552
    a6[4] = -0.00010571266915519636
    a6[5] = -0.00018205518233582244
    a6[6] = -0.0001293708359536426
    a6[7] = (a6[0] + a6[1] + a6[2] + a6[3] + a6[4] + a6[5] + a6[6]) / 7

    a12[0] = 1.3110382884146086
    a12[1] = 1.122565971648433
    a12[2] = 1.242986440727083
    a12[3] = 0.6746274887774437
    a12[4] = 1.1997586117349448
    a12[5] = 1.7077564759738402
    a12[6] = 0.8350926646710639
    a12[7] = 1.15346058170636

    a2 = -1.4232470733591636
    a3 = 0.27374638893387426
    a4 = -0.857542999517149
    a5 = -0.7670371323587445
    a9 = 0.2738160530618464
    a10 = 2.1180122269325077
    a13 = -0.14951511371661336
    a14 = -0.29870465579908617
    t11 = 0.010360299840613398
    C4 = 18.523332829838687

    tau = 0.42  # truth: 0.49006863478774887
    phi = 0.55  # truth: 0.6283312794068001

    # prediction
    intercept = a1[regid] + a4 * C1slab_C1face * Fs

    f_geom = (a2 + a14 * Fs + a3 * (M - 7.8)) * log(Rrup + C4 * exp((M - 6) * a9))

    f_attn = a6[regid] * Rrup + a10 * Fs

    f_mag = Msmall * (a4 * (M - C1) + a13 * (10 - M) ** 2) + Mbig * (a5 * (M - C1) + a13 * (10 - M) ** 2)

    f_depth = Zsmall * t11 * (ZTOR - 60) * Fs + Zbig * t11 * (100 - 60) * Fs

    siterock = (a12[regid] + b * n) * log(1100 / Vlin)
    IMrock = exp(intercept + f_mag + f_geom + f_depth + f_attn + siterock)

    f_site = (
        vssmall * a12[regid] * log(Vs30_new / Vlin)
        - b * log(IMrock + c)
        + b * log(IMrock + c * (Vs30_new / Vlin) ** n)
        + vsbig * (a12[regid] + b * n) * log(Vs30_new / Vlin)
    )

    lnpgv = intercept + f_mag + f_geom + f_depth + f_site + f_attn

    sigma = np.sqrt(tau ** 2 + phi ** 2)

    return lnpgv, tau, phi, sigma
