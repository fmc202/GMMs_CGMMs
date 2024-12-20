import numpy as np

# code verified on Jan 9th, 2022
# no bugs found


def TM_crustal(M):
    # find the Tpgv of PSA(Tpgv) which is magnitude dependent
    Tpgv = np.exp(-4.09 + 0.66 * M)

    return Tpgv


def TM_sub(M, mechanism):
    # find the Tpgv of PSA(Tpgv) which is magnitude dependent
    # mechanism: interface =0; intraslab = 1
    F = mechanism
    coef = np.array([[-4.06, 0.57], [-2.32, 0.28]])
    Tpgv = np.exp(coef[F, 0] + coef[F, 1] * M)

    return Tpgv


def con_pgv_crustal(M, Rrup, Vs30, lnSa_Tpgv, tra_sigma=None):
    # PGV output in cm/s
    # Sa_Tpgv in g
    # tra_sigma: sigma of lnSa_Tpgv
    a1 = 5.39
    a2 = 0.799
    a3 = 0.654
    a4 = 0.479
    a5 = -0.062
    a6 = -0.359
    a7 = -0.134
    a8 = 0.023
    tau = 0.16
    phi = 0.29
    sigma_con = 0.33

    # median model
    f_mag = (M < 5) * a2 + (M >= 5) * (M <= 7.5) * (a2 + (a3 - a2) * (M - 5) / 2.5) + (M > 7.5) * a3
    lnPGV = a1 + f_mag * lnSa_Tpgv + a4 * (M - 6) + a5 * (8.5 - M) ** 2 + a6 * np.log(Rrup + 5 * np.exp(0.4 * (M - 6))) + (a7 + a8 * (M - 5)) * np.log(Vs30 / 425)

    # propagation of error
    if tra_sigma is not None:
        sigma = np.sqrt(sigma_con ** 2 + (f_mag ** 2) * (tra_sigma ** 2))
    else:
        sigma = None

    return lnPGV, sigma



def con_pgv_sub(M, Rrup, Vs30, lnSa_Tpgv, mechanism, tra_sigma=None):

    F = mechanism
    a = np.array([[7.12, 0.633, 0.674, 0.66, -0.038, -0.51, 0.018, 0.093, 0.22, 0.33], [6.12, 0.76, 0.71, 0.77, 0.017, -0.57, -0.12, -0.01, 0.21, 0.31]])
    sigma_con = np.sqrt(a[F, -2] ** 2 + a[F, -1] ** 2)

    # median model
    f_mag = (M < 5) * a[F, 1] + (M >= 5) * (M <= 7.5) * (a[F, 1] + (a[F, 2] - a[F, 1]) * (M - 5.5) / 2.5) + (M > 7.5) * a[F, 2]

    M1 = np.array([8.4, 7.2])
    M2 = 10
    Vs0 = np.array([450, 460])

    lnPGV = (
        a[F, 0]
        + f_mag * lnSa_Tpgv
        + a[F, 3] * (M - M1[F])
        + a[F, 4] * (M2 - M) ** 2
        + a[F, 5] * np.log(Rrup + 5 * np.exp(0.4 * (M - 6)))
        + (a[F, 6] + a[F, 7] * (M - M1[F])) * np.log(Vs30 / Vs0[F])
    )

    # propagation of errors
    if tra_sigma is not None:
        sigma = np.sqrt(sigma_con ** 2 + (f_mag ** 2) * (tra_sigma ** 2))
    else:
        sigma = None

    return lnPGV, sigma
