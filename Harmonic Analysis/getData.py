from math import sin, cos, atan, sqrt
from ConfigParser import ConfigParser
from datetime import datetime, timedelta
from dateutil import parser
from numpy import array
from re import findall


def get_sigma_V(raw_data, sigma_const, V0_const, t):

    Doodson = [int(num) for num in findall('[-+]?\d+', raw_data)]
    sigma = Doodson[0] * sigma_const[0] + \
            Doodson[1] * sigma_const[1] + \
            Doodson[2] * sigma_const[2] + \
            Doodson[3] * sigma_const[3] - \
            Doodson[4] * sigma_const[4] + \
            Doodson[5] * sigma_const[5]
    V0 = Doodson[0] * V0_const[0] + \
         Doodson[1] * V0_const[1] + \
         Doodson[2] * V0_const[2] + \
         Doodson[3] * V0_const[3] - \
         Doodson[4] * V0_const[4] + \
         Doodson[5] * V0_const[5] + \
         Doodson[6] * 90.
    return sigma, V0 + sigma * t


def get_f(raw_data, N, p):
    conv = 3.1415926535897932 / 180.
    N *= conv
    p *= conv
    m_m = 1.0000 - 0.1300 * cos(N) + 0.0013 * cos(2*N)
    m_f = 1.0429 + 0.4135 * cos(N) - 0.0040 * cos(2*N)
    o_1 = 1.0089 + 0.1871 * cos(N) - 0.0147 * cos(2*N) + 0.0014 * cos(3*N)
    k_1 = 1.0060 + 0.1150 * cos(N) - 0.0088 * cos(2*N) + 0.0006 * cos(3*N)
    j_1 = 1.0129 + 0.1676 * cos(N) - 0.0170 * cos(2*N) + 0.0016 * cos(3*N)
    oo_1 = 1.1027 + 0.6504 * cos(N) + 0.0317 * cos(2*N) - 0.0014 * cos(3*N)
    m_2 = 1.0004 - 0.0373 * cos(N) + 0.0003 * cos(2*N)
    k_2 = 1.0241 + 0.2863 * cos(N) + 0.0083 * cos(2*N) - 0.0015 * cos(3*N)
    m_3 = -0.5 + 1.5 * m_2
    m_1 = sqrt((2*cos(p)+0.4*cos(p-N))**2+(sin(p)+0.2*sin(p-N))**2)
    l_2 = sqrt((1.0000-0.2505*cos(2*p)-0.1103*cos(2*p-N)-0.0156*cos(2*p-2*N)-0.0366*cos(N)+0.0047*cos(2*p+N)) ** 2 +
               (-0.2505*sin(2*p)-0.1103*sin(2*p-N)-0.0156*sin(2*p-2*N)-0.0366*sin(N)+0.0047*sin(2*p+N)) ** 2)
    p_1 = 1.
    return eval(raw_data)


def get_u(raw_data, N, p):
    pi = 3.1415926535897932
    conv = pi / 180.
    N *= conv
    p *= conv
    m_m = 0.
    m_f = -23.74 * sin(N) + 2.68 * sin(2*N) - 0.38 * sin(3*N)
    o_1 = 10.80 * sin(N) - 1.34 * sin(2*N) + 0.19 * sin(3*N)
    k_1 = -8.86 * sin(N) + 0.68 * sin(2*N) - 0.07 * sin(3*N)
    j_1 = -12.94 * sin(N) + 1.34 * sin(2*N) - 0.19 * sin(3*N)
    oo_1 = -36.68 * sin(N) + 4.02 * sin(2*N) - 0.57 * sin(3*N)
    m_2 = -2.14 * sin(N)
    k_2 = -17.74 * sin(N) + 0.68 * sin(2*N) - 0.04 * sin(3*N)
    m_3 = 1.5 * m_2
    tmp1 = (sin(p)+0.2*sin(p-N))
    tmp2 = (2*cos(p)+0.4*cos(p-N))
    m_1 = atan(tmp1 / tmp2) / conv
    if sin(m_1*conv) * tmp1 < 0.:
        m_1 += 180.
    tmp1 = (-0.2505*sin(2*p)-0.1103*sin(2*p-N)-0.0156*sin(2*p-2*N)-0.0366*sin(N)+0.0047*sin(2*p+N))
    tmp2 = (1.0000-0.2505*cos(2*p)-0.1103*cos(2*p-N)-0.0156*cos(2*p-2*N)-0.0366*cos(N)+0.0047*cos(2*p+N))
    l_2 = atan(tmp1 / tmp2) / conv
    if sin(l_2*conv) * tmp1 < 0.:
        l_2 += 180.
    p_1 = 0.
    return eval(raw_data)


def get_variables(date_str, predict=['m_2', 's_2', 'k_1', 'o_1', 'm_4', 'ms_4'], delta_g = 50.):
    cfg_name = 'constituent.ini'
    config = ConfigParser()

    try:
        config.read(cfg_name)
    except:
        print 'Error: config file not found or illegal format!'
        return

    if len(predict) == 0:
        predict = config.sections()

    # return values
    P = 0
    Q = 0
    sigma = []
    V = []
    u = []
    f = []
    # these are only for accompany constituents
    kappa = []
    alpha = []
    corr = []
    # this records the constituents for prediction
    mask = []
    # this records the relationship between constituents and its index
    record = {}

    # check the global config and init date
    date = parser.parse(date_str)
    D = (datetime(year=date.year, month=date.month, day=date.day) -
         datetime(year=date.year, month=1, day=1)).days
    Y = int((date.year - 1901) / 4)
    y = date.year
    t = date.hour + date.minute / 60. + date.second / 3600.

    sigma_tau = 14.49205212
    sigma_s = 0.54901653
    sigma_h = 0.04106864
    sigma_p = 0.00464183
    sigma_N = 0.00220641
    sigma_p_ = 0.00000196

    T_0 = 180.
    s_0 = (277.025 + 129.38481 * (y - 1900) + 13.17640 * (D + Y))
    h_0 = (280.190 - 0.23872 * (y - 1900) + 0.98565 * (D + Y))
    p_0 = (334.385 + 40.66249 * (y - 1900) + 0.11140 * (D + Y))
    N_0 = (259.157 - 19.32818 * (y - 1900) - 0.05295 * (D + Y))
    _p_0 = (281.221 + 0.01718 * (y - 1900) + 0.000047 * (D + Y))
    tau_0 = T_0 - s_0 + h_0

    sigma_const = [sigma_tau, sigma_s, sigma_h, sigma_p, sigma_N, sigma_p_]
    V0_const = [tau_0, s_0, h_0, p_0, N_0, _p_0]

    p = sigma_p * t + p_0
    N = sigma_N * t + N_0

    main_index = {}
    # get the main constituent
    for sec in config.sections():
        if not ('main' in config.options(sec)):
            P += 1
            record[P] = sec
            if sec in predict:
                mask.append(P)
            for opt in config.options(sec):
                raw_data = config.get(sec, opt)
                if opt == 'doodson':
                    tmp_sigma, tmp_V = get_sigma_V(raw_data, sigma_const, V0_const, t)
                    sigma.append(tmp_sigma)
                    V.append(tmp_V)
                elif opt == 'f':
                    f.append(get_f(raw_data, N, p))
                elif opt == 'u':
                    u.append(get_u(raw_data, N, p))
                else:
                    print 'Error: unknown option name!'
                    return
            main_index[sec.lower()] = P

    # get the accompany constituent
    for sec in config.sections():
        if 'main' in config.options(sec):
            Q += 1
            record[P + Q] = sec
            if sec in predict:
                mask.append(P + Q)
            for opt in config.options(sec):
                raw_data = config.get(sec, opt)
                if opt == 'doodson':
                    tmp_sigma, tmp_V = get_sigma_V(raw_data, sigma_const, V0_const, t)
                    sigma.append(tmp_sigma)
                    V.append(tmp_V)
                elif opt == 'f':
                    f.append(get_f(raw_data, N, p))
                elif opt == 'u':
                    u.append(get_u(raw_data, N, p))
                elif opt == 'kappa':
                    kappa.append(float(raw_data))
                elif opt == 'alpha':
                    alpha.append(float(raw_data))
                elif opt == 'main':
                    try:
                        corr.append(main_index[raw_data.lower()])
                    except:
                        print 'Error: main constituent not found!'
                else:
                    print 'Error: unknown option name!'
                    return

    # Convert to array of numpy
    sigma = array(sigma)
    V = array(V)
    f = array(f)
    u = array(u)
    kappa = array(kappa)
    alpha = array(alpha)
    corr = array(corr, dtype=int)
    mask = array(mask, dtype=int)
    pi = 3.141592653589793
    conv = pi / 180.
    return [P, Q, sigma * conv, (V + u) * conv, f, kappa, alpha * conv * delta_g, corr, mask, record]


if __name__ == '__main__':
    print get_variables('1997.8.1 1:00')
