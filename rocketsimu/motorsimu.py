# 火箭发动机两相流内弹道微分方程的龙格库塔解法程序 @ZQWEI Quix 2019.9.9
# MIT License, All right reserved.

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import simps
import warnings
import math

def ibc_function(t, p, e, grain, two_ph=True, Erosion_ratio=1):  # ibc系internal_ballistic_model 简写
    [Vc, Ab] = grain.burning(e)
    cp_fraction = grain.pro.cp_fraction
    density_gr = grain.pro.density_gr
    dagama = grain.pro.dagama
    c = grain.pro.c
    Ctp = grain.Ctp
    At = grain.At
    [_, _, burnrate] = grain.pro.burnrate(np.real(p))

    if two_ph:
        dy = (1 / Vc) * (((1 - cp_fraction) * density_gr) * pow(dagama, 2) * pow(c, 2) * Ab * burnrate * Erosion_ratio - pow(dagama,
                                                                                                                 2) * c * At * (
                                 Ctp * pow(p, 0.000835)) * p)
    else:
        dy = (1 / Vc) * (density_gr * pow(dagama, 2) * pow(c, 2) * Ab * burnrate * Erosion_ratio - pow(dagama,
                                                                                                          2) * c * At * (
                                 Ctp * pow(p, 0.000835)) * p)

    return np.real(dy)


def pressure_calc(x0, y0, h, N, grain, two_ph_model=True, Erosion_ratio=1, longer_time=0.1):
    if IsLogEnabled:
        print("Calculating the internal ballistic model....")
    e = 0
    n = 1
    pressure = [1]
    t0 = [0]
    k = 1
    if N == -1:
        j = 0
        while (j < (longer_time / h)):
            x1 = x0 + h
            k1 = ibc_function(x0, y0, e, grain, two_ph_model, Erosion_ratio)
            k2 = ibc_function(x0 + h / 2, y0 + h * k1 / 2, e, grain, two_ph_model, Erosion_ratio)
            k3 = ibc_function(x0 + h / 2, y0 + h * k2 / 2, e, grain, two_ph_model, Erosion_ratio)
            k4 = ibc_function(x1, y0 + h * k3, e, grain, two_ph_model, Erosion_ratio)
            y1 = y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            t0.append(x1)
            pressure.append(y1)
            a, b, r = grain.pro.burnrate(pressure[n - 1])
            e += -(x0 - x1) * r
            n = n + 1
            x0 = x1
            y0 = y1
            if y1 < 101325:
                j += 1
            else:
                j = 0
            k = __progress_calc(j, k, (longer_time / h) + 1)
    else:
        while n != N / h :
            x1 = x0 + h
            k1 = ibc_function(x0, y0, e, grain, two_ph_model, Erosion_ratio)
            k2 = ibc_function(x0 + h / 2, y0 + h * k1 / 2, e, grain, two_ph_model, Erosion_ratio)
            k3 = ibc_function(x0 + h / 2, y0 + h * k2 / 2, e, grain, two_ph_model, Erosion_ratio)
            k4 = ibc_function(x1, y0 + h * k3, e, grain, two_ph_model, Erosion_ratio)
            y1 = y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            t0.append(x1)
            pressure.append(y1)
            [_, _, r] = grain.pro.burnrate(pressure[n - 1])
            e += -(x0 - x1) * r
            n = n + 1
            x0 = x1
            y0 = y1
            k = __progress_calc(n, k, (longer_time / h) + 1)
    k = __progress_calc(j, -1, (longer_time / h) + 1)
    return pressure, t0


def thrust_calc(pressure, grain, nozzle, Pa=101325):
    if IsLogEnabled:
        print("Calculating the thrust performance....")
    length = len(pressure)
    i = 0
    j = 1
    F = []
    k = grain.pro.k_chamber
    k2 = grain.pro.k_2phase
    dagama = grain.pro.dagama
    if nozzle.nozzle_type == "straight":  # 直喷管
        while i <= length - 1:
            pe = pressure[i] / (1 + (k - 1) * 0.5) ** (k / (k - 1))
            Cf = (dagama * math.sqrt((2 * k2 / (k2 - 1)) * (1 - (pe / pressure[i]) ** ((k2 - 1) / k2))) + (
                    pe / pressure[i] - Pa / pressure[i]))
            force = Cf * nozzle.effiency * grain.At * pressure[i]
            if force <= 0:
                F.append(0)
            else:
                F.append(Cf * nozzle.effiency * grain.At * pressure[i])
            i += 1
            j = __progress_calc(i, j, length)
    elif nozzle.nozzle_type == "auto_laval":  # 自动设计最优化拉法尔喷管
        global p, optimum_arearatio, a
        optimum_arearatio = optimum_arearatio_calc(pressure,grain,nozzle, Pa)
        warnings.filterwarnings("ignore")
        # mae=math.sqrt((2/(k-1))*((pressure[i]/Pa)**((k-1)/k)-1))
        # pe=pressure[i]/(1+(k-1)*0.5*mae**2)**(k/(k-1))
        while i <= length - 1:
            def pressure_solving_function(x):
                ans = ((((k2 + 1) / 2) ** (1 / (k2 - 1))) * ((x / p) ** (1 / k2)) * (
                        (((k2 + 1) / (k2 - 1)) * (1 - ((x / p) ** ((k - 1) / k2)))) ** 0.5)) ** (-1) - optimum_arearatio
                return ans

            p = pressure[i]
            pe = max(fsolve(pressure_solving_function, [101325]))
            if pe <= Pa or pe == None:
                pe = Pa
            try:
                Cf = (dagama * math.sqrt((2 * k2 / (k2 - 1)) * (1 - (pe / pressure[i]) ** ((k2 - 1) / k2))) + (
                        pe / pressure[i] - Pa / pressure[i]))
                F.append(Cf * nozzle.effiency * grain.At * pressure[i])
            except BaseException:
                F.append(0)
            finally:
                i += 1

            j = __progress_calc(i, j, length)
    elif nozzle.nozzle_type == "customize_laval":  # 自定义拉法尔喷管
        arearatio = nozzle.aeat
        warnings.filterwarnings("ignore")
        while i <= length - 1:
            def pressure_solving_function(x):
                ans = ((((k2 + 1) / 2) ** (1 / (k2 - 1))) * ((x / p) ** (1 / k2)) * (
                        (((k2 + 1) / (k2 - 1)) * (1 - ((x / p) ** ((k2 - 1) / k2)))) ** 0.5)) ** (-1) - arearatio
                return ans

            p = pressure[i]
            pe = max(fsolve(pressure_solving_function, [101325]))
            if pe <= Pa or pe == None:
                pe = Pa
            try:
                Cf = (dagama * math.sqrt((2 * k2 / (k2 - 1)) * (1 - (pe / pressure[i]) ** ((k2 - 1) / k2))) + (
                        pe / pressure[i] - Pa / pressure[i]))
                F.append(Cf * nozzle.effiency * grain.At * pressure[i])
            except BaseException:
                F.append(0)
            finally:
                i += 1
            j = __progress_calc(i, j, length)
    return F


def optimum_arearatio_calc(pressure, grain,nozzle, pa=101325, limit=100000):
    length = len(pressure)
    i = 0
    global endp1, endp2
    while i < length - 1:
        if pressure[i] > limit:
            endp1 = i
            break
        else:
            i += 1
    i += 1
    while i < length - 1:
        if (pressure[i] < limit):
            endp2 = i
            break
        else:
            i += 1
    avg_pressure = np.mean(pressure[endp1:endp2])
    k = grain.pro.k_2phase
    optimum_arearatio = ((((k + 1) / 2) ** (1 / (k - 1))) * ((pa / avg_pressure) ** (1 / k)) * (
            (((k + 1) / (k - 1)) * (1 - ((pa / avg_pressure) ** ((k - 1) / k)))) ** 0.5)) ** -1
    nozzle.aeat=optimum_arearatio
    return optimum_arearatio

def impulse_calc(F,t):
    impulse= simps(F, t)
    return impulse

def thrust2pressure(F, grain, arearatio, xl, pa=101325):
    if IsLogEnabled:
        print("Converting the thrust into pressure data....")
    length = len(F)
    i = 0
    j = 1
    po = []
    k = grain.pro.k_chamber
    k2 = grain.pro.k_2phase
    dagama = grain.pro.dagama
    At = grain.At
    global F_temp
    while i < length:
        def mae_solving_function(x):
            ans = (1 / x) * ((1 + ((k2 - 1) / 2) * x ** 2) * (2 / (k2 + 1))) ** ((k2 + 1) / (2 * (k2 - 1))) - arearatio
            return ans

        def po_solving_function(x):
            ans = xl * (dagama * math.sqrt((2 * k2 / (k2 - 1)) * (
                    1 - ((x / (1 + (k - 1) * 0.5 * mae ** 2) ** (k / (k - 1))) / x) ** ((k2 - 1) / k2))) + arearatio * (
                                (x / (
                                        1 + (k - 1) * 0.5 * mae ** 2) ** (
                                         k / (k - 1))) / x - pa / x)) * At * x - F_temp
            return ans

        mae = fsolve(mae_solving_function, [1])
        F_temp = F[i]
        po_temp = fsolve(po_solving_function, [10e6])
        if po_temp <= pa or po_temp == None:
            po.append(pa)
        else:
            po.append(po_temp[0])
        i += 1

        j = __progress_calc(i, j, length)
    return po


def pro_gasphase_perform_calc(po, t, grain, density_gr, c=-1):
    if IsLogEnabled:
        print("Solving for the propellent performance using the pure gas model....")
    At = grain.At
    if c == -1:
        c = (At / grain.mp) * simps(po, t)
    e = 0
    i = 0
    told = 0
    burnrate = [0]
    pressure_new = [0]
    length = len(po)
    j = 1
    while i < length - 1:
        [_, Ab] = grain.burning(e)
        if (Ab == 0):
            __progress_calc(i, -1, length)
            break
        else:
            r = po[i] / (((Ab / At) + 0.0001) * density_gr * c)
        if r < 0.3 * burnrate[i - 1]:
            __progress_calc(i, -1, length)
            break
        e = e + r * (t[i + 1] - told)
        burnrate.append(r)
        pressure_new.append(po[i])
        told = t[i + 1]
        i += 1
        j = __progress_calc(i, j, length)
    return burnrate, pressure_new


def pro_2phase_perform_calc(pressure, t, grain, density_gr, density_s, k, cp_fraction, erotion_ratio, c=-1):
    if IsLogEnabled:
        print("Solving for the propellent performance using the two phase flow model....")
    At = grain.At
    if c==-1:
        c = (At/grain.mp)*simps(pressure,t)
    length = len(pressure)
    i = 0
    dagama = math.sqrt(k) * (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)))
    dc0 = 5.2e-07
    e = 0
    told = 0
    j = 1
    burnrate = [0]
    pressure_new = [0]
    while i < length - 1:
        [Vc, Ab] = grain.burning(e)
        if Ab == 0:
            __progress_calc(i, -1, length)
            break
        dpdt = ((pressure[i + 1] - pressure[i]) / (t[i + 1] - t[i]))
        if cp_fraction == 0:
            Ctp = 1
        else:
            Ctp = (0.74721 ** cp_fraction) * (0.12754 * ((cp_fraction / density_s) ** (1 / 6)) * (
                    (((1 - math.e ** (-0.0001575 * grain.Lo)) * (
                                1 + 0.001772 * grain.rt * 2)) ** 0.5) / dc0)) ** 0.00501
        r = (Vc * dpdt + ((dagama ** 2) * c * At * (Ctp * (pressure[i] ** 0.000835))) * pressure[i]) / (
                ((1 - cp_fraction) * density_gr) * (dagama ** 2) * (c ** 2) * Ab * erotion_ratio)
        if r < 0.3 * burnrate[i-1]:
            __progress_calc(i, -1, length)
            break
        e = e + r * (t[i] - told)
        told = t[i]
        burnrate.append(r)
        pressure_new.append(pressure[i])
        i += 1
        j = __progress_calc(i, j, length)
    return burnrate, pressure_new


def pro_burnrate_calc(pro, p_max=10e6, step=1000):
    pressure = []
    burnrate = []
    limit = p_max / step
    i = 0
    p = 0
    while i < limit:
        p += step
        a, b, r = pro.burnrate(p)
        burnrate.append(r)
        pressure.append(p)
        i += 1
    return burnrate, pressure


def __progress_calc(i, j, length, constant=25):
    # function for calculating the progress
    global IsLogEnabled
    if IsLogEnabled:
        if j==-1 :
            print("\rFinished %s" % str(100.00) + "%", end='\n')
        else:
            if i > j * (length - 1) / constant:
                if j / (constant * 0.01)!=100:
                    print("\rFinished %s" % str(j / (constant * 0.01)) + "%", end='')
                else:
                    print("\rFinished %s" % str(100.00) + "%", end='\n')
                j += 1
        return j

def config_log(config):
    if isinstance(config,str):
        global IsLogEnabled
        if config=="true" or config =="1" or config =="True":
            IsLogEnabled=True
        elif config=="false" or config =="0" or config =="False":
            IsLogEnabled = False
        else:
            raise Exception("Invalid config!", config)
    elif isinstance(config,bool):
            IsLogEnabled = config
    else:
        raise Exception("Invalid input type!", config)

# initiate config_log to false
config_log(False)