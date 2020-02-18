import math


# 定义装药类 ZQWEI Quix 2019.8.22
class tube_grain:  # 管状装药类
    def __init__(self, R, r, rt, end_face_num, L, pro,nozzle_type=0,aeat=1):
        self.R = R
        self.r = r
        self.rt = rt
        self.end_face_num = end_face_num
        self.L = L
        self.At = math.pi * pow(rt, 2)
        self.Vci = math.pi * pow(r, 2) * L
        self.pro = pro
        self.Lo = self.Vci / self.At
        self.nozzle_type=nozzle_type
        self.aeat=aeat
        self.mp=self.pro.density_gr * math.pi * (self.R **2 -self.r **2) * self.L

        def Ctp_calc():
            dc0 = 5.2e-07
            Ctp = (0.74721 ** self.pro.cp_fraction) * (0.12754 * ((self.pro.cp_fraction / self.pro.density_s) ** (1 / 6)) * (
                    (((1 - math.e ** (-0.0001575 * self.Lo)) * (1 + 0.001772 * self.rt * 2)) ** 0.5) / dc0)) ** 0.00501
            return Ctp

        self.Ctp = Ctp_calc()

    def burning(self, e):
        eh = self.R - self.r
        if e <= eh and e <= self.L:
            Vc = math.pi * pow(self.R, 2) * self.L - (self.L - self.end_face_num * e) * math.pi * (
                    pow(self.R, 2) - pow((self.r + e), 2))
            Ab = self.end_face_num * math.pi * (pow(self.R, 2) - pow((self.r + e), 2)) + math.pi * 2 * (
                    self.L - self.end_face_num * e) * (self.r + e)
        else:
            Vc = math.pi * pow(self.R, 2) * self.L
            Ab = 0.000
        return Vc, Ab


class endfaceonly_grain:  # 纯端面装药类
    def __init__(self, R, rt, end_face_num, L, pro, nozzle_type=1,L_intialspace=0.02):
        self.R = R
        self.rt = rt
        self.end_face_num = end_face_num
        self.L = L
        self.At = math.pi * pow(rt, 2)
        self.Vci = math.pi * R ** 2 * L_intialspace
        self.pro = pro
        self.Lo = self.Vci / self.At
        self.mp=self.pro.density_gr * math.pi * (self.R **2) * self.L
        def Ctp_calc():
            dc0 = 5.2e-07
            Ctp = (0.74721 ** self.pro.cp_fraction) * (0.12754 * ((self.pro.cp_fraction / self.pro.density_s) ** (1 / 6)) * (
                    (((1 - math.e ** (-0.0001575 * self.Lo)) * (1 + 0.001772 * self.rt * 2)) ** 0.5) / dc0)) ** 0.00501
            return Ctp

        self.Ctp = Ctp_calc()
        self.nozzle_type = nozzle_type

    def burning(self, e):
        if e <= self.L:
            Ab = self.end_face_num * math.pi * self.R ** 2
            Vc = math.pi * self.R ** 2 * e + self.Vci
        else:
            Vc = math.pi * self.R ^ 2 * self.L + self.Vci
            Ab = 0

        return Vc, Ab
