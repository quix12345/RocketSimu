import math


# 固体火箭燃料类定义 ZQWEI Quix 2019.8.22
class KNDX:
    def __init__(self, density_gr_fraction=0.97, density_gr=1800, k_chamber=1.1308, k_2phase=1.0435, c=912, Mr=42.39, Tp=1600,
                 cp_fraction=0.425,
                 density_s=2430):
        self.density_gr = density_gr_fraction * density_gr
        self.density_gr_fration = density_gr_fraction
        self.k_chamber = k_chamber
        self.k_2phase = k_2phase
        self.c = c
        self.Mr = Mr
        self.Tp = Tp
        self.dagama = math.sqrt(k_chamber) * pow((2 / (k_chamber + 1)), ((k_chamber + 1) / (2 * (k_chamber - 1))))
        self.cp_fraction = cp_fraction
        self.density_s = density_s
        self.name= 'KNDX'

    def burnrate(self, p):
        if p < 779000:
            n = 0.619
            b = 1.714621881352520e-06

        elif 779000 <= p < 2572000:
            n = -0.009
            b = 0.008553019943518

        elif 2572000 <= p < 5930000:
            n = 0.688
            b = 2.860515511772906e-07

        elif 5930000 <= p < 8502000:
            n = -0.148
            b = 0.132901060636461

        elif p <= 0:
            n = 0
            b = 0
        else:
            n = 0.442
            b = 1.064077783776570e-5
        burning_rate = b * pow(p, n)
        return b, n, burning_rate


class KNSB:
    def __init__(self, density_gr_fraction=0.97, density_gr=1800, k_chamber=1.1361, k_2phase=1.03, c=909, Mr=39.86, Tp=1600,
                 cp_fraction=0.436,
                 density_s=2430):
        self.density_gr = density_gr_fraction * density_gr
        self.density_gr_fration = density_gr_fraction
        self.k_chamber = k_chamber
        self.k_2phase = k_2phase
        self.c = c
        self.Mr = Mr
        self.Tp = Tp
        self.dagama = math.sqrt(k_chamber) * pow((2 / (k_chamber + 1)), ((k_chamber + 1) / (2 * (k_chamber - 1))))
        self.cp_fraction = cp_fraction
        self.density_s = density_s
        self.name = 'KNSB'

    def burnrate(self, p):
        if p < 807000:
            n = 0.625
            b = 1.904181592269679e-6

        elif 807000 <= p < 1503000:
            n = -0.314
            b = 0.670892306636334

        elif 1503000 <= p < 3792000:
            n = -0.013
            b = 0.009396806651824

        elif 3792000 <= p < 7033000:
            n = 0.535
            b = 2.409036672272810e-6

        else:
            n = 0.064
            b = 0.003987147536711
        burning_rate = b * pow(p, n)
        return b, n, burning_rate


class KNSU:
    def __init__(self, density_gr_fraction=0.97, density_gr=1800, k_chamber=1.133, k_2phase=1.0437, c=919, Mr=41.98, Tp=1600,
                 cp_fraction=0.424,
                 density_s=2430):
        self.density_gr = density_gr_fraction * density_gr
        self.density_gr_fration = density_gr_fraction
        self.k_chamber = k_chamber
        self.k_2phase = k_2phase
        self.c = c
        self.Mr = Mr
        self.Tp = Tp
        self.dagama = math.sqrt(k_chamber) * pow((2 / (k_chamber + 1)), ((k_chamber + 1) / (2 * (k_chamber - 1))))
        self.cp_fraction = cp_fraction
        self.density_s = density_s
        self.name = 'KNSB'

    def burnrate(self, p):
        n = 0.319
        b = 1.007251105591617e-4
        burning_rate = b * pow(p, n)
        return b, n, burning_rate


class APCP_NOAL:
    def __init__(self, density_gr_fraction=1, density_gr=1650, k_chamber=1.3, k_2phase=1.3, c=1000, Mr=25, Tp=2900, cp_fraction=0.01,
                 density_s=2430):
        self.density_gr = density_gr_fraction * density_gr
        self.density_gr_fration = density_gr_fraction
        self.k_chamber = k_chamber
        self.k_2phase = k_2phase
        self.c = c
        self.Mr = Mr
        self.Tp = Tp
        self.dagama = math.sqrt(k_chamber) * pow((2 / (k_chamber + 1)), ((k_chamber + 1) / (2 * (k_chamber - 1))))
        self.cp_fraction = cp_fraction
        self.density_s = density_s
        self.name = 'APCP(No Al)'

    def burnrate(self, p):
        if 0 <= p < 8e5:
            n = 0.969
            b = 7.588000000000000e-09

        elif p >= 8e5:
            n = 0.3203
            b = 4.982000000000001e-05

        else:
            n = 0
            b = 0
        burning_rate = b * pow(p, n)
        return b, n, burning_rate

class NONE:
    pass
