# 火箭内弹道算法调用测试
# @ZQWEI Quix 2019.8.24  MIT License

import rocketsimu as rs
import rocketsimu.motorsimu as ms
from rocketsimu import propellent,grain,nozzle
import matplotlib.pyplot as plt

# 定义药柱参数
R = 0.01  # 药柱外径
r = 0.002  # 药柱内径
rt = 0.002  # 喉管半径
L = 0.05  # 药柱总长
end_faces = 0  # 可燃端面数

total_time = 1.5  # 计算总时间
step_length = 0.0001  # 计算时间步
t_intial = 0  # 初始时间（一般为0）
p_intial = 101325  # 初始压强（一般为大气压）
two_phase_model_swtich = True  # 是否打开两相流模型
nozzle_type = 1  # 选择尾喷管模型
erosion_ratio = 1  # 默认平均侵蚀比

# ms.config_log(True)

propellent_test = rs.propellent.KNSB()  # 定义燃料
grain_test = rs.grain.tube_grain(R, r, rt, end_faces, L, propellent_test, nozzle_type)  # 实例化药柱类
nozzle_test=rs.nozzle.straight_nozzle(effiency=0.9)

P, t = ms.pressure_calc(t_intial, p_intial, step_length, -1, grain_test, two_phase_model_swtich,
                         erosion_ratio)  # 四阶龙格库塔法求解/自动时间步求解
F = ms.thrust_calc(P, grain_test, nozzle_test)
It = ms.impulse_calc(F,t)

# 展示结果
print("Impulse："+str(It))
plt.plot(t,P)
plt.xlabel("time/s")
plt.ylabel("pressure/MPa")
plt.figure()
plt.plot(t,F)
plt.xlabel("time/s")
plt.ylabel("force/N")
plt.show()