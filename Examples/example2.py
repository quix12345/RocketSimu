# 火箭内弹道算法调用测试
# @ZQWEI Quix 2020.2.18 MIT License

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

ms.config_log(True)

propellent_test = rs.propellent.KNSB()  # 定义燃料
grain_test = rs.grain.tube_grain(R, r, rt, end_faces, L, propellent_test, nozzle_type)  # 实例化药柱类
laval_nozzle_test=rs.nozzle.laval_nozzle(effiency=0.9,IsCustomized=False) # 实例化拉法尔喷管，指定为自动优化设计类型
straight_nozzle_test=rs.nozzle.straight_nozzle(effiency=0.9)

P, t = ms.pressure_calc(t_intial, p_intial, step_length, -1, grain_test, two_phase_model_swtich,
                         erosion_ratio)  # 四阶龙格库塔法求解/自动时间步求解
F_laval = ms.thrust_calc(P, grain_test, laval_nozzle_test)
F_str = ms.thrust_calc(P, grain_test, straight_nozzle_test)
It_laval = ms.impulse_calc(F_laval,t)
It_str = ms.impulse_calc(F_str,t)

# 展示结果
print("Impulse  Laval vs Straight：\r\n"+str(It_laval)+" vs "+str(It_str))
plt.plot(t,P)
plt.xlabel("time/s")
plt.ylabel("pressure/MPa")
plt.figure()
p1=plt.plot(t,F_laval)
p2=plt.plot(t,F_str)
plt.xlabel("time/s")
plt.ylabel("force/N")
plt.show()