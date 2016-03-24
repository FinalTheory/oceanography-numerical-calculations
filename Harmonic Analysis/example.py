# -*- coding: utf-8 -*-

import Tides
import dateutil
import datetime
import numpy as np
from getData import get_variables
import matplotlib.pyplot as plt

#先定义一个用来绘图的函数
def show_diff(original, predict):
	fig, axes = plt.subplots(figsize=(18,4))
	axes.plot(original, 'b', lw=2)
	axes.plot(predict, 'r', lw = 0.5)
	axes.set_title('Diff between real data and predicted')
	axes.legend(["real", "predicted"])
	plt.show()

#计算日期
date = dateutil.parser.parse('21/09/2014 1:20:30.000 PM')
middle_date = (date + datetime.timedelta(hours = 371) - datetime.timedelta(hours=8)).ctime()
#读入数据
zeta = np.loadtxt('Zeta.txt')
#单位转换成毫米
zeta = zeta * 1000
#输入中间时刻，计算天文要素
P, Q, sigma, V, f, kappa, alpha, corr, mask, record = get_variables(middle_date)
#进行调和分析，计算调和常数
H, g, S_0 = Tides.analyze(P, zeta, sigma, V, f, kappa, alpha, corr)
#带入调和常数进行回报
zeta_predict = Tides.predict(-172, 172, sigma, V, f, H, g, S_0, mask, 1.)
#回报逐分钟的潮位用以计算高低潮时和潮高
zeta_predict_minute = Tides.predict(-172*60, 172*60, sigma, V, f, H, g, S_0, mask, 1./60.)
#同图绘制回报过程曲线与实测过程曲线
show_diff(zeta, zeta_predict);
#计算回报的过程曲线与真实测量值的相关系数，若其数值相当接近1，说明回报结果正确
print '相关系数：%.5f' % plt.xcorr(zeta, zeta_predict, maxlags=1)[1][1]
print '平均海面：%.2fmm' % S_0
print '分潮及其对应的H、g：'
for k in record.keys():
    print '分潮名称: %s\t\t振幅: %.2fmm\t\t迟角: %.2frad' % (record[k], H[k-1], g[k-1])

    
def print_max_min(data, time):
    ssum1 = 0.
    ssum2 = 0.
    num1 = 0
    num2 = 0
    t = dateutil.parser.parse(time)
    for i in range(1, data.size - 1):
        if data[i] > data[i - 1] and data[i] > data[i + 1]:
            print '高潮时: %s\t潮位:%dmm' % ( (t + datetime.timedelta(minutes=i)).strftime('%c'), data[i] )
            ssum1 += data[i]
            num1 += 1
        elif data[i] < data[i - 1] and data[i] < data[i + 1]:
            print '低潮时: %s\t潮位:%dmm' % ( (t + datetime.timedelta(minutes=i)).strftime('%c'), data[i] )
            ssum2 += data[i]
            num2 += 1
    print '平均潮差：%.2fmm' % ( ssum1 / float(num1) - ssum2 / float(num2)  )

print_max_min(zeta_predict_minute, '14/09/2014 9:20:30.000 AM')
fig, axes = plt.subplots(figsize=(18,4))
axes.plot(zeta - zeta_predict);
axes.set_title('Residual water level')
plt.show();