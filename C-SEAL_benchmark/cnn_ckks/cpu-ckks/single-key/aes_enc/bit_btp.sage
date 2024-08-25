# 导入必需的 SageMath 模块
from sage.all import *

# 定义函数 f(x) = (1 - cos(2πx))/2
def f(x):
    return (1 - cos( pi * x)) / 2

def AND(x, y):
    x = (x+y)/3
    return ( 1 - 2 * sin( 2*pi*x + pi/6 ) ) / 3
def NAND(x, y):
    x = (x+y)/3
    return ( 2 * ( 1 + sin( 2*pi*x + pi/6 ) ) ) / 3
def OR(x, y):
    x = (x+y)/3
    return 2*(1 - cos(2*pi*x)) / 3
def XOR(x,y):
    x = (x+y)/3
    return ( 1 + 2*sin(2*pi*x - pi/6 )) / 3
def NOR(x,y):
    x = (x+y)/3
    return (1 + 2 * cos(2*pi*x)) / 3
def XNOR(x,y):
    x = (x+y)/3
    return  2* (1 - sin(2*pi*x - pi/6 ) ) / 3

# 指定计算的点
x_value = 351.0
# x_value = x_value/2
# 计算 f(x_value)
result = f(x_value).n()  # 或者使用 .numerical_approx() 同样效果

a = 1.0
b = 0.0
y = AND(a, b).n()
# 输出结果
print(f"AND: f({a},{b}) = {y}")
y = NAND(a, b).n()
# 输出结果
print(f"NAND: f({a},{b}) = {y}")

y = OR(a, b).n()
# 输出结果
print(f"OR: f({a},{b}) = {y}")

y = XOR(a, b).n()
# 输出结果
print(f"XOR: f({a},{b}) = {y}")

y = NOR(a, b).n()
# 输出结果
print(f"NOR: f({a},{b}) = {y}")

y = XNOR(a, b).n()
# 输出结果
print(f"XNOR: f({a},{b}) = {y}")

# 输出结果
print(f"f({x_value}) = {result}")