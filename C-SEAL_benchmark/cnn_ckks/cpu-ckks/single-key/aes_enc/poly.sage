
import order

sbox = [
	0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5,
	0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
	0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
	0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
	0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC,
	0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
	0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
	0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
	0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
	0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
	0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B,
	0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
	0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85,
	0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
	0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
	0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
	0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17,
	0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
	0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
	0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
	0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
	0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
	0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9,
	0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
	0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6,
	0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
	0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
	0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
	0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94,
	0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
	0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68,
	0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
]

bit_indexes = [[(num >> bit) & 1 for bit in range(8)] for num in range(256)]
bit_sboxes = [[(sbox[num] >> bit) & 1 for bit in range(8)] for num in range(256)]
sbox_polys = []
coefficients_dict = {}


def make_monomial_coefficient():
    R = PolynomialRing(QQ, 'x', 8)  # 生成一个多项式环, x0 到 x7
    x = R.gens()  # 获取生成元, 即变量 x0 到 x7
    for bit in range(8):
        polynomial_total = []  # 初始化结果列表
        for i in range(256):
            # 根据bit_indexes[i][0] 初始化 poly
            poly = x[0] if bit_indexes[i][0] == 1 else (1 - x[0])
            # 循环处理剩余的二进制位
            for j in range(1, 8):
                # 根据 bit_indexes[i][j] 选择 x[j] 或 (1-x[j]) 加入到当前多项式
                poly *= x[j] if bit_indexes[i][j] == 1 else (1 - x[j])
            poly *= bit_sboxes[i][bit]
            # poly *= sbox[i]
            # 将构建的多项式添加到结果列表
            polynomial_total.append(poly)
        for i in range(1, 256):
            polynomial_total[0] += polynomial_total[i]
        sbox_polys.append(polynomial_total[0])
    return x

def make_all_monomial():
    R = PolynomialRing(QQ, 'x', 8)  # 生成一个多项式环, x0 到 x7
    x = R.gens()  # 获取生成元, 即变量 x0 到 x7
    for i in range(256):
        poly = x[0] if bit_indexes[i][0] == 1 else 1
        for j in range(1, 8):
            poly *= x[j] if bit_indexes[i][j] == 1 else 1
        coefficients_dict[ poly ] = 0


# 调用函数
x = make_monomial_coefficient()
make_all_monomial()


for i in range(8):
    print("monomial number is:", len( sbox_polys[i].monomials()) )
    for monomial in sbox_polys[i].monomials():
        coeff = sbox_polys[i].monomial_coefficient(monomial)
        coefficients_dict[monomial] += coeff

        # print(f"{coeff} * {monomial}")

# 现在，我们可以找到所有系数为零的单项式
zero_coefficients_monomials = [monomial for monomial, coeff in coefficients_dict.items() if coeff == 0]

# 打印系数为零的单项式
# print("Monomials with zero coefficient:")
# for monomial in zero_coefficients_monomials:
#     print(monomial)

# 打印所有单项式的系数->C++ format
coefficient_ordered = []

file_name = 'sbox_coeff.txt'
with open(file_name, 'w') as file:
    file.write('# Sbox monomial coefficient is as follows\n')

    R = PolynomialRing(QQ, 'x', 8)  # 生成一个多项式环, x0 到 x7
    x = R.gens()  # 获取生成元, 即变量 x0 到 x7
    for i in range(8):
        c_str = f'static int sbox_{i}[255] = {{ '
        file.write(c_str)
        poly_coeff = []
        for j in range(1, 256):
            poly = x[0] if bit_indexes[j][0] == 1 else 1
            for k in range(1, 8):
                poly *= x[k] if bit_indexes[j][k] == 1 else 1
            # 获取系数
            coeff = sbox_polys[i].monomial_coefficient(poly)
            poly_coeff.append(coeff)
            # 写入系数和多项式到文件
            if j < 255:  # 最后一个元素后不加逗号
                file.write(f"\n {coeff},    //{j}: {poly}")
            else:
                file.write(f"\n{coeff}    //{j}: {poly}\n")
        coefficient_ordered.append(poly_coeff)
        file.write('};\n\n')

print("coeff ordered: ", coefficient_ordered)
# 定义变量
x = var('x0 x1 x2 x3 x4 x5 x6 x7')
# x = var('x0 x1 x2 x3')


# 定义组合函数
def combine(group1, group2):
    result = []
    for g1 in group1:
        for g2 in group2:
            result.append(g1 * g2)
    return result

# 定义基于层次的分组函数
def layered_combine(variables):
    assert len(variables) % 2 == 0, "变量数量应为偶数"

    if len(variables) == 2:
        # 返回变量本身以及它们的乘积
        return [variables[0], variables[1], variables[0]*variables[1]]
    
    mid = len(variables) // 2
    left = layered_combine(variables[:mid])
    right = layered_combine(variables[mid:])
    tmp = combine(left, right) + left + right

    # 组合左右两边的结果
    return combine(left, right) + left + right

def sub_byte( _all_monomials, _monomial_order, _assignment, bit_pos ):
    result = bit_sboxes[0][bit_pos]
    print(hex(result))
    # 使用 subs 替换变量
    for i in range( 0, len( _all_monomials ) ):
        ind = _monomial_order[i]-1
        # print("coeff: ", coefficient_ordered[bit_pos][ind])
        # print(ind)
        result += coefficient_ordered[bit_pos][ind] * _all_monomials[i].subs(_assignment) 
    print(f"final result{bit_pos}: ", hex(result))



# 调用函数生成所有单项式
all_monomials = layered_combine(x)


bitset = []
for i in range(8):
    bitset.append(1<<i)

monomial_order = order.layered_combine_bin(bitset)

assignment = {xi: value for xi, value in zip(x, [1, 1, 0, 0, 0, 0, 0, 0])}

print("bitsbox[0]: ", bit_sboxes[0])

for i in range(8):
    bit_sboxes[0][i]
    sub_byte(all_monomials, monomial_order, assignment, i)


# #打印所有单项式
# print(f"Total number of monomials: {len(all_monomials)}")
# for i in range( len(all_monomials) ):
#     print( "number:", monomial_order[i], all_monomials[i])