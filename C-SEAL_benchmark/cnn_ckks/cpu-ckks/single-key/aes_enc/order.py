# 定义组合函数
def combine_bin(group1, group2):
    result = []
    for g1 in group1:
        for g2 in group2:
            result.append( (g1^g2)&0xff )
    return result

# 定义基于层次的分组函数
def layered_combine_bin(variables):
    assert len(variables) % 2 == 0, "变量数量应为偶数"
    if len(variables) == 2:
        # 返回变量本身以及它们的乘积
        return [variables[0], variables[1], variables[0]^variables[1]]
    
    mid = len(variables) // 2
    left = layered_combine_bin(variables[:mid])
    right = layered_combine_bin(variables[mid:])

    # 组合左右两边的结果
    return combine_bin(left, right) + left + right


