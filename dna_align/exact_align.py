import numpy as np


class BMAlign:
    """
    Boyer-Moore比对算法
    """

    def __init__(self, p, t):
        """
        p: primer or pattern
        t: genome or text
        loc: aligned loc of p in t， -1 for no map
        """
        self.p = p
        self.t = t
        self._element = ['A', 'C', 'G', 'T', 'N']  # 基因组序列组成元素

    @property
    def element(self):
        return self._element
    
    
    @property
    def loc(self):
        return self.get_loc(p, t)


    # @staticmethod
    def get_nearst_index(self, substr, s):
        """
        substr中查找chr,返回最大索引，查不到则返回-1
        substr: string
        s: 单个字符
        """
        index = []
        for i in range(len(substr)):
            if substr[i] == s:
                index.append(i)

        if len(index) >= 1:
            return max(index)
        else:
            return -1

    def bc_shift_table(self, p):
        """bad character rule shift table. ie: row for genome element, col for pattern
            T C G C
        A [[1 2 3 4]
        C [1 '-' 1 '-']
        G [1 2 '-' 1]
        T ['-' 1 2 3]
        N [1 2 3 4]]

        """
        t = self._element

        # 创建空表，行为基因组元素， 列为pattern
        data = np.empty(shape=(len(t), len(p)), dtype='object')

        # 计算bc shift
        # i for row(t element), j for col(primer or pattern)
        for i in range(len(t)):
            for j in range(len(p)):
                cur_t = t[i]
                cur_p = p[j]

                if cur_t == cur_p:
                    data[i][j] = '-'
                elif j == 0 and cur_t != cur_p:
                    data[i][j] = 1
                elif j != 0 and cur_t != cur_p:
                    index = self.get_nearst_index(p[:j], t[i])
                    if index != -1:
                        data[i][j] = j - index
                    else:
                        data[i][j] = j + 1
        return data

    def get_bc_step(self, bc_table, bad_c, j):
        """获取bc step
        j表示p中错配位置索引, bad_c表示坏字符
        """
        index = self._element.index(bad_c)
        step = bc_table[index][j]
        return step

    # @staticmethod
    def max_suffix_len(self, p, j):
        """
        计算p的后缀与前缀能够比对上的最大长度。未找到则返回-1
        如：ACDBAAC,B字符串出错配，则最大后缀AC与前缀AC能够比对。
        j表示错配位置
        """
        k = 0  # 对p错配前进行迭代
        j = j + 1  # 对p错配后进行迭代
        while j < len(p):
            if p[k] == p[j]:
                k, j = k + 1, j + 1
            else:
                j = j - k + 1
                k = 0

        if k != 0:  # 找到最大后缀与前缀匹配
            return k  # 此处k表示最大后缀长度
        else:
            return -1

    def gs_shift_table(self, p):
        """生成good character shift table
        ie: CDAACD: [4, 4, 4, 4, 4, 1]
        """

        # 初始化
        plen = len(p)
        GS = [plen] * plen  # 列表中每个索引表示p该索引位置字符为错配，错配后为好字符。
        GS[plen - 1] = 1  # p最后一个字符错配，则移动一个单位

        for i in range(plen - 1):  # i表示错配位置
            index = p[:i - 1].find(p[i + 1:])
            d = self.max_suffix_len(p, i)

            if index != -1:
                # 1.搜寻p中错配前的字符串中是否还有与p[i + 1:]相等的子串
                GS[i] = i + 1 - index
            elif d != -1:
                # 2.p中i处错配后的后缀串中，找到1个与前缀想匹配
                GS[i] = plen - d
            else:
                # 3.好后缀以及子后缀均未匹配到，需要移动p整个字符串的长度
                GS[i] = plen
        return GS

    @staticmethod
    def get_gs_step(self, gs_table, j):
        """j为pattern中错配索引"""
        return gs_table[j]

    def bc_align(self, p, t):
        """
        仅仅使用bad character rule进行比对
        """
        occurence = []
        bc_table = self.bc_shift_table(p)
        # print(bc_table)

        i = 0  # i for t index, 正序
        j = len(p) - 1  # j for p index, 倒序
        while i <= len(t) - len(p) and j >= 0:
            cur_t = t[i + j]
            if t[i + j] == p[j]:  # match
                j -= 1
            else:
                step = self.get_bc_step(bc_table, cur_t, j)
                print(i, step)
                i += step
                j = len(p) - 1

            if j == -1:
                occurence.append(i)
                i += 1
                j = len(p) - 1
        return occurence

    def gs_align(self, p, t):
        """
        仅仅使用good character rule进行比对
        """
        occurence = []
        gs_table = self.gs_shift_table(p)
        # print(gs_table)
        i = 0  # i for t index, 正序
        j = len(p) - 1  # j for p index, 倒序
        while i <= len(t) - len(p) and j >= 0:
            # while j >= 0:
            if t[i + j] == p[j]:  # match
                j -= 1
            else:
                step = self.get_gs_step(gs_table, j)
                print(i, step)
                i += step
                j = len(p) - 1

            if j == -1:
                occurence.append(i)
                print(i, 1)
                i += 1
                j = len(p) - 1

        return occurence

    def bm_align(self, p, t):
        """使用bc rule 和gs rule 进行Boyer-Moore精确比对"""
        p_len = len(p)
        t_len = len(t)
        occurrence = []

        bc_table = self.bc_shift_table(p)
        gs_table = self.gs_shift_table(p)

        i, j = p_len - 1, p_len - 1  # i for t j  for p
        d = 0  # good char length
        while i < t_len and j < p_len:
            if t[i] == p[j]:
                d += 1
                i, j = i - 1, j - 1
            else:
                index = self._element.index(t[i])
                bc_step = bc_table[index][j]
                gs_step = gs_table[j]
                step = max(bc_step, gs_step)
                i = i + step + d
                j = p_len - 1
                d = 0

            if j == -1:
                occurrence.append([i + 1, i + p_len])
                i = i + 1 + p_len
                j = p_len - 1
                d = 0

        if len(occurrence) > 0:
            return occurrence
        else:
            return -1

    def get_loc(self, p, t):
        """进行bm比对，并返回比对索引"""
        return self.bm_align(p, t)


if __name__ == '__main__':
    p = 'ATACT'
    t = 'ATAGTATAATACTATACT'
    test = BMAlign(p, t)
    print(test.loc)
