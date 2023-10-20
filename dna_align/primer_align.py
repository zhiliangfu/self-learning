# 描述
# 给定1对引物（或接头）和1条fa序列，返回引物（或接头）比对位置。
# 比对默认允许有15%的错配。
import numpy as np
import math

class SeqAlign:
    """
    1条adapter或single primer序列，和1条target 序列进行比对（计算edist距离）。
    返回：最小距离，比对到的链和位置、比对到的序列
    """

    def __init__(self, p, t, end, error=0.15):
        """
        :param p: 5'->3' primer or adapter sequence
        :param t: 5'->3' target sub-sequence
        :param end: int value. 5 or 3 end
        :param error: error rate, default 0.15
        """
        self.p = p
        self.plen = len(p)
        self.t = t
        self.tlen = len(t)
        self.end = end
        self.threshold = math.ceil(len(p)*error)
        self.edist, self.strand = self.edist_matrixs(p, t)
        self.mismatch = self.get_dist(self.edist)
        self.loc = self.get_loc(self.edist, end, self.strand, self.threshold, self.mismatch)
        self.sseq = self.get_aligned(t, self.loc)  # 仅供参考，可能有1bp出入


    def reverse_complement(self, str):
        """反向互补序列"""
        completement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        s = ''
        for base in str:
            s = completement[base] + s
        return s


    def edist_matrix(self, p, t):
        """
        calculate edist distance between p and t
        p: primer or pattern
        t: target sequence
        """
        p_len = len(p)
        t_len = len(t)
        dist_matrix = np.zeros(shape=(p_len + 1, t_len + 1), dtype=int)  # p for rows, t for cols

        # i for row(primer), j for col(target sequence)
        for i in range(p_len + 1):
            dist_matrix[i][0] = i

        # edist distance
        for i in range(1, p_len + 1):
            for j in range(1, t_len + 1):
                hor = dist_matrix[i][j - 1] + 1  # horizontal
                ver = dist_matrix[i - 1][j] + 1  # vertical

                # diag
                if p[i - 1] == t[j - 1]:
                    delta = 0
                else:
                    delta = 1
                dia = dist_matrix[i - 1][j - 1] + delta

                dist_matrix[i][j] = min(dia, ver, hor)
        return dist_matrix

    def get_dist(self, dist_m):
        return min(dist_m[len(dist_m) - 1])

    def edist_matrixs(self, p, t):
        """
        p和t以及t的反向互补各比对1次
        """
        edist_matrix1 = self.edist_matrix(p, t)
        edist_matrix2 = self.edist_matrix(p, self.reverse_complement(t))
        dist1 = self.get_dist(edist_matrix1)
        dist2 = self.get_dist(edist_matrix2)
        if dist1 <= dist2:
            return edist_matrix1, '+'
        else:
            return edist_matrix2, '-'

    def get_loc(self, dist_m, end, strand, threshold, mismatch):
        """
        当有多种比对结果时: 切割序列时，尽可能多（贪婪）切割序列。
        输出p在t中的比对位置
        """
        p_len = len(dist_m)-1
        t_len = len(dist_m[0])-1

        data = dist_m[len(dist_m) - 1].tolist()
        data.pop(0)
        min_data = min(data)
        index = []
        for i in range(len(data)):
            if data[i] == min_data:
                index.append(i)

        # 返回索引
        if mismatch <= threshold and strand == '+':
            if end == 5:
                return [max(index)-p_len+1, max(index)]
            if end == 3:
                return [min(index)-p_len+1, min(index)]
        elif mismatch <= threshold and strand == '-':
            if end == 5:
                return [t_len-min(index)+p_len-2, t_len-min(index)-1]
            if end == 3:
                return [t_len-max(index)+p_len-2, t_len-max(index)-1]
        else:
            return -1

    def get_aligned(self, t, loc):
        """提取t中被匹配到的序列"""
        if loc != -1:
            return t[min(loc):max(loc)+1]
        else:
            return -1

class SeqsAlign(SeqAlign):
    """1对引物或接头进行比对，比对上则返回比对上的位置"""
    def __init__(self, p1, p2, t, error=0.15):
        self.p1 = p1   # 5' seq
        self.p2 = self.reverse_complement(p2)  # 3' seq
        self.p1len = len(p1)
        self.p2len = len(p2)
        self.t = t
        self.tlen = len(t)
        self.error = error
        self.p1_obj, self.p2_obj = self.align(p1, p2, t, error)
        self.loc = self.get_loc(self.p1_obj, self.p2_obj)
        self.sseq = self.get_aligned(self.p1_obj, self.p2_obj)

    def align(self, p1, p2, t, error):
        """双端进行比对，为提高运行速度，只提取两端部分序列进行分析"""
        ratio = 1/3
        loc1_index = math.ceil((len(t)-1) * ratio)
        loc2_index = math.ceil((len(t)-1) * (1-ratio))

        forward = SeqAlign(p1, t[:loc1_index], end=5, error = error)
        reverse = SeqAlign(p2, t[loc2_index:], end=3, error = error)

        if reverse.loc != -1:
            reverse.loc[0] += loc2_index
            reverse.loc[1] += loc2_index
        return forward, reverse

    def get_loc(self, p1_obj, p2_obj):
        """提取双端比对上的位置"""
        return [p1_obj.loc, p2_obj.loc]

    def get_aligned(self, p1_obj, p2_obj):
        return [p1_obj.sseq, p2_obj.sseq]


if __name__ == '__main__':
    # a1 = 'ACTTGCCTGTCGCTCTATCTTC'
    # a2 = 'GCAATATCAGCACCAACAGAAA'
    # fa = 'CTTGTACTTCGTTCCAGTTGCGGTATTGCTCGAACCTTTGAATCAAATTAACCTTTTCTGTTGGTGCTGATATTGCTTTATTCAGTCCTTCCACGCCTACTATAATAACTGTAATTAGTTATTTCGGTTTCTTATTGGCTTCTATCATTTTCACCTAATTCTATTCATTGGTTCGAGCAAGATACAACTTATTCAAAATGAAGGAAGGAAACCCTCTCAGTTCACTATTTCAATTTCAACTAATTTGTTTATTCCTCTATATCAATTTCTGATCATTGAGATTCATAGCAAATTAAGGTACTATCCAGTATAGCTATCATCCCTACTTATTCCGTTGGATAGAAATGATTGAGGCTTTGCCATCCGGAATTGTGTTAGGTCCCAATTCCTATAACTTTGTCGGATTATTCATGACTGCATATTCACAATACCGACGTGGTGATCAATTGGATGTTGATCCTGGACATCCTCCAATTTCATTTGGCCTCCTACTTTTCCTGGAGGTGATTCCATTGCTCGTTCAGTTGGAATGAATTATCGATAATTTTTGAGGGTTGGATTTGATAGGAACAGAATCACGCTCTGTAGGATTTGAACCTACGGCATCAGGTTTTGGAGACCGAAGATAGAGCGACAGGCAAGTAGGTTAATTCTGATTCAAAGGTTGGTTGTTCAGCGGCAATACG'
    # test = SeqsAlign(a1, a2, fa)
    # print(test.loc)
    # print(test.sseq)

    a1 = 'ATCG'
    # a2 = 'GCAATATCAGCACCAACAGAAA'
    fa = 'AAAATCATAGGAAATC'
    test = SeqAlign(a1, fa, 5)
    print(test.mismatch)
    # print(test.sseq)





