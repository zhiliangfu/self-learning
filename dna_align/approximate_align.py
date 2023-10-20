import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class GlobalAlign:
    """needleman-wunch algorithm"""

    def __init__(self, p, t, match=5, mismatch=-4, gap_open=-10, gap_extend=-0.5):
        """
        end_weight控制从score_matrix哪个位置进行回溯。
        - True表示从右下角回溯，
        - False表示从整个score_matrix进行回溯
        """
        self.p = p
        self.t = t
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.score_matrix, self.trace_matrix = self.score_matrix(p, t)
        self.top, self.middle, self.bottom = self.trace(p, t, self.trace_matrix)
        self.score = self.get_score(self.score_matrix)
        self.mismatch = self.get_mismatch(self.middle)

    @property
    def align(self):
        print('Align Seqs:')
        print('\tp:\t' + self.top)
        print('\t\t' + self.middle)
        print('\tt:\t' + self.bottom)

    def score_matrix(self, p, t):
        """计算得分矩阵
        行为primer or pattern， 列为genome or t。
        第一行和第一列均罚分-1
        """

        p_len = len(p) + 1
        t_len = len(t) + 1
        score_matrix = np.zeros(shape=(p_len, t_len))
        trace_matrix = np.zeros(shape=(p_len, t_len), dtype=int)  # 0: up, 1:diag, 2:left

        # 初始化score_matrix第一行第一列
        # i for row(primer), j for col(genome)
        for i in range(1, p_len):
            if i == 1:
                score_matrix[i][0] = self.gap_open
            else:
                score_matrix[i][0] = self.gap_open + (i - 1) * self.gap_extend
            trace_matrix[i][0] = 0

        for j in range(1, t_len):
            if j == 1:
                score_matrix[0][j] = self.gap_open
            else:
                score_matrix[0][j] = self.gap_open + (j - 1) * self.gap_extend
            trace_matrix[0][j] = 2

        for i in range(1, p_len):
            for j in range(1, t_len):
                # 每次判定是否为extend
                # 0: up
                if trace_matrix[i-1][j] == 0:
                    up = score_matrix[i-1][j] + self.gap_extend
                else:
                    up = score_matrix[i-1][j] + self.gap_open

                # 1: diag
                if p[i-1] == t[j-1]:
                    diag = score_matrix[i-1][j-1] + self.match
                else:
                    diag = score_matrix[i-1][j-1] + self.mismatch

                # 2: left
                if trace_matrix[i][j-1] == 2:
                    left = score_matrix[i][j-1] + self.gap_extend
                else:
                    left = score_matrix[i][j-1] + self.gap_open

                scores = [up, diag, left]
                score = max(scores)
                trace = scores.index(score)

                score_matrix[i][j] = score
                trace_matrix[i][j] = trace

        return score_matrix, trace_matrix

    def get_score(self, score_matrix):
        """返回比对得分"""
        return score_matrix[len(score_matrix)-1][len(score_matrix[0])-1]

    def trace(self, p, t, trace_matrix):
        """根据trace_matrix输出比对结果
        """
        p = '-' + p
        t = '-' + t
        seq_p, seq_m, seq_t = '', '', ''

        i, j = len(trace_matrix)-1, len(trace_matrix[0])-1
        while i > 0 or j > 0:
            if trace_matrix[i][j] == 0:
                cur_p = p[i]
                cur_t = '-'
                cur_m = ' '
                i -= 1
            elif trace_matrix[i][j] == 1:
                cur_p = p[i]
                cur_t = t[j]
                if cur_p == cur_t:
                    cur_m = '|'
                else:
                    cur_m = ' '
                i, j = i-1, j-1
            elif trace_matrix[i][j] == 2:
                cur_p = '-'
                cur_t = t[j]
                cur_m = ' '
                j -= 1

            seq_p = cur_p + seq_p
            seq_t = cur_t + seq_t
            seq_m = cur_m + seq_m
        return seq_p, seq_m, seq_t

    def plot_trace(self):
        """plot trace path"""
        p = self.p
        t = self.t
        trace_matrix = self.trace_matrix

        fig, ax = plt.subplots()
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        df = pd.DataFrame(trace_matrix, columns=list(" " + t))
        color = [['w'] * (len(t) + 1) for _ in range(len(p) + 1)]

        # 标记数据来源
        i, j = len(trace_matrix) - 1, len(trace_matrix[0]) - 1
        color[i][j] = 'g'    # g for green
        while i > 0 or j > 0:
            if trace_matrix[i][j] == 0:
                color[i - 1][j] = 'g'
                i -= 1
            elif trace_matrix[i][j] == 1:
                color[i - 1][j - 1] = 'g'
                i, j = i - 1, j - 1
            elif trace_matrix[i][j] == 2:
                color[i][j - 1] = 'g'
                j -= 1

        ax.table(cellText=df.values, colLabels=df.columns,
                 loc='center', rowLabels=" " + p,
                 cellColours=color, colWidths=[0.05 for x in df.columns], )
        fig.tight_layout()
        plt.show()

    def get_mismatch(self, middle):
        """计算p和t的错配数"""
        n = 0
        for i in middle:
            if i == ' ':
                n += 1
        return n


class LocalAlign:
    def __init__(self, p, t, match=5, mismatch=-4, gap_open=-10, gap_extend=-0.5):
        """
        end_weight控制从score_matrix哪个位置进行回溯。
        - True表示从右下角回溯，
        - False表示从整个score_matrix进行回溯
        """
        self.p = p
        self.t = t
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.score_matrix, self.trace_matrix = self.score_matrix(p, t)
        self.top, self.middle, self.bottom, self.p_loc, self.t_loc = self.trace(p, t, self.trace_matrix)
        self.score = self.get_score(self.score_matrix)
        self.mismatch = self.get_mismatch(self.middle, p)

    @property
    def align(self):
        print('Align Seqs:')
        print('\tp:\t' + self.top)
        print('\t\t' + self.middle)
        print('\tt:\t' + self.bottom)

    def score_matrix(self, p, t):
        """计算得分矩阵
        行为primer or pattern， 列为genome or t。
        """

        p_len = len(p) + 1
        t_len = len(t) + 1
        score_matrix = np.zeros(shape=(p_len, t_len))
        trace_matrix = np.zeros(shape=(p_len, t_len), dtype=int)  # 0: up, 1:diag, 2:left

        # 初始化score_matrix第一行第一列
        # i for row(primer), j for col(genome)
        for i in range(1, p_len):
            if i == 1:
                score_matrix[i][0] = self.gap_open
            else:
                score_matrix[i][0] = self.gap_open + (i - 1) * self.gap_extend
            trace_matrix[i][0] = 0

        for j in range(1, t_len):
            if j == 1:
                score_matrix[0][j] = self.gap_open
            else:
                score_matrix[0][j] = self.gap_open + (j - 1) * self.gap_extend
            trace_matrix[0][j] = 2

        for i in range(1, p_len):
            for j in range(1, t_len):
                # 每次判定是否为extend
                # 0: up
                if trace_matrix[i-1][j] == 0:
                    up = score_matrix[i-1][j] + self.gap_extend
                else:
                    up = score_matrix[i-1][j] + self.gap_open

                # 1: diag
                if p[i-1] == t[j-1]:
                    diag = score_matrix[i-1][j-1] + self.match
                else:
                    diag = score_matrix[i-1][j-1] + self.mismatch

                # 2: left
                if trace_matrix[i][j-1] == 2:
                    left = score_matrix[i][j-1] + self.gap_extend
                else:
                    left = score_matrix[i][j-1] + self.gap_open

                scores = [up, diag, left, 0]
                score = max(scores)
                trace = scores.index(score)  # index 3 for stop

                score_matrix[i][j] = score
                trace_matrix[i][j] = trace

        return score_matrix, trace_matrix

    def get_max_score_index(self, score_matrix):
        """返回最高分索引"""
        max_score = -1
        for i in range(len(score_matrix)):
            for j in range(len(score_matrix[0])):
                score = score_matrix[i][j]
                if score > max_score:
                    max_score, max_i, max_j = score, i, j
        return max_i, max_j

    def get_score(self, score_matrix):
        """返回最高得分"""
        i, j = self.get_max_score_index(score_matrix)
        return score_matrix[i][j]

    def trace(self, p, t, trace_matrix):
        """根据trace_matrix输出比对结果
        """

        p = '-' + p
        t = '-' + t
        seq_p, seq_m, seq_t = '', '', ''
        score_matrix = self.score_matrix
        i, j = self.get_max_score_index(score_matrix)
        p_end, t_end = i - 1, j - 1
        while i > 0 or j > 0:
            if trace_matrix[i][j] == 0:
                cur_p = p[i]
                cur_t = '-'
                cur_m = ' '
                i -= 1
            elif trace_matrix[i][j] == 1:
                cur_p = p[i]
                cur_t = t[j]
                if cur_p == cur_t:
                    cur_m = '|'
                else:
                    cur_m = ' '
                i, j = i - 1, j - 1
            elif trace_matrix[i][j] == 2:
                cur_p = '-'
                cur_t = t[j]
                cur_m = ' '
                j -= 1
            elif trace_matrix[i][j] == 3: # stop
                p_start, t_start = i , j  # 等于3时，不进行匹配
                p_loc = [p_start, p_end]
                t_loc = [t_start, t_end]
                return seq_p, seq_m, seq_t, p_loc, t_loc

            seq_p = cur_p + seq_p
            seq_t = cur_t + seq_t
            seq_m = cur_m + seq_m

    def plot_trace(self):
        """plot trace path"""
        p = self.p
        t = self.t
        trace_matrix = self.trace_matrix
        score_matrix = self.score_matrix

        fig, ax = plt.subplots()
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        df = pd.DataFrame(trace_matrix, columns=list(" " + t))
        color = [['w'] * (len(t) + 1) for _ in range(len(p) + 1)]

        # 标记数据来源
        i, j = self.get_max_score_index(score_matrix)
        color[i][j] = 'g'
        while i > 0 or j > 0:
            if trace_matrix[i][j] == 0:
                color[i - 1][j] = 'g'
                i -= 1
            elif trace_matrix[i][j] == 1:
                color[i - 1][j - 1] = 'g'
                i, j = i - 1, j - 1
            elif trace_matrix[i][j] == 2:
                color[i][j - 1] = 'g'
                j -= 1
            elif trace_matrix[i][j] == 3:
                break

        ax.table(cellText=df.values, colLabels=df.columns,
                 loc='center', rowLabels=" " + p,
                 cellColours=color, colWidths=[0.05 for x in df.columns], )
        fig.tight_layout()
        plt.show()

    def get_mismatch(self, middle, p):
        """计算p和t的错配数"""
        n = 0
        for i in middle:
            if i == '|':
                n += 1
        return len(p)-n



if __name__ == '__main__':
    # p = 'ATCGTTCC'
    # t = 'AAATCGCGTTCCGGG'
    # test = GlobalAlign(p, t)
    # print(test.align)
    # print(test.score_matrix)
    # test.plot_trace()
    # print(test.mismatch)

    p = 'ATCGTTCC'
    t = 'AAATCGCGTTCCGGG'
    test = LocalAlign(p, t)
    print(test.align)
    # print(test.score_matrix)
    print(test.score)
    # test.plot_trace()
    print(test.p_loc)
    print(test.t_loc)
    # print(test.score_matrix)
    # print(test.score_matrix[8][12])
    # test.plot_trace()
    # print(test.middle)
    # print(test.mismatch)


