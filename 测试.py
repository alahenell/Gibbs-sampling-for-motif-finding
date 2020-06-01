#!/usr/bin/python3
# !/usr/bin/python3
import random

# 读fa文件写入raw_sequences字典
raw_seqences = {}
head = ''
seq = ''
file = open('motif_crp.fa', 'r')
for line in file:
    if line[0] == '>' and seq != '':
        seq = ''
    if line[0] == '>':  # 取每条序列的注释为键
        a = line[1:]
        head = a.strip()
    else:
        seq = seq + line
    raw_seqences[head] = seq.strip()  # 每条序列为raw_sequence字典的值
file.close()

# 获得反向互补序列字典
rev_com_sequences = {}
for key in raw_seqences:
    seq = raw_seqences[key]
    regx = str.maketrans('ACGT', 'TGCA')  # 制作翻译规则
    com_seq = seq.translate(regx)  # 翻译，得到互补DNA序列
    head = 'rev_com_' + key
    rev_com_sequences[head] = com_seq[::-1]  # 序列反向
# 将raw_sequcens与rev_com_sequences合并为一个字典，便于后续操作
sequences = {**raw_seqences, **rev_com_sequences}

# 存储每条序列得分的字典
score = {}
# 存储每条序列的fragment的字典
fragment_matrix = {}
# fragment矩阵的每种核苷酸的每个位置的数量作为打分矩阵
A = [0, 0, 0, 0, 0, 0, 0, 0]
C = [0, 0, 0, 0, 0, 0, 0, 0]
G = [0, 0, 0, 0, 0, 0, 0, 0]
T = [0, 0, 0, 0, 0, 0, 0, 0]

# 得到每个序列，长度为8的随机片段
for key in sequences:
    random.seed()
    num = random.randint(0, len(sequences[key]) - 8)  # 取一个随机位点
    head = key
    seq = sequences[key]
    fragment = seq[num:num + 8]  # 取一段长度为8bp的随机segment
    fragment_matrix[head] = fragment


# 制作并更新得分矩阵
def update_matrix():
    global A, C, G, T
    for key in fragment_matrix:  # key以fragment_matrix的键循环，循环36次
        j = 0
        for i in fragment_matrix[key]:  # 以某条序列的segment的8个position依次循环
            if i == 'A':
                A[j] = A[j] + 1/36
            elif i == 'C':
                C[j] = C[j] + 1/36
            elif i == 'G':
                G[j] = G[j] + 1/36
            elif i == 'T':
                T[j] = T[j] + 1/36
            j = j + 1

    # 根据得分更新每条序列的随机片段，并记录这些片段的位置

# 更新矩阵并打分
def update_fragment():
    global score
    random.seed()
    a = random.choice(list(sequences))  # 随机选取sequences字典中键并赋予a
    seq = sequences[a]
    i = 0
    score_list = []
    while i <= len(seq) - 8:
        p = 1
        fragment = seq[i:i + 8]  # 以8bp大小的窗口滑动
        k = 0
        for j in fragment:
            if j == 'A':
                p = p * A[k]
            elif j == 'C':
                p = p * C[k]
            elif j == 'G':
                p = p * G[k]
            elif j == 'T':
                p = p * T[k]
            k = k + 1
        score_list.append(p)
        i = i + 1
    # 根据滑动窗口法，算出每个片段的得分
    score[a] = score_list  # 每条序列的所有片段的得分
    index = score_list.index(max(score_list))  # 得到最高得分的索引
    fragment = seq[index:index + 8]
    fragment_matrix[a] = fragment  # 更新每条序列的随机片段


# 主程序
update_matrix()
m = 0
while m < 1000:  # 循环1000次
    update_fragment()  # 更新片段
    A = [0, 0, 0, 0, 0, 0, 0, 0]  # 将各个碱基每个位置的得分归零
    C = [0, 0, 0, 0, 0, 0, 0, 0]
    G = [0, 0, 0, 0, 0, 0, 0, 0]
    T = [0, 0, 0, 0, 0, 0, 0, 0]
    update_matrix()  # 更新矩阵
    m += 1
# 输出最后的打分矩阵
print("pos"+"\t"+"A"+"\t"+"C"+"\t"+"G"+"\t"+"T")
print("1"+"\t"+'%.3f'% A[0]+"\t"+'%.3f'% C[0]+"\t"+'%.3f'% G[0]+"\t"+'%.3f'% T[0])
print("1"+"\t"+'%.3f'% A[1]+"\t"+'%.3f'% C[1]+"\t"+'%.3f'% G[1]+"\t"+'%.3f'% T[1])
print("1"+"\t"+'%.3f'% A[2]+"\t"+'%.3f'% C[2]+"\t"+'%.3f'% G[2]+"\t"+'%.3f'% T[2])
print("1"+"\t"+'%.3f'% A[3]+"\t"+'%.3f'% C[3]+"\t"+'%.3f'% G[3]+"\t"+'%.3f'% T[3])
print("1"+"\t"+'%.3f'% A[4]+"\t"+'%.3f'% C[4]+"\t"+'%.3f'% G[4]+"\t"+'%.3f'% T[4])
print("1"+"\t"+'%.3f'% A[5]+"\t"+'%.3f'% C[5]+"\t"+'%.3f'% G[5]+"\t"+'%.3f'% T[5])
print("1"+"\t"+'%.3f'% A[6]+"\t"+'%.3f'% C[6]+"\t"+'%.3f'% G[6]+"\t"+'%.3f'% T[6])
print("1"+"\t"+'%.3f'% A[7]+"\t"+'%.3f'% C[7]+"\t"+'%.3f'% G[7]+"\t"+'%.3f'% T[7])

# 画图
import matplotlib.pyplot as plt
def plot_histogram(seq_name):
    global score
    num_list = score[seq_name]  # 将得分储存在num_list里
    plt.bar(range(len(num_list)), num_list, color='g')  # 用range函数产生一个序列作为x轴的位置序列，num_list作为柱形图的高度
    plt.title(seq_name)


# 画出cole1和反向互补序列中分数大的直方图
a = max(score['cole1'])
b = max(score['rev_com_cole1'])
plt.subplot(211)
if a > b:
    plot_histogram('cole1')
else:
    plot_histogram('rev_com_cole1')
# 画出ecoarabop和反向互补序列中分数大的直方图
c = max(score['ecoarabop'])
d = max(score['rev_com_ecoarabop'])
plt.subplot(212)
if a > b:
    plot_histogram('ecoarabop')
else:
    plot_histogram('rev_com_ecoarabop')
plt.show()
