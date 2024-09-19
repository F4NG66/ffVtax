import csv
import pandas as pd
import sys
csv.field_size_limit(sys.maxsize)
def calculate_scores(input_csv, output_csv, jacc_threshold=0.95, qcov_threshold=0.8):
    """
    读取输入的 CSV 文件，计算 Jaccard 指数和查询覆盖率，并将结果写入新的 CSV 文件。
    同时根据用户定义的阈值进行过滤。

    :param input_csv: 输入的包含匹配结果的 CSV 文件。
    :param output_csv: 输出包含 Jacc 和 Qcov 的新 CSV 文件。
    :param jacc_threshold: Jacc 的过滤阈值，默认值为 0.95。
    :param qcov_threshold: Qcov 的过滤阈值，默认值为 0.8。
    """
    data = []

    # 读取输入的 CSV 并计算 Jacc 和 Qcov
    with open(input_csv, mode='r', newline='') as infile:
        reader = csv.reader(infile)

        # 读取表头
        header = next(reader)
        header.extend(['Jacc', 'Qcov'])  # 添加 Jacc 和 Qcov 列

        # 遍历每一行并计算 Jacc 和 Qcov
        for row in reader:
            sequence_name = row[0]
            total_kmers = int(row[1])  # 查询序列的总 k-mers 数量
            matched_kmers = int(row[3])  # 匹配到的 k-mers 数量
            reference_total_kmers = int(row[4])  # 参考序列的总 k-mers 数量

            query_kmers = set(row[6].split(','))  # 查询序列的 k-mers 集合
            reference_kmers = set(row[7].split(','))  # 参考序列的 k-mers 集合

            # 计算 Jacc 和 Qcov
            if len(query_kmers) > 0 and reference_total_kmers > 0:
                jacc = len(query_kmers & reference_kmers) / len(query_kmers | reference_kmers)  # Jaccard 指数
                qcov = len(query_kmers & reference_kmers) / len(reference_kmers) # Query Coverage
            else:
                jacc = 0
                qcov = 0

            # 将结果存储在列表中
            data.append([sequence_name, total_kmers, row[2], matched_kmers, reference_total_kmers, row[5], jacc, qcov])

    # 转换为 pandas DataFrame
    df = pd.DataFrame(data, columns=['Sequence Name', 'Total Input k-mers', 'Reference', 'Matched k-mers', 'Reference Total k-mers', 'GCA Name', 'Jacc', 'Qcov'])

    # 过滤结果，保留满足阈值的行
    filtered_results = df[(df['Jacc'] >= jacc_threshold) & (df['Qcov'] >= qcov_threshold)]

    # 按照 Sequence Name（#query），然后按 Jacc 和 Qcov 排序
    filtered_results = filtered_results.sort_values(by=['Sequence Name', 'Jacc', 'Qcov'], ascending=[True, False, False])

    # 去除重复的 Sequence Name，只保留每个 query 的最佳匹配
    filtered_result = filtered_results.drop_duplicates(subset=['Sequence Name'], keep='first')

    # 将过滤后的结果写入输出 CSV 文件
    filtered_result.to_csv(output_csv, index=False)

    print(f"Filtered results saved to {output_csv}")
