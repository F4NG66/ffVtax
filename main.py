import argparse
import os
from bloom_filter_handler import BloomFilterHandler
from sequence_matcher import SequenceMatcher
from scoring import calculate_scores
from gca_to_taxid import add_taxid_to_scored_output

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Run viral sequence matching, scoring, and taxid lookup.")
    parser.add_argument('-i', '--input', required=True, help="Input query sequence file (FASTA format)")
    parser.add_argument('-d', '--database', required=True, help="Path to the reference database directory")
    parser.add_argument('-o', '--output', required=True, help="Output directory for results")
    parser.add_argument('--kmer_size', type=int, default=21, help="k-mer size (default: 21)")
    parser.add_argument('--error_rate', type=float, default=0.001, help="Error rate for the Bloom filter (default: 0.001)")
    parser.add_argument('--jacc_threshold', type=float, default=0.95, help="Jaccard index threshold (default: 0.95)")
    parser.add_argument('--qcov_threshold', type=float, default=0.8, help="Query coverage threshold (default: 0.8)")

    args = parser.parse_args()

    # 确保输出目录存在
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # 初始化 BloomFilterHandler 并加载数据库
    print("Loading reference database...")
    bloom_handler = BloomFilterHandler(kmer_size=args.kmer_size, factor=4)
    bloom_handler.error_rate = args.error_rate  # 设置用户指定的错误率
    bloom_handler.load_database(args.database)
    print("Database loaded successfully.")

    # 初始化 SequenceMatcher 并进行序列匹配
    print(f"Matching sequences from {args.input}...")
    matcher = SequenceMatcher(bloom_handler)
    matched_references = matcher.match_sequences(args.input)
    print("Sequence matching completed.")

    # 输出匹配结果到 CSV 文件
    matching_result_file = os.path.join(args.output, "matching_results.csv")
    print(f"Writing matching results to {matching_result_file}...")
    matcher.output_to_csv(matched_references, matching_result_file)
    print("Matching results written to CSV.")

    # 计算 Jaccard 和 Qcov 分数，并输出到 scored_output.csv 文件

    scored_output_file = os.path.join(args.output, "scored_output.csv")
    print(f"Calculating Jaccard and Query Coverage scores...")
    calculate_scores(matching_result_file, scored_output_file,
                 jacc_threshold=args.jacc_threshold,
                 qcov_threshold=args.qcov_threshold)
    print(f"Scores added and saved to {scored_output_file} and filtered results saved to output directory.")

    # 将 TaxID 添加到 scored_output 文件
    taxid_map_file = "taxid.map"  # 假设 taxid.map 文件固定位置
    output_with_taxid = os.path.join(args.output, "scored_output_with_taxid.csv")
    print(f"Adding TaxID information based on GCA from {taxid_map_file}...")
    add_taxid_to_scored_output(scored_output_file, output_with_taxid, taxid_map_file)
    print(f"TaxID column added to {output_with_taxid}")



if __name__ == "__main__":
    main()
