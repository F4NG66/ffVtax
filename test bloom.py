import unittest
from bloom_filter_handler import BloomFilterHandler
from sequence_matcher import SequenceMatcher
from io import StringIO
import csv
import os
from scoring import calculate_scores  # 导入新添加的 scoring 模块

class TestBloomFilterHandler(unittest.TestCase):

    def setUp(self):
        """
        初始化 BloomFilterHandler 用于测试。
        """
        self.bloom_handler = BloomFilterHandler(kmer_size=5, factor=4)

        # 模拟 Bloom filter 和参考序列映射
        self.bloom_handler.bloom = set(['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'])
        self.bloom_handler.kmer_to_reference = {
            'AAAAA': {'reference1'},
            'CCCCC': {'reference1'},
            'GGGGG': {'reference2'},
            'TTTTT': {'reference3'}
        }

    def test_extract_kmers(self):
        sequence = 'AAAAACCCCCTTTTTGGGGG'
        expected_kmers = ['AAAAA', 'AAAAC', 'AAACC', 'AACCC', 'ACCCC', 'CCCCC', 'CCCCT', 'CCCTT', 'CCTTT', 'CTTTT', 'TTTTT', 'TTTTG', 'TTTGG', 'TTGGG', 'TGGGG', 'GGGGG']
        kmers = self.bloom_handler.extract_kmers(sequence)
        self.assertEqual(kmers, expected_kmers)

class TestSequenceMatcher(unittest.TestCase):

    def setUp(self):
        """
        初始化 SequenceMatcher 用于测试。
        """
        self.bloom_handler = BloomFilterHandler(kmer_size=5, factor=4)
        self.bloom_handler.bloom = set(['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT'])
        self.bloom_handler.kmer_to_reference = {
            'AAAAA': {'reference1'},
            'CCCCC': {'reference1'},
            'GGGGG': {'reference2'},
            'TTTTT': {'reference3'}
        }

        self.matcher = SequenceMatcher(self.bloom_handler)

    def test_match_sequence(self):
        query_sequence = 'AAAAACCCCCTTTTTGGGGG'
        matches = self.matcher._match_sequence(query_sequence)
        expected_matches = {'reference1', 'reference2', 'reference3'}
        self.assertEqual(matches, expected_matches)

    def test_match_sequences_from_file(self):
        fasta_content = """>seq1
        AAAAA
        >seq2
        CCCCC
        >seq3
        GGGGG
        """

        fasta_file = StringIO(fasta_content)

        matched_references = self.matcher.match_sequences(fasta_file)
        expected_references = {
            'seq1': {'reference1'},
            'seq2': {'reference1'},
            'seq3': {'reference2'}
        }

        self.assertEqual(matched_references, expected_references)

    def test_output_to_csv(self):
        """
        测试将匹配结果输出到 CSV 文件。
        """
        matched_references = {
            'seq1': {'total_kmers': 150, 'matched_references': {'reference1': 50, 'reference2': 30}, 'reference_gca_names': {'reference1': 'GCA_000001', 'reference2': 'GCA_000002'}},
            'seq2': {'total_kmers': 200, 'matched_references': {'reference3': 100}, 'reference_gca_names': {'reference3': 'GCA_000003'}}
        }

        output_file = 'test_output.csv'

        # 输出到 CSV 文件
        self.matcher.output_to_csv(matched_references, output_file)

        # 读取 CSV 文件并验证内容
        with open(output_file, newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)

        expected_rows = [
            ['Sequence Name', 'Total k-mers', 'Reference', 'Matched k-mers', 'Reference Total k-mers', 'GCA Name'],
            ['seq1', '150', 'reference1', '50', '200', 'GCA_000001'],
            ['seq1', '150', 'reference2', '30', '180', 'GCA_000002'],
            ['seq2', '200', 'reference3', '100', '300', 'GCA_000003']
        ]

        self.assertEqual(rows, expected_rows)

        # 清理测试生成的文件
        os.remove(output_file)

class TestScoring(unittest.TestCase):

    def test_calculate_scores(self):
        """
        测试计算 Jacc 和 Qcov 分数并将结果写入新文件。
        """
        input_content = """Sequence Name,Total k-mers,Reference,Matched k-mers,Reference Total k-mers,GCA Name
        seq1,150,reference1,50,200,GCA_000001
        seq1,150,reference2,30,180,GCA_000002
        seq2,200,reference3,100,300,GCA_000003
        """

        input_file = 'test_input.csv'
        output_file = 'test_scored_output.csv'

        # 写入模拟的输入 CSV 文件
        with open(input_file, mode='w', newline='') as f:
            f.write(input_content)

        # 调用 scoring 模块中的 calculate_scores 函数
        calculate_scores(input_file, output_file)

        # 验证输出文件的内容
        with open(output_file, newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)

        expected_rows = [
            ['Sequence Name', 'Total k-mers', 'Reference', 'Matched k-mers', 'Reference Total k-mers', 'GCA Name', 'Jacc', 'Qcov'],
            ['seq1', '150', 'reference1', '50', '200', 'GCA_000001', '0.2500', '0.3333'],
            ['seq1', '150', 'reference2', '30', '180', 'GCA_000002', '0.1538', '0.2000'],
            ['seq2', '200', 'reference3', '100', '300', 'GCA_000003', '0.3333', '0.5000']
        ]

        self.assertEqual(rows, expected_rows)

        # 清理测试生成的文件
        os.remove(input_file)
        os.remove(output_file)

if __name__ == '__main__':
    unittest.main()
