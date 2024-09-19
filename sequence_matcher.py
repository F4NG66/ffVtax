import csv
from collections import defaultdict
import gzip

class SequenceMatcher:
    def __init__(self, bloom_filter_handler):
        self.bloom_filter_handler = bloom_filter_handler

    def match_sequences(self, query_file):
        """
        从输入文件中读取序列并匹配到参考序列。
        返回字典，记录每个序列的 k-mer 数量、匹配的参考序列以及匹配的 k-mer 数量。
        """
        stats = defaultdict(dict)

        if isinstance(query_file, str):
            open_func = gzip.open if query_file.endswith('.gz') else open
            with open_func(query_file, 'rt') as file:
                for sequence_name, sequence_data in self._read_fasta_file(file):
                    stats[sequence_name] = self._match_sequence(sequence_name, sequence_data)
        else:
            for sequence_name, sequence_data in self._read_fasta_file(query_file):
                stats[sequence_name] = self._match_sequence(sequence_name, sequence_data)

        return stats

    def _read_fasta_file(self, file):
        sequence_name = None
        sequence_data = []

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    yield sequence_name, ''.join(sequence_data)
                sequence_name = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)

        if sequence_name is not None:
            yield sequence_name, ''.join(sequence_data)

    def _match_sequence(self, sequence_name, sequence_data):
        """
        计算输入序列的 k-mer 数量，匹配到的参考序列及其匹配的 k-mer 数量。
        """
        input_kmers = set(self._extract_kmers(sequence_data))  # 使用集合存储唯一的 k-mers
        reference_kmer_matches = defaultdict(set)  # 用集合存储匹配的参考 k-mers

        for kmer in input_kmers:
            if kmer in self.bloom_filter_handler.bloom:
                for ref in self.bloom_filter_handler.kmer_to_reference[kmer]:
                    reference_kmer_matches[ref].add(kmer)
        print(f"Query k-mers: {len(input_kmers)}, Reference k-mers matches: {[len(kmers) for kmers in reference_kmer_matches.values()]}")
        return input_kmers, reference_kmer_matches


    def _extract_kmers(self, sequence):
        """
        从序列中提取 k-mers。
        """
        return [sequence[i:i + self.bloom_filter_handler.kmer_size] for i in range(len(sequence) - self.bloom_filter_handler.kmer_size + 1)]

    def _calculate_total_kmers(self, ref):
        """
        根据参考序列的名称计算其 k-mers 总数。
        需要根据参考序列名从数据库中提取序列数据。
        """
        return self.bloom_filter_handler.reference_kmer_count.get(ref, 0)  # 计算 k-mers 数量


    def output_to_csv(self, matched_references, output_file):
        """
        将匹配结果和统计信息导出到 CSV 文件，包含每个参考序列的 k-mer 数量。
        """
        with open(output_file, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Sequence Name', 'Total Input k-mers', 'Reference', 'Matched k-mers', 'Reference Total k-mers', 'GCA Name', 'Query k-mers', 'Reference k-mers'])

            for sequence_name, (query_kmers, reference_kmer_matches) in matched_references.items():
                for ref, matched_kmers in reference_kmer_matches.items():
                    reference_total_kmers = self.bloom_filter_handler.reference_kmer_count.get(ref, 0)
                    gca_name = self.bloom_filter_handler.reference_gca_map.get(ref, 'Unknown')
                    reference_kmers = self.bloom_filter_handler.get_reference_kmers(ref)
                    writer.writerow([sequence_name, len(query_kmers), ref, len(matched_kmers), reference_total_kmers, gca_name, ','.join(query_kmers), ','.join(reference_kmers)])
