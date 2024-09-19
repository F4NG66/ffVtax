import csv

def load_taxid_map(taxid_map_file):
    """
    加载 taxid.map 文件并返回 GCA 到 TaxID 的映射字典。
    :param taxid_map_file: taxid.map 文件路径。
    :return: 包含 GCA 到 TaxID 映射的字典。
    """
    gca_to_taxid = {}
    with open(taxid_map_file, mode='r') as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                gca, taxid = parts
                gca_to_taxid[gca] = taxid
    return gca_to_taxid


def add_taxid_to_scored_output(scored_output_file, output_file, taxid_map_file):
    """
    读取 scored_output 文件，使用 taxid.map 文件将 GCA 转换为 TaxID，并将其添加到新文件中。
    :param scored_output_file: 输入的包含 GCA 的 scored_output 文件路径。
    :param output_file: 输出包含 taxid 的新文件路径。
    :param taxid_map_file: GCA 到 TaxID 的转换文件。
    """
    # 加载 GCA 到 TaxID 的映射表
    gca_to_taxid = load_taxid_map(taxid_map_file)

    # 打开 scored_output 文件并添加 taxid 列
    with open(scored_output_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # 读取表头并添加新的 TaxID 列
        header = next(reader)
        header.append('TaxID')
        writer.writerow(header)

        # 遍历每一行并添加对应的 TaxID
        for row in reader:
            gca = row[5]  # 假设 GCA 号在第3列 (索引2)
            taxid = gca_to_taxid.get(gca, 'N/A')  # 如果找不到对应的 TaxID，标记为 'N/A'
            row.append(taxid)
            writer.writerow(row)

    print(f"TaxID column added to {output_file}")


