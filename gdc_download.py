#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GDC TCGA manifest 下载 + 表达矩阵合并工具
用法：直接把 manifest txt 文件拖到 exe（或脚本）上即可
输出：<exe所在目录>/dataset_download/
    ├── raw/                原始 tsv 文件
    ├── expression_matrix_TPM.csv      合并后的表达矩阵（默认 TPM）
    └── failed_ids.txt                 失败清单（若有）
"""

import os
import re
import sys
import time
import hashlib
import csv
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============ 配置区（按需修改）============
GDC_DATA_URL = "https://api.gdc.cancer.gov/data/"
THREADS = 4                  # 并发下载数
VALUE_COLUMN = "tpm_unstranded"  # 矩阵取值列，可选:
                                 # "unstranded"(原始count) / "tpm_unstranded" / "fpkm_unstranded" / "fpkm_uq_unstranded"
# ==========================================


def base_dir():
    """exe 双击/拖拽运行时，输出到 exe 所在目录；脚本运行时输出到脚本目录"""
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def get_input_file():
    """从命令行参数获取拖拽的 txt；没有则提示手动输入"""
    if len(sys.argv) >= 2 and os.path.isfile(sys.argv[1]):
        return sys.argv[1]
    print("未检测到拖拽文件，请手动输入 manifest txt 路径（可直接把文件拖进本窗口后回车）：")
    path = input("> ").strip().strip('"').strip("'")
    if os.path.isfile(path):
        return path
    print("路径无效！")
    return None


def parse_manifest(path):
    entries = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [l for l in f.read().strip().splitlines() if l.strip()]
    for line in lines[1:]:  # 跳过表头
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) >= 5:
            entries.append({"id": parts[0], "filename": parts[1],
                            "md5": parts[2], "size": int(parts[3])})
    return entries


def md5_file(path, chunk=1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_file(entry, out_dir, max_retries=5):
    url = GDC_DATA_URL + entry["id"]
    dest = os.path.join(out_dir, entry["filename"])
    tmp = dest + ".part"

    # 已下载完成则跳过
    if os.path.exists(dest) and os.path.getsize(dest) == entry["size"]:
        return entry["id"], "SKIP"

    for attempt in range(1, max_retries + 1):
        try:
            done = os.path.getsize(tmp) if os.path.exists(tmp) else 0
            headers = {"Range": f"bytes={done}-"} if done else {}
            with requests.get(url, headers=headers, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(tmp, "ab" if done else "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

            if os.path.getsize(tmp) != entry["size"]:
                raise IOError("文件大小不符")

            if md5_file(tmp) != entry["md5"]:
                raise ValueError("MD5 校验失败")  # MD5 错说明数据损坏，删掉重来
            # 上面校验失败时其实会走到 except 分支删除重下，逻辑在 except 里处理

            os.rename(tmp, dest)
            return entry["id"], "OK"

        except Exception as e:
            # MD5 失败（文件损坏）时清空 .part 从头下载；网络中断则续传
            if isinstance(e, ValueError) and os.path.exists(tmp):
                os.remove(tmp)
            time.sleep(min(2 * attempt, 10))

    return entry["id"], "FAILED"


def merge_matrices(raw_dir, out_csv):
    """把所有 *.rna_seq.augmented_star_gene_counts.tsv 合并为 gene x sample 矩阵"""
    files = [f for f in os.listdir(raw_dir)
             if f.endswith(".rna_seq.augmented_star_gene_counts.tsv")]
    if not files:
        print("没有找到可合并的 tsv 文件！")
        return

    gene_names = {}          # gene_id -> gene_name
    gene_order = []          # 保持文件中出现顺序
    matrix = {}              # gene_id -> {sample: value}
    samples = []             # 样本列顺序
    tcga_re = re.compile(r"TCGA-\w\w-\w{4}")

    total = len(files)
    for idx, fname in enumerate(sorted(files), 1):
        sample_m = tcga_re.search(fname)
        sample = sample_m.group(0) if sample_m else fname[:20]
        if sample in samples:
            sample = f"{sample}_{idx}"  # 防止重名
        samples.append(sample)

        path = os.path.join(raw_dir, fname)
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            rows = csv.reader(f, delimiter="\t")
            # 跳过前4行元信息（# 开头）
            for i, row in enumerate(rows):
                if i < 4:
                    continue
                if not row or row[0].startswith("#") or row[0] == "__para":
                    continue
                gid = row[0]
                if gid not in matrix:  # PAR_Y 重复基因只保留第一次出现
                    matrix[gid] = {}
                    gene_order.append(gid)
                    gene_names[gid] = row[1]
                # 定位列索引（第一个文件时按表头找）
                col = header_idx.get(VALUE_COLUMN, None)
                if col is None:
                    print(f"错误：找不到列 {VALUE_COLUMN}")
                    return
                matrix[gid][sample] = row[col]

        print(f"\r合并进度: {idx}/{total}", end="")

    # 写出 CSV
    out_path = os.path.join(base_dir(), out_csv)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene_id", "gene_name"] + samples)
        for gid in gene_order:
            w.writerow([gid, gene_names[gid]] +
                       [matrix[gid].get(s, "") for s in samples])

    print(f"\n✓ 表达矩阵已保存: {out_path}")
    print(f"  基因数: {len(gene_order)}, 样本数: {len(samples)}")


# 表头索引缓存（GDC 流水线所有文件格式一致，只需解析一次）
header_idx = {}


def build_header_index(raw_dir):
    """从第一个文件读取列名 -> 索引映射"""
    for f in os.listdir(raw_dir):
        if f.endswith(".rna_seq.augmented_star_gene_counts.tsv"):
            path = os.path.join(raw_dir, f)
            with open(path, "r", encoding="utf-8", errors="ignore") as fp:
                for i, line in enumerate(fp):
                    if i == 4:  # 第5行是列名行
                        cols = line.strip().split("\t")
                        header_idx.update({c: j for j, c in enumerate(cols)})
                        return
    print("错误：找不到列名行！")


def main():
    print("=" * 60)
    print("   GDC manifest 下载 + 表达矩阵合并工具")
    print("=" * 60)

    manifest = get_input_file()
    if not manifest:
        input("\n按回车退出...")
        return

    out_root = os.path.join(base_dir(), "dataset_download")
    raw_dir = os.path.join(out_root, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    entries = parse_manifest(manifest)
    print(f"共解析到 {len(entries)} 个文件\n")

    # ---------- 下载 ----------
    results = {}
    with ThreadPoolExecutor(max_workers=THREADS) as pool:
        futs = {pool.submit(download_file, e, raw_dir): e["id"] for e in entries}
        for n, fut in enumerate(as_completed(futs), 1):
            fid, status = fut.result()
            results[fid] = status
            print(f"\r[{n}/{len(entries)}] {status:6s} {fid}", end="")
    print()

    failed = [k for k, v in results.items() if v == "FAILED"]
    ok = [k for k, v in results.items() if v in ("OK", "SKIP")]
    print(f"\n下载完成: {len(ok)} 成功, {len(failed)} 失败")

    if failed:
        with open(os.path.join(out_root, "failed_ids.txt"), "w") as f:
            f.write("\n".join(failed))
        print("失败清单已保存到 failed_ids.txt（重新运行程序会自动续传）")

    # ---------- 合并 ----------
    if ok:
        build_header_index(raw_dir)
        merge_matrices(raw_dir, f"expression_matrix_{VALUE_COLUMN}.csv")

    input("\n全部完成，按回车退出...")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n发生错误: {e}")
        input("按回车退出...")
