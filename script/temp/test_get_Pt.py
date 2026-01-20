"""Batch fetch sample info (alias/title) from SRA Experiment accessions (e.g., SRX...).

给定 SRX 号（SRA Experiment），到 NCBI E-utilities 查询其 SAMPLE 信息，
提取 sample alias / title，方便批量映射到 GEO 页面显示的样本名（如 Pt1_Pre_AD101148-6）。

Usage examples:
	python test_get_Pt.py SRX2405808 SRX123456
	python test_get_Pt.py --from-file ids.txt --out samples.csv
	python test_get_Pt.py SRX2405808 --delay 0.4 --api-key YOUR_NCBI_KEY

输出默认 CSV 到 stdout（或 --out 指定文件）。
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import requests


NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class SampleFetchError(RuntimeError):
	"""Raised when fetching/parsing SRA sample info fails."""


def fetch_sample_info(
	accession: str,
	api_key: Optional[str] = None,
	email: Optional[str] = None,
	timeout: float = 20.0,
) -> Dict[str, Optional[str]]:
	"""Fetch sample alias/title for an SRX accession using NCBI efetch.

	Returns dict with keys: accession, sample_accession, sample_alias, sample_title.
	"""

	params = {
		"db": "sra",
		"id": accession,
		"retmode": "xml",
	}
	if api_key:
		params["api_key"] = api_key
	if email:
		params["email"] = email

	try:
		resp = requests.get(NCBI_EFETCH_URL, params=params, timeout=timeout)
	except requests.RequestException as exc:
		raise SampleFetchError(f"Request failed for {accession}: {exc}") from exc

	if resp.status_code != 200:
		raise SampleFetchError(
			f"HTTP {resp.status_code} for {accession}: {resp.text[:200]}"
		)

	try:
		root = ET.fromstring(resp.text)
	except ET.ParseError as exc:
		raise SampleFetchError(f"XML parse error for {accession}: {exc}") from exc

	# SRA XML structure: EXPERIMENT_PACKAGE_SET / EXPERIMENT_PACKAGE / SAMPLE
	pkg = root.find(".//EXPERIMENT_PACKAGE")
	sample_elem = pkg.find("SAMPLE") if pkg is not None else None
	if sample_elem is None:
		raise SampleFetchError(f"No SAMPLE element found for {accession}")

	sample_accession = None
	ids_elem = sample_elem.find("IDENTIFIERS")
	if ids_elem is not None:
		prim = ids_elem.find("PRIMARY_ID")
		if prim is not None and prim.text:
			sample_accession = prim.text.strip()

	sample_alias = sample_elem.attrib.get("alias")

	title_elem = sample_elem.find("TITLE")
	sample_title = title_elem.text.strip() if title_elem is not None and title_elem.text else None

	return {
		"accession": accession,
		"sample_accession": sample_accession,
		"sample_alias": sample_alias,
		"sample_title": sample_title,
	}


def iter_accessions(args_ids: List[str], from_file: Optional[Path]) -> List[str]:
	ids: List[str] = []
	if args_ids:
		ids.extend(args_ids)
	if from_file:
		with from_file.open() as fh:
			for line in fh:
				cleaned = line.strip()
				if cleaned:
					ids.append(cleaned)
	return ids


def write_csv(rows: Iterable[Dict[str, Optional[str]]], out_path: Optional[Path]):
	fieldnames = ["accession", "sample_accession", "sample_alias", "sample_title"]
	if out_path:
		out_path.parent.mkdir(parents=True, exist_ok=True)
		outfile = out_path.open("w", newline="")
		close_needed = True
	else:
		outfile = sys.stdout
		close_needed = False

	writer = csv.DictWriter(outfile, fieldnames=fieldnames)
	writer.writeheader()
	for row in rows:
		writer.writerow(row)

	if close_needed:
		outfile.close()


def main(argv: Optional[List[str]] = None) -> int:
	parser = argparse.ArgumentParser(
		description="Batch fetch SRA sample alias/title from SRX accessions via NCBI efetch",
	)
	parser.add_argument("accessions", nargs="*", help="SRX accessions (e.g., SRX2405808)")
	parser.add_argument("--from-file", type=Path, help="Text file with one SRX accession per line")
	parser.add_argument("--api-key", dest="api_key", help="NCBI API key to increase rate limit")
	parser.add_argument("--email", help="Contact email for NCBI etiquette")
	parser.add_argument("--delay", type=float, default=0.35, help="Sleep seconds between requests")
	parser.add_argument("--timeout", type=float, default=20.0, help="HTTP timeout seconds")
	parser.add_argument("--out", type=Path, help="Output CSV file (default stdout)")

	args = parser.parse_args(argv)

	ids = iter_accessions(args.accessions, args.from_file)
	if not ids:
		parser.error("Please provide SRX accessions as arguments or via --from-file")

	results = []
	for i, acc in enumerate(ids, 1):
		try:
			info = fetch_sample_info(
				acc,
				api_key=args.api_key,
				email=args.email,
				timeout=args.timeout,
			)
			results.append(info)
			print(
				f"[{i}/{len(ids)}] {acc} -> alias={info['sample_alias']} title={info['sample_title']}",
				file=sys.stderr,
			)
		except SampleFetchError as exc:
			print(f"[{i}/{len(ids)}] {acc} failed: {exc}", file=sys.stderr)
			results.append(
				{
					"accession": acc,
					"sample_accession": None,
					"sample_alias": None,
					"sample_title": f"ERROR: {exc}",
				}
			)
		time.sleep(max(args.delay, 0))

	write_csv(results, args.out)
	return 0


if __name__ == "__main__":
	sys.exit(main())

