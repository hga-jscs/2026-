"""并行分块下载 AC 所需的参考数据归档，并支持断点续传与尺寸校验。

仅负责下载归档；解包和 required-file 校验由 ``prepare_references.ps1`` 完成。
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import threading
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


BASE_URL = "https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect"
REFERENCES = ("GRCh37", "hg19", "GRCh38_viral", "mm10")
REFERENCE_SIZES = {
    "GRCh37": 1_163_687_112,
    "hg19": 1_152_267_451,
    "GRCh38_viral": 1_100_983_085,
    "mm10": 1_108_257_901,
}
USER_AGENT = "task2-reference-download/1.0"
BLOCK_SIZE = 8 * 1024 * 1024


def remote_size(url: str) -> int:
    request = Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    for attempt in range(1, 13):
        try:
            with urlopen(request, timeout=120) as response:
                length = response.headers.get("Content-Length")
                if not length or not length.isdigit():
                    raise ValueError(f"missing Content-Length for {url}")
                return int(length)
        except (HTTPError, URLError, TimeoutError, OSError):
            if attempt == 12:
                raise
            time.sleep(min(30, attempt * 2))
    raise RuntimeError("unreachable")


def download_range(
    url: str,
    target: Path,
    start: int,
    end: int,
    label: str,
    print_lock: threading.Lock,
) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = end - start + 1
    for attempt in range(1, 31):
        existing = target.stat().st_size if target.exists() else 0
        if existing == expected:
            return
        if existing > expected:
            target.unlink()
            existing = 0
        request_start = start + existing
        request = Request(
            url,
            headers={
                "User-Agent": USER_AGENT,
                "Range": f"bytes={request_start}-{end}",
            },
        )
        try:
            with urlopen(request, timeout=45) as response:
                if response.status != 206:
                    raise OSError(f"range request returned HTTP {response.status}")
                with target.open("ab") as handle:
                    downloaded = existing
                    next_report = downloaded + 128 * 1024 * 1024
                    while True:
                        block = response.read(BLOCK_SIZE)
                        if not block:
                            break
                        handle.write(block)
                        downloaded += len(block)
                        if downloaded >= next_report:
                            with print_lock:
                                print(
                                    f"REF_PROGRESS ref={label} chunk={target.name} "
                                    f"bytes={downloaded}/{expected}",
                                    flush=True,
                                )
                            next_report = downloaded + 128 * 1024 * 1024
            if target.stat().st_size == expected:
                return
            raise OSError(f"incomplete range: {target.stat().st_size}/{expected}")
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if attempt == 30:
                raise
            with print_lock:
                print(
                    f"REF_RETRY ref={label} chunk={target.name} attempt={attempt} "
                    f"error={type(exc).__name__}",
                    flush=True,
                )
            time.sleep(min(30, attempt * 2))


def download_reference(reference: str, archive_root: Path, workers: int) -> None:
    url = f"{BASE_URL}/{reference}.tar.gz"
    target = archive_root / f"{reference}.tar.gz"
    total = REFERENCE_SIZES.get(reference) or remote_size(url)
    if target.exists() and target.stat().st_size == total:
        print(f"REF_ARCHIVE_REUSED ref={reference} bytes={total}", flush=True)
        return

    chunk_root = archive_root / ".chunks" / reference
    chunk_size = (total + workers - 1) // workers
    parts: list[tuple[Path, int, int]] = []
    for index in range(workers):
        start = index * chunk_size
        if start >= total:
            break
        end = min(total - 1, start + chunk_size - 1)
        parts.append((chunk_root / f"part{index:02d}", start, end))

    if target.exists() and 0 < target.stat().st_size <= parts[0][2] - parts[0][1] + 1:
        parts[0][0].parent.mkdir(parents=True, exist_ok=True)
        if not parts[0][0].exists():
            os.replace(target, parts[0][0])
        else:
            target.unlink()

    print(f"REF_DOWNLOAD ref={reference} bytes={total} chunks={len(parts)}", flush=True)
    print_lock = threading.Lock()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(download_range, url, path, start, end, reference, print_lock)
            for path, start, end in parts
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    assembled = target.with_suffix(target.suffix + ".assembling")
    with assembled.open("wb") as output:
        for part, start, end in parts:
            expected = end - start + 1
            if part.stat().st_size != expected:
                raise OSError(f"bad chunk size for {part}: {part.stat().st_size}/{expected}")
            with part.open("rb") as source:
                while True:
                    block = source.read(BLOCK_SIZE)
                    if not block:
                        break
                    output.write(block)
    if assembled.stat().st_size != total:
        raise OSError(f"assembled size mismatch for {reference}")
    os.replace(assembled, target)
    for part, _, _ in parts:
        part.unlink()
    try:
        chunk_root.rmdir()
    except OSError:
        pass
    print(f"REF_ARCHIVE_READY ref={reference} bytes={total}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Parallel range downloader for official AC reference archives.")
    parser.add_argument("--archive-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.archive_root.mkdir(parents=True, exist_ok=True)
    for reference in REFERENCES:
        download_reference(reference, args.archive_root, max(1, args.workers))
    print(f"REF_DOWNLOADS_COMPLETE count={len(REFERENCES)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
