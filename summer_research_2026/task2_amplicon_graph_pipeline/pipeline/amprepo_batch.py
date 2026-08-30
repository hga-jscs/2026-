"""下载、校验并安全解包 AmpliconRepository 全部公开项目。

逐步逻辑：读取 REST 项目清单；按项目断点续传归档；校验字节数和 SHA-256；
拒绝绝对路径/``..`` 成员；只提取主 graph/cycles/summary；最后写数据集清单。
脚本不会把 BFBArchitect 派生输出误当作主输入。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


API = "https://ampliconrepository.org/api/v1"
USER_AGENT = "task2-amprepo-run/1.0"
CHUNK_SIZE = 8 * 1024 * 1024
PROGRESS_BYTES = 128 * 1024 * 1024


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def request_json(url: str) -> object:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urlopen(request, timeout=120) as response:
        return json.load(response)


def wait_with_progress(seconds: int, reason: str) -> None:
    remaining = max(1, seconds)
    while remaining > 0:
        step = min(30, remaining)
        print(f"WAIT seconds={remaining} reason={reason}", flush=True)
        time.sleep(step)
        remaining -= step


def content_total(response, existing: int) -> int | None:
    content_range = response.headers.get("Content-Range", "")
    match = re.search(r"/(\d+)$", content_range)
    if match:
        return int(match.group(1))
    length = response.headers.get("Content-Length")
    if length and length.isdigit():
        return existing + int(length)
    return None


def download_archive(project_id: str, target: Path) -> tuple[int, str]:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, 13):
        existing = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": USER_AGENT}
        if existing:
            headers["Range"] = f"bytes={existing}-"
        request = Request(f"{API}/projects/{project_id}/download/", headers=headers)
        try:
            with urlopen(request, timeout=180) as response:
                status = getattr(response, "status", 200)
                if existing and status != 206:
                    partial.unlink(missing_ok=True)
                    existing = 0
                mode = "ab" if existing and status == 206 else "wb"
                total = content_total(response, existing)
                written = existing
                next_report = written + PROGRESS_BYTES
                with partial.open(mode) as handle:
                    while True:
                        block = response.read(CHUNK_SIZE)
                        if not block:
                            break
                        handle.write(block)
                        written += len(block)
                        if written >= next_report:
                            total_text = str(total) if total is not None else "unknown"
                            print(
                                f"DOWNLOAD_PROGRESS id={project_id} bytes={written} total={total_text}",
                                flush=True,
                            )
                            next_report = written + PROGRESS_BYTES
                if total is not None and written != total:
                    raise OSError(f"incomplete download: {written} of {total} bytes")
            os.replace(partial, target)
            digest = hashlib.sha256()
            with target.open("rb") as handle:
                for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
                    digest.update(block)
            return target.stat().st_size, digest.hexdigest()
        except HTTPError as exc:
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            if exc.code == 429 or 500 <= exc.code < 600:
                delay = int(retry_after) + 2 if retry_after and retry_after.isdigit() else min(60, 5 * attempt)
                wait_with_progress(delay, f"HTTP_{exc.code}_id_{project_id}")
                continue
            raise
        except (URLError, TimeoutError, OSError) as exc:
            if attempt == 12:
                raise
            wait_with_progress(min(60, 5 * attempt), f"{type(exc).__name__}_id_{project_id}")
    raise RuntimeError(f"download attempts exhausted for {project_id}")


def safe_member_path(member_name: str) -> PurePosixPath:
    path = PurePosixPath(member_name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe archive member: {member_name}")
    return path


def wanted_member(member_name: str) -> bool:
    path = PurePosixPath(member_name)
    lowered = str(path).lower()
    lowered_parts = tuple(part.lower() for part in path.parts)
    if any(part == "bfbarchitect_outputs" or part.endswith("_bfbarchitect_outputs") for part in lowered_parts):
        return False
    if any(part.startswith("._") for part in path.parts):
        return False
    if "bpg_converted" in lowered or "_classification/files/" in lowered:
        return False
    name = path.name
    if name.endswith("_annotated_cycles.txt"):
        return False
    if name.endswith("_features_to_graph.txt") or name.endswith("_feature_to_graph.txt"):
        return False
    return name.endswith("_graph.txt") or name.endswith("_cycles.txt")


def filesystem_path(path: Path) -> Path:
    """Use the Windows extended-length form for long sample names."""

    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return Path("\\\\?\\" + resolved)
    return Path(resolved)


def extract_inputs(archive: Path, destination: Path) -> tuple[int, int, int]:
    destination.mkdir(parents=True, exist_ok=True)
    graph_prefixes: set[str] = set()
    cycle_prefixes: set[str] = set()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            if not member.isfile() or not wanted_member(member.name):
                continue
            relative = safe_member_path(member.name)
            output = filesystem_path(destination.joinpath(*relative.parts))
            output.parent.mkdir(parents=True, exist_ok=True)
            source = bundle.extractfile(member)
            if source is None:
                raise OSError(f"could not read archive member: {member.name}")
            with source, output.open("wb") as sink:
                shutil.copyfileobj(source, sink, CHUNK_SIZE)
            normalized = str(relative).replace("\\", "/")
            if normalized.endswith("_graph.txt"):
                graph_prefixes.add(normalized[: -len("_graph.txt")])
            elif normalized.endswith("_cycles.txt"):
                cycle_prefixes.add(normalized[: -len("_cycles.txt")])
    return len(graph_prefixes), len(cycle_prefixes), len(graph_prefixes & cycle_prefixes)


def write_manifests(rows: list[dict[str, object]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "dataset_manifest.json"
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    columns = [
        "project_id",
        "project_name",
        "reference_genome",
        "reported_sample_count",
        "archive_bytes",
        "sha256",
        "graph_count",
        "cycles_count",
        "paired_amplicon_count",
        "status",
        "checked_at_utc",
        "error",
    ]
    with (output_dir / "dataset_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def existing_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Download and audit every public AmpliconRepository archive.")
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--manifest-dir", required=True, type=Path)
    args = parser.parse_args()

    archives_dir = args.data_root / "archives"
    extracted_dir = args.data_root / "extracted"
    projects = request_json(f"{API}/projects/")
    if not isinstance(projects, list):
        raise TypeError("projects endpoint did not return a list")
    projects.sort(key=lambda item: str(item.get("project_name", "")).casefold())
    rows: list[dict[str, object]] = []
    print(f"PUBLIC_PROJECTS count={len(projects)}", flush=True)

    for index, project in enumerate(projects, start=1):
        project_id = str(project["id"])
        archive = archives_dir / f"{project_id}.tar.gz"
        destination = extracted_dir / project_id
        row: dict[str, object] = {
            "project_id": project_id,
            "project_name": project.get("project_name", ""),
            "reference_genome": project.get("reference_genome", ""),
            "reported_sample_count": project.get("sample_count", 0),
            "archive_bytes": "",
            "sha256": "",
            "graph_count": 0,
            "cycles_count": 0,
            "paired_amplicon_count": 0,
            "status": "PENDING",
            "checked_at_utc": utc_now(),
            "error": "",
        }
        print(
            f"PROJECT {index}/{len(projects)} id={project_id} name={project.get('project_name', '')}",
            flush=True,
        )
        try:
            if archive.exists() and tarfile.is_tarfile(archive):
                size = archive.stat().st_size
                digest = existing_digest(archive)
                print(f"DOWNLOAD_REUSED id={project_id} bytes={size}", flush=True)
            else:
                archive.unlink(missing_ok=True)
                size, digest = download_archive(project_id, archive)
                print(f"DOWNLOAD_DONE id={project_id} bytes={size}", flush=True)
            graph_count, cycles_count, pair_count = extract_inputs(archive, destination)
            row.update(
                {
                    "archive_bytes": size,
                    "sha256": digest,
                    "graph_count": graph_count,
                    "cycles_count": cycles_count,
                    "paired_amplicon_count": pair_count,
                    "status": "READY" if pair_count else "NO_COMPATIBLE_PAIR",
                    "checked_at_utc": utc_now(),
                }
            )
            print(
                f"INPUTS id={project_id} graphs={graph_count} cycles={cycles_count} pairs={pair_count}",
                flush=True,
            )
        except Exception as exc:  # pragma: no cover - boundary audit
            row["status"] = "ERROR"
            row["error"] = f"{type(exc).__name__}: {exc}"
            row["checked_at_utc"] = utc_now()
            print(f"PROJECT_ERROR id={project_id} error={row['error']}", file=sys.stderr, flush=True)
        rows.append(row)
        write_manifests(rows, args.manifest_dir)

    ready = sum(row["status"] == "READY" for row in rows)
    errors = sum(row["status"] == "ERROR" for row in rows)
    pairs = sum(int(row["paired_amplicon_count"]) for row in rows)
    print(f"AUDIT_COMPLETE projects={len(rows)} ready={ready} errors={errors} pairs={pairs}", flush=True)
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
