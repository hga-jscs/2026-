"""Deterministic file discovery and SHA-256 provenance for the task2 pipeline."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable


def filesystem_path(path: Path) -> Path:
    """Return an extended-length Windows path while preserving path identity."""

    if os.name != "nt":
        return path
    absolute = str(path.resolve())
    if absolute.startswith("\\\\?\\"):
        return Path(absolute)
    if absolute.startswith("\\\\"):
        return Path("\\\\?\\UNC\\" + absolute.lstrip("\\"))
    return Path("\\\\?\\" + absolute)


def sha256(path: Path) -> str:
    """Hash one file without loading it completely into memory."""

    digest = hashlib.sha256()
    with filesystem_path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(payload: object) -> str:
    """Hash a JSON-serializable value using a stable representation."""

    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def graph_identity(graph_path: Path) -> tuple[str, str]:
    """Read ``sample`` and ``ampliconN`` from an AA/CoRAL graph filename."""

    suffix = "_graph.txt"
    if not graph_path.name.endswith(suffix):
        raise ValueError(f"not a graph filename: {graph_path.name}")
    prefix = graph_path.name[: -len(suffix)]
    marker = prefix.rfind("_amplicon")
    if marker < 1 or not prefix[marker + 1 :].removeprefix("amplicon").isdigit():
        raise ValueError(f"cannot parse sample/amplicon from {graph_path.name}")
    return prefix[:marker], prefix[marker + 1 :]


def discover_paired_graphs(project_input: Path) -> list[Path]:
    """Discover true graph/cycles pairs and exclude AC-derived graph-like files."""

    graphs: list[Path] = []
    for graph in filesystem_path(project_input).rglob("*_graph.txt"):
        lowered_parts = tuple(part.lower() for part in graph.parts)
        if any(
            part == "bfbarchitect_outputs" or part.endswith("_bfbarchitect_outputs")
            for part in lowered_parts
        ):
            continue
        lowered = graph.as_posix().lower()
        if "bpg_converted" in lowered or "/_classification/files/" in lowered:
            continue
        if graph.name.endswith(("_features_to_graph.txt", "_feature_to_graph.txt")):
            continue
        cycles = graph.with_name(graph.name[: -len("_graph.txt")] + "_cycles.txt")
        if cycles.is_file():
            graphs.append(graph)
    return sorted(graphs, key=lambda path: path.as_posix())


def relative_path(path: Path, root: Path) -> str:
    """Return a stable POSIX-style relative path, including on long Windows paths."""

    return path.relative_to(filesystem_path(root)).as_posix()


def paired_input_snapshot(project_input: Path) -> dict[str, object]:
    """Hash every selected graph/cycles pair and the ordered project snapshot."""

    pairs: list[dict[str, str]] = []
    for graph in discover_paired_graphs(project_input):
        cycles = graph.with_name(graph.name[: -len("_graph.txt")] + "_cycles.txt")
        sample, amplicon = graph_identity(graph)
        pairs.append(
            {
                "sample": sample,
                "amplicon": amplicon,
                "graph_relative_path": relative_path(graph, project_input),
                "graph_sha256": sha256(graph),
                "cycles_relative_path": relative_path(cycles, project_input),
                "cycles_sha256": sha256(cycles),
            }
        )
    compact = [
        {
            "graph_relative_path": row["graph_relative_path"],
            "graph_sha256": row["graph_sha256"],
            "cycles_relative_path": row["cycles_relative_path"],
            "cycles_sha256": row["cycles_sha256"],
        }
        for row in pairs
    ]
    return {
        "pair_count": len(pairs),
        "snapshot_sha256": canonical_sha256(compact),
        "pairs": pairs,
    }


def tree_sha256(root: Path, relative_files: Iterable[str] | None = None) -> str:
    """Hash file names and contents in a source tree, ignoring build residue."""

    root = root.resolve()
    if relative_files is None:
        paths = [
            path
            for path in root.rglob("*")
            if path.is_file()
            and not any(part in {".git", "__pycache__", ".pytest_cache"} for part in path.parts)
            and path.suffix.lower() not in {".pyc", ".pyo"}
        ]
    else:
        paths = [root / relative for relative in relative_files]
    records = []
    for path in sorted(paths, key=lambda item: item.relative_to(root).as_posix()):
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256(path),
            }
        )
    return canonical_sha256(records)


def ac_runtime_tree_sha256(ac_source: Path) -> str:
    """Hash the AC runtime code, configuration and bundled resources."""

    top_level = {
        "amplicon_classifier.py",
        "feature_similarity.py",
        "make_results_table.py",
        "make_input.sh",
    }
    files = [path for path in ac_source.iterdir() if path.is_file() and path.name in top_level]
    files.extend(
        path
        for path in (ac_source / "ampclasslib").rglob("*")
        if path.is_file() and "__pycache__" not in path.parts and path.suffix.lower() != ".pyc"
    )
    records = [
        {
            "path": path.relative_to(ac_source).as_posix(),
            "sha256": sha256(path),
        }
        for path in sorted(files, key=lambda item: item.relative_to(ac_source).as_posix())
    ]
    return canonical_sha256(records)
