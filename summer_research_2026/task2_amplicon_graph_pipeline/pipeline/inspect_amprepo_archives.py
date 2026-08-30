"""只读审计 AmpliconRepository 公开项目及下载归档的 HTTP 元数据。"""

from __future__ import annotations

import concurrent.futures
import json
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import Request, urlopen


API = "https://ampliconrepository.org/api/v1"


def inspect(project: dict[str, object]) -> dict[str, object]:
    project_id = str(project["id"])
    url = f"{API}/projects/{project_id}/download/"
    result = {
        "id": project_id,
        "project_name": project.get("project_name", ""),
        "sample_count": project.get("sample_count", 0),
        "reference_genome": project.get("reference_genome", ""),
        "reconstruction_tools": project.get("reconstruction_tools", ""),
        "status": None,
        "archive_bytes": None,
        "archive_host": None,
        "error": "",
    }
    try:
        request = Request(url, method="HEAD", headers={"User-Agent": "task2-data-audit/1.0"})
        with urlopen(request, timeout=90) as response:
            result["status"] = response.status
            length = response.headers.get("Content-Length")
            result["archive_host"] = urlparse(response.url).hostname
        result["archive_bytes"] = int(length) if length and length.isdigit() else None
    except HTTPError as exc:
        result["status"] = exc.code
        result["error"] = exc.reason
    except Exception as exc:  # pragma: no cover - network boundary
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    request = Request(
        f"{API}/projects/",
        headers={"Accept": "application/json", "User-Agent": "task2-data-audit/1.0"},
    )
    with urlopen(request, timeout=90) as response:
        projects = json.load(response)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(inspect, projects))
    rows.sort(key=lambda row: str(row["project_name"]).casefold())
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
