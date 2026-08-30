"""允许使用 ``python -m original_graph_lwcn`` 启动命令行。"""

from .command_line import main


if __name__ == "__main__":
    raise SystemExit(main())

