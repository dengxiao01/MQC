import sys
from pathlib import Path
# main_directory = Path(__file__).resolve().parent.parent
# sys.path.append(str(main_directory))
FILE = Path(__file__).resolve()
ROOT = FILE.parents[0]
print('FILE:',FILE,'ROOT:', ROOT)
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT)) # add ROOT to PATH

from old_main import run

def main():
    # 解析命令行参数
    if len(sys.argv) != 3:
        print("Usage: mqc <model_path> <output_path>")
        sys.exit(1)

    arg1 = sys.argv[1]
    arg2 = sys.argv[2]

    # 调用 main.py 中的函数，执行程序逻辑
    run(arg1, arg2)

if __name__ == "__main__":
    main()
