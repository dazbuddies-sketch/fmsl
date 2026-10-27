import sys
from .compiler import run


def main():
    if len(sys.argv) != 2:
        print("usage: fmsl PROGRAM.fms")
        sys.exit(2)
    with open(sys.argv[1]) as f:
        run(f.read())


if __name__ == "__main__":
    main()
