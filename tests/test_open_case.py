"""examples/open_case.py: curve zeros of the denominator, checked against exact circle ranges."""
import os, sys, unittest, contextlib, io

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples"))


class TestOpenCase(unittest.TestCase):
    def test_all_agree(self):
        import open_case
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(open_case.main())


if __name__ == "__main__":
    unittest.main()
