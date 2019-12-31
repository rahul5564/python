from __future__ import division, absolute_import, print_function

import pytest

import numpy as np
from numpy.testing import assert_array_equal, assert_equal
from numpy.f2py.crackfortran import markinnerspaces
from . import util


class TestNoSpace(util.F2PyTest):
    # issue gh-15035: add handling for endsubroutine, endfunction with no space
    # between "end" and the block name
    code = """
        subroutine subb(k)
          real(8), intent(inout) :: k(:)
          k=k+1
        endsubroutine

        subroutine subc(w,k)
          real(8), intent(in) :: w(:)
          real(8), intent(out) :: k(size(w))
          k=w+1
        endsubroutine

        function t0(value)
          character value
          character t0
          t0 = value
        endfunction
    """

    def test_module(self):
        k = np.array([1, 2, 3], dtype=np.float64)
        w = np.array([1, 2, 3], dtype=np.float64)
        self.module.subb(k)
        assert_array_equal(k, w + 1)
        self.module.subc([w, k])
        assert_array_equal(k, w + 1)
        assert self.module.t0(23) == b'2'


class TestMarkinnerspaces():
    # issue #14118: markinnerspaces does not handle multiple quotations

    def test_do_not_touch_normal_spaces(self):
        test_list = ["a ", " a", "a b c", "'abcdefghij'"]
        for i in test_list:
            assert_equal(markinnerspaces(i), i)

    def test_one_relevant_space(self):
        assert_equal(markinnerspaces("a 'b c' \\\' \\\'"), "a 'b@_@c' \\' \\'")

    def test_multiple_relevant_spaces(self):
        assert_equal(markinnerspaces("a 'b c' 'd e'"), "a 'b@_@c' 'd@_@e'")
