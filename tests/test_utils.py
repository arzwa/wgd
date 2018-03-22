#!/usr/bin/python3.5
"""
Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be
"""
import unittest
from wgd.utils import translate_cds


class TestFastaMethods(unittest.TestCase):

    def setUp(self):
        self.cds = {'g1': "ATGACAAAGGCTTATTCAACACGTGTCCTCACGTTTCTGATATTGATCTAA",
                    'g2': "ATGGCATTGGCTTATTCAACACGTGTCCTCACGTGTCTGATATTGATCTAA",
                    'g3': "ATGACAAAGGCTTATTCAACACGTGTCCTCACGTTTCTGATATTGTAATAA"}
        self.prt = {'g1': "MTKAYSTRVLTFLILI",
                    'g2': 'MALAYSTRVLTCLILI'}

    def test_translate(self):
        self.assertDictEqual(translate_cds(self.cds), self.prt)


if __name__ == '__main__':
    unittest.main()

