__author__ = 'lustelzl'

from kinetics.kinetics import *
from numpy.testing import assert_array_equal


class TestKineticsRemd:

    def __init__(self):
        self.state_ar = np.array([[1, 0, 0], [1, 0, 1], [0, 0, 1]])
        self.temp_ar = np.array([[0, 1, 2], [2, 0, 1], [2, 1, 0]])
        self.time_temp_ar = np.column_stack(([0, 1, 2], self.temp_ar))

    def loop_reorder_indices(self, state_ar, temp_ar):
        temp_dict = {}
        for t in temp_ar[0, :]:
            out_l = []
            for i, row in enumerate(state_ar):
                value = row[temp_ar[i, :] == t]
                out_l.append(list(value)[0])
            temp_dict[t] = out_l
        return temp_dict

    def test_reorder_indices(self):
        re_order_dict = self.loop_reorder_indices(self.state_ar, self.temp_ar)
        re_order = reorder_temperature_indices(self.state_ar, self.time_temp_ar, arrange_up_to=3)

        print "state array"
        print self.state_ar
        print "temperature array"
        print self.temp_ar

        print "re ordering indices"
        print re_order
        print "re ordering via dictionary"
        print re_order_dict
        assert_array_equal(re_order, np.vstack(re_order_dict.values()).T)