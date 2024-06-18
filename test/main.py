from src.MatDisp import *
import numpy

n1 = Node(0, 0)
n2 = Node(12, 0)
n3 = Node(0, 6)
n4 = Node(12, 6)
s = Structure()
e1 = s.link(n1, n2, A=0.63, I=1 / 12, E=6.94 * 10 ** -3 * 12 * 12)
e2 = s.link(n1, n3, A=0.5, I=1 / 24, E=6.94 * 10 ** -3 * 24 * 6)
e3 = s.link(n2, n4, A=0.5, I=1 / 24, E=6.94 * 10 ** -3 * 24 * 6)
e2.add_restrain([False, False, False, True, True, True])
e3.add_restrain([False, False, False, True, True, True])
e2.set_load([2, 1])

s.size_of_K = s.get_freedom()
s.load_process()
s.get_entire_k()
s.freedom_process()
print("整体刚度矩阵为:\n", s.K)
print("荷载列阵为:", s.load)
s.result = numpy.linalg.solve(s.K, s.load)
print("结果为：", s.result)
s.get_internal_force()
for element in s.elements:
    print(f"{element.index}号单元杆端内力为：{element.force}")
