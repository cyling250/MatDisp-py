import numpy as np
import sympy as sym


class Element:
    def __init__(self, node1=None, node2=None, link_way1=False, link_way2=False, E=sym.Symbol("E"), I=sym.Symbol("I"),
                 A=sym.Symbol("A")):
        super(Element, self).__init__()
        self.nodes = [[node1, link_way1], [node2, link_way2]]  # 单元所包含的节点编号
        self.freedom = [0, 0, 0, 0, 0, 0]  # 定位向量及自由度
        self.E = E  # 单元弹性模量初始化1
        self.I = I  # 单元截面惯性矩初始化1
        self.L = 0  # 单元长度
        self.A = A  # 单元截面面积
        self.k_e = []  # 单元局部刚度矩阵
        self.T = []  # 坐标变换矩阵
        self.restrain = [0, 0, 0, 0, 0, 0]  # 约束条件
        self.load = np.array([0, 0, 0, 0, 0, 0])  # 局部坐标荷载列阵
        self.delta = np.zeros(6)  # 杆端位移
        self.force = []  # 杆端内力

    def add_restrain(self, restrain):
        self.restrain = restrain

    def set_load(self, load):
        # 转化等效节点力
        # Load包含了各种荷载形式，荷载形式可以自行定义，只需要最终作用到self.Load变量上即可
        L = self.L
        if load[0] == 0:
            # 0表示杆上集中荷载，杆上集中荷载的形式[0,Fp,a],Fp是力的大小，a是集中力到杆件正方向左端的距离
            b = self.L - load[2]
            a = load[2]
            fp = load[1]
            self.load = self.load + np.array(
                [0, -fp * b ** 2 * (1 + 2 * a / self.L) / self.L ** 2, -fp * a * b ** 2 / self.L ** 2,
                 0, -fp * a ** 2 * (1 + 2 * b / self.L) / self.L ** 2, fp * a * b ** 2 / self.L ** 2])
        elif load[0] == 1:
            # 1表示杆端弯矩，Load[1]存放杆左端弯矩,Load[2]存放杆右端弯矩
            self.load[2] -= load[1]
            self.load[5] -= load[2]
        elif load[0] == 2:
            # 2表示均布荷载，Load[1]存放均布荷载的大小
            q = load[1]
            self.load = self.load + np.array([0., -q * L / 2., -q * L ** 2 / 12.,
                                              0., -q * L / 2., q * L ** 2 / 12.])
            # print("最初加上的荷载列阵为:", self.Load)
        elif load[0] == 3:
            # 3表示具有等斜度的均布荷载
            q = load[1]
            cosa = self.T[0][0]
            self.load = self.load + np.matmul(self.T, np.array([0., -q * L / 2. * cosa, -q * L ** 2 / 12. * cosa ** 2,
                                                                0., -q * L / 2. * cosa, q * L ** 2 / 12. * cosa ** 2]))

        return
